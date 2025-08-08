import argparse
import os
import sys
import numpy as np
import hashlib
import time
import json
import shutil
from datetime import datetime
from Bio.PDB import PDBParser, PDBIO, Select
from Bio.PDB.Atom import Atom
from typing import Dict, List, Tuple
import logging

class NumpyEncoder(json.JSONEncoder):
    """Numpy 타입을 JSON 직렬화 가능하게 변환"""
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        elif isinstance(obj, np.floating):
            return float(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        return super(NumpyEncoder, self).default(obj)

def log_distance_analysis(initial_distance, final_distance, improvement_threshold=0.1, logger=None):
    """거리 변화 상세 분석 및 로깅"""
    if logger is None:
        return
    
    distance_change = initial_distance - final_distance
    improvement_percentage = (distance_change / initial_distance) * 100 if initial_distance > 0 else 0
    
    logger.info(f"📊 거리 변화 분석:")
    logger.info(f"   초기 거리: {initial_distance:.3f} Å")
    logger.info(f"   최종 거리: {final_distance:.3f} Å")
    logger.info(f"   절대 변화: {distance_change:.3f} Å")
    logger.info(f"   상대 변화: {improvement_percentage:.2f}%")
    
    if distance_change > improvement_threshold:
        logger.info(f"   ✅ 개선됨 ({distance_change:.3f}Å > {improvement_threshold}Å)")
        return True
    elif abs(distance_change) <= 0.01:
        logger.info(f"   ⚪ 거의 변화 없음 (±0.01Å 이내)")
        return False
    else:
        logger.info(f"   ❌ 개선 부족 ({distance_change:.3f}Å ≤ {improvement_threshold}Å)")
        return False

def log_retry_status(iteration, retry_count, max_retries, logger=None):
    """재시도 상태 로깅"""
    if logger is None:
        return
    
    progress = (retry_count / max_retries) * 100
    remaining = max_retries - retry_count
    
    logger.info(f"🔄 재시도 상태:")
    logger.info(f"   현재 iteration: {iteration}")
    logger.info(f"   재시도 횟수: {retry_count}/{max_retries} ({progress:.1f}%)")
    logger.info(f"   남은 시도: {remaining}회")
    
    if progress > 80:
        logger.warning(f"   ⚠️  재시도 한계 근접 ({progress:.1f}%)")
    elif progress > 50:
        logger.info(f"   🟡 재시도 중간 단계 ({progress:.1f}%)")

def create_iteration_summary(iteration, retry_count, initial_distance, final_distance, 
                           rmsd_vs_golden, simulation_time_ns, samples_info, logger=None):
    """Iteration 요약 정보 생성 및 로깅"""
    distance_improvement = initial_distance - final_distance
    
    summary = {
        "iteration": iteration,
        "retry_count": retry_count,
        "distances": {
            "initial": round(initial_distance, 3),
            "final": round(final_distance, 3),
            "improvement": round(distance_improvement, 3),
            "improvement_percentage": round((distance_improvement/initial_distance)*100, 2) if initial_distance > 0 else 0
        },
        "rmsd_vs_golden": round(rmsd_vs_golden, 3),
        "simulation_info": {
            "time_ns": simulation_time_ns,
            "time_ps": simulation_time_ns * 1000,
            "successful_samples": samples_info.get("successful", 0),
            "failed_samples": samples_info.get("failed", 0),
            "total_samples": samples_info.get("total", 0)
        },
        "status": "success" if distance_improvement > 0.1 else "retry_needed"
    }
    
    if logger:
        logger.info(f"📋 Iteration {iteration} 요약 (재시도 {retry_count}):")
        logger.info(f"   거리: {initial_distance:.3f}Å → {final_distance:.3f}Å ({distance_improvement:+.3f}Å)")
        logger.info(f"   RMSD vs Golden: {rmsd_vs_golden:.3f}Å")
        logger.info(f"   샘플 성공률: {samples_info.get('successful', 0)}/{samples_info.get('total', 0)}")
        logger.info(f"   시뮬레이션: {simulation_time_ns}ns ({simulation_time_ns*1000:.0f}ps)")
        logger.info(f"   상태: {'✅ 채택' if summary['status'] == 'success' else '🔄 재시도 필요'}")
    
    return summary

def save_retry_log(job_output_dir, iteration_summaries, logger=None):
    """재시도 로그 파일 저장"""
    try:
        retry_log_file = os.path.join(job_output_dir, "retry_analysis.json")
        
        retry_analysis = {
            "total_iterations": len([s for s in iteration_summaries if s.get("retry_count", 0) == 0]),
            "total_retries": sum([s.get("retry_count", 0) for s in iteration_summaries]),
            "retry_statistics": {},
            "iteration_details": iteration_summaries
        }
        
        # 재시도 통계 계산
        retry_counts = [s.get("retry_count", 0) for s in iteration_summaries]
        if retry_counts:
            retry_analysis["retry_statistics"] = {
                "min_retries": min(retry_counts),
                "max_retries": max(retry_counts),
                "avg_retries": sum(retry_counts) / len(retry_counts),
                "iterations_with_retries": len([r for r in retry_counts if r > 0])
            }
        
        with open(retry_log_file, 'w', encoding='utf-8') as f:
            json.dump(retry_analysis, f, indent=2, ensure_ascii=False)
            
        if logger:
            logger.info(f"재시도 분석 로그 저장: {retry_log_file}")
            logger.info(f"총 재시도: {retry_analysis['total_retries']}회")
            
    except Exception as e:
        if logger:
            logger.error(f"재시도 로그 저장 실패: {e}")

# robust_sumd_master.py의 메인 루프에서 사용할 개선된 거리 확인 함수
def check_distance_improvement(current_distance, final_distance, improvement_threshold=0.1, logger=None):
    """거리 개선 여부를 체크하고 상세 정보를 로깅"""
    distance_change = current_distance - final_distance
    
    # 로깅
    is_improved = log_distance_analysis(current_distance, final_distance, improvement_threshold, logger)
    
    # 추가 컨텍스트 정보
    if logger:
        if distance_change < -0.05:  # 거리가 더 멀어진 경우
            logger.warning(f"   ⚠️  거리가 오히려 증가했습니다 ({distance_change:.3f}Å)")
        elif 0 <= distance_change <= improvement_threshold:
            logger.info(f"   🟡 소폭 개선되었으나 임계값 미달")
    
    return is_improved, distance_change

# 사용 예시를 위한 함수
def print_usage_example():
    """새로운 기능의 사용 예시"""
    print("""
=== 새로운 거리 기반 반복 SuMD 사용법 ===

기본 사용법:
./robust_sumd_master.sh input.pdb golden_standard.pdb A B

300ps 고정 시뮬레이션:
./robust_sumd_master.sh input.pdb golden_standard.pdb A B 0.3

재시도 횟수 조정:
./robust_sumd_master.sh input.pdb golden_standard.pdb A B 0.3 5.0 1.5 10 5 /output "" false 50

주요 개선사항:
✓ 300ps 고정으로 빠른 반복
✓ 거리 개선 기반 진행 (0.1Å 임계값)
✓ 자동 재시도 (최대 100회)
✓ 재시도 한계 초과시 원본으로 리셋
✓ 상세한 거리 변화 모니터링
✓ Golden Standard 기반 절대적 수렴 판정

로그 파일:
- 메인 로그: sumd_[JOB_ID].log  
- 재시도 분석: retry_analysis.json
- 결과 요약: simulation_summary.txt
""")