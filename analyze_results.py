#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Simple SuMD 결과 분석 유틸리티
"""

import os
import sys
import json
import glob,re
import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime

def load_results(output_dir):
    """결과 파일들 로드 - final_results.json이 없어도 진행 중인 결과 분석 가능"""
    final_results_file = os.path.join(output_dir, "final_results.json")
    
    # final_results.json 로드 시도
    final_results = None
    if os.path.exists(final_results_file):
        try:
            with open(final_results_file) as f:
                final_results = json.load(f)
            print("완료된 시뮬레이션 결과를 분석합니다.")
        except Exception as e:
            print(f"final_results.json 읽기 실패: {e}")
    else:
        print("진행 중인 시뮬레이션 결과를 분석합니다.")
    
    # Iteration별 상세 결과 로드
    iteration_files = glob.glob(os.path.join(output_dir, "iteration_*_summary.json"))

    def extract_iteration_number(filename):
        """파일명에서 iteration 번호 추출"""
        match = re.search(r'iteration_(\d+)_summary\.json', os.path.basename(filename))
        return int(match.group(1)) if match else 0
    
    # 숫자 순서로 정렬
    iteration_files_sorted = sorted(iteration_files, key=extract_iteration_number)

    iteration_results = []
    
    for file in iteration_files_sorted:
        try:
            with open(file) as f:
                iteration_results.append(json.load(f))
        except Exception as e:
            print(f"iteration 파일 읽기 실패 {file}: {e}")
    
    # Attempt별 상세 결과 로드
    attempt_files=[]
    for i in range(1, 9999):
        iteration_attempts = os.path.join(output_dir, f"iteration_{i}", f"iteration_{i}_attempt_1.json")
        if os.path.exists(iteration_attempts):
            attempt_files.extend(glob.glob(os.path.join(output_dir, f"iteration_{i}", f"iteration_{i}_attempt_*.json")))
    
    attempt_results = []
    for file in sorted(attempt_files):
        try:
            with open(file) as f:
                attempt_results.append(json.load(f))
        except Exception as e:
            print(f"attempt 파일 읽기 실패 {file}: {e}")
    
    # final_results가 없으면 iteration 결과로부터 생성
    if final_results is None:
        final_results = create_final_results_from_iterations(output_dir, iteration_results)
    
    return final_results, iteration_results, attempt_results

def create_final_results_from_iterations(output_dir, iteration_results):
    """iteration 결과들로부터 final_results 정보 생성"""
    if not iteration_results:
        return {
            "input_pdb": "unknown",
            "chains": ["unknown", "unknown"],
            "start_time": "unknown",
            "end_time": "unknown",
            "total_iterations": 0,
            "successful_iterations": 0,
            "status": "진행 중"
        }
    
    # 첫 번째와 마지막 iteration으로부터 정보 추출
    first_iter = iteration_results[0]
    last_iter = iteration_results[-1]
    
    successful_iterations = len([r for r in iteration_results if r.get('success', False)])
    
    # 기본 정보 생성
    final_results = {
        "input_pdb": "unknown",
        "chains": ["unknown", "unknown"],
        "start_time": "unknown",
        "end_time": "진행 중",
        "total_iterations": len(iteration_results),
        "successful_iterations": successful_iterations,
        "status": "진행 중",
        "long_md_executed": any(r.get('long_md', False) for r in iteration_results)
    }
    
    # 시작 시간 추정 (첫 번째 iteration 디렉토리의 생성 시간)
    try:
        first_iter_dir = os.path.join(output_dir, f"iteration_{first_iter['iteration']}")
        if os.path.exists(first_iter_dir):
            start_time = os.path.getctime(first_iter_dir)
            final_results["start_time"] = datetime.fromtimestamp(start_time).isoformat()

        last_iter_dir = os.path.join(output_dir, f"iteration_{last_iter['iteration']}")
        if os.path.exists(last_iter_dir):
            last_time = os.path.getctime(last_iter_dir)
            final_results["end_time"] = datetime.fromtimestamp(last_time).isoformat()
    except:
        pass
    
    return final_results

def calculate_additional_stats(final_results, attempt_results):
    """추가 통계 계산"""
    stats = {}
    
    # 총 attempt 수 계산
    stats['total_attempts'] = len(attempt_results)
    
    # 평균 attempt 수 (iteration당)
    if final_results['total_iterations'] > 0:
        stats['avg_attempts_per_iteration'] = stats['total_attempts'] / final_results['total_iterations']
    else:
        stats['avg_attempts_per_iteration'] = 0
    
    # 평균 attempt 시간 계산
    stats['avg_attempt_time'] = "계산 불가"
    if final_results['start_time'] != "unknown" and stats['total_attempts'] > 0:
        try:
            start_time = datetime.fromisoformat(final_results['start_time'])
            end_time = datetime.fromisoformat(final_results['end_time'])
            total_duration = end_time - start_time
            total_seconds = total_duration.total_seconds()
            avg_seconds = total_seconds / stats['total_attempts']
            stats['avg_attempt_time'] = f"{avg_seconds:.1f}초"
        except:
            stats['avg_attempt_time'] = "계산 실패"
    
    return stats


def analyze_top_structures(attempt_results):
    """Top 구조들 분석 - RMSD 정보 포함"""
    info = {
        'has_top_structures': False,
        'attempts_with_top': 0,
        'total_top_structures': 0,
        'best_distances': [],
        'best_rmsd': [],
        'best_combined_scores': [],
        'all_top_structures': [],
        'has_rmsd_data': False
    }
    
    for attempt in attempt_results:
        top_structures = attempt.get('top_structures', [])
        if top_structures:
            info['has_top_structures'] = True
            info['attempts_with_top'] += 1
            info['total_top_structures'] += len(top_structures)
            info['all_top_structures'].extend(top_structures)
            
            # 각 attempt의 최고값들 수집
            if top_structures:
                best_distance = min([s['distance'] for s in top_structures])
                info['best_distances'].append(best_distance)
                
                # RMSD 데이터 확인 및 수집
                rmsd_values = [s.get('rmsd') for s in top_structures if s.get('rmsd') is not None]
                if rmsd_values:
                    info['has_rmsd_data'] = True
                    best_rmsd = min(rmsd_values)
                    info['best_rmsd'].append(best_rmsd)
                
                # 조합점수 수집
                combined_scores = [s.get('combined_score') for s in top_structures if s.get('combined_score') is not None]
                if combined_scores:
                    best_combined = min(combined_scores)
                    info['best_combined_scores'].append(best_combined)
    
    return info

def print_top_structures_details(attempt_results):
    """Top 구조 상세 정보 출력 - RMSD 및 조합점수 포함"""
    print("\n" + "=" * 60)
    print("Top 구조 상세 분석 (거리 + RMSD 기반)")
    print("=" * 60)
    
    top_info = analyze_top_structures(attempt_results)
    
    if not top_info['has_top_structures']:
        print("저장된 Top 구조가 없습니다.")
        return
    
    print(f"Top 구조 보유 attempts: {top_info['attempts_with_top']}개")
    print(f"총 Top 구조 수: {top_info['total_top_structures']}개")
    print(f"RMSD 데이터 포함: {'예' if top_info['has_rmsd_data'] else '아니오'}")
    
    # 전체 Top 구조들을 조합점수 또는 거리 순으로 정렬
    all_tops = top_info['all_top_structures']
    
    # 조합점수가 있으면 조합점수로, 없으면 거리로 정렬
    if any('combined_score' in s for s in all_tops):
        all_tops_sorted = sorted([s for s in all_tops if 'combined_score' in s], 
                               key=lambda x: x['combined_score'])
        sort_key = "조합점수"
    else:
        all_tops_sorted = sorted(all_tops, key=lambda x: x['distance'])
        sort_key = "거리"
    
    print(f"\n전체 Top 구조 중 최고 상위 10개 ({sort_key} 기준):")
    print("-" * 80)
    if top_info['has_rmsd_data']:
        print(f"{'순위':<4} {'Frame':<6} {'거리(Å)':<8} {'RMSD(Å)':<9} {'조합점수':<10} {'파일명'}")
    else:
        print(f"{'순위':<4} {'Frame':<6} {'거리(Å)':<8} {'파일명'}")
    print("-" * 80)
    
    for i, struct in enumerate(all_tops_sorted[:10]):
        if top_info['has_rmsd_data']:
            rmsd_str = f"{struct.get('rmsd', 0):<9.2f}" if 'rmsd' in struct else f"{'N/A':<9}"
            score_str = f"{struct.get('combined_score', 0):<10.4f}" if 'combined_score' in struct else f"{'N/A':<10}"
            print(f"{i+1:<4} {struct['frame']:<6} {struct['distance']:<8.2f} "
                  f"{rmsd_str} {score_str} {struct.get('filename', 'N/A')}")
        else:
            print(f"{i+1:<4} {struct['frame']:<6} {struct['distance']:<8.2f} {struct.get('filename', 'N/A')}")
    
    if len(all_tops_sorted) > 10:
        print(f"... 외 {len(all_tops_sorted)-10}개 구조")
    
    # 통계 정보
    distances = [s['distance'] for s in all_tops]
    print(f"\n통계 요약:")
    print(f"거리 통계:")
    print(f"  최소: {min(distances):.2f} Å")
    print(f"  최대: {max(distances):.2f} Å") 
    print(f"  평균: {np.mean(distances):.2f} Å")
    print(f"  표준편차: {np.std(distances):.2f} Å")
    
    if top_info['has_rmsd_data']:
        rmsd_values = [s['rmsd'] for s in all_tops if 'rmsd' in s and s['rmsd'] is not None]
        if rmsd_values:
            print(f"RMSD 통계:")
            print(f"  최소: {min(rmsd_values):.2f} Å")
            print(f"  최대: {max(rmsd_values):.2f} Å")
            print(f"  평균: {np.mean(rmsd_values):.2f} Å")
            print(f"  표준편차: {np.std(rmsd_values):.2f} Å")
    
    if top_info['best_combined_scores']:
        print(f"조합점수 통계:")
        print(f"  최소: {min(top_info['best_combined_scores']):.4f}")
        print(f"  최대: {max(top_info['best_combined_scores']):.4f}")
        print(f"  평균: {np.mean(top_info['best_combined_scores']):.4f}")
    
    # Frame 분포
    frames = [s['frame'] for s in all_tops]
    print(f"\nFrame 분포:")
    print(f"  최소 Frame: {min(frames)}")
    print(f"  최대 Frame: {max(frames)}")
    print(f"  평균 Frame: {np.mean(frames):.1f}")

def print_summary(final_results, iteration_results, attempt_results):
    """결과 요약 출력 - Top 구조 정보 포함"""
    print("=" * 50)
    if final_results.get("status") == "진행 중":
        print("진행 중인 Simple SuMD 결과 분석")
    else:
        print("Simple SuMD 결과 요약")
    print("=" * 50)
    
    # 추가 통계 계산
    stats = calculate_additional_stats(final_results, attempt_results)
    
    # 기본 정보
    print(f"입력 PDB: {final_results['input_pdb']}")
    print(f"체인: {' - '.join(final_results['chains'])}")
    print(f"시작 시간: {final_results['start_time']}")
    
    if final_results.get("status") == "진행 중":
        print(f"상태: 진행 중")
        print(f"현재까지 iterations: {final_results['total_iterations']}")
    else:
        print(f"종료 시간: {final_results['end_time']}")
        # 실행 시간 계산
        if final_results['start_time'] != "unknown" and final_results['end_time'] != "진행 중":
            try:
                start_time = datetime.fromisoformat(final_results['start_time'])
                end_time = datetime.fromisoformat(final_results['end_time'])
                duration = end_time - start_time
                print(f"실행 시간: {duration}")
            except:
                pass
        print(f"총 iterations: {final_results['total_iterations']}")
    
    print(f"성공한 iterations: {final_results['successful_iterations']}")
    if final_results['total_iterations'] > 0:
        print(f"성공률: {final_results['successful_iterations']/final_results['total_iterations']*100:.1f}%")
    
    # 추가 통계 출력
    print(f"\nAttempt 통계:")
    print(f"총 attempts: {stats['total_attempts']}")
    print(f"평균 attempts/iteration: {stats['avg_attempts_per_iteration']:.1f}")
    print(f"평균 attempt 시간: {stats['avg_attempt_time']}")
    
    # 긴 MD 통계 추가
    long_md_attempts = [a for a in attempt_results if a.get('long_md_executed', False)]
    close_contact_attempts = [a for a in attempt_results if a.get('close_contact_detected', False)]
    
    print(f"\n긴 MD 통계:")
    print(f"근접 접촉 감지: {len(close_contact_attempts)}회")
    print(f"긴 MD 실행: {len(long_md_attempts)}회")
    if close_contact_attempts:
        min_distances = [a.get('min_distance', float('inf')) for a in close_contact_attempts if 'min_distance' in a]
        if min_distances:
            print(f"최소 거리: {min(min_distances):.2f} Å")
    
    # Top 구조 정보 추가
    top_structures_info = analyze_top_structures(attempt_results)
    if top_structures_info['has_top_structures']:
        print(f"\n최종 Top 구조 정보:")
        print(f"Top 구조가 있는 attempts: {top_structures_info['attempts_with_top']}회")
        print(f"총 Top 구조 수: {top_structures_info['total_top_structures']}개")
        if top_structures_info['best_distances']:
            print(f"최고 거리 (Top 1들 중): {min(top_structures_info['best_distances']):.2f} Å")
            print(f"평균 최고 거리: {np.mean(top_structures_info['best_distances']):.2f} Å")
    
    print("\nIteration별 결과:")
    print("-" * 30)
    
    for iter_result in iteration_results:
        status = "성공" if iter_result.get('success', False) else "실패"
        attempts = iter_result.get('attempts_used', 'unknown')
        print(f"Iteration {iter_result.get('iteration', 'unknown')}: {status} (Attempts: {attempts})")
        
        if iter_result.get('success', False) and iter_result.get('final_result'):
            result = iter_result['final_result']
            print(f"  - 기울기: {result.get('slope', 'N/A'):.6f}")
            print(f"  - 초기 거리: {result.get('initial_distance', 'N/A'):.2f} Å")
            print(f"  - 최종 거리: {result.get('final_distance', 'N/A'):.2f} Å")
            print(f"  - 최소 거리: {result.get('min_distance', 'N/A'):.2f} Å")
            if result.get('long_md_executed', False):
                print(f"  - 긴 MD 실행됨 ✓")
                # Top 구조 정보 표시
                top_structures = result.get('top_structures', [])
                if top_structures:
                    print(f"  - Top {len(top_structures)} 구조 저장됨:")
                    for i, struct in enumerate(top_structures[:3]):  # 상위 3개만 표시
                        print(f"    {struct['rank']}위: Frame {struct['frame']}, {struct['distance']:.2f}Å")
                    if len(top_structures) > 3:
                        print(f"    ... 외 {len(top_structures)-3}개")

def plot_distance_evolution(iteration_results, attempt_results, output_dir):
    """거리 변화 그래프 생성 - RMSD 및 조합점수 정보 포함"""
    try:
        top_info = analyze_top_structures(attempt_results)
        
        # 기존 그래프들 + RMSD 관련 그래프 추가
        plt.figure(figsize=(20, 15))
        
        # Subplot 1: Iteration별 최종 거리
        plt.subplot(4, 3, 1)
        iterations = []
        final_distances = []
        min_distances = []
        
        for iter_result in iteration_results:
            if iter_result['success'] and iter_result['final_result']:
                iterations.append(int(iter_result['iteration']))
                final_distances.append(iter_result['final_result']['final_distance'])
                min_distances.append(iter_result['final_result'].get('min_distance', iter_result['final_result']['final_distance']))
        
        if final_distances:
            plt.plot(iterations, final_distances, 'bo-', linewidth=2, markersize=8, label='Final Distance')
            plt.plot(iterations, min_distances, 'ro-', linewidth=2, markersize=6, alpha=0.7, label='Min Distance')
            plt.axhline(y=10.0, color='k', linestyle='--', alpha=0.5, label='Long MD Threshold (10Å)')
            plt.xlabel('Iteration')
            plt.ylabel('Final Distance (Å)')
            plt.title('Final Distance per Iteration')
            plt.xticks(iterations)
            plt.legend()
            plt.grid(True, alpha=0.3)
        
        # Subplot 2: Iteration별 기울기
        plt.subplot(4, 3, 2)
        slopes = []
        
        for iter_result in iteration_results:
            if iter_result['success'] and iter_result['final_result']:
                slopes.append(iter_result['final_result']['slope'])
        
        if slopes:
            plt.plot(iterations, slopes, 'ro-', linewidth=2, markersize=4)
            plt.axhline(y=-0.001, color='k', linestyle='--', alpha=0.5, label='Threshold')
            plt.xlabel('Iteration')
            plt.ylabel('Linear Regression Slope')
            plt.title('Linear Regression Slope per Iteration')
            plt.xticks(iterations)
            plt.legend()
            plt.grid(True, alpha=0.3)
        
        # Subplot 3: Iteration별 Attempts 수
        plt.subplot(4, 3, 3)
        attempts_used = [iter_result['attempts_used'] for iter_result in iteration_results]
        all_iterations = [int(iter_result['iteration']) for iter_result in iteration_results]
        
        plt.bar(all_iterations, attempts_used, alpha=0.7, color='green')
        plt.xlabel('Iteration')
        plt.ylabel('Attempts Used')
        plt.title('Attempts Used per Iteration')
        plt.xticks(all_iterations)
        plt.grid(True, alpha=0.3)
        
        # Subplot 4: 긴 MD 통계
        plt.subplot(4, 3, 4)
        long_md_counts = {'Close Contact': 0, 'Long MD Executed': 0, 'Success after Long MD': 0}
        
        for attempt in attempt_results:
            if attempt.get('close_contact_detected', False):
                long_md_counts['Close Contact'] += 1
            if attempt.get('long_md_executed', False):
                long_md_counts['Long MD Executed'] += 1
                if attempt.get('success', False):
                    long_md_counts['Success after Long MD'] += 1
        
        categories = list(long_md_counts.keys())
        counts = list(long_md_counts.values())
        colors = ['orange', 'blue', 'green']
        
        plt.bar(categories, counts, color=colors, alpha=0.7)
        plt.ylabel('Count')
        plt.title('Long MD Statistics')
        plt.xticks(rotation=45)
        plt.grid(True, alpha=0.3)
        
        # Subplot 5: 거리 분포
        plt.subplot(4, 3, 5)
        all_min_distances = [a.get('min_distance', float('inf')) for a in attempt_results if 'min_distance' in a and a['min_distance'] != float('inf')]
        
        if all_min_distances:
            plt.hist(all_min_distances, bins=20, alpha=0.7, color='skyblue', edgecolor='black')
            plt.axvline(x=10.0, color='red', linestyle='--', linewidth=2, label='Long MD Threshold')
            plt.xlabel('Minimum Distance (Å)')
            plt.ylabel('Frequency')
            plt.title('Minimum Distance Distribution')
            plt.legend()
            plt.grid(True, alpha=0.3)
        
        # Subplot 6: 성공률 비교
        plt.subplot(4, 3, 6)
        long_md_attempts = [a for a in attempt_results if a.get('long_md_executed', False)]
        regular_attempts = [a for a in attempt_results if not a.get('long_md_executed', False)]
        
        categories = []
        success_rates = []
        
        if regular_attempts:
            regular_success_rate = len([a for a in regular_attempts if a['success']]) / len(regular_attempts) * 100
            categories.append('Regular MD')
            success_rates.append(regular_success_rate)
        
        if long_md_attempts:
            long_md_success_rate = len([a for a in long_md_attempts if a['success']]) / len(long_md_attempts) * 100
            categories.append('Long MD')
            success_rates.append(long_md_success_rate)
        
        if categories:
            plt.bar(categories, success_rates, color=['lightblue', 'lightgreen'], alpha=0.7)
            plt.ylabel('Success Rate (%)')
            plt.title('Success Rate by MD Type')
            plt.grid(True, alpha=0.3)
        
        # Subplot 7: RMSD 분포 (기존 Top 구조 거리 분포를 대체)
        plt.subplot(4, 3, 7)
        if top_info['has_rmsd_data'] and top_info['all_top_structures']:
            rmsd_values = [s.get('rmsd') for s in top_info['all_top_structures'] 
                          if s.get('rmsd') is not None]
            if rmsd_values:
                plt.hist(rmsd_values, bins=15, alpha=0.7, color='lightcoral', edgecolor='black')
                plt.xlabel('RMSD (Å)')
                plt.ylabel('Frequency')
                plt.title('Top Structure RMSD Distribution')
                plt.grid(True, alpha=0.3)
        
        # Subplot 8: 조합점수 분포
        plt.subplot(4, 3, 8)
        if top_info['best_combined_scores']:
            plt.hist(top_info['best_combined_scores'], bins=15, alpha=0.7, 
                    color='lightgreen', edgecolor='black')
            plt.xlabel('Combined Score')
            plt.ylabel('Frequency')
            plt.title('Combined Score Distribution')
            plt.grid(True, alpha=0.3)
        
        # Subplot 9: 거리 vs RMSD 산점도
        plt.subplot(4, 3, 9)
        if top_info['has_rmsd_data'] and top_info['all_top_structures']:
            distances = []
            rmsds = []
            for struct in top_info['all_top_structures']:
                if 'rmsd' in struct and struct['rmsd'] is not None:
                    distances.append(struct['distance'])
                    rmsds.append(struct['rmsd'])
            
            if distances and rmsds:
                plt.scatter(distances, rmsds, alpha=0.6, c='blue', s=50)
                plt.xlabel('Distance (Å)')
                plt.ylabel('RMSD (Å)')
                plt.title('Distance vs RMSD')
                plt.grid(True, alpha=0.3)
        
        # Subplot 10: 조합점수 vs 거리
        plt.subplot(4, 3, 10)
        if top_info['best_combined_scores'] and top_info['all_top_structures']:
            distances = []
            scores = []
            for struct in top_info['all_top_structures']:
                if 'combined_score' in struct and struct['combined_score'] is not None:
                    distances.append(struct['distance'])
                    scores.append(struct['combined_score'])
            
            if distances and scores:
                plt.scatter(distances, scores, alpha=0.6, c='green', s=50)
                plt.xlabel('Distance (Å)')
                plt.ylabel('Combined Score')
                plt.title('Distance vs Combined Score')
                plt.grid(True, alpha=0.3)
        
        # Subplot 11: 조합점수 vs RMSD
        plt.subplot(4, 3, 11)
        if top_info['has_rmsd_data'] and top_info['best_combined_scores']:
            rmsds = []
            scores = []
            for struct in top_info['all_top_structures']:
                if ('rmsd' in struct and struct['rmsd'] is not None and 
                    'combined_score' in struct and struct['combined_score'] is not None):
                    rmsds.append(struct['rmsd'])
                    scores.append(struct['combined_score'])
            
            if rmsds and scores:
                plt.scatter(rmsds, scores, alpha=0.6, c='red', s=50)
                plt.xlabel('RMSD (Å)')
                plt.ylabel('Combined Score')
                plt.title('RMSD vs Combined Score')
                plt.grid(True, alpha=0.3)
        
        plt.tight_layout()
        
        # 그래프 저장
        plot_file = os.path.join(output_dir, "analysis_plots.png")
        plt.savefig(plot_file, dpi=300, bbox_inches='tight')
        print(f"\n강화 그래프 저장됨: {plot_file}")
        
    except ImportError:
        print("matplotlib가 설치되지 않아 그래프를 생성할 수 없습니다.")
        print("설치: pip install matplotlib")
    except Exception as e:
        print(f"그래프 생성 중 오류: {e}")

def analyze_attempts(attempt_results):
    """Attempt 분석 - 긴 MD 정보 포함"""
    print("\n" + "=" * 50)
    print("Attempt 분석")
    print("=" * 50)
    
    total_attempts = len(attempt_results)
    successful_attempts = len([a for a in attempt_results if a['success']])
    
    print(f"총 attempts: {total_attempts}")
    print(f"성공한 attempts: {successful_attempts}")
    print(f"성공률: {successful_attempts/total_attempts*100:.1f}%")
    
    # 긴 MD 분석
    long_md_attempts = [a for a in attempt_results if a.get('long_md_executed', False)]
    close_contact_attempts = [a for a in attempt_results if a.get('close_contact_detected', False)]
    
    print(f"\n긴 MD 분석:")
    print(f"근접 접촉 감지: {len(close_contact_attempts)}/{total_attempts} ({len(close_contact_attempts)/total_attempts*100:.1f}%)")
    print(f"긴 MD 실행: {len(long_md_attempts)}/{total_attempts} ({len(long_md_attempts)/total_attempts*100:.1f}%)")
    
    if close_contact_attempts:
        min_distances = [a.get('min_distance', float('inf')) for a in close_contact_attempts if 'min_distance' in a]
        if min_distances:
            print(f"근접 접촉 시 최소 거리 통계:")
            print(f"  평균: {np.mean(min_distances):.2f} Å")
            print(f"  표준편차: {np.std(min_distances):.2f} Å")
            print(f"  최소값: {np.min(min_distances):.2f} Å")
            print(f"  최대값: {np.max(min_distances):.2f} Å")
    
    # 긴 MD vs 기본 MD 성공률 비교
    if long_md_attempts:
        long_md_success = len([a for a in long_md_attempts if a['success']])
        regular_attempts = [a for a in attempt_results if not a.get('long_md_executed', False)]
        regular_success = len([a for a in regular_attempts if a['success']])
        
        print(f"\n성공률 비교:")
        print(f"긴 MD 실행 시: {long_md_success}/{len(long_md_attempts)} ({long_md_success/len(long_md_attempts)*100:.1f}%)")
        if regular_attempts:
            print(f"기본 MD만: {regular_success}/{len(regular_attempts)} ({regular_success/len(regular_attempts)*100:.1f}%)")
    
    # 기울기 분포
    slopes = [a['slope'] for a in attempt_results if 'slope' in a]
    if slopes:
        print(f"\n기울기 통계:")
        print(f"  평균: {np.mean(slopes):.6f}")
        print(f"  표준편차: {np.std(slopes):.6f}")
        print(f"  최소값: {np.min(slopes):.6f}")
        print(f"  최대값: {np.max(slopes):.6f}")
    
    # 단계별 성공률
    stage_stats = {}
    for attempt in attempt_results:
        if 'stages' in attempt:
            for stage in attempt['stages']:
                stage_name = stage['stage']
                if stage_name not in stage_stats:
                    stage_stats[stage_name] = {'total': 0, 'success': 0}
                
                stage_stats[stage_name]['total'] += 1
                if stage['success']:
                    stage_stats[stage_name]['success'] += 1
    
    if stage_stats:
        print(f"\n단계별 성공률:")
        for stage, stats in stage_stats.items():
            success_rate = stats['success'] / stats['total'] * 100
            print(f"  {stage}: {success_rate:.1f}% ({stats['success']}/{stats['total']})")


def save_detailed_report(output_dir, final_results, iteration_results, attempt_results):
    """상세 보고서 저장"""
    report_file = os.path.join(output_dir, "detailed_report.txt")
    
    with open(report_file, 'w') as f:
        f.write("Simple SuMD 상세 분석 보고서\n")
        f.write("=" * 50 + "\n\n")
        
        # 기본 정보
        f.write(f"입력 PDB: {final_results['input_pdb']}\n")
        f.write(f"체인: {' - '.join(final_results['chains'])}\n")
        f.write(f"시작 시간: {final_results['start_time']}\n")
        f.write(f"종료 시간: {final_results['end_time']}\n\n")
        
        # Iteration별 상세 정보
        f.write("Iteration별 상세 정보:\n")
        f.write("-" * 30 + "\n")
        
        for iter_result in iteration_results:
            f.write(f"\nIteration {iter_result['iteration']}:\n")
            f.write(f"  성공 여부: {iter_result['success']}\n")
            f.write(f"  사용된 Attempts: {iter_result['attempts_used']}\n")
            
            if iter_result['success'] and iter_result['final_result']:
                result = iter_result['final_result']
                f.write(f"  기울기: {result.get('slope', 'N/A')}\n")
                f.write(f"  초기 거리: {result.get('initial_distance', 'N/A'):.2f} Å\n")
                f.write(f"  최종 거리: {result.get('final_distance', 'N/A'):.2f} Å\n")
                
                if 'distances' in result:
                    distances = result['distances']
                    f.write(f"  거리 통계:\n")
                    f.write(f"    평균: {np.mean(distances):.2f} Å\n")
                    f.write(f"    표준편차: {np.std(distances):.2f} Å\n")
                    f.write(f"    최소: {np.min(distances):.2f} Å\n")
                    f.write(f"    최대: {np.max(distances):.2f} Å\n")
    
    print(f"상세 보고서 저장됨: {report_file}")



def main():
    if len(sys.argv) != 2:
        print("사용법: python3 analyze_results.py <output_directory>")
        sys.exit(1)
    
    output_dir = sys.argv[1]
    
    if not os.path.exists(output_dir):
        print(f"오류: 출력 디렉토리 '{output_dir}'를 찾을 수 없습니다.")
        sys.exit(1)
    
    # 결과 로드 - final_results.json이 없어도 진행
    try:
        final_results, iteration_results, attempt_results = load_results(output_dir)
        
        if not iteration_results and not attempt_results:
            print("분석할 수 있는 결과 파일을 찾을 수 없습니다.")
            print("확인해야 할 파일들:")
            print("  - iteration_*_summary.json")
            print("  - iteration_*/iteration_*_attempt_*.json")
            sys.exit(1)
        
        # 분석 실행
        print_summary(final_results, iteration_results, attempt_results)
        analyze_attempts(attempt_results)
        print_top_structures_details(attempt_results)  # Top 구조 상세 분석 추가
        plot_distance_evolution(iteration_results, attempt_results, output_dir)
        save_detailed_report(output_dir, final_results, iteration_results, attempt_results)
        
        if final_results.get("status") == "진행 중":
            print(f"\n진행 중인 시뮬레이션 분석 완료! 결과는 {output_dir}에 저장되었습니다.")
            print("시뮬레이션이 완료되면 다시 분석을 실행하세요.")
        else:
            print(f"\n분석 완료! 결과는 {output_dir}에 저장되었습니다.")
            
    except Exception as e:
        print(f"분석 중 오류 발생: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()
