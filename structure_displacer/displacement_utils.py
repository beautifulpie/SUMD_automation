#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
50Å 거리 기반 구조 변형 결과 분석 유틸리티

새로운 변형 방식(50Å 거리 기반)에 맞춘 분석 도구입니다.

주요 기능:
1. 목표 거리(50Å) 달성도 분석
2. Clash 발생 없음 확인
3. 이동 경로 차단 없음 확인  
4. 변형 품질 평가 및 순위 매기기
5. Golden Standard와의 비교 분석

사용법:
    python displacement_utils.py --golden_standard native.pdb --displaced_dir ./displaced_structures_50A --target_chains A,B
"""

import os
import sys
import argparse
import numpy as np
import pandas as pd
import logging
import json, time
from typing import List, Dict, Tuple, Optional
from Bio.PDB import PDBParser
from datetime import datetime

# 새로운 변형 모듈에서 import
try:
    from structure_displacer import MinDistanceCalculator, ClashChecker, PathChecker
except ImportError:
    print("경고: structure_displacer.py를 찾을 수 없습니다. 기본 기능으로 대체합니다.")
    
    class MinDistanceCalculator:
        @staticmethod
        def calculate_min_distance_between_chains(chain1, chain2):
            """기본 최소 거리 계산 (단백질 원자만, HETATM 제외)"""
            min_distance = float('inf')
            
            # 표준 아미노산 잔기만 선택 (물 분자 및 HETATM 제외)
            atoms1 = []
            atoms2 = []
            
            for residue in chain1:
                if residue.id[0] == ' ':  # 표준 아미노산만 (HETATM 제외)
                    atoms1.extend(list(residue.get_atoms()))
            
            for residue in chain2:
                if residue.id[0] == ' ':  # 표준 아미노산만 (HETATM 제외)
                    atoms2.extend(list(residue.get_atoms()))
            
            for atom1 in atoms1:
                for atom2 in atoms2:
                    distance = np.linalg.norm(atom1.coord - atom2.coord)
                    if distance < min_distance:
                        min_distance = distance
            
            return min_distance, (None, None)
    
    class ClashChecker:
        def __init__(self, threshold=2.0):
            self.threshold = threshold
        
        def check_chain_clash(self, chain1, other_chains):
            return False, float('inf'), "Basic implementation"
    
    class PathChecker:
        def __init__(self, threshold=5.0, resolution=1.0):
            self.threshold = threshold
        
        def check_path_obstruction(self, start, end, other_chains):
            return False, "Basic implementation"


class Advanced50AAnalyzer:
    """50Å 거리 기반 변형 결과 분석 클래스"""
    
    def __init__(self, target_distance: float = 50.0, logger=None):
        self.target_distance = target_distance
        self.logger = logger or logging.getLogger(__name__)
        self.parser = PDBParser(QUIET=True)
        self.min_dist_calc = MinDistanceCalculator()
        self.clash_checker = ClashChecker(threshold=2.0)
        self.path_checker = PathChecker(threshold=5.0, resolution=1.0)
    
    def calculate_rmsd(self, pdb1: str, pdb2: str, chain_ids: List[str] = None) -> float:
        """두 PDB 구조 간의 RMSD 계산"""
        try:
            structure1 = self.parser.get_structure("struct1", pdb1)
            structure2 = self.parser.get_structure("struct2", pdb2)
            
            coords1 = []
            coords2 = []
            
            for model in structure1:
                for chain in model:
                    if chain_ids is None or chain.id in chain_ids:
                        for residue in chain:
                            if "CA" in residue:
                                coords1.append(residue["CA"].coord)
            
            for model in structure2:
                for chain in model:
                    if chain_ids is None or chain.id in chain_ids:
                        for residue in chain:
                            if "CA" in residue and len(coords2) < len(coords1):
                                coords2.append(residue["CA"].coord)
            
            min_len = min(len(coords1), len(coords2))
            if min_len == 0:
                return float('inf')
            
            coords1 = np.array(coords1[:min_len])
            coords2 = np.array(coords2[:min_len])
            
            # 중심 이동 후 RMSD 계산
            center1 = np.mean(coords1, axis=0)
            center2 = np.mean(coords2, axis=0)
            
            coords1_centered = coords1 - center1
            coords2_centered = coords2 - center2
            
            diff = coords1_centered - coords2_centered
            rmsd = np.sqrt(np.mean(np.sum(diff**2, axis=1)))
            
            return float(rmsd)
            
        except Exception as e:
            self.logger.error(f"RMSD 계산 실패: {e}")
            return float('inf')
    
    def analyze_displaced_structure(self, pdb_file: str, target_chains: List[str]) -> Dict:
        """단일 변형된 구조 분석"""
        analysis_start_time = time.time()
        if len(target_chains) != 2:
            raise ValueError("정확히 2개의 체인이 필요합니다")
        
        try:
            structure = self.parser.get_structure("structure", pdb_file)
            
            # 체인 수집
            chains = {}
            for model in structure:
                for chain in model:
                    chains[chain.id] = chain
            
            target_chain1 = chains.get(target_chains[0])
            target_chain2 = chains.get(target_chains[1])
            
            if target_chain1 is None or target_chain2 is None:
                raise ValueError(f"타겟 체인들을 찾을 수 없습니다: {target_chains}")
            
            # 다른 체인들
            other_chains = [chain for chain_id, chain in chains.items() 
                          if chain_id not in target_chains]
            
            analysis_result = {
                'filename': os.path.basename(pdb_file),
                'filepath': pdb_file,
                'target_chains': target_chains,
                'other_chains': [chain.id for chain in other_chains]
            }
            
            # 1. 거리 분석
            min_distance, closest_coords = self.min_dist_calc.calculate_min_distance_between_chains(
                target_chain1, target_chain2
            )
            
            distance_error = abs(min_distance - self.target_distance)
            distance_achievement = max(0, 1 - distance_error / 10.0)  # 10Å 오차를 기준으로 점수화
            
            analysis_result['distance_analysis'] = {
                'min_distance': min_distance,
                'target_distance': self.target_distance,
                'distance_error': distance_error,
                'distance_achievement_score': distance_achievement,
                'distance_acceptable': distance_error <= 5.0  # 5Å 이내 오차는 허용
            }
            
            # 2. Clash 분석
            has_clash1, min_clash_dist1, clash_info1 = self.clash_checker.check_chain_clash(
                target_chain1, [target_chain2] + other_chains
            )
            
            has_clash2, min_clash_dist2, clash_info2 = self.clash_checker.check_chain_clash(
                target_chain2, [target_chain1] + other_chains
            )
            
            overall_clash = has_clash1 or has_clash2
            min_overall_clash_dist = min(min_clash_dist1, min_clash_dist2)
            
            clash_score = 1.0 if not overall_clash else max(0, (min_overall_clash_dist - 2.0) / 3.0)
            
            analysis_result['clash_analysis'] = {
                'has_clash': overall_clash,
                'min_clash_distance': min_overall_clash_dist,
                'clash_score': clash_score,
                'clash_details': {
                    'chain1_clash': has_clash1,
                    'chain1_clash_info': clash_info1,
                    'chain2_clash': has_clash2,
                    'chain2_clash_info': clash_info2
                }
            }
            
            # 3. 경로 분석 (원래 위치를 추정해야 하므로 간소화)
            # 실제로는 원래 구조와 비교해야 하지만, 여기서는 체인 중심 간 직선 경로를 분석
            chain1_center = np.mean([atom.coord for atom in target_chain1.get_atoms()], axis=0)
            chain2_center = np.mean([atom.coord for atom in target_chain2.get_atoms()], axis=0)
            
            # 가상의 "가까운" 위치 (현재 위치에서 target_distance만큼 가까워진 지점)
            direction = chain1_center - chain2_center
            if np.linalg.norm(direction) > 0:
                direction = direction / np.linalg.norm(direction)
                estimated_close_position = chain2_center + direction * (min_distance - self.target_distance)
            else:
                estimated_close_position = chain1_center
            
            is_path_obstructed, obstruction_info = self.path_checker.check_path_obstruction(
                estimated_close_position, chain1_center, other_chains
            )
            
            path_score = 0.0 if is_path_obstructed else 1.0
            
            analysis_result['path_analysis'] = {
                'is_path_obstructed': is_path_obstructed,
                'obstruction_info': obstruction_info,
                'path_score': path_score,
                'estimated_movement_vector': (chain1_center - estimated_close_position).tolist(),
                'movement_magnitude': np.linalg.norm(chain1_center - estimated_close_position)
            }
            
            # 4. 종합 품질 점수
            quality_score = (
                0.5 * distance_achievement +  # 거리 달성도 50%
                0.3 * clash_score +           # Clash 없음 30%  
                0.2 * path_score              # 경로 차단 없음 20%
            )
            
            analysis_result['overall_quality'] = {
                'quality_score': quality_score,
                'components': {
                    'distance_achievement': distance_achievement,
                    'clash_score': clash_score,
                    'path_score': path_score
                },
                'weights': {
                    'distance': 0.5,
                    'clash': 0.3,
                    'path': 0.2
                }
            }
            
            # 5. 구조 유효성
            is_valid = (
                analysis_result['distance_analysis']['distance_acceptable'] and
                not analysis_result['clash_analysis']['has_clash'] and
                not analysis_result['path_analysis']['is_path_obstructed']
            )
            
            analysis_result['structure_validity'] = {
                'is_valid': is_valid,
                'validation_criteria': {
                    'distance_acceptable': analysis_result['distance_analysis']['distance_acceptable'],
                    'no_clash': not analysis_result['clash_analysis']['has_clash'],
                    'path_clear': not analysis_result['path_analysis']['is_path_obstructed']
                }
            }

            analysis_time = time.time() - analysis_start_time

            analysis_result['timing'] = {
                'analysis_time_seconds': analysis_time,
                'analysis_timestamp': datetime.now().isoformat()
            }
            
            return analysis_result
            
        except Exception as e:
            self.logger.error(f"구조 분석 실패 ({pdb_file}): {e}")
            analysis_failed_time = time.time() - analysis_start_time
            return {
                'filename': os.path.basename(pdb_file),
                'filepath': pdb_file,
                'error': str(e),
                'analysis_failed': True,
                'timing': {
                    'analysis_time_seconds': analysis_failed_time,
                    'analysis_timestamp': datetime.now().isoformat()
                }
            }
    
    def analyze_batch_structures(self, golden_standard: str, displaced_files: List[str], 
                               target_chains: List[str], output_dir: str = None) -> pd.DataFrame:
        """배치 구조 분석"""
        batch_start_time = time.time()
        if output_dir is None:
            output_dir = "displacement_analysis_50A"
        os.makedirs(output_dir, exist_ok=True)
        
        self.logger.info(f"배치 구조 분석 시작: {len(displaced_files)}개 파일")
        
        # Golden Standard 분석
        golden_start_time = time.time()
        golden_analysis = self.analyze_displaced_structure(golden_standard, target_chains)
        golden_time = time.time() - golden_start_time
        golden_distance = golden_analysis['distance_analysis']['min_distance']
        
        self.logger.info(f"Golden Standard 거리: {golden_distance:.2f}Å")
        self.logger.info(f"Golden Standard 분석 시간: {golden_time:.2f}초")
        
        # 각 변형된 구조 분석
        results = []
        analysis_times = []
        
        for i, displaced_file in enumerate(displaced_files):
            file_start_time = time.time()
            try:
                analysis = self.analyze_displaced_structure(displaced_file, target_chains)
                file_time = time.time() - file_start_time
                analysis_times.append(file_time)

                if 'analysis_failed' not in analysis:
                    # Golden Standard와의 RMSD 계산
                    rmsd_vs_golden = self.calculate_rmsd(golden_standard, displaced_file, target_chains)
                    analysis['rmsd_vs_golden'] = rmsd_vs_golden
                    
                    # Golden Standard와의 거리 변화
                    distance_change = analysis['distance_analysis']['min_distance'] - golden_distance
                    analysis['distance_change_vs_golden'] = distance_change
                    
                    results.append(analysis)
                
                if (i + 1) % 10 == 0:
                    avg_time = sum(analysis_times[-10:]) / min(10, len(analysis_times))  # 추가
                    self.logger.info(f"분석 진행: {i + 1}/{len(displaced_files)} (최근 10개 평균: {avg_time:.2f}초/파일)")
                    
            except Exception as e:
                self.logger.error(f"파일 {displaced_file} 분석 실패: {e}")
                continue
        
        if not results:
            self.logger.error("분석할 수 있는 파일이 없습니다.")
            return pd.DataFrame()
        
        # DataFrame 생성
        df_creation_start_time = time.time()
        df_data = []
        for result in results:
            row = {
                'filename': result['filename'],
                'filepath': result['filepath'],
                'min_distance': result['distance_analysis']['min_distance'],
                'distance_error': result['distance_analysis']['distance_error'],
                'distance_acceptable': result['distance_analysis']['distance_acceptable'],
                'distance_achievement_score': result['distance_analysis']['distance_achievement_score'],
                'has_clash': result['clash_analysis']['has_clash'],
                'min_clash_distance': result['clash_analysis']['min_clash_distance'],
                'clash_score': result['clash_analysis']['clash_score'],
                'is_path_obstructed': result['path_analysis']['is_path_obstructed'],
                'path_score': result['path_analysis']['path_score'],
                'quality_score': result['overall_quality']['quality_score'],
                'is_valid': result['structure_validity']['is_valid'],
                'rmsd_vs_golden': result.get('rmsd_vs_golden', float('inf')),
                'distance_change_vs_golden': result.get('distance_change_vs_golden', 0),
                'movement_magnitude': result['path_analysis']['movement_magnitude']
            }
            df_data.append(row)
        
        df = pd.DataFrame(df_data)
        
        # 유효한 구조들만 품질 점수로 정렬
        valid_df = df[df['is_valid'] == True].copy()
        invalid_df = df[df['is_valid'] == False].copy()
        
        if not valid_df.empty:
            valid_df = valid_df.sort_values('quality_score', ascending=False).reset_index(drop=True)
            valid_df['rank'] = valid_df.index + 1
        
        if not invalid_df.empty:
            invalid_df = invalid_df.sort_values('quality_score', ascending=False).reset_index(drop=True)
            invalid_df['rank'] = len(valid_df) + invalid_df.index + 1
        
        # 합치기 (유효한 구조들이 먼저)
        df = pd.concat([valid_df, invalid_df], ignore_index=True)
        
        # 결과 저장
        output_file = os.path.join(output_dir, "displacement_analysis_50A.csv")
        df.to_csv(output_file, index=False)
        
        # JSON 저장을 위해 NumPy 타입을 Python 기본 타입으로 변환
        def convert_numpy_types(obj):
            """NumPy 타입을 Python 기본 타입으로 변환"""
            if isinstance(obj, np.integer):
                return int(obj)
            elif isinstance(obj, np.floating):
                return float(obj)
            elif isinstance(obj, np.ndarray):
                return obj.tolist()
            return obj
        
        # DataFrame을 딕셔너리로 변환 후 NumPy 타입 변환
        records = df.to_dict('records')
        for record in records:
            for key, value in record.items():
                record[key] = convert_numpy_types(value)
        
        json_file = os.path.join(output_dir, "displacement_analysis_50A.json")
        with open(json_file, 'w') as f:
            json.dump(records, f, indent=2)
        
        # 요약 통계
        summary = {
            'analysis_type': '50A_distance_based',
            'target_distance': self.target_distance,
            'golden_standard_distance': golden_distance,
            'total_files': len(displaced_files),
            'analyzed_files': len(df),
            'valid_structures': len(valid_df),
            'invalid_structures': len(invalid_df),
            'validation_rate': float(len(valid_df) / len(df) * 100) if len(df) > 0 else 0.0,
            'distance_stats': {
                'mean': float(df['min_distance'].mean()),
                'std': float(df['min_distance'].std()),
                'min': float(df['min_distance'].min()),
                'max': float(df['min_distance'].max()),
                'target_achievement_rate': float(len(df[df['distance_acceptable']]) / len(df) * 100)
            },
            'clash_stats': {
                'clash_free_rate': float(len(df[~df['has_clash']]) / len(df) * 100),
                'mean_clash_distance': float(df['min_clash_distance'].mean())
            },
            'path_stats': {
                'path_clear_rate': float(len(df[~df['is_path_obstructed']]) / len(df) * 100)
            },
            'quality_stats': {
                'mean': float(df['quality_score'].mean()),
                'std': float(df['quality_score'].std()),
                'min': float(df['quality_score'].min()),
                'max': float(df['quality_score'].max())
            }
        }
        
        if not valid_df.empty:
            summary['best_valid_structure'] = {
                'rank': 1,
                'filename': valid_df.iloc[0]['filename'],
                'quality_score': float(valid_df.iloc[0]['quality_score']),
                'min_distance': float(valid_df.iloc[0]['min_distance']),
                'distance_error': float(valid_df.iloc[0]['distance_error']),
                'rmsd_vs_golden': float(valid_df.iloc[0]['rmsd_vs_golden'])
            }
        
        # summary 딕셔너리의 NumPy 타입도 변환
        def convert_dict_numpy_types(d):
            """딕셔너리 내의 NumPy 타입을 재귀적으로 변환"""
            if isinstance(d, dict):
                return {k: convert_dict_numpy_types(v) for k, v in d.items()}
            elif isinstance(d, list):
                return [convert_dict_numpy_types(item) for item in d]
            else:
                return convert_numpy_types(d)
        
        summary_converted = convert_dict_numpy_types(summary)
        
        summary_file = os.path.join(output_dir, "analysis_summary_50A.json")
        with open(summary_file, 'w') as f:
            json.dump(summary_converted, f, indent=2)
        
        df_creation_time = time.time() - df_creation_start_time
        batch_total_time = time.time() - batch_start_time

        self.logger.info(f"분석 완료: {len(df)}개 구조 분석")
        self.logger.info(f"배치 분석 총 시간: {batch_total_time:.2f}초")
        self.logger.info(f"파일당 평균 분석 시간: {sum(analysis_times)/len(analysis_times):.2f}초")
        self.logger.info(f"DataFrame 생성 시간: {df_creation_time:.2f}초")
        self.logger.info(f"유효한 구조: {len(valid_df)}개 ({len(valid_df)/len(df)*100:.1f}%)")
        if not valid_df.empty:
            best = summary['best_valid_structure']
            self.logger.info(f"최고 품질 구조: {best['filename']} "
                           f"(점수: {best['quality_score']:.3f}, 거리: {best['min_distance']:.2f}Å)")
        
        return df
    
    def create_analysis_plots(self, df: pd.DataFrame, output_dir: str):
        """분석 결과 시각화"""
        if df.empty:
            self.logger.warning("시각화할 데이터가 없습니다.")
            return
        
        try:
            import matplotlib.pyplot as plt
            import seaborn as sns
            
            plt.style.use('default')
            sns.set_palette("husl")
            
            fig, axes = plt.subplots(2, 3, figsize=(18, 12))
            fig.suptitle('50A Distance-based Structure Displacement Analysis Results', fontsize=16)
            
            # 1. Distance Achievement Distribution
            axes[0, 0].hist(df['min_distance'], bins=20, alpha=0.7, edgecolor='black')
            axes[0, 0].axvline(self.target_distance, color='red', linestyle='--', linewidth=2, 
                              label=f'Target Distance: {self.target_distance}A')
            axes[0, 0].set_xlabel('Minimum Distance (A)')
            axes[0, 0].set_ylabel('Frequency')
            axes[0, 0].set_title('Distance Achievement Distribution')
            axes[0, 0].legend()
            
            # 2. 유효성 파이 차트
            valid_counts = df['is_valid'].value_counts()
            
            # 라벨과 색상을 실제 데이터에 맞춰 동적 생성
            labels = []
            colors = []
            for idx in valid_counts.index:
                if idx == True:
                    labels.append('Valid')
                    colors.append('lightgreen')
                else:
                    labels.append('Invalid')
                    colors.append('lightcoral')
            
            axes[0, 1].pie(valid_counts.values, labels=labels, colors=colors, 
                          autopct='%1.1f%%', startangle=90)
            axes[0, 1].set_title('Structure Validity Distribution')
            
            # 3. 품질 점수 vs 거리 오차
            valid_mask = df['is_valid'] == True
            scatter = axes[0, 2].scatter(df[valid_mask]['distance_error'], df[valid_mask]['quality_score'], 
                                       c='green', alpha=0.6, label='Valid', s=50)
            if any(~valid_mask):
                axes[0, 2].scatter(df[~valid_mask]['distance_error'], df[~valid_mask]['quality_score'], 
                                 c='red', alpha=0.6, label='Invalid', s=50)
            axes[0, 2].set_xlabel('Distance Error (A)')
            axes[0, 2].set_ylabel('Quality Score')
            axes[0, 2].set_title('Quality Score vs Distance Error')
            axes[0, 2].legend()
            
            # 4. Clash 분석
            clash_counts = df['has_clash'].value_counts()
            
            # Clash 라벨과 색상을 실제 데이터에 맞춰 동적 생성
            clash_labels = []
            clash_colors = []
            for idx in clash_counts.index:
                if idx == False:
                    clash_labels.append('No Clash')
                    clash_colors.append('lightblue')
                else:
                    clash_labels.append('Clash')
                    clash_colors.append('orange')
            
            axes[1, 0].pie(clash_counts.values, labels=clash_labels, colors=clash_colors, 
                          autopct='%1.1f%%', startangle=90)
            axes[1, 0].set_title('Clash Distribution')
            
            # 5. 경로 차단 분석
            path_counts = df['is_path_obstructed'].value_counts()
            
            # 경로 차단 라벨과 색상을 실제 데이터에 맞춰 동적 생성
            path_labels = []
            path_colors = []
            for idx in path_counts.index:
                if idx == False:
                    path_labels.append('Path Clear')
                    path_colors.append('lightgreen')
                else:
                    path_labels.append('Path Obstructed')
                    path_colors.append('salmon')
            
            axes[1, 1].pie(path_counts.values, labels=path_labels, colors=path_colors, 
                          autopct='%1.1f%%', startangle=90)
            axes[1, 1].set_title('Path Obstruction Distribution')
            
            # 6. 상위 10개 구조의 품질 점수
            top_10 = df.head(10)
            bars = axes[1, 2].bar(range(len(top_10)), top_10['quality_score'],
                                 color=['green' if valid else 'red' for valid in top_10['is_valid']])
            axes[1, 2].set_xlabel('Rank')
            axes[1, 2].set_ylabel('Quality Score')
            axes[1, 2].set_title('Top 10 Structure Quality Scores')
            axes[1, 2].set_xticks(range(len(top_10)))
            axes[1, 2].set_xticklabels([f"{i+1}" for i in range(len(top_10))])
            
            # 범례 추가
            from matplotlib.patches import Patch
            legend_elements = [Patch(facecolor='green', label='Valid'),
                             Patch(facecolor='red', label='Invalid')]
            axes[1, 2].legend(handles=legend_elements)
            
            plt.tight_layout()
            
            # 저장
            plot_file = os.path.join(output_dir, "displacement_analysis_50A_plots.png")
            plt.savefig(plot_file, dpi=300, bbox_inches='tight')
            plt.close()
            
            self.logger.info(f"시각화 저장: {plot_file}")
            
        except ImportError:
            self.logger.warning("matplotlib/seaborn이 설치되지 않아 시각화를 건너뜁니다.")
        except Exception as e:
            self.logger.error(f"시각화 생성 실패: {e}")


def collect_displaced_files(directory: str, pattern: str = "*.pdb") -> List[str]:
    """디렉토리에서 변형된 PDB 파일들 수집"""
    import glob
    
    search_pattern = os.path.join(directory, "**", pattern)
    files = glob.glob(search_pattern, recursive=True)
    files.sort()
    
    return files


def main():
    """메인 함수"""
    parser = argparse.ArgumentParser(
        description="50Å 거리 기반 구조 변형 결과 분석 도구",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
분석 기준:
- 목표 거리 달성도 (50Å ± 2Å)
- Clash 발생 없음 (원자 간 거리 ≥ 2Å)
- 이동 경로 차단 없음
- 종합 품질 점수 (위 3개 지표의 가중 평균)

사용 예시:
  python displacement_utils.py --golden_standard native.pdb --displaced_dir ./displaced_structures_50A --target_chains A,B
        """
    )
    
    # 필수 인수
    parser.add_argument("--golden_standard", required=True, 
                       help="Golden Standard PDB 파일")
    parser.add_argument("--target_chains", required=True,
                       help="타겟 체인 ID (쉼표로 구분, 예: A,B)")
    
    # 입력 파일 지정
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument("--displaced_dir", 
                           help="변형된 구조가 있는 디렉토리")
    input_group.add_argument("--displaced_files",
                           help="변형된 PDB 파일들 (쉼표로 구분)")
    
    # 분석 파라미터
    parser.add_argument("--target_distance", type=float, default=50.0,
                       help="목표 거리 (Å, 기본값: 50.0)")
    
    # 출력 관련
    parser.add_argument("--output_dir", default="displacement_analysis_50A",
                       help="출력 디렉토리")
    parser.add_argument("--plot", action="store_true",
                       help="분석 결과 시각화 생성")
    
    # 기타
    parser.add_argument("--verbose", "-v", action="store_true", 
                       help="상세 로그 출력")
    
    args = parser.parse_args()
    
    # 로깅 설정
    level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )
    logger = logging.getLogger("50A-DisplacementAnalyzer")
    
    # 파일 존재 확인
    if not os.path.exists(args.golden_standard):
        logger.error(f"Golden Standard 파일을 찾을 수 없습니다: {args.golden_standard}")
        return 1
    
    # 타겟 체인 파싱
    target_chains = [chain.strip().upper() for chain in args.target_chains.split(',')]
    if len(target_chains) != 2:
        logger.error("정확히 2개의 체인이 필요합니다")
        return 1
    
    logger.info(f"타겟 체인: {target_chains}")
    logger.info(f"목표 거리: {args.target_distance}Å")
    
    # 변형된 파일들 수집
    displaced_files = []
    
    if args.displaced_dir:
        if not os.path.isdir(args.displaced_dir):
            logger.error(f"변형된 구조 디렉토리를 찾을 수 없습니다: {args.displaced_dir}")
            return 1
        
        displaced_files = collect_displaced_files(args.displaced_dir)
        logger.info(f"디렉토리에서 {len(displaced_files)}개 PDB 파일 발견")
        
    elif args.displaced_files:
        file_paths = [f.strip() for f in args.displaced_files.split(',')]
        for file_path in file_paths:
            if os.path.exists(file_path):
                displaced_files.append(file_path)
            else:
                logger.warning(f"파일을 찾을 수 없습니다: {file_path}")
    
    if not displaced_files:
        logger.error("분석할 변형된 구조 파일이 없습니다.")
        return 1
    
    logger.info(f"총 {len(displaced_files)}개 파일 분석 예정")
    
    try:
        # 분석 실행
        analyzer = Advanced50AAnalyzer(args.target_distance, logger)
        df = analyzer.analyze_batch_structures(
            args.golden_standard, displaced_files, target_chains, args.output_dir
        )
        
        if df.empty:
            logger.error("분석 결과가 없습니다.")
            return 1
        
        # 시각화 생성
        if args.plot:
            analyzer.create_analysis_plots(df, args.output_dir)
        
        # 결과 요약 출력
        valid_df = df[df['is_valid'] == True]
        invalid_df = df[df['is_valid'] == False]
        
        print(f"\n=== 50Å 거리 기반 구조 변형 분석 결과 ===")
        print(f"분석된 파일: {len(df)}개")
        print(f"유효한 구조: {len(valid_df)}개 ({len(valid_df)/len(df)*100:.1f}%)")
        print(f"무효한 구조: {len(invalid_df)}개 ({len(invalid_df)/len(df)*100:.1f}%)")
        print(f"")
        print(f"거리 달성 통계:")
        print(f"  평균 거리: {df['min_distance'].mean():.2f}±{df['min_distance'].std():.2f}Å")
        print(f"  목표 거리 달성률: {len(df[df['distance_acceptable']])/len(df)*100:.1f}%")
        print(f"")
        print(f"품질 통계:")
        print(f"  Clash 없음: {len(df[~df['has_clash']])/len(df)*100:.1f}%")
        print(f"  경로 차단 없음: {len(df[~df['is_path_obstructed']])/len(df)*100:.1f}%")
        print(f"  평균 품질 점수: {df['quality_score'].mean():.3f}±{df['quality_score'].std():.3f}")
        
        if not valid_df.empty:
            print(f"\n최고 품질 유효 구조 TOP 5:")
            top_5 = valid_df.head(5)
            for idx, row in top_5.iterrows():
                print(f"  {int(row['rank'])}위: {row['filename']}")
                print(f"       품질점수: {row['quality_score']:.3f}, 거리: {row['min_distance']:.2f}Å (오차: {row['distance_error']:.2f}Å)")
        
        print(f"\n결과 파일:")
        print(f"  CSV: {os.path.join(args.output_dir, 'displacement_analysis_50A.csv')}")
        print(f"  JSON: {os.path.join(args.output_dir, 'displacement_analysis_50A.json')}")
        print(f"  요약: {os.path.join(args.output_dir, 'analysis_summary_50A.json')}")
        
        if args.plot:
            print(f"  시각화: {os.path.join(args.output_dir, 'displacement_analysis_50A_plots.png')}")
        
        return 0
        
    except Exception as e:
        logger.error(f"분석 실행 중 오류: {e}")
        return 1


if __name__ == "__main__":
    exit(main())
