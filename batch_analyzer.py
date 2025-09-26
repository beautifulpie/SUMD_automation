#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Simple SuMD 배치 실행 결과 통합 분석 스크립트 (개선된 버전)
현재 simple_sumd.py 코드 구조에 맞게 설계됨
"""

import os
import sys
import json
import glob
import re
from datetime import datetime, timedelta
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.ticker import MaxNLocator
from collections import defaultdict
import pandas as pd

class BatchAnalyzer:
    def __init__(self, output_dir, long_md_threshold=10.0):
        """배치 분석기 초기화"""
        self.output_dir = output_dir
        self.long_md_threshold = long_md_threshold  # longMD 기준값 설정 가능
        self.data = {
            'batch_info': {},
            'pdb_results': [],
            'structure_results': [],
            'iteration_results': [],
            'attempt_results': [],
            'execution_log': {}
        }
    
    def extract_pdb_info_from_folder_name(self, folder_name):
        """폴더명에서 PDB 정보 추출 (개선된 버전)"""
        try:
            # 패턴: (pdb_code)_(chain1)_(chain2) 또는 (pdb_code)_(receptor)_(ligand)
            parts = folder_name.split('_')
            if len(parts) >= 3:
                return {
                    'pdb_code': parts[0],
                    'chain1': parts[1],
                    'chain2': parts[2],
                    'receptor_chain': parts[1],
                    'ligand_chain': parts[2]
                }
            elif len(parts) == 2:
                # 구조명만 있는 경우
                return {
                    'pdb_code': parts[0],
                    'chain1': 'A',
                    'chain2': 'B',
                    'receptor_chain': 'A',
                    'ligand_chain': 'B'
                }
            else:
                return {
                    'pdb_code': folder_name,
                    'chain1': 'A',
                    'chain2': 'B',
                    'receptor_chain': 'A',
                    'ligand_chain': 'B'
                }
        except:
            return {
                'pdb_code': folder_name,
                'chain1': 'A',
                'chain2': 'B',
                'receptor_chain': 'A',
                'ligand_chain': 'B'
            }

    def reconstruct_pdb_results_from_structures(self):
        """구조 결과에서 PDB 결과 재구성"""
        print("구조 결과에서 PDB 정보 재구성 중...")
        
        pdb_data = defaultdict(lambda: {
            'pdb_code': '',
            'chain1': '',
            'chain2': '',
            'receptor_chain': '',
            'ligand_chain': '',
            'total_structures': 0,
            'successful_structures': 0,
            'max_successful_iterations': 0,
            'best_structure_name': '',
            'has_successful_result': False,
            'structures': []
        })
        
        # 구조 결과에서 PDB 정보 집계
        for structure in self.data['structure_results']:
            pdb_name = structure.get('pdb_name', '')
            pdb_code = structure.get('pdb_code', '')
            
            # PDB 코드가 없으면 pdb_name에서 추출
            if not pdb_code and pdb_name:
                pdb_info = self.extract_pdb_info_from_folder_name(pdb_name)
                pdb_code = pdb_info['pdb_code']
                structure['pdb_code'] = pdb_code
                structure['chain1'] = pdb_info['chain1']
                structure['chain2'] = pdb_info['chain2']
            
            if not pdb_code:
                continue
            
            pdb_entry = pdb_data[pdb_code]
            
            # 기본 정보 설정 (첫 번째 구조에서)
            if not pdb_entry['pdb_code']:
                pdb_entry['pdb_code'] = pdb_code
                pdb_entry['chain1'] = structure.get('chain1', 'A')
                pdb_entry['chain2'] = structure.get('chain2', 'B')
                pdb_entry['receptor_chain'] = structure.get('receptor_chain', structure.get('chain1', 'A'))
                pdb_entry['ligand_chain'] = structure.get('ligand_chain', structure.get('chain2', 'B'))
            
            # 통계 집계
            pdb_entry['total_structures'] += 1
            pdb_entry['structures'].append(structure)
            
            success_count = structure.get('successful_iterations', 0)
            if success_count > 0:
                pdb_entry['successful_structures'] += 1
                pdb_entry['has_successful_result'] = True
                
                if success_count > pdb_entry['max_successful_iterations']:
                    pdb_entry['max_successful_iterations'] = success_count
                    pdb_entry['best_structure_name'] = structure.get('structure_name', '')
        
        # pdb_results 업데이트
        self.data['pdb_results'] = list(pdb_data.values())
        
        print(f"PDB 정보 재구성 완료: {len(self.data['pdb_results'])}개 PDB")

    def load_detailed_results(self):
        """상세 결과 파일들 로드 (개선된 버전)"""
        print("상세 결과 파일들 로드 중...")
        
        # 각 PDB 폴더에서 상세 결과 수집
        pdb_folders = [d for d in glob.glob(os.path.join(self.output_dir, "*")) 
                       if os.path.isdir(d) and not d.endswith('temp_batch_input')]
        
        iteration_count = 0
        attempt_count = 0
        
        for pdb_folder in pdb_folders:
            pdb_name = os.path.basename(pdb_folder)
            
            # 구조별 폴더 찾기
            structure_folders = glob.glob(os.path.join(pdb_folder, "structure_*"))
            
            for structure_folder in structure_folders:
                structure_name = os.path.basename(structure_folder).replace("structure_", "")
                
                # structure_results.json 로드
                struct_result_file = os.path.join(structure_folder, "structure_results.json")
                if os.path.exists(struct_result_file):
                    try:
                        with open(struct_result_file) as f:
                            struct_data = json.load(f)
                            struct_data['pdb_name'] = pdb_name
                            struct_data['folder_path'] = structure_folder
                            
                            # PDB 정보가 없으면 폴더명에서 추출
                            if 'pdb_code' not in struct_data:
                                pdb_info = self.extract_pdb_info_from_folder_name(pdb_name)
                                struct_data.update(pdb_info)
                            
                            self.data['structure_results'].append(struct_data)
                    except Exception as e:
                        print(f"구조 결과 파일 로드 실패 {struct_result_file}: {e}")
                
                # iteration 요약 파일들 로드
                iter_files = glob.glob(os.path.join(structure_folder, "iteration_*_summary.json"))
                for iter_file in iter_files:
                    try:
                        with open(iter_file) as f:
                            iter_data = json.load(f)
                            iter_data['pdb_name'] = pdb_name
                            iter_data['structure_name'] = structure_name
                            iter_data['file_path'] = iter_file
                            
                            # PDB 정보 추가
                            pdb_info = self.extract_pdb_info_from_folder_name(pdb_name)
                            iter_data.update(pdb_info)
                            
                            self.data['iteration_results'].append(iter_data)
                            iteration_count += 1
                    except Exception as e:
                        print(f"iteration 파일 로드 실패 {iter_file}: {e}")
                
                # attempt 상세 파일들 로드
                attempt_files = glob.glob(os.path.join(structure_folder, "iteration_*", "iteration_*_attempt_*.json"))
                for attempt_file in attempt_files:
                    try:
                        with open(attempt_file) as f:
                            attempt_data = json.load(f)
                            attempt_data['pdb_name'] = pdb_name
                            attempt_data['structure_name'] = structure_name
                            attempt_data['file_path'] = attempt_file
                            
                            # 파일 시간 정보 추가
                            file_stat = os.stat(attempt_file)
                            attempt_data['file_time'] = datetime.fromtimestamp(file_stat.st_mtime)
                            
                            # PDB 정보 추가
                            pdb_info = self.extract_pdb_info_from_folder_name(pdb_name)
                            attempt_data.update(pdb_info)
                            
                            self.data['attempt_results'].append(attempt_data)
                            attempt_count += 1
                    except Exception as e:
                        print(f"attempt 파일 로드 실패 {attempt_file}: {e}")
        
        print(f"상세 결과 로드 완료: {iteration_count}개 iteration, {attempt_count}개 attempt")
        
        # PDB 결과가 비어있으면 구조 결과에서 재구성
        if not self.data['pdb_results'] and self.data['structure_results']:
            self.reconstruct_pdb_results_from_structures()

    def plot_comprehensive_analysis(self):
        """종합 분석 그래프 생성 - longMD 기준값 적용"""
        try:
            # 구조별 거리 데이터 수집 및 그래프 개수 계산
            distance_graph_count = self._calculate_distance_graph_count()
            total_graphs = 4 + distance_graph_count  # 기본 4개 + 거리 그래프들
            
            # 동적 레이아웃 계산 (최대 3열)
            cols = 3
            rows = (total_graphs + cols - 1) // cols
            
            fig = plt.figure(figsize=(18, 6 * rows))
            
            # 영어 폰트 설정
            plt.rcParams['font.family'] = 'DejaVu Sans'
            plt.rcParams['axes.unicode_minus'] = False
            
            current_subplot = 1
            
            # 1. PDB별 성공률
            ax1 = plt.subplot(rows, cols, current_subplot)
            current_subplot += 1
            
            if self.data['pdb_results']:
                pdb_names = [p['pdb_code'] for p in self.data['pdb_results'][:15]]
                pdb_success_rates = [p.get('max_successful_iterations', 0) for p in self.data['pdb_results'][:15]]
                
                bars = ax1.bar(range(len(pdb_names)), pdb_success_rates, alpha=0.7, color='skyblue')
                ax1.set_xlabel('PDB Code')
                ax1.set_ylabel('Max Successful Iterations')
                ax1.set_title('Max Successful Iterations by PDB')
                ax1.set_xticks(range(len(pdb_names)))
                ax1.set_xticklabels(pdb_names, rotation=45)
                ax1.yaxis.set_major_locator(MAXNLocator(integer=True))  # 정수 축
                ax1.grid(True, alpha=0.3)
            else:
                ax1.text(0.5, 0.5, 'No PDB data available', ha='center', va='center', transform=ax1.transAxes)
                ax1.set_title('Max Successful Iterations by PDB')
            
            # 2-N. 구조별 iteration 순차적 거리 변화 (동적 그래프 생성)
            current_subplot = self._plot_distance_graphs(fig, rows, cols, current_subplot)
            
            # N+1. Iteration당 Attempt 수 분포
            ax_attempts = plt.subplot(rows, cols, current_subplot)
            current_subplot += 1
            
            if self.data['iteration_results']:
                attempts_per_iteration = [i.get('attempts_used', 0) for i in self.data['iteration_results'] if i.get('attempts_used', 0) > 0]
                
                if attempts_per_iteration:
                    max_attempts = max(attempts_per_iteration)
                    bins = range(1, max_attempts + 2)  # 1부터 max_attempts+1까지
                    
                    ax_attempts.hist(attempts_per_iteration, bins=bins, alpha=0.7, color='lightcoral', edgecolor='black', align='left')
                    ax_attempts.set_xlabel('Attempts per Iteration')
                    ax_attempts.set_ylabel('Frequency')
                    ax_attempts.set_title('Attempts per Iteration Distribution')
                    ax_attempts.set_xticks(range(1, max_attempts + 1))
                    ax_attempts.yaxis.set_major_locator(MAXNLocator(integer=True))  # 정수 축
                    ax_attempts.grid(True, alpha=0.3)
                else:
                    ax_attempts.text(0.5, 0.5, 'No iteration data', ha='center', va='center', transform=ax_attempts.transAxes)
                    ax_attempts.set_title('Attempts per Iteration Distribution')
            else:
                ax_attempts.text(0.5, 0.5, 'No iteration data', ha='center', va='center', transform=ax_attempts.transAxes)
                ax_attempts.set_title('Attempts per Iteration Distribution')
            
            # N+2. 거리 통계 (성공한 attempt들) - 동적 threshold 적용
            ax_dist = plt.subplot(rows, cols, current_subplot)
            current_subplot += 1
            
            if self.data['attempt_results']:
                fin_distances = []
                
                for attempt in self.data['attempt_results']:
                    if attempt.get('success', False) and 'final_distance' in attempt:
                        fin_distances.append(attempt['final_distance'])
                
                if fin_distances:
                    ax_dist.hist(fin_distances, bins=20, alpha=0.7, label='Min Distance', color='skyblue')
                    # 동적 longMD threshold 적용
                    ax_dist.axvline(x=self.long_md_threshold, color='red', linestyle='--', 
                                   label=f'Long MD Threshold ({self.long_md_threshold}Å)')
                    ax_dist.set_xlabel('Distance (Å)')
                    ax_dist.set_ylabel('Frequency')
                    ax_dist.set_title('Minimum Distance Distribution (Successful)')
                    ax_dist.legend()
                    ax_dist.yaxis.set_major_locator(MAXNLocator(integer=True))  # 정수 축
                    ax_dist.grid(True, alpha=0.3)
                else:
                    ax_dist.text(0.5, 0.5, 'No distance data', ha='center', va='center', transform=ax_dist.transAxes)
                    ax_dist.set_title('Minimum Distance Distribution')
            else:
                ax_dist.text(0.5, 0.5, 'No distance data', ha='center', va='center', transform=ax_dist.transAxes)
                ax_dist.set_title('Minimum Distance Distribution')
            
            # N+3. Iteration별 성공/실패 횟수
            ax_success = plt.subplot(rows, cols, current_subplot)
            self._plot_iteration_success_failure_counts(ax_success)
            
            plt.tight_layout()
            
            # 그래프 저장
            plot_path = os.path.join(self.output_dir, "comprehensive_analysis.png")
            plt.savefig(plot_path, dpi=300, bbox_inches='tight')
            plt.close()
            
            print(f"\n📊 Comprehensive analysis graph saved: {plot_path}")
            
        except Exception as e:
            print(f"그래프 생성 중 오류: {e}")
            import traceback
            traceback.print_exc()

    def _calculate_distance_graph_count(self):
        """거리 그래프 개수 계산"""
        try:
            # 구조별 거리 데이터 수집
            structure_distance_data = self._collect_structure_distance_data()
            
            if not structure_distance_data:
                return 1  # 데이터가 없어도 빈 그래프 1개는 표시
            
            # PDB별로 그룹화
            pdb_groups = defaultdict(list)
            for structure_name, distance_data in structure_distance_data.items():
                pdb_code = structure_name.split('_')[0] if '_' in structure_name else structure_name
                pdb_groups[pdb_code].append((structure_name, distance_data))
            
            graph_count = 0
            for pdb_code, structures in pdb_groups.items():
                # 같은 PDB에서 나온 구조가 5개를 넘어도 1개 그래프에 모두 표현
                graph_count += 1
            
            # 전체 구조를 5개씩 나누어서 추가 그래프 계산
            remaining_structures = []
            for pdb_code, structures in pdb_groups.items():
                if len(structures) <= 5:
                    continue  # 이미 위에서 계산됨
                remaining_structures.extend(structures)
            
            # 나머지 구조들을 5개씩 그룹화
            if len(remaining_structures) > 0:
                additional_graphs = max(0, (len(remaining_structures) - 1) // 5)
                # 이미 PDB별로 1개씩은 계산했으므로, 추가분만 더함
                pass  # PDB별로 1개 그래프에 모든 구조 표현하므로 추가 그래프 없음
            
            return max(1, graph_count)
            
        except Exception as e:
            print(f"거리 그래프 개수 계산 오류: {e}")
            return 1

    def _collect_structure_distance_data(self):
        """구조별 거리 데이터 수집 (iteration 번호 수정)"""
        structure_distance_data = defaultdict(list)
        
        try:
            for attempt in self.data['attempt_results']:
                if not attempt.get('success', False) or 'min_distance' not in attempt:
                    continue
                
                structure_name = attempt.get('structure_name', 'Unknown')
                pdb_code = attempt.get('pdb_code', 'Unknown')
                # min_distance = attempt['min_distance']
                final_distance = attempt.get('final_distance', float('inf'))  # 추가
                
                # iteration 번호 추출 - 여러 방법으로 시도
                iteration = self._extract_iteration_number(attempt)
                
                if iteration is None or iteration < 1:
                    continue  # 유효하지 않은 iteration은 건너뜀
                
                # 구조 이름에 PDB 코드 포함
                full_structure_name = f"{pdb_code}_{structure_name}" if pdb_code != 'Unknown' else structure_name
                
                structure_distance_data[full_structure_name].append((iteration, final_distance))
            
            # 각 구조별 데이터를 iteration 순으로 정렬 및 중복 제거
            for structure_name in structure_distance_data:
                # iteration 순으로 정렬
                structure_distance_data[structure_name].sort(key=lambda x: x[0])
                
                # 같은 iteration의 중복 데이터가 있다면 평균값 사용
                unique_data = {}
                for iteration, distance in structure_distance_data[structure_name]:
                    if iteration in unique_data:
                        # 평균값 계산
                        unique_data[iteration] = (unique_data[iteration] + distance) / 2
                    else:
                        unique_data[iteration] = distance
                
                # 다시 리스트로 변환
                structure_distance_data[structure_name] = [(iter_num, dist) for iter_num, dist in sorted(unique_data.items())]
            
            # 최소 2개 이상의 데이터 포인트가 있는 구조만 반환
            filtered_data = {name: data for name, data in structure_distance_data.items() if len(data) >= 2}
            
            if filtered_data:
                print(f"거리 데이터 수집 완료: {len(filtered_data)}개 구조")
                for name, data in list(filtered_data.items())[:3]:  # 처음 3개만 출력
                    iterations = [d[0] for d in data]
                    print(f"  {name}: iterations {min(iterations)}-{max(iterations)} ({len(data)} points)")
            
            return filtered_data
            
        except Exception as e:
            print(f"구조별 거리 데이터 수집 오류: {e}")
            import traceback
            traceback.print_exc()
            return {}
        
    def _extract_iteration_number(self, attempt):
        """attempt 데이터에서 iteration 번호 추출"""
        try:
            # 방법 1: attempt 데이터에 직접 기록된 iteration 번호
            if 'iteration' in attempt and attempt['iteration'] is not None:
                iteration = int(attempt['iteration'])
                if iteration >= 1:
                    return iteration
            
            # 방법 2: 파일 경로에서 iteration 추출
            file_path = attempt.get('file_path', '')
            if file_path:
                # 경로에서 iteration_N 패턴 찾기
                import re
                match = re.search(r'iteration_(\d+)', file_path)
                if match:
                    iteration = int(match.group(1))
                    if iteration >= 1:
                        return iteration
            
            # 방법 3: attempt 번호를 기반으로 iteration 추정
            attempt_num = attempt.get('attempt', 0)
            if attempt_num >= 1:
                # attempt 번호가 있다면 이를 기반으로 iteration 추정
                # 보통 iteration은 1부터 시작
                return attempt_num  # attempt가 곧 iteration일 수 있음
            
            return None
            
        except Exception as e:
            print(f"Iteration 번호 추출 오류: {e}")
            return None

    def _plot_distance_graphs(self, fig, rows, cols, start_subplot):
        """구조별 iteration 순차적 거리 변화 그래프들 생성 - 순서 개선"""
        try:
            structure_distance_data = self._collect_structure_distance_data()
            
            if not structure_distance_data:
                # 데이터가 없는 경우 빈 그래프 1개 생성
                ax = plt.subplot(rows, cols, start_subplot)
                ax.text(0.5, 0.5, 'No distance data available', ha='center', va='center', transform=ax.transAxes)
                ax.set_title('Distance Progress by Structure')
                return start_subplot + 1
            
            # PDB별로 그룹화 및 구조 정렬 개선
            pdb_groups = defaultdict(list)
            for structure_name, distance_data in structure_distance_data.items():
                pdb_code = structure_name.split('_')[0] if '_' in structure_name else structure_name
                pdb_groups[pdb_code].append((structure_name, distance_data))
            
            # PDB별로 구조들을 일관된 순서로 정렬
            for pdb_code in pdb_groups:
                pdb_groups[pdb_code].sort(key=lambda x: x[0])  # 구조명 기준 정렬
            
            # 미리 정의된 색상 팔레트 - 더 구분되는 색상들
            colors = [
                '#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd',
                '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf',
                '#aec7e8', '#ffbb78', '#98df8a', '#ff9896', '#c5b0d5',
                '#c49c94', '#f7b6d3', '#c7c7c7', '#dbdb8d', '#9edae5'
            ]
            
            current_subplot = start_subplot
            
            # PDB별로 그래프 생성 (PDB 코드 순서대로 정렬)
            for pdb_idx, (pdb_code, structures) in enumerate(sorted(pdb_groups.items())):
                ax = plt.subplot(rows, cols, current_subplot)
                current_subplot += 1
                
                plotted_count = 0
                
                # 같은 PDB의 모든 구조를 하나의 그래프에 표시 (정렬된 순서로)
                for struct_idx, (structure_name, distance_data) in enumerate(structures):
                    if len(distance_data) >= 2:  # 최소 2개 데이터 포인트 필요
                        iterations, distances = zip(*distance_data)
                        
                        # 구조 이름에서 PDB 코드 제거 (표시용)
                        display_name = structure_name.replace(f"{pdb_code}_", "") if structure_name.startswith(f"{pdb_code}_") else structure_name
                        
                        # 구조별로 일관된 색상 할당 (struct_idx 기준)
                        color = colors[struct_idx % len(colors)]
                        
                        ax.plot(iterations, distances, 'o-', label=display_name, 
                               color=color, linewidth=2, markersize=4)
                        plotted_count += 1
                
                if plotted_count > 0:
                    # iteration 축 범위 설정
                    all_iterations = []
                    for structure_name, distance_data in structures:
                        if len(distance_data) >= 2:
                            iterations, _ = zip(*distance_data)
                            all_iterations.extend(iterations)
                    
                    if all_iterations:
                        min_iter = max(1, min(all_iterations))
                        max_iter = max(all_iterations)
                        ax.set_xlim(min_iter - 0.5, max_iter + 0.5)
                        # x축 틱을 정수로만 설정
                        ax.set_xticks(range(int(min_iter), int(max_iter) + 1))
                    
                    ax.set_xlabel('Iteration')
                    ax.set_ylabel('Distance (Å)')
                    ax.set_title(f'Distance Progress - {pdb_code} ({plotted_count} structures)')
                    ax.axhline(y=self.long_md_threshold, color='red', linestyle='--', 
                                   label=f'Long MD Threshold ({self.long_md_threshold}Å)')
                    ax.xaxis.set_major_locator(MAXNLocator(integer=True))  # 정수 축
                    ax.grid(True, alpha=0.3)
                    
                    # 범례 설정 (구조가 많으면 작게)
                    if plotted_count <= 5:
                        ax.legend()
                    else:
                        ax.legend(fontsize='small', ncol=2)
                else:
                    ax.text(0.5, 0.5, f'No valid data for {pdb_code}', ha='center', va='center', transform=ax.transAxes)
                    ax.set_title(f'Distance Progress - {pdb_code}')
            
            return current_subplot
            
        except Exception as e:
            print(f"거리 그래프 생성 오류: {e}")
            # 오류 발생 시 빈 그래프 1개 생성
            ax = plt.subplot(rows, cols, start_subplot)
            ax.text(0.5, 0.5, f'Error: {str(e)}', ha='center', va='center', transform=ax.transAxes)
            ax.set_title('Distance Progress by Structure')
            return start_subplot + 1

    def _plot_iteration_success_failure_counts(self, ax):
        """Iteration별 성공/실패 횟수 그래프"""
        try:
            if not self.data['iteration_results']:
                ax.text(0.5, 0.5, 'No iteration data', ha='center', va='center', transform=ax.transAxes)
                ax.set_title('Success/Failure Count by Iteration')
                return
            
            # iteration별 성공/실패 횟수 계산
            iteration_counts = defaultdict(lambda: {'success': 0, 'failure': 0})
            
            for iteration in self.data['iteration_results']:
                iter_num = iteration.get('iteration', 0)
                if iteration.get('success', False):
                    iteration_counts[iter_num]['success'] += 1
                else:
                    iteration_counts[iter_num]['failure'] += 1
            
            if iteration_counts:
                iteration_nums = sorted(iteration_counts.keys())
                success_counts = [iteration_counts[i]['success'] for i in iteration_nums]
                failure_counts = [iteration_counts[i]['failure'] for i in iteration_nums]
                
                width = 0.35
                x_pos = np.arange(len(iteration_nums))
                
                ax.bar(x_pos - width/2, success_counts, width, alpha=0.7, label='Successful', color='green')
                ax.bar(x_pos + width/2, failure_counts, width, alpha=0.7, label='Failed', color='red')
                
                ax.set_xlabel('Iteration Number')
                ax.set_ylabel('Count')
                ax.set_title('Success/Failure Count by Iteration')
                ax.set_xticks(x_pos)
                ax.set_xticklabels(iteration_nums)
                ax.yaxis.set_major_locator(MaxNLocator(integer=True))  # 정수 축
                ax.yaxis.set_major_locator(MAXNLocator(integer=True))  # 정수 축
                ax.legend()
                ax.grid(True, alpha=0.3)
            else:
                ax.text(0.5, 0.5, 'No data available', ha='center', va='center', transform=ax.transAxes)
                ax.set_title('Success/Failure Count by Iteration')
        except Exception as e:
            ax.text(0.5, 0.5, f'Error: {str(e)}', ha='center', va='center', transform=ax.transAxes)
            ax.set_title('Success/Failure Count by Iteration')

    def print_comprehensive_summary(self):
        """종합 결과 요약 출력 - longMD 기준값 표시"""
        print("\n" + "="*80)
        print("Simple SuMD Batch Execution Comprehensive Analysis Results")
        print("="*80)
        
        # 설정 정보 표시
        print(f"Analysis Configuration:")
        print(f"  • Long MD threshold: {self.long_md_threshold} Å")
        print()
        
        # 전체 통계
        total_pdbs = len(self.data['pdb_results'])
        successful_pdbs = len([p for p in self.data['pdb_results'] if p.get('has_successful_result', False)])
        total_structures = len(self.data['structure_results'])
        successful_structures = len([s for s in self.data['structure_results'] if s.get('successful_iterations', 0) > 0])
        total_iterations = len(self.data['iteration_results'])
        successful_iterations = len([i for i in self.data['iteration_results'] if i.get('success', False)])
        total_attempts = len(self.data['attempt_results'])
        successful_attempts = len([a for a in self.data['attempt_results'] if a.get('success', False)])
        
        print(f"📊 Overall Statistics:")
        print(f"  • PDB files: {successful_pdbs}/{total_pdbs} successful ({successful_pdbs/max(1,total_pdbs)*100:.1f}%)")
        print(f"  • Structures: {successful_structures}/{total_structures} successful ({successful_structures/max(1,total_structures)*100:.1f}%)")
        print(f"  • Iterations: {successful_iterations}/{total_iterations} successful ({successful_iterations/max(1,total_iterations)*100:.1f}%)")
        print(f"  • Attempts: {successful_attempts}/{total_attempts} successful ({successful_attempts/max(1,total_attempts)*100:.1f}%)")
        
        # PDB 정보가 없는 경우 알림
        if total_pdbs == 0 and total_structures > 0:
            print("\n⚠️  Note: PDB-level statistics are unavailable.")
            print("   This may happen when analyzing ongoing runs or missing batch_summary.json")
            
            # 폴더에서 추정된 PDB 수 계산
            unique_pdbs = set()
            for structure in self.data['structure_results']:
                pdb_code = structure.get('pdb_code', structure.get('pdb_name', ''))
                if pdb_code:
                    unique_pdbs.add(pdb_code)
            
            if unique_pdbs:
                print(f"   Estimated PDBs from folder structure: {len(unique_pdbs)} unique PDBs")
                print(f"   PDB codes: {', '.join(sorted(unique_pdbs))}")
        
        # 시간 통계
        if self.data['attempt_results']:
            attempt_times = [a.get('estimated_duration', 0) for a in self.data['attempt_results']]
            print(f"\n⏱️  Timing Statistics:")
            print(f"  • Average attempt time: {np.mean(attempt_times):.1f}s")
            print(f"  • Median attempt time: {np.median(attempt_times):.1f}s")
            print(f"  • Min/Max attempt time: {np.min(attempt_times):.1f}s / {np.max(attempt_times):.1f}s")
            
            total_estimated_time = sum(attempt_times)
            print(f"  • Total estimated time: {total_estimated_time/3600:.2f} hours")
        
        # 긴 MD 통계 - 동적 threshold 적용
        long_md_attempts = [a for a in self.data['attempt_results'] if a.get('long_md_executed', False)]
        close_contact_attempts = [a for a in self.data['attempt_results'] 
                                if a.get('min_distance', float('inf')) <= self.long_md_threshold]
        
        print(f"\n🔬 Long MD Statistics (Threshold: {self.long_md_threshold}Å):")
        print(f"  • Close contact detected: {len(close_contact_attempts)} times ({len(close_contact_attempts)/max(1,total_attempts)*100:.1f}%)")
        print(f"  • Long MD executed: {len(long_md_attempts)} times ({len(long_md_attempts)/max(1,total_attempts)*100:.1f}%)")


    # 다른 함수들도 동일하게 영어로 변경 (기존 코드 유지하되 출력 메시지만 영어로 변경)
    def run_complete_analysis(self):
        """전체 분석 실행 (개선된 버전)"""
        print("=== Simple SuMD Batch Execution Comprehensive Analysis Started ===")
        
        if not os.path.exists(self.output_dir):
            raise FileNotFoundError(f"Output directory not found: {self.output_dir}")
        
        # 1단계: 통계 파일 우선 로드
        print("Step 1: Loading summary files...")
        batch_loaded = self.load_batch_summary()
        log_loaded = self.load_execution_log()
        
        statistics_available = batch_loaded and len(self.data['structure_results']) > 0
        
        if statistics_available:
            print("✓ Data loaded from summary files")
        else:
            print("! Summary files incomplete - supplementing with folder scan")
        
        # 2단계: 상세 파일 및 폴더 스캔
        if not statistics_available or not self.data['structure_results']:
            print("Step 2: Direct folder structure scan...")
            self.load_detailed_results()
        else:
            print("Step 2: Supplementing with folder timing information...")
            self._supplement_with_folder_timing()
        
        # 데이터 유효성 검증
        if not self.data['structure_results'] and not self.data['iteration_results'] and not self.data['attempt_results']:
            raise ValueError("No analyzable data found. Please check the output directory.")
        
        # 3단계: 소요시간 통계 계산
        print("Step 3: Calculating timing statistics...")
        self.calculate_timing_statistics()
        
        # 4단계: 결과 출력
        print("Step 4: Outputting analysis results...")
        self.print_comprehensive_summary()
        self.analyze_timing_patterns()
        self.analyze_failure_patterns()
        
        # 5단계: 그래프 및 보고서 생성
        print("Step 5: Generating visualization and reports...")
        try:
            self.plot_comprehensive_analysis()
        except Exception as e:
            print(f"Graph generation failed: {e}")
        
        self.save_detailed_report()
        
        print(f"\n🎉 Comprehensive analysis complete! Results saved in {self.output_dir}")
        print(f"   📊 Graph: comprehensive_analysis.png")
        print(f"   📄 Report: comprehensive_analysis_report.txt")

    # 기존의 다른 함수들은 그대로 유지...
    def _calculate_single_attempt_duration(self, attempt_dir):
        """단일 attempt 디렉토리의 소요시간 계산 - 명확한 규칙 적용"""
        try:
            # attempt 시작시간: input.pdb 파일 시간
            input_pdb = os.path.join(attempt_dir, 'input.pdb')
            start_time = None
            
            if os.path.exists(input_pdb):
                stat = os.stat(input_pdb)
                start_time = datetime.fromtimestamp(stat.st_ctime)
            
            # attempt 종료시간: md.gro 파일의 시간
            md_gro = os.path.join(attempt_dir, 'md.gro')
            end_time = None
            
            if os.path.exists(md_gro):
                stat = os.stat(md_gro)
                end_time = datetime.fromtimestamp(stat.st_mtime)
            
            # 둘 다 있으면 정확한 시간 계산
            if start_time and end_time:
                duration = (end_time - start_time).total_seconds()
                return max(0, duration)
            
            # 백업: 기존 방식 사용
            return self._fallback_attempt_duration_calculation(attempt_dir)
            
        except Exception as e:
            print(f"Attempt duration 계산 실패 {attempt_dir}: {e}")
            return 0

    def _fallback_attempt_duration_calculation(self, attempt_dir):
        """백업용 attempt 시간 계산"""
        try:
            key_files = ['md.gro', 'md.xtc', 'md.log', 'npt.gro', 'nvt.gro', 'em.gro']
            timestamps = []
            
            for filename in key_files:
                filepath = os.path.join(attempt_dir, filename)
                if os.path.exists(filepath):
                    stat = os.stat(filepath)
                    timestamps.append(datetime.fromtimestamp(stat.st_mtime))
            
            if len(timestamps) >= 2:
                duration = (max(timestamps) - min(timestamps)).total_seconds()
                return max(0, duration)
            
            return 0
        except:
            return 0

    def load_batch_summary(self):
        """batch_summary.json 로드"""
        batch_file = os.path.join(self.output_dir, "batch_summary.json")
        if os.path.exists(batch_file):
            try:
                with open(batch_file, encoding='utf-8') as f:
                    summary = json.load(f)
                    self.data['batch_info'] = summary.get('batch_execution_summary', {})
                    self.data['pdb_results'] = summary.get('pdb_results', [])
                    self.data['structure_results'] = summary.get('detailed_structure_results', [])
                print(f"배치 요약 로드 완료: {len(self.data['pdb_results'])}개 PDB, {len(self.data['structure_results'])}개 구조")
                return True
            except Exception as e:
                print(f"batch_summary.json 로드 실패: {e}")
        return False
    
    def load_execution_log(self):
        """execution_log.json 로드"""
        log_file = os.path.join(self.output_dir, "execution_log.json")
        if os.path.exists(log_file):
            try:
                with open(log_file, encoding='utf-8') as f:
                    self.data['execution_log'] = json.load(f)
                print("실행 로그 로드 완료")
                return True
            except Exception as e:
                print(f"execution_log.json 로드 실패: {e}")
        return False

    def calculate_timing_statistics(self):
        """소요시간 통계 계산"""
        print("소요시간 통계 계산 중...")
        
        # Attempt별 소요시간 추정 (GROMACS 단계별 기본 추정치 사용)
        stage_estimates = {
            'pdb2gmx': 5, 'editconf': 2, 'solvate': 10,
            'ions_grompp': 3, 'genion': 5, 'em': 60,
            'nvt': 120, 'npt': 120, 'npt2': 120, 'MD': 300, '긴 MD': 1800
        }
        
        for attempt in self.data['attempt_results']:
            total_duration = 0
            stages = attempt.get('stages', [])
            
            for stage in stages:
                stage_name = stage.get('stage', '')
                if stage.get('success', False):
                    estimate = stage_estimates.get(stage_name, 30)
                    total_duration += estimate
            
            # 긴 MD 추가 시간
            if attempt.get('long_md_executed', False):
                total_duration += stage_estimates['긴 MD']
            
            attempt['estimated_duration'] = total_duration
        
        # Iteration별 소요시간 계산
        for iteration in self.data['iteration_results']:
            related_attempts = [a for a in self.data['attempt_results'] 
                              if (a['pdb_name'] == iteration['pdb_name'] and 
                                  a['structure_name'] == iteration['structure_name'] and 
                                  a.get('iteration', 0) == iteration.get('iteration', 0))]
            
            total_time = sum(a.get('estimated_duration', 0) for a in related_attempts)
            iteration['estimated_duration'] = total_time
            iteration['attempt_count'] = len(related_attempts)
        
        # 구조별 소요시간 계산
        for structure in self.data['structure_results']:
            related_iterations = [i for i in self.data['iteration_results']
                                if (i['pdb_name'] == structure['pdb_name'] and 
                                    i['structure_name'] == structure['structure_name'])]
            
            total_time = sum(i.get('estimated_duration', 0) for i in related_iterations)
            structure['estimated_total_duration'] = total_time
            structure['total_iteration_count'] = len(related_iterations)
            structure['total_attempt_count'] = sum(i.get('attempt_count', 0) for i in related_iterations)

    def _supplement_with_folder_timing(self):
        """통계 파일 기반 데이터에 폴더 타이밍 정보 보완 - 개선된 버전"""
        print("폴더 타이밍 정보 보완 중 (명확한 규칙 적용)...")
        
        # 구조별 실제 시간 정보 추가
        for structure in self.data['structure_results']:
            if 'folder_path' in structure:
                folder_path = structure['folder_path']
            else:
                # 경로 추정
                pdb_name = structure.get('pdb_code', '')
                structure_name = structure.get('structure_name', '')
                folder_path = os.path.join(self.output_dir, f"{pdb_name}*", f"structure_{structure_name}")
                matching_folders = glob.glob(folder_path)
                folder_path = matching_folders[0] if matching_folders else None
            
            if folder_path and os.path.exists(folder_path):
                actual_duration, iter_durations = self.calculate_actual_duration_from_folders(folder_path)
                structure['actual_total_duration'] = actual_duration
                structure['detailed_iteration_timing'] = iter_durations
                
                # 추가 타이밍 정보
                timing_info = self.extract_folder_timing_info(folder_path)
                structure['folder_timing'] = timing_info

    def calculate_actual_duration_from_folders(self, folder_path):
        """폴더 구조에서 실제 소요시간 계산 - 개선된 버전"""
        try:
            # 구조 레벨 시간 계산 (요구사항에 따른 명확한 규칙 적용)
            total_duration, start_time, end_time = self.calculate_structure_duration_from_folder(folder_path)
            
            # iteration별 상세 시간도 함께 계산
            iteration_durations = []
            iteration_folders = sorted([d for d in glob.glob(os.path.join(folder_path, "iteration_*")) if os.path.isdir(d)])
            
            for iter_folder in iteration_folders:
                iter_duration, iter_start, iter_end = self.calculate_iteration_duration_from_folder(iter_folder)
                iteration_durations.append({
                    'folder': os.path.basename(iter_folder),
                    'duration': iter_duration,
                    'start_time': iter_start,
                    'end_time': iter_end
                })
            
            return total_duration, iteration_durations
            
        except Exception as e:
            print(f"폴더 duration 계산 실패 {folder_path}: {e}")
            return 0, []

    def calculate_structure_duration_from_folder(self, structure_folder):
        """구조 폴더에서 실제 소요시간 계산 - 명확한 규칙 적용"""
        try:
            if not os.path.exists(structure_folder):
                return 0, None, None
            
            # PDB에서 파생된 구조의 시작시간: iteration_1/attempt_1/input.pdb 파일 시간
            iter1_attempt1_dir = os.path.join(structure_folder, "iteration_1", "attempt_1")
            start_time = None
            
            if os.path.exists(iter1_attempt1_dir):
                input_pdb = os.path.join(iter1_attempt1_dir, 'input.pdb')
                if os.path.exists(input_pdb):
                    stat = os.stat(input_pdb)
                    start_time = datetime.fromtimestamp(stat.st_ctime)
            
            # 구조 종료시간: final_structure.pdb로 하되 없다면 마지막 iteration 폴더의 시간
            end_time = None
            
            # 1순위: final_structure.pdb
            final_structure = os.path.join(structure_folder, "final_structure.pdb")
            if os.path.exists(final_structure):
                stat = os.stat(final_structure)
                end_time = datetime.fromtimestamp(stat.st_mtime)
            else:
                # 2순위: 마지막 iteration 폴더의 시간
                iteration_folders = sorted([d for d in glob.glob(os.path.join(structure_folder, "iteration_*")) if os.path.isdir(d)])
                if iteration_folders:
                    last_iter_folder = iteration_folders[-1]
                    stat = os.stat(last_iter_folder)
                    end_time = datetime.fromtimestamp(stat.st_mtime)
            
            # 정확한 시간 계산
            if start_time and end_time:
                duration = (end_time - start_time).total_seconds()
                return max(0, duration), start_time, end_time
            
            # 백업: iteration들의 이합으로 계산
            iteration_folders = [d for d in glob.glob(os.path.join(structure_folder, "iteration_*")) if os.path.isdir(d)]
            total_duration = 0
            
            for iter_folder in iteration_folders:
                duration, _, _ = self.calculate_iteration_duration_from_folder(iter_folder)
                total_duration += duration
            
            return total_duration, start_time, end_time
            
        except Exception as e:
            print(f"Structure duration 계산 실패 {structure_folder}: {e}")
            return 0, None, None

    def calculate_iteration_duration_from_folder(self, iteration_folder):
        """iteration 폴더에서 실제 소요시간 계산 - 명확한 규칙 적용"""
        try:
            if not os.path.exists(iteration_folder):
                return 0, None, None
            
            # iteration 시작시간: attempt_1의 input.pdb 파일 시간
            attempt1_dir = os.path.join(iteration_folder, "attempt_1")
            start_time = None
            
            if os.path.exists(attempt1_dir):
                input_pdb = os.path.join(attempt1_dir, 'input.pdb')
                if os.path.exists(input_pdb):
                    stat = os.stat(input_pdb)
                    start_time = datetime.fromtimestamp(stat.st_ctime)
            
            # iteration 종료시간: next_structure.pdb의 파일 시간
            next_structure = os.path.join(os.path.dirname(iteration_folder), "next_structure.pdb")
            end_time = None
            
            if os.path.exists(next_structure):
                stat = os.stat(next_structure)
                end_time = datetime.fromtimestamp(stat.st_mtime)
            
            # 정확한 시간 계산
            if start_time and end_time:
                duration = (end_time - start_time).total_seconds()
                return max(0, duration), start_time, end_time
            
            # 백업: attempt들의 이합으로 계산
            attempt_dirs = [d for d in glob.glob(os.path.join(iteration_folder, "attempt_*")) if os.path.isdir(d)]
            total_duration = 0
            
            for attempt_dir in attempt_dirs:
                duration = self._calculate_single_attempt_duration(attempt_dir)
                total_duration += duration
            
            return total_duration, start_time, end_time
            
        except Exception as e:
            print(f"Iteration duration 계산 실패 {iteration_folder}: {e}")
            return 0, None, None

    def extract_folder_timing_info(self, folder_path):
        """폴더에서 타이밍 정보 추출"""
        timing_info = {
            'start_time': None,
            'end_time': None,
            'duration_seconds': 0,
            'iteration_count': 0,
            'file_timestamps': []
        }
        
        try:
            if not os.path.exists(folder_path):
                return timing_info
            
            # 모든 파일의 타임스탬프 수집
            all_files = []
            for root, dirs, files in os.walk(folder_path):
                for file in files:
                    file_path = os.path.join(root, file)
                    try:
                        stat = os.stat(file_path)
                        all_files.append({
                            'path': file_path,
                            'created': datetime.fromtimestamp(stat.st_ctime),
                            'modified': datetime.fromtimestamp(stat.st_mtime),
                            'size': stat.st_size
                        })
                    except:
                        continue
            
            if all_files:
                # 시작/종료 시간 계산
                all_files.sort(key=lambda x: x['created'])
                timing_info['start_time'] = all_files[0]['created']
                timing_info['end_time'] = max(f['modified'] for f in all_files)
                timing_info['duration_seconds'] = (timing_info['end_time'] - timing_info['start_time']).total_seconds()
                timing_info['file_timestamps'] = [(f['path'], f['created'], f['modified']) for f in all_files]
            
            # iteration 개수 계산
            iteration_folders = [d for d in glob.glob(os.path.join(folder_path, "iteration_*")) if os.path.isdir(d)]
            timing_info['iteration_count'] = len(iteration_folders)
            
            return timing_info
            
        except Exception as e:
            print(f"타이밍 정보 추출 실패 {folder_path}: {e}")
            return timing_info

    def analyze_failure_patterns(self):
        """실패 패턴 분석"""
        print(f"\n🔍 Failure Pattern Analysis:")
        print("-" * 50)
        
        # GROMACS 단계별 실패 분석
        stage_stats = defaultdict(lambda: {'total': 0, 'success': 0})
        
        for attempt in self.data['attempt_results']:
            stages = attempt.get('stages', [])
            for stage in stages:
                stage_name = stage.get('stage', 'unknown')
                stage_stats[stage_name]['total'] += 1
                if stage.get('success', False):
                    stage_stats[stage_name]['success'] += 1
        
        print("GROMACS Stage Success Rates:")
        for stage, stats in sorted(stage_stats.items()):
            if stats['total'] > 0:
                success_rate = stats['success'] / stats['total'] * 100
                print(f"  • {stage}: {success_rate:.1f}% ({stats['success']}/{stats['total']})")
        
        # 실패 이유 분석
        failure_reasons = defaultdict(int)
        for attempt in self.data['attempt_results']:
            if not attempt.get('success', False):
                reason = attempt.get('reason', 'unknown')
                failure_reasons[reason] += 1
        
        if failure_reasons:
            print(f"\nFailure Reason Distribution:")
            total_failures = sum(failure_reasons.values())
            for reason, count in sorted(failure_reasons.items(), key=lambda x: x[1], reverse=True):
                print(f"  • {reason}: {count} times ({count/total_failures*100:.1f}%)")

    def analyze_timing_patterns(self):
        """시간 패턴 분석"""
        print(f"\n⏱️  Detailed Timing Pattern Analysis:")
        print("-" * 50)
        
        if not self.data['attempt_results']:
            print("No attempt data available for analysis.")
            return
        
        # 성공/실패별 시간 분석
        successful_attempts = [a for a in self.data['attempt_results'] if a.get('success', False)]
        failed_attempts = [a for a in self.data['attempt_results'] if not a.get('success', False)]
        
        if successful_attempts:
            success_times = [a.get('estimated_duration', 0) for a in successful_attempts]
            print(f"Successful attempt timing statistics:")
            print(f"  • Average: {np.mean(success_times):.1f}s")
            print(f"  • Median: {np.median(success_times):.1f}s")
            print(f"  • Std Dev: {np.std(success_times):.1f}s")
        
        if failed_attempts:
            fail_times = [a.get('estimated_duration', 0) for a in failed_attempts]
            print(f"Failed attempt timing statistics:")
            print(f"  • Average: {np.mean(fail_times):.1f}s")
            print(f"  • Median: {np.median(fail_times):.1f}s")
            print(f"  • Std Dev: {np.std(fail_times):.1f}s")

    def save_detailed_report(self):
        """상세 분석 보고서 저장 (한글 버전)"""
        report_path = os.path.join(self.output_dir, "comprehensive_analysis_report.txt")
        
        with open(report_path, 'w', encoding='utf-8') as f:
            f.write("Simple SuMD 배치 실행 종합 분석 보고서\n")
            f.write("=" * 80 + "\n\n")
            
            # 기본 통계
            f.write("1. 전체 실행 통계\n")
            f.write("-" * 40 + "\n")
            
            total_pdbs = len(self.data['pdb_results'])
            successful_pdbs = len([p for p in self.data['pdb_results'] if p.get('has_successful_result', False)])
            total_structures = len(self.data['structure_results'])
            successful_structures = len([s for s in self.data['structure_results'] if s.get('successful_iterations', 0) > 0])
            total_iterations = len(self.data['iteration_results'])
            successful_iterations = len([i for i in self.data['iteration_results'] if i.get('success', False)])
            total_attempts = len(self.data['attempt_results'])
            successful_attempts = len([a for a in self.data['attempt_results'] if a.get('success', False)])
            
            f.write(f"PDB 파일 수: {total_pdbs}개 (성공: {successful_pdbs}개, {successful_pdbs/max(1,total_pdbs)*100:.1f}%)\n")
            f.write(f"생성된 구조 수: {total_structures}개 (성공: {successful_structures}개, {successful_structures/max(1,total_structures)*100:.1f}%)\n")
            f.write(f"총 iteration 수: {total_iterations}개 (성공: {successful_iterations}개, {successful_iterations/max(1,total_iterations)*100:.1f}%)\n")
            f.write(f"총 attempt 수: {total_attempts}개 (성공: {successful_attempts}개, {successful_attempts/max(1,total_attempts)*100:.1f}%)\n\n")
            
            # 시간 통계
            if self.data['attempt_results']:
                attempt_times = [a.get('estimated_duration', 0) for a in self.data['attempt_results']]
                f.write("2. 소요시간 통계\n")
                f.write("-" * 40 + "\n")
                f.write(f"평균 attempt 시간: {np.mean(attempt_times):.1f}초\n")
                f.write(f"중간값 attempt 시간: {np.median(attempt_times):.1f}초\n")
                f.write(f"최단 attempt 시간: {np.min(attempt_times):.1f}초\n")
                f.write(f"최장 attempt 시간: {np.max(attempt_times):.1f}초\n")
                f.write(f"시간 표준편차: {np.std(attempt_times):.1f}초\n")
                f.write(f"전체 예상 소요시간: {sum(attempt_times)/3600:.2f}시간\n\n")
            
            # PDB별 상세 결과
            if self.data['pdb_results']:
                f.write("3. PDB별 상세 결과\n")
                f.write("-" * 40 + "\n")
                
                for pdb_data in sorted(self.data['pdb_results'], key=lambda x: x.get('max_successful_iterations', 0), reverse=True):
                    f.write(f"PDB: {pdb_data['pdb_code']} (체인: {pdb_data.get('chain1', 'A')}-{pdb_data.get('chain2', 'B')})\n")
                    f.write(f"  구조 수: {pdb_data['total_structures']}개 (성공: {pdb_data['successful_structures']}개)\n")
                    f.write(f"  최대 성공 iteration: {pdb_data['max_successful_iterations']}회\n")
                    f.write(f"  최고 구조: {pdb_data.get('best_structure_name', 'N/A')}\n")
                    f.write("\n")
            
            # 실패 패턴 분석
            f.write("4. 실패 패턴 분석\n")
            f.write("-" * 40 + "\n")
            
            stage_stats = defaultdict(lambda: {'total': 0, 'success': 0})
            for attempt in self.data['attempt_results']:
                stages = attempt.get('stages', [])
                for stage in stages:
                    stage_name = stage.get('stage', 'unknown')
                    stage_stats[stage_name]['total'] += 1
                    if stage.get('success', False):
                        stage_stats[stage_name]['success'] += 1
            
            f.write("GROMACS 단계별 성공률:\n")
            for stage, stats in sorted(stage_stats.items()):
                if stats['total'] > 0:
                    success_rate = stats['success'] / stats['total'] * 100
                    f.write(f"  {stage}: {success_rate:.1f}% ({stats['success']}/{stats['total']})\n")
            
            f.write("\n실패 이유 분포:\n")
            failure_reasons = defaultdict(int)
            for attempt in self.data['attempt_results']:
                if not attempt.get('success', False):
                    reason = attempt.get('reason', 'unknown')
                    failure_reasons[reason] += 1
            
            total_failures = sum(failure_reasons.values())
            for reason, count in sorted(failure_reasons.items(), key=lambda x: x[1], reverse=True):
                f.write(f"  {reason}: {count}회 ({count/total_failures*100:.1f}%)\n")
            
            # 긴 MD 통계
            f.write("\n5. 긴 MD 통계\n")
            f.write("-" * 40 + "\n")
            
            long_md_attempts = [a for a in self.data['attempt_results'] if a.get('long_md_executed', False)]
            close_contact_attempts = [a for a in self.data['attempt_results'] if a.get('close_contact_detected', False)]
            
            f.write(f"근접 접촉 감지: {len(close_contact_attempts)}회 ({len(close_contact_attempts)/max(1,total_attempts)*100:.1f}%)\n")
            f.write(f"긴 MD 실행: {len(long_md_attempts)}회 ({len(long_md_attempts)/max(1,total_attempts)*100:.1f}%)\n")
            
            if long_md_attempts:
                long_md_success = len([a for a in long_md_attempts if a.get('success', False)])
                regular_attempts = [a for a in self.data['attempt_results'] if not a.get('long_md_executed', False)]
                regular_success = len([a for a in regular_attempts if a.get('success', False)])
                
                f.write(f"긴 MD 성공률: {long_md_success}/{len(long_md_attempts)} ({long_md_success/len(long_md_attempts)*100:.1f}%)\n")
                if regular_attempts:
                    f.write(f"기본 MD 성공률: {regular_success}/{len(regular_attempts)} ({regular_success/len(regular_attempts)*100:.1f}%)\n")
        
        print(f"📄 상세 분석 보고서 저장: {report_path}")

def main():
    """메인 함수 - longMD 기준값 설정 가능"""
    if len(sys.argv) < 2 or len(sys.argv) > 3:
        print("Usage: python3 batch_analyzer.py <output_directory> [long_md_threshold]")
        print("Example: python3 batch_analyzer.py batch_sumd_output")
        print("Example: python3 batch_analyzer.py batch_sumd_output 8.0")
        print()
        print("Parameters:")
        print("  output_directory: Directory containing SuMD results")
        print("  long_md_threshold: Distance threshold for long MD (default: 10.0 Å)")
        sys.exit(1)
    
    output_dir = sys.argv[1]
    
    # longMD 기준값 설정 (기본값: 10.0 Å)
    long_md_threshold = 10.0
    if len(sys.argv) == 3:
        try:
            long_md_threshold = float(sys.argv[2])
            if long_md_threshold <= 0:
                raise ValueError("Threshold must be positive")
            print(f"Using custom Long MD threshold: {long_md_threshold} Å")
        except ValueError as e:
            print(f"Invalid threshold value: {sys.argv[2]} ({e})")
            print("Using default threshold: 10.0 Å")
            long_md_threshold = 10.0
    
    try:
        # 분석기 초기화 및 실행
        print("Initializing Batch Analyzer...")
        analyzer = BatchAnalyzer(output_dir, long_md_threshold=long_md_threshold)

        print("Running comprehensive analysis...")
        analyzer.run_complete_analysis()
        
    except Exception as e:
        print(f"Analysis error occurred: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)