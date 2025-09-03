#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Simple SuMD 성능 분석 및 시각화 모듈
각 attempt, iteration, 구조별 소요시간과 성능 지표를 분석하고 시각화합니다.
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
from collections import defaultdict
import pandas as pd

class PerformanceAnalyzer:
    def __init__(self, output_dir):
        """성능 분석기 초기화"""
        self.output_dir = output_dir
        self.performance_data = {
            'attempts': [],
            'iterations': [],
            'structures': [],
            'timeline': []
        }
        
    def parse_datetime(self, datetime_str):
        """다양한 형식의 datetime 문자열 파싱"""
        try:
            if 'T' in datetime_str:
                return datetime.fromisoformat(datetime_str.replace('Z', '+00:00'))
            else:
                return datetime.strptime(datetime_str, '%Y-%m-%d %H:%M:%S')
        except:
            return None
    
    def extract_attempt_data(self):
        """Attempt별 성능 데이터 추출"""
        print("Attempt 데이터 추출 중...")
        
        attempt_files = glob.glob(os.path.join(self.output_dir, "**", "iteration_*_attempt_*.json"), recursive=True)
        
        for attempt_file in attempt_files:
            try:
                with open(attempt_file) as f:
                    data = json.load(f)
                
                # 파일명에서 정보 추출
                basename = os.path.basename(attempt_file)
                match = re.search(r'iteration_(\d+)_attempt_(\d+)\.json', basename)
                if not match:
                    continue
                
                iteration_num = int(match.group(1))
                attempt_num = int(match.group(2))
                
                # 구조 정보 추출 (경로에서)
                structure_match = re.search(r'structure_([^/]+)', attempt_file)
                structure_name = structure_match.group(1) if structure_match else "unknown"
                
                # 시간 정보는 stages에서 추출
                stages = data.get('stages', [])
                total_duration = 0
                stage_durations = {}
                
                # 파일 생성/수정 시간으로 추정
                file_stats = os.stat(attempt_file)
                estimated_duration = 0
                
                # 기본 추정치 (단계별 평균 시간)
                stage_estimates = {
                    'pdb2gmx': 5, 'editconf': 2, 'solvate': 10,
                    'ions_grompp': 3, 'genion': 5, 'em': 60,
                    'nvt': 120, 'npt': 120, 'MD': 300, '긴 MD': 1800
                }
                
                for stage in stages:
                    stage_name = stage.get('stage', 'unknown')
                    if stage.get('success', False):
                        estimate = stage_estimates.get(stage_name, 30)
                        stage_durations[stage_name] = estimate
                        total_duration += estimate
                
                # 긴 MD 여부 확인
                is_long_md = data.get('long_md_executed', False)
                if is_long_md:
                    total_duration += stage_estimates['긴 MD']
                
                attempt_data = {
                    'structure': structure_name,
                    'iteration': iteration_num,
                    'attempt': attempt_num,
                    'success': data.get('success', False),
                    'duration_seconds': total_duration,
                    'stage_durations': stage_durations,
                    'slope': data.get('slope', None),
                    'min_distance': data.get('min_distance', None),
                    'long_md_executed': is_long_md,
                    'close_contact_detected': data.get('close_contact_detected', False),
                    'file_path': attempt_file,
                    'file_time': datetime.fromtimestamp(file_stats.st_mtime)
                }
                
                self.performance_data['attempts'].append(attempt_data)
                
            except Exception as e:
                print(f"Attempt 파일 처리 실패 {attempt_file}: {e}")
        
        print(f"총 {len(self.performance_data['attempts'])}개 attempt 데이터 추출 완료")
    
    def extract_iteration_data(self):
        """Iteration별 성능 데이터 추출"""
        print("Iteration 데이터 추출 중...")
        
        iteration_files = glob.glob(os.path.join(self.output_dir, "**", "iteration_*_summary.json"), recursive=True)
        
        for iter_file in iteration_files:
            try:
                with open(iter_file) as f:
                    data = json.load(f)
                
                # 파일명에서 정보 추출
                basename = os.path.basename(iter_file)
                match = re.search(r'iteration_(\d+)_summary\.json', basename)
                if not match:
                    continue
                
                iteration_num = int(match.group(1))
                
                # 구조 정보 추출
                structure_match = re.search(r'structure_([^/]+)', iter_file)
                structure_name = structure_match.group(1) if structure_match else "unknown"
                
                # 해당 iteration의 모든 attempt 찾기
                attempts_in_iter = [a for a in self.performance_data['attempts'] 
                                  if a['structure'] == structure_name and a['iteration'] == iteration_num]
                
                total_duration = sum(a['duration_seconds'] for a in attempts_in_iter)
                
                iteration_data = {
                    'structure': structure_name,
                    'iteration': iteration_num,
                    'success': data.get('success', False),
                    'attempts_used': data.get('attempts_used', 0),
                    'duration_seconds': total_duration,
                    'avg_attempt_duration': total_duration / len(attempts_in_iter) if attempts_in_iter else 0,
                    'long_md': data.get('long_md', False),
                    'final_slope': data.get('final_result', {}).get('slope', None) if data.get('final_result') else None,
                    'final_min_distance': data.get('final_result', {}).get('min_distance', None) if data.get('final_result') else None
                }
                
                self.performance_data['iterations'].append(iteration_data)
                
            except Exception as e:
                print(f"Iteration 파일 처리 실패 {iter_file}: {e}")
        
        print(f"총 {len(self.performance_data['iterations'])}개 iteration 데이터 추출 완료")
    
    def extract_structure_data(self):
        """구조별 성능 데이터 추출"""
        print("구조별 데이터 집계 중...")
        
        structure_files = glob.glob(os.path.join(self.output_dir, "structure_*", "structure_results.json"))
        
        for struct_file in structure_files:
            try:
                with open(struct_file) as f:
                    data = json.load(f)
                
                structure_name = data.get('structure_name', 'unknown')
                
                # 시간 정보 파싱
                start_time = self.parse_datetime(data.get('start_time', ''))
                end_time = self.parse_datetime(data.get('end_time', ''))
                
                total_duration = 0
                if start_time and end_time:
                    total_duration = (end_time - start_time).total_seconds()
                else:
                    # attempt들의 시간 합산으로 추정
                    structure_attempts = [a for a in self.performance_data['attempts'] if a['structure'] == structure_name]
                    total_duration = sum(a['duration_seconds'] for a in structure_attempts)
                
                # 구조별 통계 계산
                structure_iterations = [i for i in self.performance_data['iterations'] if i['structure'] == structure_name]
                structure_attempts = [a for a in self.performance_data['attempts'] if a['structure'] == structure_name]
                
                structure_data = {
                    'structure_name': structure_name,
                    'assigned_gpu': data.get('assigned_gpu', 'unknown'),
                    'total_iterations': data.get('total_iterations', 0),
                    'successful_iterations': data.get('successful_iterations', 0),
                    'total_attempts': len(structure_attempts),
                    'successful_attempts': len([a for a in structure_attempts if a['success']]),
                    'total_duration_seconds': total_duration,
                    'avg_iteration_duration': total_duration / max(1, len(structure_iterations)),
                    'avg_attempt_duration': np.mean([a['duration_seconds'] for a in structure_attempts]) if structure_attempts else 0,
                    'success_rate_iterations': data.get('successful_iterations', 0) / max(1, data.get('total_iterations', 1)),
                    'success_rate_attempts': len([a for a in structure_attempts if a['success']]) / max(1, len(structure_attempts)),
                    'long_md_executed': data.get('long_md_executed', False),
                    'start_time': start_time,
                    'end_time': end_time
                }
                
                self.performance_data['structures'].append(structure_data)
                
            except Exception as e:
                print(f"구조 파일 처리 실패 {struct_file}: {e}")
        
        print(f"총 {len(self.performance_data['structures'])}개 구조 데이터 집계 완료")
    
    def create_timeline_data(self):
        """시간대별 성능 데이터 생성"""
        print("타임라인 데이터 생성 중...")
        
        # 모든 attempt를 시간순으로 정렬
        attempts_with_time = []
        for attempt in self.performance_data['attempts']:
            if attempt.get('file_time'):
                attempts_with_time.append(attempt)
        
        attempts_with_time.sort(key=lambda x: x['file_time'])
        
        # 시간별 통계 생성
        current_time = None
        time_window_minutes = 10  # 10분 단위 집계
        
        for attempt in attempts_with_time:
            attempt_time = attempt['file_time']
            
            # 10분 단위로 반올림
            rounded_time = attempt_time.replace(second=0, microsecond=0)
            rounded_time = rounded_time.replace(minute=(rounded_time.minute // time_window_minutes) * time_window_minutes)
            
            timeline_entry = {
                'time': rounded_time,
                'structure': attempt['structure'],
                'success': attempt['success'],
                'duration': attempt['duration_seconds'],
                'long_md': attempt['long_md_executed']
            }
            
            self.performance_data['timeline'].append(timeline_entry)
        
        print(f"타임라인 데이터 {len(self.performance_data['timeline'])}개 엔트리 생성 완료")
    
    def analyze_all_data(self):
        """모든 데이터 추출 및 분석"""
        print("=== Simple SuMD 성능 분석 시작 ===")
        
        if not os.path.exists(self.output_dir):
            raise FileNotFoundError(f"출력 디렉토리를 찾을 수 없습니다: {self.output_dir}")
        
        self.extract_attempt_data()
        self.extract_iteration_data()  
        self.extract_structure_data()
        self.create_timeline_data()
        
        print("=== 데이터 추출 완료 ===")
    
    def print_performance_summary(self):
        """성능 요약 정보 출력"""
        print("\n" + "="*60)
        print("성능 분석 요약")
        print("="*60)
        
        # 전체 통계
        total_attempts = len(self.performance_data['attempts'])
        successful_attempts = len([a for a in self.performance_data['attempts'] if a['success']])
        total_iterations = len(self.performance_data['iterations'])
        successful_iterations = len([i for i in self.performance_data['iterations'] if i['success']])
        total_structures = len(self.performance_data['structures'])
        
        print(f"📊 전체 통계:")
        print(f"  • 총 구조 수: {total_structures}")
        print(f"  • 총 iteration 수: {total_iterations}")
        print(f"  • 총 attempt 수: {total_attempts}")
        print(f"  • Iteration 성공률: {successful_iterations}/{total_iterations} ({successful_iterations/max(1,total_iterations)*100:.1f}%)")
        print(f"  • Attempt 성공률: {successful_attempts}/{total_attempts} ({successful_attempts/max(1,total_attempts)*100:.1f}%)")
        
        # 시간 통계
        if self.performance_data['attempts']:
            attempt_durations = [a['duration_seconds'] for a in self.performance_data['attempts']]
            print(f"\n⏱️  시간 통계:")
            print(f"  • 평균 attempt 시간: {np.mean(attempt_durations):.1f}초")
            print(f"  • 최단 attempt 시간: {np.min(attempt_durations):.1f}초")
            print(f"  • 최장 attempt 시간: {np.max(attempt_durations):.1f}초")
            print(f"  • 시간 표준편차: {np.std(attempt_durations):.1f}초")
        
        # 구조별 성능
        if self.performance_data['structures']:
            print(f"\n🏗️  구조별 성능:")
            for struct in sorted(self.performance_data['structures'], key=lambda x: x['success_rate_iterations'], reverse=True):
                print(f"  • {struct['structure_name']}: "
                      f"{struct['successful_iterations']}/{struct['total_iterations']} iteration 성공 "
                      f"({struct['success_rate_iterations']*100:.1f}%), "
                      f"평균 {struct['avg_iteration_duration']:.1f}초/iteration, "
                      f"GPU {struct['assigned_gpu']}")
        
        # 긴 MD 통계
        long_md_attempts = len([a for a in self.performance_data['attempts'] if a['long_md_executed']])
        close_contact_attempts = len([a for a in self.performance_data['attempts'] if a['close_contact_detected']])
        
        print(f"\n🔬 긴 MD 통계:")
        print(f"  • 근접 접촉 감지: {close_contact_attempts}회")
        print(f"  • 긴 MD 실행: {long_md_attempts}회")
        if long_md_attempts > 0:
            long_md_success = len([a for a in self.performance_data['attempts'] if a['long_md_executed'] and a['success']])
            print(f"  • 긴 MD 성공률: {long_md_success}/{long_md_attempts} ({long_md_success/long_md_attempts*100:.1f}%)")
    
    def plot_attempt_duration_analysis(self):
        """Attempt 소요시간 분석 그래프"""
        fig, axes = plt.subplots(2, 2, figsize=(15, 12))
        
        attempts = self.performance_data['attempts']
        if not attempts:
            return fig
        
        # 1. Attempt 소요시간 히스토그램
        durations = [a['duration_seconds'] for a in attempts]
        successful_durations = [a['duration_seconds'] for a in attempts if a['success']]
        failed_durations = [a['duration_seconds'] for a in attempts if not a['success']]
        
        axes[0,0].hist(successful_durations, bins=20, alpha=0.7, label='성공', color='green')
        axes[0,0].hist(failed_durations, bins=20, alpha=0.7, label='실패', color='red')
        axes[0,0].set_xlabel('소요시간 (초)')
        axes[0,0].set_ylabel('빈도')
        axes[0,0].set_title('Attempt 소요시간 분포')
        axes[0,0].legend()
        axes[0,0].grid(True, alpha=0.3)
        
        # 2. 구조별 평균 attempt 시간
        struct_durations = defaultdict(list)
        for attempt in attempts:
            struct_durations[attempt['structure']].append(attempt['duration_seconds'])
        
        struct_names = list(struct_durations.keys())
        struct_avg_durations = [np.mean(struct_durations[name]) for name in struct_names]
        struct_std_durations = [np.std(struct_durations[name]) for name in struct_names]
        
        x_pos = range(len(struct_names))
        axes[0,1].bar(x_pos, struct_avg_durations, yerr=struct_std_durations, 
                     capsize=5, alpha=0.7, color='skyblue')
        axes[0,1].set_xlabel('구조')
        axes[0,1].set_ylabel('평균 소요시간 (초)')
        axes[0,1].set_title('구조별 평균 Attempt 시간')
        axes[0,1].set_xticks(x_pos)
        axes[0,1].set_xticklabels(struct_names, rotation=45)
        axes[0,1].grid(True, alpha=0.3)
        
        # 3. 시간에 따른 성능 변화 (러닝 평균)
        attempts_by_time = sorted(attempts, key=lambda x: x.get('file_time', datetime.min))
        window_size = min(10, len(attempts_by_time)//4)  # 적응적 윈도우 크기
        
        if window_size > 0:
            running_avg = []
            running_success_rate = []
            time_points = []
            
            for i in range(window_size, len(attempts_by_time)):
                window_attempts = attempts_by_time[i-window_size:i]
                avg_duration = np.mean([a['duration_seconds'] for a in window_attempts])
                success_rate = np.mean([a['success'] for a in window_attempts])
                
                running_avg.append(avg_duration)
                running_success_rate.append(success_rate * 100)  # 퍼센트로 변환
                time_points.append(i)
            
            axes[1,0].plot(time_points, running_avg, 'b-', label='평균 시간')
            axes[1,0].set_xlabel('Attempt 순서')
            axes[1,0].set_ylabel('러닝 평균 시간 (초)', color='b')
            axes[1,0].tick_params(axis='y', labelcolor='b')
            
            ax2 = axes[1,0].twinx()
            ax2.plot(time_points, running_success_rate, 'r-', label='성공률')
            ax2.set_ylabel('러닝 평균 성공률 (%)', color='r')
            ax2.tick_params(axis='y', labelcolor='r')
            
            axes[1,0].set_title(f'시간 경과에 따른 성능 변화 (윈도우: {window_size})')
            axes[1,0].grid(True, alpha=0.3)
        
        # 4. 긴 MD vs 일반 MD 비교
        regular_md = [a for a in attempts if not a['long_md_executed']]
        long_md = [a for a in attempts if a['long_md_executed']]
        
        categories = []
        durations_box = []
        
        if regular_md:
            categories.append('일반 MD')
            durations_box.append([a['duration_seconds'] for a in regular_md])
        
        if long_md:
            categories.append('긴 MD')
            durations_box.append([a['duration_seconds'] for a in long_md])
        
        if durations_box:
            axes[1,1].boxplot(durations_box, labels=categories)
            axes[1,1].set_ylabel('소요시간 (초)')
            axes[1,1].set_title('MD 타입별 소요시간 비교')
            axes[1,1].grid(True, alpha=0.3)
        
        plt.tight_layout()
        return fig
    
    def plot_iteration_analysis(self):
        """Iteration 분석 그래프"""
        fig, axes = plt.subplots(2, 2, figsize=(15, 12))
        
        iterations = self.performance_data['iterations']
        if not iterations:
            return fig
        
        # 1. Iteration당 attempt 수 분포
        attempts_counts = [i['attempts_used'] for i in iterations]
        successful_attempts_counts = [i['attempts_used'] for i in iterations if i['success']]
        failed_attempts_counts = [i['attempts_used'] for i in iterations if not i['success']]
        
        bins = range(1, max(attempts_counts) + 2)
        axes[0,0].hist(successful_attempts_counts, bins=bins, alpha=0.7, label='성공', color='green')
        axes[0,0].hist(failed_attempts_counts, bins=bins, alpha=0.7, label='실패', color='red')
        axes[0,0].set_xlabel('Iteration당 Attempt 수')
        axes[0,0].set_ylabel('빈도')
        axes[0,0].set_title('Iteration당 Attempt 수 분포')
        axes[0,0].legend()
        axes[0,0].grid(True, alpha=0.3)
        
        # 2. 구조별 평균 iteration 시간
        struct_iter_durations = defaultdict(list)
        for iteration in iterations:
            struct_iter_durations[iteration['structure']].append(iteration['duration_seconds'])
        
        struct_names = list(struct_iter_durations.keys())
        struct_avg_iter_durations = [np.mean(struct_iter_durations[name]) for name in struct_names]
        
        axes[0,1].bar(struct_names, struct_avg_iter_durations, alpha=0.7, color='lightcoral')
        axes[0,1].set_xlabel('구조')
        axes[0,1].set_ylabel('평균 소요시간 (초)')
        axes[0,1].set_title('구조별 평균 Iteration 시간')
        axes[0,1].tick_params(axis='x', rotation=45)
        axes[0,1].grid(True, alpha=0.3)
        
        # 3. 성공/실패 iteration의 attempt 수 비교
        success_attempts = [i['attempts_used'] for i in iterations if i['success']]
        fail_attempts = [i['attempts_used'] for i in iterations if not i['success']]
        
        data_to_plot = []
        labels = []
        if success_attempts:
            data_to_plot.append(success_attempts)
            labels.append('성공')
        if fail_attempts:
            data_to_plot.append(fail_attempts)
            labels.append('실패')
        
        if data_to_plot:
            axes[1,0].boxplot(data_to_plot, labels=labels)
            axes[1,0].set_ylabel('Attempt 수')
            axes[1,0].set_title('성공/실패별 Attempt 수 비교')
            axes[1,0].grid(True, alpha=0.3)
        
        # 4. Iteration 번호별 성공률
        iter_numbers = defaultdict(list)
        for iteration in iterations:
            iter_numbers[iteration['iteration']].append(iteration['success'])
        
        iter_nums = sorted(iter_numbers.keys())
        success_rates = [np.mean(iter_numbers[num]) * 100 for num in iter_nums]
        
        axes[1,1].bar(iter_nums, success_rates, alpha=0.7, color='gold')
        axes[1,1].set_xlabel('Iteration 번호')
        axes[1,1].set_ylabel('성공률 (%)')
        axes[1,1].set_title('Iteration 번호별 성공률')
        axes[1,1].grid(True, alpha=0.3)
        
        plt.tight_layout()
        return fig
    
    def plot_structure_comparison(self):
        """구조별 성능 비교 그래프"""
        fig, axes = plt.subplots(2, 2, figsize=(15, 12))
        
        structures = self.performance_data['structures']
        if not structures:
            return fig
        
        # 구조를 성공률로 정렬
        structures_sorted = sorted(structures, key=lambda x: x['success_rate_iterations'], reverse=True)
        
        struct_names = [s['structure_name'] for s in structures_sorted]
        
        # 1. 구조별 성공률 비교
        iter_success_rates = [s['success_rate_iterations'] * 100 for s in structures_sorted]
        attempt_success_rates = [s['success_rate_attempts'] * 100 for s in structures_sorted]
        
        x_pos = np.arange(len(struct_names))
        width = 0.35
        
        axes[0,0].bar(x_pos - width/2, iter_success_rates, width, label='Iteration 성공률', alpha=0.7)
        axes[0,0].bar(x_pos + width/2, attempt_success_rates, width, label='Attempt 성공률', alpha=0.7)
        axes[0,0].set_xlabel('구조')
        axes[0,0].set_ylabel('성공률 (%)')
        axes[0,0].set_title('구조별 성공률 비교')
        axes[0,0].set_xticks(x_pos)
        axes[0,0].set_xticklabels(struct_names, rotation=45)
        axes[0,0].legend()
        axes[0,0].grid(True, alpha=0.3)
        
        # 2. 구조별 총 소요시간
        total_durations = [s['total_duration_seconds'] for s in structures_sorted]
        colors = plt.cm.viridis(np.linspace(0, 1, len(struct_names)))
        
        axes[0,1].bar(struct_names, total_durations, color=colors, alpha=0.7)
        axes[0,1].set_xlabel('구조')
        axes[0,1].set_ylabel('총 소요시간 (초)')
        axes[0,1].set_title('구조별 총 소요시간')
        axes[0,1].tick_params(axis='x', rotation=45)
        axes[0,1].grid(True, alpha=0.3)
        
        # 3. 효율성 분석 (성공률 vs 시간)
        axes[1,0].scatter(total_durations, iter_success_rates, alpha=0.7, s=100)
        for i, name in enumerate(struct_names):
            axes[1,0].annotate(name, (total_durations[i], iter_success_rates[i]), 
                             xytext=(5, 5), textcoords='offset points', fontsize=8)
        axes[1,0].set_xlabel('총 소요시간 (초)')
        axes[1,0].set_ylabel('Iteration 성공률 (%)')
        axes[1,0].set_title('효율성 분석: 시간 vs 성공률')
        axes[1,0].grid(True, alpha=0.3)
        
        # 4. GPU별 성능 비교
        gpu_performance = defaultdict(list)
        for struct in structures:
            gpu_id = struct.get('assigned_gpu', 'unknown')
            gpu_performance[gpu_id].append(struct['success_rate_iterations'])
        
        gpu_ids = list(gpu_performance.keys())
        gpu_success_rates = [np.mean(gpu_performance[gpu_id]) * 100 for gpu_id in gpu_ids]
        
        if len(gpu_ids) > 1:
            axes[1,1].bar(gpu_ids, gpu_success_rates, alpha=0.7, color='lightgreen')
            axes[1,1].set_xlabel('GPU ID')
            axes[1,1].set_ylabel('평균 성공률 (%)')
            axes[1,1].set_title('GPU별 평균 성능')
            axes[1,1].grid(True, alpha=0.3)
        else:
            axes[1,1].text(0.5, 0.5, 'GPU 1개만 사용됨', 
                          horizontalalignment='center', verticalalignment='center',
                          transform=axes[1,1].transAxes, fontsize=12)
            axes[1,1].set_title('GPU 사용 현황')
        
        plt.tight_layout()
        return fig
    
    def plot_timeline_analysis(self):
        """시간대별 성능 분석 그래프"""
        fig, axes = plt.subplots(2, 1, figsize=(15, 10))
        
        timeline = self.performance_data['timeline']
        if not timeline:
            return fig
        
        # 시간별 데이터 집계
        time_stats = defaultdict(lambda: {'total': 0, 'success': 0, 'durations': []})
        
        for entry in timeline:
            time_key = entry['time']
            time_stats[time_key]['total'] += 1
            if entry['success']:
                time_stats[time_key]['success'] += 1
            time_stats[time_key]['durations'].append(entry['duration'])
        
        if not time_stats:
            return fig
        
        times = sorted(time_stats.keys())
        success_rates = [time_stats[t]['success'] / time_stats[t]['total'] * 100 for t in times]
        avg_durations = [np.mean(time_stats[t]['durations']) for t in times]
        attempt_counts = [time_stats[t]['total'] for t in times]
        
        # 1. 시간대별 성공률과 attempt 수
        axes[0].plot(times, success_rates, 'g-o', label='성공률 (%)', markersize=4)
        axes[0].set_ylabel('성공률 (%)', color='g')
        axes[0].tick_params(axis='y', labelcolor='g')
        
        ax2 = axes[0].twinx()
        ax2.bar(times, attempt_counts, alpha=0.3, color='blue', label='Attempt 수', width=0.003)
        ax2.set_ylabel('Attempt 수', color='b')
        ax2.tick_params(axis='y', labelcolor='b')
        
        axes[0].set_title('시간대별 성공률 및 Attempt 수')
        axes[0].grid(True, alpha=0.3)
        
        # X축 포맷팅
        axes[0].xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
        axes[0].xaxis.set_major_locator(mdates.MinuteLocator(interval=30))
        fig.autofmt_xdate()
        
        # 2. 시간대별 평균 소요시간
        axes[1].plot(times, avg_durations, 'r-o', label='평균 소요시간', markersize=4)
        axes[1].set_xlabel('시간')
        axes[1].set_ylabel('평균 소요시간 (초)')
        axes[1].set_title('시간대별 평균 소요시간')
        axes[1].grid(True, alpha=0.3)
        
        # X축 포맷팅
        axes[1].xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
        axes[1].xaxis.set_major_locator(mdates.MinuteLocator(interval=30))
        
        plt.tight_layout()
        return fig
    
    def save_performance_report(self):
        """성능 분석 보고서 저장"""
        report_path = os.path.join(self.output_dir, "performance_report.txt")
        
        with open(report_path, 'w', encoding='utf-8') as f:
            f.write("Simple SuMD 성능 분석 보고서\n")
            f.write("=" * 50 + "\n\n")
            
            # 기본 통계
            total_attempts = len(self.performance_data['attempts'])
            successful_attempts = len([a for a in self.performance_data['attempts'] if a['success']])
            total_iterations = len(self.performance_data['iterations'])
            successful_iterations = len([i for i in self.performance_data['iterations'] if i['success']])
            
            f.write("전체 통계:\n")
            f.write(f"  총 구조 수: {len(self.performance_data['structures'])}\n")
            f.write(f"  총 iteration 수: {total_iterations}\n")
            f.write(f"  총 attempt 수: {total_attempts}\n")
            f.write(f"  Iteration 성공률: {successful_iterations}/{total_iterations} ({successful_iterations/max(1,total_iterations)*100:.1f}%)\n")
            f.write(f"  Attempt 성공률: {successful_attempts}/{total_attempts} ({successful_attempts/max(1,total_attempts)*100:.1f}%)\n\n")
            
            # 시간 통계
            if self.performance_data['attempts']:
                durations = [a['duration_seconds'] for a in self.performance_data['attempts']]
                f.write("시간 통계:\n")
                f.write(f"  평균 attempt 시간: {np.mean(durations):.1f}초\n")
                f.write(f"  최단/최장 attempt 시간: {np.min(durations):.1f}초 / {np.max(durations):.1f}초\n")
                f.write(f"  시간 표준편차: {np.std(durations):.1f}초\n\n")
            
            # 구조별 상세 정보
            f.write("구조별 상세 성능:\n")
            for struct in sorted(self.performance_data['structures'], 
                               key=lambda x: x['success_rate_iterations'], reverse=True):
                f.write(f"  {struct['structure_name']}:\n")
                f.write(f"    - Iteration 성공률: {struct['success_rate_iterations']*100:.1f}%\n")
                f.write(f"    - Attempt 성공률: {struct['success_rate_attempts']*100:.1f}%\n")
                f.write(f"    - 총 소요시간: {struct['total_duration_seconds']:.1f}초\n")
                f.write(f"    - 평균 iteration 시간: {struct['avg_iteration_duration']:.1f}초\n")
                f.write(f"    - 사용 GPU: {struct['assigned_gpu']}\n")
                f.write(f"    - 긴 MD 실행: {'예' if struct['long_md_executed'] else '아니오'}\n\n")
            
            # 단계별 시간 분석
            stage_times = defaultdict(list)
            for attempt in self.performance_data['attempts']:
                for stage, duration in attempt.get('stage_durations', {}).items():
                    stage_times[stage].append(duration)
            
            if stage_times:
                f.write("단계별 평균 소요시간:\n")
                for stage, times in stage_times.items():
                    f.write(f"  {stage}: {np.mean(times):.1f}초 (σ={np.std(times):.1f})\n")
        
        print(f"성능 분석 보고서 저장: {report_path}")
    
    def generate_all_plots(self):
        """모든 성능 분석 그래프 생성 및 저장"""
        print("\n성능 분석 그래프 생성 중...")
        
        # 한국어 폰트 설정 (Linux 환경 고려)
        plt.rcParams['font.family'] = 'DejaVu Sans'
        plt.rcParams['axes.unicode_minus'] = False
        
        plots = [
            ("attempt_analysis", "Attempt 소요시간 분석", self.plot_attempt_duration_analysis),
            ("iteration_analysis", "Iteration 분석", self.plot_iteration_analysis),
            ("structure_comparison", "구조별 성능 비교", self.plot_structure_comparison),
            ("timeline_analysis", "시간대별 성능 분석", self.plot_timeline_analysis)
        ]
        
        for filename, title, plot_func in plots:
            try:
                fig = plot_func()
                plot_path = os.path.join(self.output_dir, f"performance_{filename}.png")
                fig.savefig(plot_path, dpi=300, bbox_inches='tight')
                plt.close(fig)
                print(f"  ✓ {title}: {plot_path}")
            except Exception as e:
                print(f"  ✗ {title} 생성 실패: {e}")
        
        print("성능 분석 그래프 생성 완료!")

def main():
    """메인 함수 - 독립 실행용"""
    if len(sys.argv) != 2:
        print("사용법: python3 performance_analyzer.py <output_directory>")
        print("예시: python3 performance_analyzer.py sumd_output")
        sys.exit(1)
    
    output_dir = sys.argv[1]
    
    try:
        # 성능 분석기 초기화 및 실행
        analyzer = PerformanceAnalyzer(output_dir)
        
        # 모든 데이터 분석
        analyzer.analyze_all_data()
        
        # 결과 출력
        analyzer.print_performance_summary()
        
        # 그래프 생성
        analyzer.generate_all_plots()
        
        # 보고서 저장
        analyzer.save_performance_report()
        
        print(f"\n🎉 성능 분석 완료! 결과는 {output_dir} 디렉토리에 저장되었습니다.")
        print(f"   📊 그래프: performance_*.png")
        print(f"   📄 보고서: performance_report.txt")
        
    except Exception as e:
        print(f"❌ 성능 분석 중 오류 발생: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()