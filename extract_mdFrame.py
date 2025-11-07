#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import shutil
import glob
import json
import subprocess
import argparse
from datetime import datetime

def run_gromacs_command(cmd, cwd=None):
    """GROMACS 명령어 실행"""
    try:
        result = subprocess.run(
            cmd, shell=True, cwd=cwd,
            capture_output=True, text=True, timeout=3600
        )
        if result.returncode != 0:
            print(f"GROMACS 명령어 실행 실패: {cmd}")
            print(f"Error: {result.stderr}")
            return False
        return True
    except Exception as e:
        print(f"GROMACS 명령어 실행 중 오류: {e}")
        return False

def extract_trajectory_frames(tpr_file, xtc_file, output_dir, prefix="frame"):
    """궤적에서 모든 프레임을 개별 PDB로 추출"""
    print(f"  궤적 추출 중: {os.path.basename(xtc_file)}")
    
    # trjconv로 protein만 선택하여 각 프레임을 개별 PDB로 추출
    cmd = f"echo 'Protein' | gmx trjconv -s {tpr_file} -f {xtc_file} -o {prefix}.pdb -sep"
    
    if not run_gromacs_command(cmd, cwd=output_dir):
        print(f"  ⚠ 궤적 추출 실패: {os.path.basename(xtc_file)}")
        return []
    
    # 생성된 프레임 파일들 찾기
    frame_files = []
    frame_num = 0
    while True:
        frame_file = os.path.join(output_dir, f"{prefix}{frame_num}.pdb")
        if os.path.exists(frame_file) and os.path.getsize(frame_file) > 0:
            frame_files.append(frame_file)
            frame_num += 1
        else:
            break
    
    print(f"  ✓ 추출된 프레임 수: {len(frame_files)}")
    return frame_files

def find_successful_attempt(iteration_dir, iteration_summary):
    """성공한 attempt 디렉토리 찾기"""
    if not iteration_summary.get("success", False):
        return None
    
    attempts_used = iteration_summary.get("attempts_used", 1)
    
    # 마지막 성공한 attempt 디렉토리 확인
    for attempt_num in range(attempts_used, 0, -1):
        attempt_dir = os.path.join(iteration_dir, f"attempt_{attempt_num}")
        if os.path.exists(attempt_dir):
            # md.tpr과 md.xtc 파일이 있는지 확인
            tpr_file = os.path.join(attempt_dir, "md.tpr")
            xtc_file = os.path.join(attempt_dir, "md.xtc")
            
            if os.path.exists(tpr_file) and os.path.exists(xtc_file) and os.path.getsize(xtc_file) > 0:
                return attempt_dir, attempt_num
    
    return None

def get_trajectory_time_info(tpr_file, xtc_file, work_dir):
    """XTC 파일의 시간 정보 추출"""
    try:
        cmd = f"gmx check -f {xtc_file} 2>&1"
        result = subprocess.run(
            cmd, shell=True, cwd=work_dir,
            capture_output=True, text=True, timeout=30
        )
        
        # Reading frame 출력에서 시작/끝 시간 추출
        output = result.stdout + result.stderr
        lines = output.split('\n')
        
        start_time = 0.0
        end_time = 0.0
        
        for line in lines:
            if 'time' in line.lower() and 'last frame' in line.lower():
                # "time 300.000 ps" 같은 형식에서 시간 추출
                import re
                times = re.findall(r'(\d+\.?\d*)\s*', line)
                if times:
                    end_time = float(times[-1])
        
        # 또는 gmx dump로 정확한 정보 확인
        if end_time == 0.0:
            cmd = f"gmx dump -f {xtc_file} 2>&1 | grep 'time' | tail -1"
            result = subprocess.run(
                cmd, shell=True, cwd=work_dir,
                capture_output=True, text=True, timeout=30
            )
            import re
            times = re.findall(r'(\d+\.?\d*)', result.stdout)
            if times:
                end_time = float(times[-1])
        
        return start_time, end_time
        
    except Exception as e:
        print(f"  ⚠ 시간 정보 추출 실패: {e}")
        # 기본값: 설정 파일의 시뮬레이션 시간 사용
        try:
            from simple_config import SIMULATION_TIME_NS, LONG_MD_TIME_NS
            return 0.0, SIMULATION_TIME_NS * 1000  # ns -> ps
        except:
            return 0.0, 300.0  # 기본 300ps

def adjust_trajectory_time(input_xtc, output_xtc, start_time_offset, tpr_file, work_dir):
    """XTC 파일의 시간을 오프셋만큼 조정"""
    try:
        # 폴백: Protein만 (리간드가 Protein 그룹에 포함된 경우도 있음)
        cmd = f"echo 'Protein' | gmx trjconv -s {tpr_file} -f {input_xtc} -o {output_xtc} -t0 {start_time_offset}"
    
        result = subprocess.run(
            cmd, shell=True, cwd=work_dir,
            capture_output=True, text=True, timeout=3600
        )
        
        if result.returncode != 0:
            print(f"  ⚠ 시간 조정 실패: {result.stderr}")
            return False
            
        return True
        
    except Exception as e:
        print(f"  ⚠ 시간 조정 중 오류: {e}")
        return False
    
def concatenate_trajectories_to_xtc(trajectory_info_list, output_xtc, temp_dir):
    """
    여러 XTC 파일들을 시간 보정하여 하나의 XTC로 연결
    
    Args:
        trajectory_info_list: [(xtc_file, tpr_file, is_long_md, iter_num), ...] 형식의 리스트
        output_xtc: 출력 XTC 파일 경로
        temp_dir: 임시 작업 디렉토리
    """
    if not trajectory_info_list:
        return False, {}
    
    print(f"  XTC 파일 시간 보정 및 연결 중: {len(trajectory_info_list)}개 궤적")
    
    adjusted_files = []
    time_info = {
        "iterations": [],
        "total_time_ps": 0.0
    }
    
    cumulative_time = 0.0  # 누적 시간 (ps)
    
    # 1단계: 각 XTC의 시간 보정
    for idx, (src_xtc, src_tpr, is_long_md, iter_num) in enumerate(trajectory_info_list):
        print(f"  처리 중: Iteration {iter_num} (시작 시간: {cumulative_time:.1f} ps)")
        
        # 원본 XTC의 시간 정보 추출
        start_time, end_time = get_trajectory_time_info(src_tpr, src_xtc, temp_dir)
        duration = end_time - start_time
        
        # 시간 보정된 XTC 파일 생성
        adjusted_xtc = os.path.join(temp_dir, f"adjusted_{idx:03d}.xtc")
        
        # trjconv로 시작 시간 조정
        cmd = f"echo 'System' | gmx trjconv -s {src_tpr} -f {src_xtc} -o {adjusted_xtc} -t0 {cumulative_time}"
        
        result = subprocess.run(
            cmd, shell=True, cwd=temp_dir,
            capture_output=True, text=True, timeout=3600
        )
        
        if result.returncode != 0:
            print(f"  ⚠ 시간 조정 실패: {result.stderr}")
            return False, {}
        
        adjusted_files.append(adjusted_xtc)
        
        # 시간 정보 기록
        iter_time_info = {
            "iteration": iter_num,
            "start_time_ps": cumulative_time,
            "end_time_ps": cumulative_time + duration,
            "duration_ps": duration,
            "is_long_md": is_long_md
        }
        time_info["iterations"].append(iter_time_info)
        
        cumulative_time += duration
        print(f"  ✓ 조정 완료: {cumulative_time:.1f} ps까지")
    
    time_info["total_time_ps"] = cumulative_time
    time_info["total_time_ns"] = cumulative_time / 1000.0
    
    # 2단계: 모든 조정된 XTC를 하나로 연결
    try:
        print(f"\n  최종 XTC 파일 생성 중...")
        
        # trjcat으로 연결
        input_files = " ".join(adjusted_files)
        cmd = f"gmx trjcat -f {input_files} -o {output_xtc} -cat"
        
        result = subprocess.run(
            cmd, shell=True, cwd=temp_dir,
            capture_output=True, text=True, timeout=7200
        )
        
        if result.returncode != 0:
            print(f"  ⚠ XTC 연결 실패: {result.stderr}")
            return False, {}
        
        print(f"  ✓ 최종 XTC 생성 완료")
        return True, time_info
        
    except Exception as e:
        print(f"  ⚠ XTC 연결 중 오류: {e}")
        return False, {}

def find_all_structure_directories(base_dir):
    """
    주어진 경로에서 structure_{PDBID}_struct_{번호} 디렉토리들을 재귀적으로 탐색
    structure_pool 같은 다른 structure_ 디렉토리는 제외
    
    Returns:
        structure_dirs: [(structure_path, pdb_info), ...] 형식의 리스트
    """
    import re
    
    structure_dirs = []
    
    # structure_{PDBID}_struct_{번호} 패턴
    # 예: structure_1ABI_struct_000, structure_7KI0_struct_001
    structure_pattern = re.compile(r'^structure_([A-Z0-9]{4})_struct_(\d+)$')
    
    # 직접 iteration이 있는 경우 (structure_{PDBID}_struct_{번호} 레벨)
    iteration_dirs = [d for d in os.listdir(base_dir) 
                     if d.startswith("iteration_") and os.path.isdir(os.path.join(base_dir, d))]
    
    if iteration_dirs:
        # 현재 디렉토리가 structure 디렉토리인지 확인
        dir_name = os.path.basename(base_dir)
        if structure_pattern.match(dir_name):
            parent_name = os.path.basename(os.path.dirname(base_dir))
            structure_dirs.append((base_dir, {
                'structure_name': dir_name,
                'pdb_name': parent_name
            }))
            return structure_dirs
    
    # 하위 디렉토리 탐색
    try:
        for item in os.listdir(base_dir):
            item_path = os.path.join(base_dir, item)
            if os.path.isdir(item_path):
                # structure_{PDBID}_struct_{번호} 패턴에 맞는 디렉토리만
                if structure_pattern.match(item):
                    sub_iterations = [d for d in os.listdir(item_path)
                                     if d.startswith("iteration_") and 
                                     os.path.isdir(os.path.join(item_path, d))]
                    if sub_iterations:
                        parent_name = os.path.basename(base_dir)
                        structure_dirs.append((item_path, {
                            'structure_name': item,
                            'pdb_name': parent_name
                        }))
                # structure_로 시작하지 않는 디렉토리는 재귀 탐색 (1ABI_H_A 같은 폴더)
                elif not item.startswith("structure_"):
                    structure_dirs.extend(find_all_structure_directories(item_path))
    except PermissionError:
        pass
    
    return structure_dirs
    
def collect_trajectories_as_xtc(sumd_output_dir, output_dir="collected_trajectories"):
    """
    SuMD 시뮬레이션의 모든 성공한 attempt 궤적을 시간 보정하여 XTC 파일로 수집
    
    Args:
        sumd_output_dir: SuMD 출력 디렉토리 경로 (job_output 또는 structure_XXXX 레벨 모두 가능)
        output_dir: 수집된 궤적을 저장할 디렉토리 이름
    """
    print(f"SuMD 궤적 XTC 수집 시작 (시간 보정): {sumd_output_dir}")
    
    # 모든 structure 디렉토리 찾기
    structure_dirs = find_all_structure_directories(sumd_output_dir)
    
    if not structure_dirs:
        print("❌ 처리할 structure 디렉토리를 찾을 수 없습니다.")
        return None, 0
    
    print(f"발견된 structure 디렉토리 수: {len(structure_dirs)}")
    
    total_processed = 0
    
    # 각 structure 디렉토리 처리
    for structure_path, structure_info in structure_dirs:
        print(f"\n{'='*60}")
        print(f"처리 중: {structure_info['pdb_name']}/{structure_info['structure_name']}")
        print(f"{'='*60}")
        
        # 출력 디렉토리 생성
        full_output_dir = os.path.join(structure_path, output_dir)
        if os.path.exists(full_output_dir):
            shutil.rmtree(full_output_dir)
        os.makedirs(full_output_dir)
        
        # 임시 작업 디렉토리
        temp_dir = os.path.join(full_output_dir, "temp_xtc")
        os.makedirs(temp_dir)
        
        collected_info = {
            "collection_time": datetime.now().isoformat(),
            "source_directory": structure_path,
            "pdb_name": structure_info['pdb_name'],
            "structure_name": structure_info['structure_name'],
            "format": "xtc",
            "time_corrected": True,
            "iterations": []
        }
        
        trajectory_info_list = []  # [(xtc, tpr, is_long_md, iter_num), ...]
        
        # Iteration별 성공한 attempt의 궤적 수집
        iteration_dirs = []
        for item in os.listdir(structure_path):
            if item.startswith("iteration_") and os.path.isdir(os.path.join(structure_path, item)):
                try:
                    iter_num_str = item.replace("iteration_", "")
                    if "_" in iter_num_str:
                        iter_num = int(iter_num_str.split("_")[0])
                    else:
                        iter_num = int(iter_num_str)
                    iteration_dirs.append((iter_num, item))
                except ValueError:
                    continue
        
        # 번호 순으로 정렬
        iteration_dirs.sort(key=lambda x: x[0])
        
        print(f"발견된 iteration 디렉토리 수: {len(iteration_dirs)}")
        
        for iter_num, iter_dir in iteration_dirs:
            print(f"\n=== Iteration {iter_num} 정보 수집 중 ===")
            
            iter_path = os.path.join(structure_path, iter_dir)
            summary_file = os.path.join(structure_path, f"{iter_dir}_summary.json")
            
            # summary 파일에서 성공 정보 읽기
            iteration_summary = {}
            if os.path.exists(summary_file):
                try:
                    with open(summary_file, 'r') as f:
                        iteration_summary = json.load(f)
                except:
                    print(f"  ⚠ {iter_dir}_summary.json 읽기 실패")
                    continue
            else:
                print(f"  ⚠ {iter_dir}_summary.json 파일이 없습니다")
                continue
            
            # 성공한 attempt 찾기
            attempt_info = find_successful_attempt(iter_path, iteration_summary)
            if not attempt_info:
                print(f"  ⚠ Iteration {iter_num}: 성공한 attempt를 찾을 수 없습니다")
                continue
            
            attempt_dir, attempt_num = attempt_info
            print(f"  ✓ 성공한 attempt: {attempt_num}")
            
            # XTC/TPR 파일 확인
            src_xtc = os.path.join(attempt_dir, "md.xtc")
            src_tpr = os.path.join(attempt_dir, "md.tpr")
            
            if not os.path.exists(src_xtc) or not os.path.exists(src_tpr):
                print(f"  ⚠ Iteration {iter_num}: XTC 또는 TPR 파일이 없습니다")
                continue
            
            is_long_md = iteration_summary.get("long_md", False)
            
            # 궤적 정보 리스트에 추가
            trajectory_info_list.append((src_xtc, src_tpr, is_long_md, iter_num))
            
            # iteration 기본 정보 기록
            iteration_info = {
                "type": "md_trajectory",
                "iteration_number": iter_num,
                "attempt_number": attempt_num,
                "source": f"{iter_dir}/attempt_{attempt_num}/md.xtc",
                "description": f"Iteration {iter_num} MD 궤적 (Attempt {attempt_num})",
                "is_long_md": is_long_md
            }
            
            # summary에서 추가 정보
            if iteration_summary.get("final_result"):
                result = iteration_summary["final_result"]
                iteration_info.update({
                    "slope": result.get("slope"),
                    "min_distance": result.get("min_distance"),
                    "initial_distance": result.get("initial_distance"),
                    "final_distance": result.get("final_distance")
                })
            
            collected_info["iterations"].append(iteration_info)
            
            long_md_mark = "(Long MD)" if is_long_md else ""
            print(f"  ✓ 궤적 추가 대기: iter_{iter_num:03d} {long_md_mark}")
        
        # 모든 궤적을 시간 보정하여 하나의 XTC로 연결
        if trajectory_info_list:
            final_xtc = os.path.join(full_output_dir, "complete_trajectory.xtc")
            print(f"\n{'='*60}")
            print(f"시간 보정 및 연속 궤적 생성 시작")
            print(f"{'='*60}")
            
            success, time_info = concatenate_trajectories_to_xtc(
                trajectory_info_list, final_xtc, temp_dir
            )
            
            if success:
                print(f"\n✓ 연속 궤적 생성 완료: complete_trajectory.xtc")
                print(f"  총 시뮬레이션 시간: {time_info['total_time_ns']:.2f} ns ({time_info['total_time_ps']:.1f} ps)")
                
                collected_info["trajectory_file"] = "complete_trajectory.xtc"
                collected_info["total_iterations"] = len(trajectory_info_list)
                collected_info["time_info"] = time_info
                
                # 각 iteration 정보에 시간 정보 추가
                for i, iter_info in enumerate(collected_info["iterations"]):
                    if i < len(time_info["iterations"]):
                        iter_info.update(time_info["iterations"][i])
                
                # 첫 번째 TPR 파일도 복사 (시각화용)
                first_tpr = trajectory_info_list[0][1]
                ref_tpr = os.path.join(full_output_dir, "reference.tpr")
                shutil.copy(first_tpr, ref_tpr)
                collected_info["reference_tpr"] = "reference.tpr"
                print(f"✓ 참조 TPR 저장: reference.tpr")
                
                total_processed += 1
            else:
                print(f"⚠ XTC 연결 실패")
        
        # 임시 디렉토리 정리
        if os.path.exists(temp_dir):
            shutil.rmtree(temp_dir)
        
        # 수집 정보 저장
        info_file = os.path.join(full_output_dir, "collection_info.json")
        with open(info_file, 'w') as f:
            json.dump(collected_info, f, indent=2, ensure_ascii=False)
        
        # 요약 텍스트 파일 생성
        summary_file = os.path.join(full_output_dir, "trajectory_summary.txt")
        with open(summary_file, 'w', encoding='utf-8') as f:
            f.write(f"SuMD 궤적 요약 (XTC 형식, 시간 보정 적용)\n")
            f.write(f"PDB: {collected_info['pdb_name']}\n")
            f.write(f"Structure: {collected_info['structure_name']}\n")
            f.write(f"수집 시간: {collected_info['collection_time']}\n")
            f.write(f"처리된 iteration 수: {collected_info.get('total_iterations', 0)}\n")
            f.write(f"총 시뮬레이션 시간: {collected_info.get('time_info', {}).get('total_time_ns', 0):.2f} ns\n")
            f.write(f"출력 파일: complete_trajectory.xtc\n")
            f.write("="*60 + "\n\n")
            
            f.write("시간별 Iteration 정보:\n")
            f.write("-"*60 + "\n")
            for info in collected_info["iterations"]:
                f.write(f"\nIteration {info['iteration_number']}: {info['description']}\n")
                if 'start_time_ps' in info:
                    f.write(f"  - 시간 범위: {info['start_time_ps']:.1f} - {info['end_time_ps']:.1f} ps\n")
                    f.write(f"  - 지속 시간: {info['duration_ps']:.1f} ps ({info['duration_ps']/1000:.2f} ns)\n")
                if 'min_distance' in info:
                    f.write(f"  - 최소거리: {info['min_distance']:.2f}Å\n")
                if 'slope' in info:
                    f.write(f"  - 기울기: {info['slope']:.6f}\n")
                if info.get("is_long_md"):
                    f.write(f"  - Long MD 궤적\n")
        
        print(f"\n{'='*60}")
        print(f"Structure 처리 완료: {structure_info['structure_name']}")
        print(f"{'='*60}")
        print(f"처리된 iteration 수: {collected_info.get('total_iterations', 0)}")
        print(f"총 시뮬레이션 시간: {collected_info.get('time_info', {}).get('total_time_ns', 0):.2f} ns")
        print(f"저장 위치: {full_output_dir}")
    
    print(f"\n{'='*60}")
    print(f"전체 수집 완료")
    print(f"{'='*60}")
    print(f"처리된 structure 수: {total_processed}/{len(structure_dirs)}")
    
    return sumd_output_dir, total_processed

def collect_sumd_trajectories(sumd_output_dir, output_dir="collected_trajectories"):
    """
    SuMD 시뮬레이션의 모든 성공한 attempt 궤적을 프레임별로 수집
    
    Args:
        sumd_output_dir: SuMD 출력 디렉토리 경로 (job_output 또는 structure_XXXX 레벨 모두 가능)
        output_dir: 수집된 구조들을 저장할 디렉토리 이름
    """
    print(f"SuMD 궤적 프레임 수집 시작: {sumd_output_dir}")
    
    # 모든 structure 디렉토리 찾기
    structure_dirs = find_all_structure_directories(sumd_output_dir)
    
    if not structure_dirs:
        print("❌ 처리할 structure 디렉토리를 찾을 수 없습니다.")
        return None, 0
    
    print(f"발견된 structure 디렉토리 수: {len(structure_dirs)}")
    
    total_frames_all = 0
    
    # 각 structure 디렉토리 처리
    for structure_path, structure_info in structure_dirs:
        print(f"\n{'='*60}")
        print(f"처리 중: {structure_info['pdb_name']}/{structure_info['structure_name']}")
        print(f"{'='*60}")
        
        # 출력 디렉토리 생성
        full_output_dir = os.path.join(structure_path, output_dir)
        if os.path.exists(full_output_dir):
            shutil.rmtree(full_output_dir)
        os.makedirs(full_output_dir)
        
        # 임시 작업 디렉토리
        temp_dir = os.path.join(full_output_dir, "temp_extraction")
        os.makedirs(temp_dir)
        
        structure_count = 0
        collected_info = {
            "collection_time": datetime.now().isoformat(),
            "source_directory": structure_path,
            "pdb_name": structure_info['pdb_name'],
            "structure_name": structure_info['structure_name'],
            "total_frames": 0,
            "iterations": []
        }
        
        # Iteration별 성공한 attempt의 궤적 수집
        iteration_dirs = []
        for item in os.listdir(structure_path):
            if item.startswith("iteration_") and os.path.isdir(os.path.join(structure_path, item)):
                try:
                    iter_num_str = item.replace("iteration_", "")
                    if "_" in iter_num_str:
                        iter_num = int(iter_num_str.split("_")[0])
                    else:
                        iter_num = int(iter_num_str)
                    iteration_dirs.append((iter_num, item))
                except ValueError:
                    continue
        
        # 번호 순으로 정렬
        iteration_dirs.sort(key=lambda x: x[0])
        
        print(f"발견된 iteration 디렉토리 수: {len(iteration_dirs)}")
        
        for iter_num, iter_dir in iteration_dirs:
            print(f"\n=== Iteration {iter_num} 처리 중 ===")
            
            iter_path = os.path.join(structure_path, iter_dir)
            summary_file = os.path.join(structure_path, f"{iter_dir}_summary.json")
            
            # summary 파일에서 성공 정보 읽기
            iteration_summary = {}
            if os.path.exists(summary_file):
                try:
                    with open(summary_file, 'r') as f:
                        iteration_summary = json.load(f)
                except:
                    print(f"  ⚠ {iter_dir}_summary.json 읽기 실패")
                    continue
            else:
                print(f"  ⚠ {iter_dir}_summary.json 파일이 없습니다")
                continue
            
            # 성공한 attempt 찾기
            attempt_info = find_successful_attempt(iter_path, iteration_summary)
            if not attempt_info:
                print(f"  ⚠ Iteration {iter_num}: 성공한 attempt를 찾을 수 없습니다")
                continue
            
            attempt_dir, attempt_num = attempt_info
            print(f"  ✓ 성공한 attempt: {attempt_num}")
            
            # 궤적 추출
            tpr_file = os.path.join(attempt_dir, "md.tpr")
            xtc_file = os.path.join(attempt_dir, "md.xtc")
            
            # 임시 디렉토리에서 프레임 추출
            temp_prefix = f"iter_{iter_num:03d}"
            frame_files = extract_trajectory_frames(tpr_file, xtc_file, temp_dir, temp_prefix)
            
            if not frame_files:
                print(f"  ⚠ Iteration {iter_num}: 프레임 추출 실패")
                continue
            
            # 프레임들을 최종 위치로 복사하고 번호 매기기
            frame_start = structure_count
            for i, frame_file in enumerate(frame_files):
                dest_path = os.path.join(full_output_dir, f"frame_{structure_count:06d}.pdb")
                shutil.copy(frame_file, dest_path)
                os.remove(frame_file)  # 임시 파일 삭제
                structure_count += 1
            
            # iteration 정보 기록
            iteration_info = {
                "type": "md_trajectory",
                "iteration_number": iter_num,
                "attempt_number": attempt_num,
                "source": f"{iter_dir}/attempt_{attempt_num}/md.xtc",
                "frame_start": frame_start,
                "frame_count": len(frame_files),
                "description": f"Iteration {iter_num} MD 궤적 (Attempt {attempt_num})",
                "is_long_md": iteration_summary.get("long_md", False)
            }
            
            # summary에서 추가 정보
            if iteration_summary.get("final_result"):
                result = iteration_summary["final_result"]
                iteration_info.update({
                    "slope": result.get("slope"),
                    "min_distance": result.get("min_distance"),
                    "initial_distance": result.get("initial_distance"),
                    "final_distance": result.get("final_distance")
                })
            
            collected_info["iterations"].append(iteration_info)
            
            long_md_mark = "(Long MD)" if iteration_info.get("is_long_md") else ""
            print(f"  ✓ Frames {frame_start:06d}-{structure_count-1:06d}: {len(frame_files)}개 프레임 수집 {long_md_mark}")
        
        # 임시 디렉토리 정리
        if os.path.exists(temp_dir):
            shutil.rmtree(temp_dir)
        
        # 수집 정보 업데이트 및 저장
        collected_info["total_frames"] = structure_count
        collected_info["total_iterations"] = len([info for info in collected_info["iterations"] 
                                                   if info["type"] == "md_trajectory"])
        
        info_file = os.path.join(full_output_dir, "collection_info.json")
        with open(info_file, 'w') as f:
            json.dump(collected_info, f, indent=2, ensure_ascii=False)
        
        # 요약 텍스트 파일 생성
        summary_file = os.path.join(full_output_dir, "frames_summary.txt")
        with open(summary_file, 'w', encoding='utf-8') as f:
            f.write(f"SuMD 궤적 프레임 요약\n")
            f.write(f"PDB: {collected_info['pdb_name']}\n")
            f.write(f"Structure: {collected_info['structure_name']}\n")
            f.write(f"수집 시간: {collected_info['collection_time']}\n")
            f.write(f"총 프레임 수: {structure_count}\n")
            f.write(f"처리된 iteration 수: {collected_info['total_iterations']}\n")
            f.write("="*50 + "\n\n")
            
            for info in collected_info["iterations"]:
                f.write(f"Frames {info['frame_start']:06d}-{info['frame_start']+info['frame_count']-1:06d}: {info['description']}\n")
                if 'min_distance' in info:
                    f.write(f"  - 최소거리: {info['min_distance']:.2f}Å\n")
                if 'slope' in info:
                    f.write(f"  - 기울기: {info['slope']:.6f}\n")
                if info.get("is_long_md"):
                    f.write(f"  - Long MD 궤적\n")
                f.write("\n")
        
        print(f"\n=== Structure 수집 완료 ===")
        print(f"수집된 프레임 수: {structure_count}")
        print(f"처리된 iteration 수: {collected_info['total_iterations']}")
        print(f"저장 위치: {full_output_dir}")
        
        total_frames_all += structure_count
    
    print(f"\n{'='*60}")
    print(f"전체 수집 완료")
    print(f"{'='*60}")
    print(f"처리된 structure 수: {len(structure_dirs)}")
    print(f"총 프레임 수: {total_frames_all}")
    
    return sumd_output_dir, total_frames_all

def main():
    """메인 함수 - argparse로 명령행 인자 처리"""
    parser = argparse.ArgumentParser(
        description="SuMD 궤적 프레임 수집 - 성공한 각 attempt의 전체 궤적을 프레임별로 추출",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
사용 예시:
  python %(prog)s                                    # 현재 디렉토리에서 실행 (PDB 형식)
  python %(prog)s -i /path/to/sumd/output           # 특정 경로 지정
  python %(prog)s -f xtc                             # XTC 형식으로 출력
  python %(prog)s -i ./results -o my_frames -f pdb  # PDB 형식 명시적 지정
  python %(prog)s --format xtc --verbose             # XTC 형식, 상세 로그
        """
    )
    
    parser.add_argument(
        '-i', '--input',
        type=str,
        default='.',
        help='SuMD 출력 디렉토리 경로 (기본값: 현재 디렉토리)'
    )
    
    parser.add_argument(
        '-o', '--output',
        type=str,
        default='collected_trajectories',
        help='수집된 프레임을 저장할 디렉토리 이름 (기본값: collected_trajectories)'
    )
    
    parser.add_argument(
        '-f', '--format',
        type=str,
        choices=['pdb', 'xtc'],
        default='pdb',
        help='출력 형식: pdb (개별 프레임) 또는 xtc (연속 궤적) (기본값: pdb)'
    )
    
    parser.add_argument(
        '-v', '--verbose',
        action='store_true',
        help='상세한 로그 출력'
    )
    
    parser.add_argument(
        '--check-only',
        action='store_true',
        help='디렉토리 구조만 확인하고 실제 추출은 하지 않음'
    )
    
    args = parser.parse_args()
    
    # 입력 디렉토리 절대 경로로 변환
    input_dir = os.path.abspath(args.input)
    
    if not os.path.exists(input_dir):
        print(f"❌ 입력 디렉토리가 존재하지 않습니다: {input_dir}")
        return 1
    
    if not os.path.isdir(input_dir):
        print(f"❌ 입력 경로가 디렉토리가 아닙니다: {input_dir}")
        return 1
    
    # SuMD 출력 디렉토리인지 확인
    iteration_dirs = find_all_structure_directories(input_dir)
    
    if not iteration_dirs:
        print(f"❌ SuMD 출력 디렉토리가 아닌 것 같습니다: {input_dir}")
        print("   target_chains.pdb 또는 iteration_ 디렉토리들이 없습니다.")
        return 1
    
    print(f"📁 입력 디렉토리: {input_dir}")
    print(f"📂 Iteration 디렉토리 수: {len(iteration_dirs)}")
    print(f"📊 출력 형식: {args.format.upper()}")
    
    if args.check_only:
        print("\n🔍 디렉토리 구조 확인 완료 (--check-only 모드)")
        return 0
    
    # GROMACS 설치 확인
    if args.verbose:
        print("\n🔧 GROMACS 확인 중...")
    
    try:
        result = subprocess.run("gmx --version", shell=True, capture_output=True, timeout=10)
        if result.returncode != 0:
            print("❌ GROMACS가 설치되어 있지 않거나 PATH에 없습니다.")
            print("   gmx 명령어를 사용할 수 있는지 확인해주세요.")
            return 1
        elif args.verbose:
            print("✓ GROMACS 확인 완료")
    except subprocess.TimeoutExpired:
        print("❌ GROMACS 확인 시간 초과")
        return 1
    except Exception as e:
        print(f"❌ GROMACS 확인 중 오류: {e}")
        return 1
    
    # 궤적 수집 실행 - 형식에 따라 분기
    print(f"\n🚀 궤적 수집 시작 ({args.format.upper()} 형식)...")
    try:
        if args.format == 'xtc':
            output_dir, count = collect_trajectories_as_xtc(input_dir, args.output)
            print(f"\n✅ 수집 완료!")
            print(f"📊 연결된 궤적 수: {count:,}")
            print(f"📁 저장 위치: {output_dir}")
            print(f"\n💡 사용법:")
            print(f"   PyMOL 확인: pymol {os.path.basename(output_dir)}/initial_structure.pdb")
            print(f"   궤적 로드: PyMOL에서 'load {os.path.basename(output_dir)}/complete_trajectory.xtc'")
            print(f"   GROMACS 분석: gmx trjconv -s reference.tpr -f complete_trajectory.xtc ...")
        else:  # pdb
            output_dir, count = collect_sumd_trajectories(input_dir, args.output)
            print(f"\n✅ 수집 완료!")
            print(f"📊 총 프레임 수: {count:,}")
            print(f"📁 저장 위치: {output_dir}")
            print(f"\n💡 사용법:")
            print(f"   PyMOL 확인: pymol {os.path.basename(output_dir)}/frame_*.pdb")
            print(f"   애니메이션: PyMOL에서 'mplay' 명령")
            print(f"   압축: tar -czvf trajectory_frames.tar.gz {os.path.basename(output_dir)}/")
        
        return 0
        
    except KeyboardInterrupt:
        print("\n⚠️  사용자에 의해 중단됨")
        return 1
    except Exception as e:
        print(f"\n❌ 수집 중 오류 발생: {e}")
        if args.verbose:
            import traceback
            traceback.print_exc()
        return 1

if __name__ == "__main__":
    import sys
    sys.exit(main())