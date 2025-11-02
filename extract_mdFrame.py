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
def concatenate_trajectories_to_xtc(trajectory_files, output_xtc, temp_dir):
    """여러 XTC 파일들을 하나의 XTC로 연결"""
    if not trajectory_files:
        return False
    
    print(f"  XTC 파일 연결 중: {len(trajectory_files)}개 궤적")
    
    # trjcat으로 여러 XTC를 하나로 연결
    traj_list = " ".join(trajectory_files)
    cmd = f"gmx trjcat -f {traj_list} -o {output_xtc} -settime"
    
    # settime 입력: 모든 프레임을 연속적으로 (0 입력)
    input_text = "0\n" * len(trajectory_files)
    
    try:
        result = subprocess.run(
            cmd, shell=True, cwd=temp_dir,
            input=input_text.encode(), 
            capture_output=True, text=True
        )
        if result.returncode != 0:
            print(f"  ⚠ XTC 연결 실패: {result.stderr}")
            return False
        return True
    except Exception as e:
        print(f"  ⚠ XTC 연결 중 오류: {e}")
        return False

def collect_trajectories_as_xtc(sumd_output_dir, output_dir="collected_trajectories"):
    """
    SuMD 시뮬레이션의 모든 성공한 attempt 궤적을 XTC 파일로 수집
    
    Args:
        sumd_output_dir: SuMD 출력 디렉토리 경로
        output_dir: 수집된 궤적을 저장할 디렉토리 이름
    """
    print(f"SuMD 궤적 XTC 수집 시작: {sumd_output_dir}")
    
    # 출력 디렉토리 생성
    full_output_dir = os.path.join(sumd_output_dir, output_dir)
    if os.path.exists(full_output_dir):
        shutil.rmtree(full_output_dir)
    os.makedirs(full_output_dir)
    
    # 임시 작업 디렉토리
    temp_dir = os.path.join(full_output_dir, "temp_xtc")
    os.makedirs(temp_dir)
    
    collected_info = {
        "collection_time": datetime.now().isoformat(),
        "source_directory": sumd_output_dir,
        "format": "xtc",
        "iterations": []
    }
    
    trajectory_files = []  # 연결할 XTC 파일들
    
    # 1. 초기 구조 (target_chains.pdb) -> initial_structure.pdb
    target_chains_path = os.path.join(sumd_output_dir, "target_chains.pdb")
    if os.path.exists(target_chains_path):
        initial_pdb = os.path.join(full_output_dir, "initial_structure.pdb")
        shutil.copy(target_chains_path, initial_pdb)
        
        collected_info["initial_structure"] = "initial_structure.pdb"
        print(f"✓ 초기 구조 저장: initial_structure.pdb")
    else:
        print("⚠ target_chains.pdb를 찾을 수 없습니다.")
    
    # 2. Iteration별 성공한 attempt의 궤적 수집
    iteration_dirs = []
    for item in os.listdir(sumd_output_dir):
        if item.startswith("iteration_") and os.path.isdir(os.path.join(sumd_output_dir, item)):
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
        
        iter_path = os.path.join(sumd_output_dir, iter_dir)
        summary_file = os.path.join(sumd_output_dir, f"{iter_dir}_summary.json")
        
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
        
        # XTC 파일 직접 복사
        src_xtc = os.path.join(attempt_dir, "md.xtc")
        src_tpr = os.path.join(attempt_dir, "md.tpr")
        
        if not os.path.exists(src_xtc) or not os.path.exists(src_tpr):
            print(f"  ⚠ Iteration {iter_num}: XTC 또는 TPR 파일이 없습니다")
            continue
        
        # 임시 디렉토리에 복사
        temp_xtc = os.path.join(temp_dir, f"iter_{iter_num:03d}.xtc")
        shutil.copy(src_xtc, temp_xtc)
        trajectory_files.append(temp_xtc)
        
        # iteration 정보 기록
        iteration_info = {
            "type": "md_trajectory",
            "iteration_number": iter_num,
            "attempt_number": attempt_num,
            "source": f"{iter_dir}/attempt_{attempt_num}/md.xtc",
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
        print(f"  ✓ XTC 파일 추가: iter_{iter_num:03d}.xtc {long_md_mark}")
    
    # 모든 궤적을 하나의 XTC로 연결
    if trajectory_files:
        final_xtc = os.path.join(full_output_dir, "complete_trajectory.xtc")
        print(f"\n연속 궤적 생성 중: {len(trajectory_files)}개 XTC 연결...")
        
        if concatenate_trajectories_to_xtc(trajectory_files, final_xtc, temp_dir):
            print(f"✓ 연속 궤적 생성 완료: complete_trajectory.xtc")
            collected_info["trajectory_file"] = "complete_trajectory.xtc"
            collected_info["total_iterations"] = len(trajectory_files)
            
            # 첫 번째 TPR 파일도 복사 (시각화용)
            first_iter_dir = os.path.join(sumd_output_dir, iteration_dirs[0][1])
            first_attempt_info = find_successful_attempt(
                first_iter_dir, 
                json.load(open(os.path.join(sumd_output_dir, f"{iteration_dirs[0][1]}_summary.json")))
            )
            if first_attempt_info:
                first_tpr = os.path.join(first_attempt_info[0], "md.tpr")
                ref_tpr = os.path.join(full_output_dir, "reference.tpr")
                shutil.copy(first_tpr, ref_tpr)
                collected_info["reference_tpr"] = "reference.tpr"
                print(f"✓ 참조 TPR 저장: reference.tpr")
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
        f.write(f"SuMD 궤적 요약 (XTC 형식)\n")
        f.write(f"수집 시간: {collected_info['collection_time']}\n")
        f.write(f"처리된 iteration 수: {collected_info.get('total_iterations', 0)}\n")
        f.write(f"출력 파일: complete_trajectory.xtc\n")
        f.write("="*50 + "\n\n")
        
        for info in collected_info["iterations"]:
            f.write(f"Iteration {info['iteration_number']}: {info['description']}\n")
            if 'min_distance' in info:
                f.write(f"  - 최소거리: {info['min_distance']:.2f}Å\n")
            if 'slope' in info:
                f.write(f"  - 기울기: {info['slope']:.6f}\n")
            if info.get("is_long_md"):
                f.write(f"  - Long MD 궤적\n")
            f.write("\n")
    
    print(f"\n=== 수집 완료 ===")
    print(f"처리된 iteration 수: {collected_info.get('total_iterations', 0)}")
    print(f"저장 위치: {full_output_dir}")
    print(f"궤적 파일: complete_trajectory.xtc")
    print(f"상세 정보: {info_file}")
    print(f"요약 파일: {summary_file}")
    
    return full_output_dir, collected_info.get('total_iterations', 0)

def collect_sumd_trajectories(sumd_output_dir, output_dir="collected_trajectories"):
    """
    SuMD 시뮬레이션의 모든 성공한 attempt 궤적을 프레임별로 수집    
    Args:
        sumd_output_dir: SuMD 출력 디렉토리 경로
        output_dir: 수집된 구조들을 저장할 디렉토리 이름
    """
    print(f"SuMD 궤적 프레임 수집 시작: {sumd_output_dir}")
    
    # 출력 디렉토리 생성
    full_output_dir = os.path.join(sumd_output_dir, output_dir)
    if os.path.exists(full_output_dir):
        shutil.rmtree(full_output_dir)
    os.makedirs(full_output_dir)
    
    # 임시 작업 디렉토리
    temp_dir = os.path.join(full_output_dir, "temp_extraction")
    os.makedirs(temp_dir)
    
    structure_count = 0
    collected_info = {
        "collection_time": datetime.now().isoformat(),
        "source_directory": sumd_output_dir,
        "total_frames": 0,
        "iterations": []
    }
    
    # 1. 초기 구조 (target_chains.pdb) -> frame_000000.pdb
    target_chains_path = os.path.join(sumd_output_dir, "target_chains.pdb")
    if os.path.exists(target_chains_path):
        dest_path = os.path.join(full_output_dir, f"frame_{structure_count:06d}.pdb")
        shutil.copy(target_chains_path, dest_path)
        
        collected_info["iterations"].append({
            "type": "initial_structure",
            "source": "target_chains.pdb",
            "frame_start": structure_count,
            "frame_count": 1,
            "description": "초기 타겟 체인 구조"
        })
        
        print(f"✓ Frame {structure_count:06d}: 초기 구조 (target_chains.pdb)")
        structure_count += 1
    else:
        print("⚠ target_chains.pdb를 찾을 수 없습니다.")
    
    # 2. Iteration별 성공한 attempt의 궤적 수집
    iteration_dirs = []
    for item in os.listdir(sumd_output_dir):
        if item.startswith("iteration_") and os.path.isdir(os.path.join(sumd_output_dir, item)):
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
        
        iter_path = os.path.join(sumd_output_dir, iter_dir)
        summary_file = os.path.join(sumd_output_dir, f"{iter_dir}_summary.json")
        
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
    collected_info["total_iterations"] = len([info for info in collected_info["iterations"] if info["type"] == "md_trajectory"])
    
    info_file = os.path.join(full_output_dir, "collection_info.json")
    with open(info_file, 'w') as f:
        json.dump(collected_info, f, indent=2, ensure_ascii=False)
    
    # 요약 텍스트 파일 생성
    summary_file = os.path.join(full_output_dir, "frames_summary.txt")
    with open(summary_file, 'w', encoding='utf-8') as f:
        f.write(f"SuMD 궤적 프레임 요약\n")
        f.write(f"수집 시간: {collected_info['collection_time']}\n")
        f.write(f"총 프레임 수: {structure_count}\n")
        f.write(f"처리된 iteration 수: {collected_info['total_iterations']}\n")
        f.write("="*50 + "\n\n")
        
        for info in collected_info["iterations"]:
            if info["type"] == "initial_structure":
                f.write(f"Frame {info['frame_start']:06d}: {info['description']}\n")
            else:
                f.write(f"Frames {info['frame_start']:06d}-{info['frame_start']+info['frame_count']-1:06d}: {info['description']}\n")
                if 'min_distance' in info:
                    f.write(f"  - 최소거리: {info['min_distance']:.2f}Å\n")
                if 'slope' in info:
                    f.write(f"  - 기울기: {info['slope']:.6f}\n")
                if info.get("is_long_md"):
                    f.write(f"  - Long MD 궤적\n")
            f.write("\n")
    
    print(f"\n=== 수집 완료 ===")
    print(f"총 수집된 프레임 수: {structure_count}")
    print(f"처리된 iteration 수: {collected_info['total_iterations']}")
    print(f"저장 위치: {full_output_dir}")
    print(f"상세 정보: {info_file}")
    print(f"요약 파일: {summary_file}")
    
    return full_output_dir, structure_count

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
    target_chains = os.path.join(input_dir, "target_chains.pdb")
    iteration_dirs = [f for f in os.listdir(input_dir) 
                     if f.startswith("iteration_") and os.path.isdir(os.path.join(input_dir, f))]
    
    if not os.path.exists(target_chains) and not iteration_dirs:
        print(f"❌ SuMD 출력 디렉토리가 아닌 것 같습니다: {input_dir}")
        print("   target_chains.pdb 또는 iteration_ 디렉토리들이 없습니다.")
        return 1
    
    print(f"📁 입력 디렉토리: {input_dir}")
    print(f"📄 Target chains: {'✓' if os.path.exists(target_chains) else '✗'}")
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