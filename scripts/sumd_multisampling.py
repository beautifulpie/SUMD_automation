#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import os
import subprocess
import numpy as np
from Bio.PDB import PDBParser
import logging
import sys
import time
from datetime import datetime
import shutil
import concurrent.futures
from multiprocessing import cpu_count, Process, Queue, Manager
import multiprocessing as mp
import random

# 로깅 설정
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("gromacs_simulation.log"),
        logging.StreamHandler(sys.stdout)
    ]
)
logger = logging.getLogger("GROMACS-MultiSample")

def run_command(cmd, shell=False, cwd=None, check=True, input_str=None, timeout=None):
    """명령어 실행 함수 - 개선된 버전"""
    cmd_str = cmd if isinstance(cmd, str) else ' '.join(cmd)
    logger.debug(f"실행 명령어: {cmd_str}")
    
    try:
        result = subprocess.run(
            cmd, 
            shell=shell, 
            cwd=cwd, 
            check=check, 
            input=input_str.encode() if input_str else None,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            timeout=timeout
        )
        
        if result.returncode == 0:
            logger.debug("명령어 성공적으로 실행됨")
        else:
            logger.error(f"명령어 실행 실패 (코드: {result.returncode})")
            stderr = result.stderr.decode() if result.stderr else ""
            if stderr:
                logger.error(f"오류 출력: {stderr}")
        
        return result
    except subprocess.TimeoutExpired:
        logger.error(f"명령어 실행 시간 초과: {cmd_str}")
        raise
    except Exception as e:
        logger.error(f"명령어 실행 중 오류 발생: {e}")
        raise

def calculate_distance_from_structure(structure_file, chain1, chain2, output_dir):
    """구조 파일(PDB/GRO)에서 두 체인 간의 거리 계산"""
    logger.debug(f"구조 파일에서 체인 {chain1}와 {chain2} 사이의 거리를 계산합니다.")
    
    file_ext = os.path.splitext(structure_file)[1].lower()
    
    if file_ext == '.pdb':
        # PDB 파일일 경우 BioPython 사용
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("structure", structure_file)
        
        model = structure[0]
        chain1_atoms = []
        chain2_atoms = []
        
        for chain in model:
            chain_id = chain.id
            
            if chain_id == chain1:
                for residue in chain:
                    for atom in residue:
                        chain1_atoms.append(atom.coord)
            elif chain_id == chain2:
                for residue in chain:
                    for atom in residue:
                        chain2_atoms.append(atom.coord)
        
        if not chain1_atoms or not chain2_atoms:
            logger.error(f"체인 {chain1} 또는 {chain2}의 원자를 찾을 수 없습니다.")
            return float('inf')
        
        # 무게중심 계산
        chain1_center = np.mean(chain1_atoms, axis=0)
        chain2_center = np.mean(chain2_atoms, axis=0)
        
        # 거리 계산 (nm 단위로 변환)
        distance = np.linalg.norm(chain1_center - chain2_center) / 10  # Å에서 nm로 변환
        logger.debug(f"체인 {chain1}와 {chain2} 사이의 무게중심 거리: {distance:.3f} nm")
        
        return distance
    
    elif file_ext == '.gro':
        # GRO 파일일 경우 GROMACS 도구 사용
        current_dir = os.getcwd()
        os.chdir(output_dir)
        
        try:
            # 체인 기반 인덱스 생성 및 거리 계산
            # 간단한 방법으로 처리 - 실제로는 더 정교한 체인 인식이 필요
            with open("distance.mdp", "w") as f:
                f.write("""
integrator          = steep
nsteps              = 0
                """)
            
            tpr_file = "distance.tpr"
            run_command([
                "gmx", "grompp",
                "-f", "distance.mdp",
                "-c", structure_file,
                "-o", tpr_file,
                "-maxwarn", "10"
            ], check=False)
            
            # 기본 거리 계산 (체인별 상세 분석 필요시 별도 구현)
            return 1.0  # 임시값 - 실제 구현에서는 proper chain analysis 필요
            
        except Exception as e:
            logger.error(f"GROMACS 거리 계산 중 오류 발생: {e}")
            return float('inf')
        finally:
            os.chdir(current_dir)
    
    else:
        logger.error(f"지원되지 않는 파일 형식: {file_ext}")
        return float('inf')

def prepare_pdb_system(pdb_file, output_dir, chain1, chain2, force_field="charmm36-jul2022"):
    """PDB 파일에서 시스템 준비 (토폴로지 생성, 박스 설정, 솔베이션, 이온화)"""
    logger.info(f"PDB 파일에서 시스템 준비 시작: 체인 {chain1}와 {chain2}를 사용합니다.")
    
    os.makedirs(output_dir, exist_ok=True)
    
    # 입력 PDB 파일 복사
    base_name = os.path.basename(pdb_file)
    pdb_path = os.path.join(output_dir, base_name)
    
    if pdb_file != pdb_path:
        run_command(f"cp {pdb_file} {pdb_path}", shell=True)
    
    current_dir = os.getcwd()
    os.chdir(output_dir)
    
    try:
        # 포스필드 리스트 - 순서대로 시도
        force_fields = [force_field, "amber99sb-ildn", "gromos54a7", "oplsaa"]
        
        success = False
        for ff in force_fields:
            try:
                logger.info(f"{ff} 포스필드로 시도합니다...")
                
                # pdb2gmx 실행: PDB에서 GRO 생성
                cmd_result = run_command([
                    "gmx", "pdb2gmx", 
                    "-f", base_name,
                    "-o", "complex.gro", 
                    "-p", "topol.top", 
                    "-water", "tip3p", 
                    "-ff", ff,
                    "-ignh"
                ], check=False)
                
                # 파일이 생성되었는지 확인
                if os.path.exists("complex.gro") and os.path.exists("topol.top"):
                    logger.info(f"{ff} 포스필드로 성공적으로 실행됨")
                    success = True
                    break
                else:
                    logger.warning(f"{ff} 포스필드로 실행했지만 필요한 파일이 생성되지 않았습니다.")
            except Exception as e:
                logger.warning(f"{ff} 포스필드로 시도 중 오류 발생: {e}")
        
        if not success:
            raise Exception(f"시스템 준비 실패: 체인 {chain1}와 {chain2}에 대한 모든 포스필드 시도가 실패했습니다.")
        
        # editconf 실행: 박스 설정
        run_command([
            "gmx", "editconf",
            "-f", "complex.gro",
            "-o", "box.gro",
            "-c",
            "-d", "1.0",
            "-bt", "cubic"
        ])
        
        # solvate 실행: 솔베이션
        run_command([
            "gmx", "solvate",
            "-cp", "box.gro",
            "-cs", "spc216.gro",
            "-o", "solv.gro",
            "-p", "topol.top"
        ])
        
        # 이온화를 위한 MDP 파일 생성
        with open("ions.mdp", "w") as f:
            f.write("""
; 이온화 설정
integrator          = steep
nsteps              = 1000
emtol               = 1000.0
emstep              = 0.02
nstlist             = 1
cutoff-scheme       = Verlet
ns_type             = grid
coulombtype         = cutoff
rcoulomb            = 1.0
rvdw                = 1.0
pbc                 = xyz
            """)
        
        # grompp 실행: 이온화 준비
        run_command([
            "gmx", "grompp",
            "-f", "ions.mdp",
            "-c", "solv.gro",
            "-p", "topol.top",
            "-o", "ions.tpr",
            "-maxwarn", "10"
        ])
        
        # genion 실행: 이온 추가
        run_command([
            "gmx", "genion",
            "-s", "ions.tpr",
            "-o", "solv_ions.gro",
            "-p", "topol.top",
            "-pname", "NA",
            "-nname", "CL",
            "-neutral"
        ], input_str="SOL")
        
        return "solv_ions.gro"
    
    finally:
        os.chdir(current_dir)

def run_single_sample_worker(args):
    """멀티프로세싱을 위한 워커 함수"""
    sample_id, input_gro, output_dir, simulation_time, chain1, chain2, distance_threshold_angstrom = args
    
    # 프로세스별 로거 설정
    process_logger = logging.getLogger(f"Sample-{sample_id}")
    
    sample_dir = os.path.join(output_dir, f"sample_{sample_id}")
    os.makedirs(sample_dir, exist_ok=True)
    
    process_logger.info(f"샘플 {sample_id} 시뮬레이션 시작 (PID: {os.getpid()})")
    
    current_dir = os.getcwd()
    
    try:
        # 필요한 파일들을 샘플 디렉토리로 복사
        shutil.copy2(input_gro, sample_dir)
        shutil.copy2(os.path.join(output_dir, "topol.top"), sample_dir)
        if os.path.exists(os.path.join(output_dir, "posre.itp")):
            shutil.copy2(os.path.join(output_dir, "posre.itp"), sample_dir)
        
        os.chdir(sample_dir)
        
        # 초기 거리 계산
        initial_distance = calculate_distance_from_structure(os.path.basename(input_gro), chain1, chain2, sample_dir)
        initial_distance_angstrom = initial_distance * 10  # nm를 Å로 변환
        
        process_logger.info(f"샘플 {sample_id} 초기 거리: {initial_distance_angstrom:.1f}Å ({initial_distance:.3f}nm)")
        
        # 에너지 최소화 (항상 수행)
        em_gro, em_tpr = run_energy_minimization_process(sample_id, process_logger)
        em_distance = calculate_distance_from_structure(em_gro, chain1, chain2, sample_dir)
        em_distance_angstrom = em_distance * 10
        
        # 거리에 따라 MD 시뮬레이션 수행 여부 결정
        if initial_distance_angstrom <= distance_threshold_angstrom:
            # 5Å 이하: MD 시뮬레이션까지 수행
            process_logger.info(f"샘플 {sample_id}: 초기 거리 {initial_distance_angstrom:.1f}Å ≤ {distance_threshold_angstrom}Å → MD 시뮬레이션 수행")
            
            md_gro, md_xtc, md_tpr = run_md_simulation_process(em_gro, sample_id, simulation_time, process_logger)
            md_distance = calculate_distance_from_structure(md_gro, chain1, chain2, sample_dir)
            md_distance_angstrom = md_distance * 10
            
            process_logger.info(f"샘플 {sample_id} 완료: 초기={initial_distance_angstrom:.1f}Å, EM={em_distance_angstrom:.1f}Å, MD={md_distance_angstrom:.1f}Å")
            
            return {
                'sample_id': sample_id,
                'initial_distance': initial_distance,
                'em_distance': em_distance,
                'md_distance': md_distance,
                'final_distance': md_distance,
                'final_structure': os.path.join(sample_dir, md_gro),
                'trajectory': os.path.join(sample_dir, md_xtc),
                'sample_dir': sample_dir,
                'simulation_type': 'EM+MD'
            }
        else:
            # 5Å 초과: EM만 수행
            process_logger.info(f"샘플 {sample_id}: 초기 거리 {initial_distance_angstrom:.1f}Å > {distance_threshold_angstrom}Å → EM만 수행")
            
            process_logger.info(f"샘플 {sample_id} 완료: 초기={initial_distance_angstrom:.1f}Å, EM={em_distance_angstrom:.1f}Å (MD 건너뜀)")
            
            return {
                'sample_id': sample_id,
                'initial_distance': initial_distance,
                'em_distance': em_distance,
                'md_distance': None,
                'final_distance': em_distance,
                'final_structure': os.path.join(sample_dir, em_gro),
                'trajectory': None,
                'sample_dir': sample_dir,
                'simulation_type': 'EM_only'
            }
    
    except Exception as e:
        process_logger.error(f"샘플 {sample_id} 실행 중 오류: {e}")
        return {
            'sample_id': sample_id,
            'initial_distance': float('inf'),
            'em_distance': float('inf'),
            'md_distance': float('inf'),
            'final_distance': float('inf'),
            'final_structure': None,
            'trajectory': None,
            'sample_dir': sample_dir,
            'simulation_type': 'failed',
            'error': str(e)
        }
    
    finally:
        os.chdir(current_dir)

def run_energy_minimization_process(sample_id, process_logger):
    """에너지 최소화 실행 - 프로세스별 실행"""
    process_logger.debug(f"샘플 {sample_id} 에너지 최소화 실행")
    
    # 에너지 최소화를 위한 MDP 파일 생성
    with open("em.mdp", "w") as f:
        f.write("""
; 에너지 최소화 설정
integrator          = steep
nsteps              = 5000
emtol               = 1000.0
emstep              = 0.01
nstlist             = 1
cutoff-scheme       = Verlet
ns_type             = grid
coulombtype         = PME
rcoulomb            = 1.0
rvdw                = 1.0
pbc                 = xyz
        """)
    
    input_gro = [f for f in os.listdir('.') if f.endswith('.gro')][0]
    
    # grompp 실행: 에너지 최소화 준비
    run_command([
        "gmx", "grompp",
        "-f", "em.mdp",
        "-c", input_gro,
        "-p", "topol.top",
        "-o", "em.tpr",
        "-maxwarn", "10"
    ])
    
    # mdrun 실행: 에너지 최소화 (단일 스레드 강제)
    run_command([
        "gmx", "mdrun",
        "-v",
        "-deffnm", "em",
        "-ntmpi", "1",
        "-ntomp", "1"
    ])
    
    return "em.gro", "em.tpr"

def run_md_simulation_process(input_gro, sample_id, simulation_time, process_logger):
    """MD 시뮬레이션 실행 - 프로세스별 실행"""
    process_logger.debug(f"샘플 {sample_id} MD 시뮬레이션 실행 ({simulation_time} ns)")
    
    # 시뮬레이션 시간에 따른 스텝 수 계산
    dt = 0.002  # 2 fs 타임스텝
    nsteps = int(simulation_time * 1000 / dt)
    
    # 프로세스별로 다른 초기 속도를 위한 랜덤 시드
    gen_seed = random.randint(1, 999999) + sample_id * 1000
    
    # MD 시뮬레이션을 위한 MDP 파일 생성
    with open("md.mdp", "w") as f:
        f.write(f"""
; SuMD 다중 샘플 시뮬레이션 설정
integrator          = md
dt                  = {dt}
nsteps              = {nsteps}
nstenergy           = 1000
nstlog              = 1000
nstxout-compressed  = 1000
continuation        = no
constraint_algorithm = lincs
constraints         = all-bonds
lincs_iter          = 1
lincs_order         = 4

; 속도 생성 (프로세스별 다른 시드)
gen_vel             = yes
gen_temp            = 300
gen_seed            = {gen_seed}

; 온도 커플링
tcoupl              = V-rescale
tc-grps             = System
tau_t               = 0.1
ref_t               = 300

; 압력 커플링
pcoupl              = Parrinello-Rahman
pcoupltype          = isotropic
tau_p               = 2.0
ref_p               = 1.0
compressibility     = 4.5e-5

; 비결합 상호작용
cutoff-scheme       = Verlet
nstlist             = 10
ns_type             = grid
coulombtype         = PME
rcoulomb            = 1.0
rvdw                = 1.0
DispCorr            = EnerPres
pbc                 = xyz
        """)
    
    # grompp 실행: MD 준비
    run_command([
        "gmx", "grompp",
        "-f", "md.mdp",
        "-c", input_gro,
        "-p", "topol.top",
        "-o", "md.tpr",
        "-maxwarn", "10"
    ])
    
    # mdrun 실행: MD 시뮬레이션 (단일 스레드 강제)
    run_command([
        "gmx", "mdrun",
        "-v",
        "-deffnm", "md",
        "-ntmpi", "1",
        "-ntomp", "1"  # 프로세스별 단일 스레드
    ])
    
    return "md.gro", "md.xtc", "md.tpr"

def run_multi_sample_iteration(input_file, output_dir, chain1, chain2, iteration, num_samples=5, simulation_time=2, distance_threshold_angstrom=5.0):
    """다중 샘플 반복 실행 - 프로세스 풀 사용"""
    logger.info(f"반복 {iteration}: {num_samples}개 샘플 실행 시작 (거리 임계값: {distance_threshold_angstrom}Å)")
    
    # 현재 반복의 디렉토리 생성
    iter_dir = os.path.join(output_dir, f"iteration_{iteration}")
    os.makedirs(iter_dir, exist_ok=True)
    
    # 입력 파일이 PDB인 경우 시스템 준비
    file_ext = os.path.splitext(input_file)[1].lower()
    if file_ext == '.pdb':
        prepared_gro = prepare_pdb_system(input_file, iter_dir, chain1, chain2)
        input_gro = os.path.join(iter_dir, prepared_gro)
    else:
        # GRO 파일인 경우 복사
        input_gro = os.path.join(iter_dir, os.path.basename(input_file))
        shutil.copy2(input_file, input_gro)
        # 관련 파일들도 복사
        input_dir = os.path.dirname(input_file)
        for file_name in ["topol.top", "posre.itp"]:
            src_file = os.path.join(input_dir, file_name)
            if os.path.exists(src_file):
                shutil.copy2(src_file, iter_dir)
    
    # CPU 개수에 따라 프로세스 수 결정
    max_workers = min(num_samples, cpu_count())
    
    # 프로세스 풀을 사용한 다중 샘플 실행
    sample_results = []
    
    # 작업 인자 준비
    work_args = []
    for sample_id in range(1, num_samples + 1):
        args = (sample_id, input_gro, iter_dir, simulation_time, chain1, chain2, distance_threshold_angstrom)
        work_args.append(args)
    
    # ProcessPoolExecutor 사용
    with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as executor:
        logger.info(f"프로세스 풀 시작: {max_workers}개 워커 프로세스")
        
        # 모든 작업 제출
        futures = [executor.submit(run_single_sample_worker, args) for args in work_args]
        
        # 결과 수집 (완료 순서대로)
        for future in concurrent.futures.as_completed(futures):
            try:
                result = future.result(timeout=3600)  # 1시간 타임아웃
                sample_results.append(result)
                logger.info(f"샘플 {result['sample_id']} 완료 ({result['simulation_type']})")
            except concurrent.futures.TimeoutError:
                logger.error("샘플 실행이 시간 초과되었습니다.")
            except Exception as e:
                logger.error(f"샘플 실행 중 오류: {e}")
    
    # 유효한 결과만 필터링
    valid_results = [r for r in sample_results if r['final_structure'] is not None]
    
    if not valid_results:
        logger.error("모든 샘플이 실패했습니다.")
        return None, float('inf')
    
    # 최적의 샘플 선택 (final_distance가 가장 짧은 것)
    best_result = min(valid_results, key=lambda x: x['final_distance'])
    
    # 시뮬레이션 유형별 통계
    em_only_count = len([r for r in valid_results if r['simulation_type'] == 'EM_only'])
    em_md_count = len([r for r in valid_results if r['simulation_type'] == 'EM+MD'])
    failed_count = len([r for r in sample_results if r['simulation_type'] == 'failed'])
    
    logger.info(f"반복 {iteration} 완료:")
    logger.info(f"  총 {len(sample_results)}개 샘플 중 {len(valid_results)}개 성공, {failed_count}개 실패")
    logger.info(f"  EM만 수행: {em_only_count}개, EM+MD 수행: {em_md_count}개")
    logger.info(f"  최적 샘플: {best_result['sample_id']} ({best_result['simulation_type']}, 거리: {best_result['final_distance']*10:.1f}Å)")
    
    # 결과 요약 저장
    with open(os.path.join(iter_dir, "iteration_summary.txt"), "w") as f:
        f.write(f"=== 반복 {iteration} 요약 ===\n\n")
        f.write(f"거리 임계값: {distance_threshold_angstrom}Å\n")
        f.write(f"실행된 샘플 수: {len(sample_results)}\n")
        f.write(f"성공한 샘플 수: {len(valid_results)}\n")
        f.write(f"실패한 샘플 수: {failed_count}\n")
        f.write(f"EM만 수행 샘플 수: {em_only_count}\n")
        f.write(f"EM+MD 수행 샘플 수: {em_md_count}\n")
        f.write(f"최적 샘플 ID: {best_result['sample_id']}\n")
        f.write(f"최적 샘플 시뮬레이션 유형: {best_result['simulation_type']}\n")
        f.write(f"최적 거리: {best_result['final_distance']*10:.1f}Å ({best_result['final_distance']:.3f}nm)\n\n")
        
        f.write("모든 샘플 결과:\n")
        for result in sample_results:
            if 'error' in result:
                f.write(f"  샘플 {result['sample_id']}: 실패 ({result['error']})\n")
            else:
                initial_ang = result['initial_distance'] * 10
                final_ang = result['final_distance'] * 10
                f.write(f"  샘플 {result['sample_id']}: {result['simulation_type']} - {initial_ang:.1f}Å → {final_ang:.1f}Å\n")
    
    return best_result['final_structure'], best_result['final_distance']_workers = min(num_samples, cpu_count())
    
    # 다중 샘플 병렬 실행
    sample_results = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
        # 모든 샘플 작업 제출
        futures = []
        for sample_id in range(1, num_samples + 1):
            future = executor.submit(
                run_single_sample, 
                sample_id, 
                input_gro, 
                iter_dir, 
                simulation_time, 
                chain1, 
                chain2,
                distance_threshold_angstrom
            )
            futures.append(future)
        
        # 결과 수집
        for future in concurrent.futures.as_completed(futures):
            try:
                result = future.result()
                sample_results.append(result)
            except Exception as e:
                logger.error(f"샘플 실행 중 오류: {e}")
    
    # 유효한 결과만 필터링
    valid_results = [r for r in sample_results if r['final_structure'] is not None]
    
    if not valid_results:
        logger.error("모든 샘플이 실패했습니다.")
        return None, float('inf')
    
    # 최적의 샘플 선택 (final_distance가 가장 짧은 것)
    best_result = min(valid_results, key=lambda x: x['final_distance'])
    
    # 시뮬레이션 유형별 통계
    em_only_count = len([r for r in valid_results if r['simulation_type'] == 'EM_only'])
    em_md_count = len([r for r in valid_results if r['simulation_type'] == 'EM+MD'])
    
    logger.info(f"반복 {iteration} 완료:")
    logger.info(f"  총 {len(sample_results)}개 샘플 중 {len(valid_results)}개 성공")
    logger.info(f"  EM만 수행: {em_only_count}개, EM+MD 수행: {em_md_count}개")
    logger.info(f"  최적 샘플: {best_result['sample_id']} ({best_result['simulation_type']}, 거리: {best_result['final_distance']*10:.1f}Å)")
    
    # 결과 요약 저장
    with open(os.path.join(iter_dir, "iteration_summary.txt"), "w") as f:
        f.write(f"=== 반복 {iteration} 요약 ===\n\n")
        f.write(f"거리 임계값: {distance_threshold_angstrom}Å\n")
        f.write(f"실행된 샘플 수: {len(sample_results)}\n")
        f.write(f"성공한 샘플 수: {len(valid_results)}\n")
        f.write(f"EM만 수행 샘플 수: {em_only_count}\n")
        f.write(f"EM+MD 수행 샘플 수: {em_md_count}\n")
        f.write(f"최적 샘플 ID: {best_result['sample_id']}\n")
        f.write(f"최적 샘플 시뮬레이션 유형: {best_result['simulation_type']}\n")
        f.write(f"최적 거리: {best_result['final_distance']*10:.1f}Å ({best_result['final_distance']:.3f}nm)\n\n")
        
        f.write("모든 샘플 결과:\n")
        for result in sample_results:
            if 'error' in result:
                f.write(f"  샘플 {result['sample_id']}: 실패 ({result['error']})\n")
            else:
                initial_ang = result['initial_distance'] * 10
                final_ang = result['final_distance'] * 10
                f.write(f"  샘플 {result['sample_id']}: {result['simulation_type']} - {initial_ang:.1f}Å → {final_ang:.1f}Å\n")
    
    return best_result['final_structure'], best_result['final_distance']

def main():
    # 멀티프로세싱 시작 방법 설정
    mp.set_start_method('spawn', force=True)
    
    parser = argparse.ArgumentParser(description="거리 기반 다중 샘플 SuMD Gromacs 시뮬레이션 (프로세스 기반)")
    parser.add_argument("--input", required=True, help="입력 구조 파일 (PDB 또는 GRO)")
    parser.add_argument("--chain1", required=True, help="첫 번째 체인 ID (예: A)")
    parser.add_argument("--chain2", required=True, help="두 번째 체인 ID (예: B)")
    parser.add_argument("--output_dir", default="./sumd_output", help="출력 디렉토리")
    parser.add_argument("--distance_threshold", type=float, default=0.5, help="수렴 거리 임계값 (nm)")
    parser.add_argument("--simulation_threshold", type=float, default=5.0, help="MD 실행 거리 임계값 (Å)")
    parser.add_argument("--max_iterations", type=int, default=10, help="최대 반복 횟수")
    parser.add_argument("--num_samples", type=int, default=5, help="각 반복당 샘플 수")
    parser.add_argument("--simulation_time", type=float, default=2, help="각 샘플의 시뮬레이션 시간 (ns)")
    parser.add_argument("--force_field", default="charmm36-jul2022", help="사용할 포스필드")
    parser.add_argument("--max_workers", type=int, default=None, help="최대 워커 프로세스 수 (기본: CPU 개수)")
    
    args = parser.parse_args()
    
    # 최대 워커 수 설정
    if args.max_workers is None:
        max_workers_limit = min(args.num_samples, cpu_count())
    else:
        max_workers_limit = min(args.max_workers, args.num_samples, cpu_count())
    
    # 출력 디렉토리 생성
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_dir = f"{args.output_dir}_{timestamp}"
    os.makedirs(output_dir, exist_ok=True)
    
    logger.info("거리 기반 다중 샘플 SuMD 시뮬레이션 시작 (프로세스 기반)")
    logger.info(f"입력 파일: {args.input}")
    logger.info(f"체인: {args.chain1} - {args.chain2}")
    logger.info(f"수렴 거리 임계값: {args.distance_threshold} nm")
    logger.info(f"MD 실행 거리 임계값: {args.simulation_threshold} Å")
    logger.info(f"최대 반복: {args.max_iterations}")
    logger.info(f"샘플 수/반복: {args.num_samples}")
    logger.info(f"시뮬레이션 시간/샘플: {args.simulation_time} ns")
    logger.info(f"최대 워커 프로세스: {max_workers_limit}")
    logger.info(f"사용 가능 CPU: {cpu_count()}")
    
    # 반복 시뮬레이션 실행
    current_structure = args.input
    convergence = False
    
    for iteration in range(1, args.max_iterations + 1):
        logger.info(f"\n=== 반복 {iteration}/{args.max_iterations} 시작 ===")
        
        best_structure, best_distance = run_multi_sample_iteration(
            current_structure,
            output_dir,
            args.chain1,
            args.chain2,
            iteration,
            args.num_samples,
            args.simulation_time,
            args.simulation_threshold
        )
        
        if best_structure is None:
            logger.error(f"반복 {iteration}에서 모든 샘플이 실패했습니다. 시뮬레이션을 종료합니다.")
            break
        
        best_distance_angstrom = best_distance * 10
        convergence_threshold_angstrom = args.distance_threshold * 10
        
        logger.info(f"반복 {iteration} 최적 거리: {best_distance_angstrom:.1f}Å ({best_distance:.3f}nm)")
        
        # 수렴 확인
        if best_distance <= args.distance_threshold:
            logger.info(f"수렴 달성! (거리: {best_distance_angstrom:.1f}Å ≤ {convergence_threshold_angstrom:.1f}Å)")
            convergence = True
            break
        
        # 다음 반복을 위한 구조 업데이트
        current_structure = best_structure
    
    # 최종 결과 출력
    logger.info(f"\n=== 거리 기반 SuMD 시뮬레이션 완료 ===")
    if convergence:
        logger.info(f"상태: 수렴 완료 (반복 {iteration})")
        logger.info(f"최종 거리: {best_distance_angstrom:.1f}Å ({best_distance:.3f}nm)")
        logger.info(f"최종 구조: {best_structure}")
        print(f"SUMD_RESULT:{best_structure}:{best_distance}")
    else:
        logger.info(f"상태: 최대 반복 횟수 도달 (수렴하지 않음)")
        logger.info(f"최종 거리: {best_distance_angstrom:.1f}Å ({best_distance:.3f}nm)")
        print(f"SUMD_RESULT:{best_structure}:{best_distance}")
    
    logger.info(f"결과 디렉토리: {output_dir}")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())