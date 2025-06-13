#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import os
import subprocess
import numpy as np
from Bio.PDB import PDBParser
import time
import logging
import sys

# 로깅 설정
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("sumd_gromacs.log"),
        logging.StreamHandler(sys.stdout)
    ]
)
logger = logging.getLogger("SuMD-Gromacs")

def run_command(cmd, shell=False, cwd=None, check=True, input_str=None, timeout=None):
    """명령어 실행 함수"""
    logger.info(f"실행 명령어: {cmd if isinstance(cmd, str) else ' '.join(cmd)}")
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
            logger.info("명령어 성공적으로 실행됨")
        else:
            logger.error(f"명령어 실행 실패: {result.stderr.decode()}")
            logger.error(f"오류 세부 정보: {result.stdout.decode()}")
        return result
    except subprocess.TimeoutExpired:
        logger.error(f"명령어 실행 시간 초과: {cmd}")
        raise
    except Exception as e:
        logger.error(f"명령어 실행 중 오류 발생: {e}")
        raise

def calculate_distance(pdb_file, peptide_res_list, protein_res_list):
    """PDB 파일에서 두 그룹 간의 무게중심 거리 계산"""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("structure", pdb_file)
    
    # 첫 번째 모델 선택
    model = structure[0]
    
    peptide_atoms = []
    protein_atoms = []
    
    # 원자 그룹화
    for chain in model:
        for residue in chain:
            res_id = str(residue.id[1])
            if res_id in peptide_res_list:
                for atom in residue:
                    peptide_atoms.append(atom.coord)
            elif res_id in protein_res_list:
                for atom in residue:
                    protein_atoms.append(atom.coord)
    
    if not peptide_atoms or not protein_atoms:
        logger.error("펩타이드 또는 단백질 원자를 찾을 수 없습니다.")
        return float('inf')
    
    # 무게중심 계산
    peptide_center = np.mean(peptide_atoms, axis=0)
    protein_center = np.mean(protein_atoms, axis=0)
    
    # 거리 계산 (nm 단위로 변환)
    distance = np.linalg.norm(peptide_center - protein_center) / 10  # Å에서 nm로 변환
    logger.info(f"계산된 무게중심 거리: {distance:.3f} nm")
    
    return distance

def create_index_file(output_file, peptide_res_list, protein_res_list):
    """인덱스 파일 생성 (Gromacs make_ndx 대체)"""
    with open(output_file, 'w') as f:
        f.write("[ Peptide ]\n")
        f.write(" ".join(peptide_res_list) + "\n")
        f.write("[ Protein ]\n")
        f.write(" ".join(protein_res_list) + "\n")

def get_chain_info(pdb_file):
    """PDB 파일에서 체인 정보 추출"""
    try:
        # gmx pdb2gmx로 체인 정보 확인 (실제 변환하지 않고 정보만 추출)
        result = subprocess.run(
            ["gmx", "pdb2gmx", "-f", pdb_file, "-q"],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False
        )
        
        output = result.stdout.decode() + result.stderr.decode()
        
        # 체인 정보 파싱
        chains = []
        start_parsing = False
        for line in output.split('\n'):
            if "chain  #res #atoms" in line:
                start_parsing = True
                continue
            if start_parsing and line.strip() and "chain" not in line:
                parts = line.split()
                if len(parts) >= 3 and "'" in line:
                    chain_id = line.split("'")[1].strip()
                    if chain_id and not chain_id.isspace() and "(only water)" not in line:
                        chains.append((len(chains) + 1, chain_id))
                        
            if start_parsing and "water" in line and len(chains) > 0:
                break
        
        return chains
    except Exception as e:
        logger.error(f"체인 정보 추출 중 오류 발생: {e}")
        return [(1, 'A')]  # 기본값 반환

def prepare_system(pdb_file, output_dir, force_field="charmm36-jul2022", selected_chain=None):
    """시스템 준비 (토폴로지 생성, 박스 설정, 솔베이션, 이온화)"""
    os.makedirs(output_dir, exist_ok=True)
    
    # 파일 복사
    base_name = os.path.basename(pdb_file)
    pdb_path = os.path.join(output_dir, base_name)
    
    if pdb_file != pdb_path:
        run_command(f"cp {pdb_file} {pdb_path}", shell=True)
    
    # 작업 디렉토리 변경
    os.chdir(output_dir)
    
    # 체인 정보 확인
    chains = get_chain_info(pdb_path)
    
    if not chains:
        logger.warning("체인 정보를 가져올 수 없습니다. 기본 체인을 사용합니다.")
        chains = [(1, 'A')]
    
    logger.info(f"사용 가능한 체인: {', '.join([f'{idx}({chain})' for idx, chain in chains])}")
    
    # 선택된 체인 처리
    chain_to_use = None
    chain_idx = None
    
    if selected_chain:
        # 사용자가 지정한 체인 ID 찾기
        for idx, chain_id in chains:
            if chain_id == selected_chain:
                chain_to_use = chain_id
                chain_idx = idx
                break
        
        if not chain_to_use:
            logger.warning(f"지정한 체인 ID '{selected_chain}'를 찾을 수 없습니다. 첫 번째 체인을 사용합니다.")
            chain_to_use = chains[0][1]
            chain_idx = chains[0][0]
    else:
        # 기본적으로 첫 번째 체인 사용
        chain_to_use = chains[0][1]
        chain_idx = chains[0][0]
    
    logger.info(f"선택된 체인: {chain_idx}({chain_to_use})")
    
    # 포스필드 리스트 - 순서대로 시도
    force_fields = [force_field, "amber99sb-ildn", "gromos54a7", "oplsaa"]
    
    success = False
    for ff in force_fields:
        try:
            logger.info(f"{ff} 포스필드로 시도합니다...")
            
            # pdb2gmx 실행: 선택된 체인만 사용하여 토폴로지 생성
            cmd_result = subprocess.run(
                f"echo -e '{chain_idx}\nq\n' | gmx pdb2gmx -f {base_name} -o complex.gro -p topol.top -water tip3p -ff {ff} -ignh",
                shell=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                check=False
            )
            
            # 명령어 실행 결과 로깅
            logger.info(f"{ff} 포스필드 실행 결과: {cmd_result.returncode}")
            if cmd_result.returncode != 0:
                logger.info(f"stdout: {cmd_result.stdout.decode()}")
                logger.info(f"stderr: {cmd_result.stderr.decode()}")
            
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
        raise Exception(f"시스템 준비 실패: 체인 {chain_to_use}에 대한 모든 포스필드 시도가 실패했습니다.")
    
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
    
    # grompp 실행: 이온화 준비
    run_command([
        "gmx", "grompp",
        "-f", "/app/mdp_templates/ions.mdp",
        "-c", "solv.gro",
        "-p", "topol.top",
        "-o", "ions.tpr",
        "-maxwarn", "10"
    ])
    
    # genion 실행: 이온 추가
    run_command([
        "echo", "SOL", "|", "gmx", "genion",
        "-s", "ions.tpr",
        "-o", "solv_ions.gro",
        "-p", "topol.top",
        "-pname", "NA",
        "-nname", "CL",
        "-neutral"
    ], shell=True)
    
    return "solv_ions.gro"

def run_energy_minimization(input_gro, output_prefix, output_dir):
    """에너지 최소화 실행"""
    os.chdir(output_dir)
    
    # grompp 실행: 에너지 최소화 준비
    run_command([
        "gmx", "grompp",
        "-f", "/app/mdp_templates/em.mdp",
        "-c", input_gro,
        "-p", "topol.top",
        "-o", f"{output_prefix}.tpr",
        "-maxwarn", "10"
    ])
    
    # mdrun 실행: 에너지 최소화
    run_command([
        "gmx", "mdrun",
        "-v",
        "-deffnm", output_prefix
    ])
    
    return f"{output_prefix}.gro"

def run_md_simulation(input_gro, output_prefix, output_dir, peptide_res, protein_res):
    """MD 시뮬레이션 실행"""
    os.chdir(output_dir)
    
    # grompp 실행: MD 준비
    run_command([
        "gmx", "grompp",
        "-f", "/app/mdp_templates/md.mdp",
        "-c", input_gro,
        "-p", "topol.top",
        "-o", f"{output_prefix}.tpr",
        "-maxwarn", "10"
    ])
    
    # mdrun 실행: MD 시뮬레이션
    run_command([
        "gmx", "mdrun",
        "-v",
        "-deffnm", output_prefix
    ])
    
    # 에너지 그룹 인덱스 파일 생성 (선택사항)
    create_index_file("index.ndx", peptide_res, protein_res)
    
    # 트라젝토리 변환 (선택사항)
    try:
        run_command([
            "echo", "Protein", "|", "gmx", "trjconv",
            "-s", f"{output_prefix}.tpr",
            "-f", f"{output_prefix}.xtc",
            "-o", "trajectory.pdb",
            "-pbc", "mol",
            "-center"
        ], shell=True)
    except Exception as e:
        logger.warning(f"트라젝토리 변환 중 오류 발생: {e}")
        logger.warning("트라젝토리 변환을 건너뜁니다.")
    
    return f"{output_prefix}.gro", f"{output_prefix}.xtc"

def main():
    parser = argparse.ArgumentParser(description="Gromacs를 이용한 SuMD 시뮬레이션 자동화")
    parser.add_argument("--pdb", required=True, help="입력 PDB 파일")
    parser.add_argument("--peptide_res", nargs="+", required=True, help="펩타이드 잔기 번호 리스트")
    parser.add_argument("--protein_res", nargs="+", required=True, help="단백질 잔기 번호 리스트")
    parser.add_argument("--output_dir", default="./output", help="출력 디렉토리")
    parser.add_argument("--distance_threshold", type=float, default=0.5, help="시뮬레이션 컷오프 거리 (nm, 기본값: 0.5)")
    parser.add_argument("--force_field", default="charmm36-jul2022", help="사용할 포스필드 (기본값: charmm36-jul2022)")
    parser.add_argument("--chain", help="사용할 체인 ID (예: A, B, C...)")
    
    args = parser.parse_args()
    
    # 출력 디렉토리 생성
    os.makedirs(args.output_dir, exist_ok=True)
    
    # 시스템 준비
    logger.info("시스템 준비 시작")
    try:
        solv_ions_gro = prepare_system(args.pdb, args.output_dir, args.force_field, args.chain)
    except Exception as e:
        logger.error(f"시스템 준비 중 오류 발생: {e}")
        return 1
    
    # 초기 무게중심 거리 계산
    initial_distance = calculate_distance(args.pdb, args.peptide_res, args.protein_res)
    logger.info(f"초기 무게중심 거리: {initial_distance:.3f} nm")
    
    # 에너지 최소화 실행
    logger.info("에너지 최소화 실행")
    try:
        em_gro = run_energy_minimization(solv_ions_gro, "em", args.output_dir)
    except Exception as e:
        logger.error(f"에너지 최소화 중 오류 발생: {e}")
        return 1
    
    # 에너지 최소화 후 거리 계산
    em_distance = calculate_distance(os.path.join(args.output_dir, em_gro), args.peptide_res, args.protein_res)
    logger.info(f"에너지 최소화 후 무게중심 거리: {em_distance:.3f} nm")
    
    # 거리에 따라 MD 시뮬레이션 실행 여부 결정
    if em_distance <= args.distance_threshold:
        logger.info("거리가 임계값보다 작음. MD 시뮬레이션 실행")
        try:
            md_gro, md_xtc = run_md_simulation(em_gro, "md", args.output_dir, args.peptide_res, args.protein_res)
        except Exception as e:
            logger.error(f"MD 시뮬레이션 중 오류 발생: {e}")
            return 1
        
        # MD 후 거리 계산
        md_distance = calculate_distance(os.path.join(args.output_dir, md_gro), args.peptide_res, args.protein_res)
        logger.info(f"MD 시뮬레이션 후 무게중심 거리: {md_distance:.3f} nm")
        
        # 결과 요약
        logger.info(f"초기 거리: {initial_distance:.3f} nm")
        logger.info(f"에너지 최소화 후 거리: {em_distance:.3f} nm")
        logger.info(f"MD 시뮬레이션 후 거리: {md_distance:.3f} nm")
        logger.info("시뮬레이션 완료")
    else:
        logger.info(f"거리가 임계값({args.distance_threshold} nm)보다 큼. MD 시뮬레이션 건너뜀")
        logger.info(f"초기 거리: {initial_distance:.3f} nm")
        logger.info(f"에너지 최소화 후 거리: {em_distance:.3f} nm")
        logger.info("에너지 최소화만 완료")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())
