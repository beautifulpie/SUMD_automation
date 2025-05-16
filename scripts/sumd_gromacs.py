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

def prepare_system(pdb_file, output_dir, force_field="charmm36-jul2022"):
    """시스템 준비 (토폴로지 생성, 박스 설정, 솔베이션, 이온화)"""
    os.makedirs(output_dir, exist_ok=True)
    
    # 파일 복사
    base_name = os.path.basename(pdb_file)
    pdb_path = os.path.join(output_dir, base_name)
    
    if pdb_file != pdb_path:
        run_command(f"cp {pdb_file} {pdb_path}", shell=True)
    
    # 작업 디렉토리 변경
    os.chdir(output_dir)
    
    # pdb2gmx 실행: 토폴로지 생성
    run_command([
        "gmx", "pdb2gmx", 
        "-f", base_name,
        "-o", "complex.gro",
        "-p", "topol.top",
        "-water", "tip3p",
        "-ff", force_field,
        "-ignh"
    ])
    
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
        "-o", "ions.tpr"
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
        "-o", f"{output_prefix}.tpr"
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
        "-o", f"{output_prefix}.tpr"
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
    run_command([
        "echo", "non-Water", "|", "gmx", "trjconv",
        "-s", f"{output_prefix}.tpr",
        "-f", f"{output_prefix}.xtc",
        "-o", "trajectory.pdb"
    ], shell=True)
    
    return f"{output_prefix}.gro", f"{output_prefix}.xtc"

def main():
    parser = argparse.ArgumentParser(description="Gromacs를 이용한 SuMD 시뮬레이션 자동화")
    parser.add_argument("--pdb", required=True, help="입력 PDB 파일")
    parser.add_argument("--peptide_res", nargs="+", required=True, help="펩타이드 잔기 번호 리스트")
    parser.add_argument("--protein_res", nargs="+", required=True, help="단백질 잔기 번호 리스트")
    parser.add_argument("--output_dir", default="./output", help="출력 디렉토리")
    parser.add_argument("--distance_threshold", type=float, default=0.5, help="시뮬레이션 컷오프 거리 (nm, 기본값: 0.5)")
    parser.add_argument("--force_field", default="charmm36-jul2022", help="사용할 포스필드 (기본값: charmm36-jul2022)")
    
    args = parser.parse_args()
    
    # 출력 디렉토리 생성
    os.makedirs(args.output_dir, exist_ok=True)
    
    # 시스템 준비
    logger.info("시스템 준비 시작")
    solv_ions_gro = prepare_system(args.pdb, args.output_dir, args.force_field)
    
    # 초기 무게중심 거리 계산
    initial_distance = calculate_distance(args.pdb, args.peptide_res, args.protein_res)
    logger.info(f"초기 무게중심 거리: {initial_distance:.3f} nm")
    
    # 에너지 최소화 실행
    logger.info("에너지 최소화 실행")
    em_gro = run_energy_minimization(solv_ions_gro, "em", args.output_dir)
    
    # 에너지 최소화 후 거리 계산
    em_distance = calculate_distance(os.path.join(args.output_dir, em_gro), args.peptide_res, args.protein_res)
    logger.info(f"에너지 최소화 후 무게중심 거리: {em_distance:.3f} nm")
    
    # 거리에 따라 MD 시뮬레이션 실행 여부 결정
    if em_distance <= args.distance_threshold:
        logger.info("거리가 임계값보다 작음. MD 시뮬레이션 실행")
        md_gro, md_xtc = run_md_simulation(em_gro, "md", args.output_dir, args.peptide_res, args.protein_res)
        
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

if __name__ == "__main__":
    main()
