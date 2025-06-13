#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import os
import subprocess
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import logging
import sys

# 로깅 설정
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("energy_analysis.log"),
        logging.StreamHandler(sys.stdout)
    ]
)
logger = logging.getLogger("Energy-Analysis")

def run_command(cmd, shell=False, cwd=None, check=True, input_str=None):
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
            stderr=subprocess.PIPE
        )
        if result.returncode == 0:
            logger.info("명령어 성공적으로 실행됨")
        else:
            logger.error(f"명령어 실행 실패: {result.stderr.decode()}")
        return result
    except Exception as e:
        logger.error(f"명령어 실행 중 오류 발생: {e}")
        raise

def create_index_file(output_dir, output_file, peptide_res_list, protein_res_list):
    """인덱스 파일 생성 (make_ndx 대체)"""
    with open(os.path.join(output_dir, output_file), 'w') as f:
        f.write("[ Peptide ]\n")
        f.write(" ".join(peptide_res_list) + "\n")
        f.write("[ Protein ]\n")
        f.write(" ".join(protein_res_list) + "\n")
        f.write("[ Peptide_Protein ]\n")
        f.write(" ".join(peptide_res_list) + " " + " ".join(protein_res_list) + "\n")

def calculate_interaction_energy(output_dir, tpr_file, trajectory_file):
    """상호작용 에너지 계산"""
    os.chdir(output_dir)
    
    # 에너지 재계산을 위한 tpr 파일 준비
    run_command([
        "gmx", "grompp",
        "-f", "/app/mdp_templates/md.mdp",
        "-c", "md.gro",
        "-p", "topol.top",
        "-n", "index.ndx",
        "-o", "energy_calc.tpr",
        "-maxwarn", "10"
    ])
    
    # 트라젝토리에서 에너지 재계산
    run_command([
        "gmx", "mdrun",
        "-rerun", trajectory_file,
        "-s", "energy_calc.tpr",
        "-deffnm", "energy_calc"
    ])
    
    # Coulomb 에너지 추출
    coulomb_result = run_command([
        "echo", "Coulomb-(SR):Peptide-Protein", "|", "gmx", "energy",
        "-f", "energy_calc.edr",
        "-s", "energy_calc.tpr",
        "-o", "coulomb.xvg"
    ], shell=True)
    
    # LJ 에너지 추출
    lj_result = run_command([
        "echo", "LJ-(SR):Peptide-Protein", "|", "gmx", "energy",
        "-f", "energy_calc.edr",
        "-s", "energy_calc.tpr",
        "-o", "lj.xvg"
    ], shell=True)
    
    return "coulomb.xvg", "lj.xvg"

def parse_xvg(file_path):
    """XVG 파일 파싱"""
    data = []
    with open(file_path, 'r') as f:
        for line in f:
            if not line.startswith(("#", "@")):
                values = line.strip().split()
                if len(values) >= 2:
                    data.append([float(values[0]), float(values[1])])
    return np.array(data)

def analyze_energy(output_dir, coulomb_file, lj_file):
    """에너지 분석 및 플롯 생성"""
    # XVG 파일 파싱
    coulomb_data = parse_xvg(os.path.join(output_dir, coulomb_file))
    lj_data = parse_xvg(os.path.join(output_dir, lj_file))
    
    # 데이터 프레임 생성
    time = coulomb_data[:, 0]
    coulomb = coulomb_data[:, 1]
    lj = lj_data[:, 1]
    total = coulomb + lj
    
    df = pd.DataFrame({
        'Time (ps)': time,
        'Coulomb (kJ/mol)': coulomb,
        'LJ (kJ/mol)': lj,
        'Total (kJ/mol)': total
    })
    
    # 통계 계산
    coulomb_avg = np.mean(coulomb)
    lj_avg = np.mean(lj)
    total_avg = np.mean(total)
    
    # 결과 출력
    logger.info(f"평균 Coulomb 에너지: {coulomb_avg:.2f} kJ/mol")
    logger.info(f"평균 LJ 에너지: {lj_avg:.2f} kJ/mol")
    logger.info(f"평균 총 상호작용 에너지: {total_avg:.2f} kJ/mol")
    
    # 플롯 생성
    plt.figure(figsize=(12, 6))
    plt.plot(time, coulomb, label='Coulomb')
    plt.plot(time, lj, label='Lennard-Jones')
    plt.plot(time, total, label='Total')
    plt.xlabel('Time (ps)')
    plt.ylabel('Energy (kJ/mol)')
    plt.title('Protein-Peptide Interaction Energy')
    plt.legend()
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.savefig(os.path.join(output_dir, 'interaction_energy.png'), dpi=300)
    plt.close()
    
    # CSV 저장
    df.to_csv(os.path.join(output_dir, 'interaction_energy.csv'), index=False)
    
    # 요약 파일 생성
    with open(os.path.join(output_dir, 'energy_summary.txt'), 'w') as f:
        f.write(f"평균 Coulomb 에너지: {coulomb_avg:.2f} kJ/mol\n")
        f.write(f"평균 LJ 에너지: {lj_avg:.2f} kJ/mol\n")
        f.write(f"평균 총 상호작용 에너지: {total_avg:.2f} kJ/mol\n")
        
        # 상호작용 강도 평가
        if total_avg < -50:
            strength = "강한 결합"
        elif total_avg < -20:
            strength = "중간 결합"
        else:
            strength = "약한 결합"
        
        f.write(f"결합 강도 평가: {strength} ({total_avg:.2f} kJ/mol)\n")
    
    return df

def main():
    parser = argparse.ArgumentParser(description="Gromacs MD 시뮬레이션 결과의 상호작용 에너지 분석")
    parser.add_argument("--output_dir", required=True, help="시뮬레이션 결과 디렉토리")
    parser.add_argument("--peptide_res", nargs="+", required=True, help="펩타이드 잔기 번호 리스트")
    parser.add_argument("--protein_res", nargs="+", required=True, help="단백질 잔기 번호 리스트")
    
    args = parser.parse_args()
    
    # 인덱스 파일 생성
    create_index_file(args.output_dir, "index.ndx", args.peptide_res, args.protein_res)
    
    # 상호작용 에너지 계산
    coulomb_file, lj_file = calculate_interaction_energy(
        args.output_dir, 
        "md.tpr", 
        "md.xtc"
    )
    
    # 에너지 분석
    energy_df = analyze_energy(args.output_dir, coulomb_file, lj_file)
    
    logger.info("에너지 분석 완료")
    logger.info(f"결과 파일: {os.path.join(args.output_dir, 'interaction_energy.png')}")
    logger.info(f"결과 데이터: {os.path.join(args.output_dir, 'interaction_energy.csv')}")
    logger.info(f"요약 파일: {os.path.join(args.output_dir, 'energy_summary.txt')}")

if __name__ == "__main__":
    main()
