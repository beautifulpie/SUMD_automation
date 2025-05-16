#!/bin/bash
# test_environment.sh - SuMD Gromacs 테스트 환경 설정 및 실행 스크립트

set -e  # 오류 발생 시 스크립트 중단

# 필요한 디렉토리 생성
mkdir -p scripts mdp_templates test_data output

# MDP 템플릿 파일 생성
echo "MDP 템플릿 파일 생성 중..."
cat > mdp_templates/ions.mdp << EOF
; LINES STARTING WITH ';' ARE COMMENTS
title               = Ion insertion    ; Title of run

; Parameters describing what to do, when to stop and what to save
integrator          = steep        ; Algorithm (steep = steepest descent minimization)
emtol               = 1000.0       ; Stop minimization when the maximum force < 10.0 kJ/mol
emstep              = 0.01         ; Energy step size
nsteps              = 50000        ; Maximum number of (minimization) steps to perform

; Parameters describing how to find the neighbors of each atom and how to calculate the interactions
nstlist             = 1            ; Frequency to update the neighbor list and long range forces
cutoff-scheme       = Verlet
ns_type             = grid         ; Method to determine neighbor list (simple, grid)
rlist               = 1.0          ; Cut-off for making neighbor list (short range forces)
coulombtype         = cutoff       ; Treatment of long range electrostatic interactions
rcoulomb            = 1.0          ; long range electrostatic cut-off
rvdw                = 1.0          ; long range Van der Waals cut-off
pbc                 = xyz          ; Periodic Boundary Conditions
EOF

cat > mdp_templates/em.mdp << EOF
; 에너지 최소화 설정
integrator          = steep
nsteps              = 5000     ; 적은 스텝으로 빠른 최소화 진행
emtol               = 1000.0   ; 큰 값으로 설정하여 빠른 수렴
emstep              = 0.01
nstlist             = 1
cutoff-scheme       = Verlet
ns_type             = grid
coulombtype         = PME
rcoulomb            = 1.0
rvdw                = 1.0
pbc                 = xyz
EOF

cat > mdp_templates/md.mdp << EOF
; 분자 시뮬레이션 설정
integrator              = md
dt                      = 0.004     ; 4 fs (더 큰 타임스텝)
nsteps                  = 2500000   ; 10 ns 시뮬레이션
nstenergy               = 5000      ; 20 ps마다 에너지 저장
nstlog                  = 5000      ; 20 ps마다 로그 파일 업데이트
nstxout-compressed      = 5000      ; 20 ps마다 압축된 좌표 저장

; 온도 커플링
tcoupl                  = V-rescale
tc-grps                 = Protein Water_and_ions
tau_t                   = 0.1 0.1
ref_t                   = 300 300

; 압력 커플링
pcoupl                  = Parrinello-Rahman
pcoupltype              = isotropic
tau_p                   = 2.0
ref_p                   = 1.0
compressibility         = 4.5e-5

; 결합 제약
constraints             = h-bonds
constraint_algorithm    = LINCS
lincs_iter              = 2
lincs_order             = 6

; 비결합 상호작용
cutoff-scheme           = Verlet
nstlist                 = 20
ns_type                 = grid
coulombtype             = PME
rcoulomb                = 1.0
rvdw                    = 1.0
DispCorr                = EnerPres
EOF

# 스크립트 파일 생성
echo "스크립트 파일 생성 중..."
cat > scripts/process_cystein.py << EOF
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import os
import sys
import numpy as np
from Bio.PDB import PDBParser, PDBIO, Atom

def process_cystein(input_pdb, output_pdb):
    """
    PDB 파일의 시스테인(CYS) 잔기에 수소 원자(HG)를 추가합니다.
    
    Args:
        input_pdb: 입력 PDB 파일 경로
        output_pdb: 출력 PDB 파일 경로
    
    Returns:
        처리된 CYS 잔기 수
    """
    print(f"시스테인 잔기 처리 중: {input_pdb} -> {output_pdb}")
    
    # PDB 파일 파싱
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("structure", input_pdb)
    
    # 처리된 CYS 잔기 수 카운트
    processed_count = 0
    
    # 모델, 체인, 잔기를 순회하면서 CYS 잔기 찾기
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.get_resname() == "CYS":
                    # 이미 HG가 존재하면, 건너뜀
                    if "HG" not in residue:
                        # SG와 CA 원자가 있는지 확인
                        if "SG" in residue and "CA" in residue:
                            sg_atom = residue["SG"]
                            ca_atom = residue["CA"]
                            
                            # SG와 CA 사이의 벡터 계산 및 단위 벡터 산출
                            vector = sg_atom.coord - ca_atom.coord
                            norm = np.linalg.norm(vector)
                            if norm == 0:
                                continue  # 0 나누기 방지
                            unit_vector = vector / norm
                            
                            # S-H 결합 길이 (Å): 일반적으로 약 1.33 Å
                            bond_length = 1.33
                            hg_coord = sg_atom.coord + unit_vector * bond_length
                            
                            # 새로운 HG 원자 생성
                            hg_atom = Atom(name="HG",
                                       coord=hg_coord,
                                       bfactor=sg_atom.get_bfactor(),
                                       occupancy=sg_atom.get_occupancy(),
                                       altloc=sg_atom.get_altloc(),
                                       fullname=" HG ",
                                       serial_number=0,
                                       element="H")
                            
                            # HG 원자를 잔기에 추가
                            residue.add(hg_atom)
                            processed_count += 1
                        else:
                            print(f"경고: 잔기 {residue.get_id()}에 SG 또는 CA 원자가 없어 HG를 추가하지 않습니다.")
    
    # 변경된 구조를 새로운 PDB 파일로 저장
    io = PDBIO()
    io.set_structure(structure)
    io.save(output_pdb)
    
    print(f"처리 완료: {processed_count}개의 시스테인 잔기에 HG 원자 추가")
    return processed_count

def main():
    parser = argparse.ArgumentParser(description="PDB 파일의 시스테인 잔기에 HG 원자 추가")
    parser.add_argument("--input", "-i", required=True, help="입력 PDB 파일 경로")
    parser.add_argument("--output", "-o", required=True, help="출력 PDB 파일 경로")
    
    args = parser.parse_args()
    
    # 입력 파일 확인
    if not os.path.exists(args.input):
        print(f"오류: 입력 파일 {args.input}을 찾을 수 없습니다.")
        return 1
    
    # 출력 디렉토리 확인
    output_dir = os.path.dirname(args.output)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)
    
    # 시스테인 잔기 처리
    processed_count = process_cystein(args.input, args.output)
    
    print(f"총 {processed_count}개의 시스테인 잔기가 처리되었습니다.")
    return 0

if __name__ == "__main__":
    sys.exit(main())
EOF

# sumd_gromacs.py와 energy_analysis.py 파일은 이전과 동일하게 생성...

# 실행 스크립트 생성
cat > scripts/run_sumd_gromacs.sh << EOF
#!/bin/bash
# run_sumd_gromacs.sh - SuMD Gromacs 직접 실행 스크립트

set -e  # 오류 발생 시 스크립트 중단

# 기본 설정
PDB_FILE=\$1
PEPTIDE_RES=\$2
PROTEIN_RES=\$3
OUTPUT_DIR=\${4:-"./output"}
DISTANCE_THRESHOLD=\${5:-0.5}
PROCESS_CYS=\${6:-"true"}  # 시스테인 처리 여부 (기본값: true)

# 인자 확인
if [ -z "\$PDB_FILE" ] || [ -z "\$PEPTIDE_RES" ] || [ -z "\$PROTEIN_RES" ]; then
    echo "사용법: \$0 <PDB_파일> <펩타이드_잔기_번호> <단백질_잔기_번호> [출력_디렉토리] [거리_임계값] [시스테인_처리]"
    echo "예시: \$0 complex.pdb \"45 46 47\" \"123 124 125\" ./results 0.5 true"
    exit 1
fi

# 필요한 디렉토리 생성
mkdir -p "\$OUTPUT_DIR"
mkdir -p "/app/mdp_templates"

# MDP 파일이 없는 경우 기본 템플릿 생성
# (MDP 파일 생성 코드는 동일하게 유지...)

echo "SuMD Gromacs 시뮬레이션 시작"
echo "PDB 파일: \$PDB_FILE"
echo "펩타이드 잔기: \$PEPTIDE_RES"
echo "단백질 잔기: \$PROTEIN_RES"
echo "출력 디렉토리: \$OUTPUT_DIR"
echo "거리 임계값: \$DISTANCE_THRESHOLD"
echo "시스테인 처리: \$PROCESS_CYS"

# 시스테인 잔기 처리
PROCESSED_PDB="\$PDB_FILE"
if [ "\$PROCESS_CYS" = "true" ]; then
    echo "시스테인 잔기 처리 중..."
    PROCESSED_PDB="\$OUTPUT_DIR/processed_cys.pdb"
    python3 /app/scripts/process_cystein.py --input "\$PDB_FILE" --output "\$PROCESSED_PDB"
fi

# SuMD Gromacs 스크립트 실행
python3 /app/scripts/sumd_gromacs.py \\
    --pdb "\$PROCESSED_PDB" \\
    --peptide_res \$PEPTIDE_RES \\
    --protein_res \$PROTEIN_RES \\
    --output_dir "\$OUTPUT_DIR" \\
    --distance_threshold "\$DISTANCE_THRESHOLD"

# MD 시뮬레이션 결과가 있는지 확인하고 에너지 분석 실행
if [ -f "\$OUTPUT_DIR/md.xtc" ]; then
    echo "MD 시뮬레이션 결과 발견. 에너지 분석 시작..."
    
    python3 /app/scripts/energy_analysis.py \\
        --output_dir "\$OUTPUT_DIR" \\
        --peptide_res \$PEPTIDE_RES \\
        --protein_res \$PROTEIN_RES
        
    echo "에너지 분석 완료"
    echo "결과 파일: \$OUTPUT_DIR/interaction_energy.png"
    echo "요약 파일: \$OUTPUT_DIR/energy_summary.txt"
else
    echo "MD 시뮬레이션 결과가 없습니다. 에너지 분석을 건너뜁니다."
    echo "무게중심 거리가 임계값 \$DISTANCE_THRESHOLD nm보다 큽니다."
    echo "에너지 최소화 결과는 \$OUTPUT_DIR/em.gro 파일에서 확인할 수 있습니다."
fi

echo "모든 작업이 완료되었습니다."
EOF

# Dockerfile 생성
echo "Dockerfile 생성 중..."
cat > Dockerfile << EOF
FROM kcaladram/gromacs-new:latest

# 기본 패키지 설치
RUN apt-get update && apt-get upgrade -y && \\
    apt-get install -y nano python3-pip git wget && \\
    pip3 install numpy biopython pandas matplotlib tqdm

# 작업 디렉토리 생성
WORKDIR /app

# 필요한 디렉토리 생성
RUN mkdir -p /app/mdp_templates /app/scripts /app/test_data /app/output

# MDP 템플릿 복사
COPY mdp_templates/ /app/mdp_templates/

# 스크립트 복사
COPY scripts/ /app/scripts/

# 실행 권한 부여
RUN chmod +x /app/scripts/*.py /app/scripts/*.sh

# 환경 변수 설정
ENV PATH="/app/scripts:\${PATH}"

# 기본 명령어 설정
CMD ["/bin/bash"]
EOF

# docker-compose.yml 생성
echo "docker-compose.yml 생성 중..."
cat > docker-compose.yml << EOF
version: '3.8'

services:
  sumd-gromacs:
    build:
      context: .
      dockerfile: Dockerfile
    image: sumd-gromacs:latest
    container_name: sumd-gromacs
    volumes:
      - ./scripts:/app/scripts
      - ./mdp_templates:/app/mdp_templates
      - ./test_data:/app/test_data
      - ./output:/app/output
    deploy:
      resources:
        reservations:
          devices:
            - driver: nvidia
              count: all
              capabilities: [gpu]
    environment:
      - NVIDIA_VISIBLE_DEVICES=all
      - NVIDIA_DRIVER_CAPABILITIES=compute,utility,graphics
    # 대화형 사용을 위해 컨테이너 유지
    tty: true
    stdin_open: true
EOF

# 실행 권한 부여
chmod +x scripts/*.py scripts/*.sh

# 테스트 데이터 다운로드
echo "테스트 데이터 다운로드 중..."
mkdir -p test_data
wget -O test_data/example.pdb "https://files.rcsb.org/download/1BRS.pdb" || echo "테스트 PDB 파일 다운로드 실패. 직접 PDB 파일을 test_data 디렉토리에 추가하세요."

# Docker 이미지 빌드 및 실행
echo "Docker 이미지 빌드 및 실행 중..."
docker-compose up -d --build

echo "테스트 환경 구축 완료!"
echo "컨테이너에 접속하려면 다음 명령어를 사용하세요:"
echo "docker exec -it sumd-gromacs bash"
echo ""
echo "시스테인 처리가 포함된 테스트 실행 예시:"
echo "docker exec -it sumd-gromacs /app/scripts/run_sumd_gromacs.sh /app/test_data/example.pdb \"45 46 47\" \"123 124 125\" /app/output 0.5 true"
