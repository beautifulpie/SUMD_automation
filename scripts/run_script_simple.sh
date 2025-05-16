#!/bin/bash
# run_sumd_gromacs.sh - SuMD Gromacs 직접 실행 스크립트

set -e  # 오류 발생 시 스크립트 중단

# 기본 설정
PDB_FILE=$1
PEPTIDE_RES=$2
PROTEIN_RES=$3
OUTPUT_DIR=${4:-"./output"}
DISTANCE_THRESHOLD=${5:-0.5}
PROCESS_CYS=${6:-"true"}  # 시스테인 처리 여부 (기본값: true)

# 인자 확인
if [ -z "$PDB_FILE" ] || [ -z "$PEPTIDE_RES" ] || [ -z "$PROTEIN_RES" ]; then
    echo "사용법: $0 <PDB_파일> <펩타이드_잔기_번호> <단백질_잔기_번호> [출력_디렉토리] [거리_임계값] [시스테인_처리]"
    echo "예시: $0 complex.pdb \"45 46 47\" \"123 124 125\" ./results 0.5 true"
    exit 1
fi

# 필요한 디렉토리 생성
mkdir -p "$OUTPUT_DIR"
mkdir -p "/app/mdp_templates"

# MDP 파일이 없는 경우 기본 템플릿 생성
if [ ! -f "/app/mdp_templates/ions.mdp" ]; then
    cat > "/app/mdp_templates/ions.mdp" << EOF
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
fi

if [ ! -f "/app/mdp_templates/em.mdp" ]; then
    cat > "/app/mdp_templates/em.mdp" << EOF
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
fi

if [ ! -f "/app/mdp_templates/md.mdp" ]; then
    cat > "/app/mdp_templates/md.mdp" << EOF
; 분자 시뮬레이션 설정
integrator              = md
; dt=0.004로 증가하고 nsteps=2500000으로 감소하여 10ns 시뮬레이션을 더 빠르게 실행
; (0.004 fs * 2500000 = 10000 ps = 10 ns)
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
lincs_iter              = 2        ; 더 큰 타임스텝을 위해 증가
lincs_order             = 6        ; 더 큰 타임스텝을 위해 증가

; 비결합 상호작용
cutoff-scheme           = Verlet
nstlist                 = 20
ns_type                 = grid
coulombtype             = PME
rcoulomb                = 1.0
rvdw                    = 1.0
DispCorr                = EnerPres
EOF
fi

echo "SuMD Gromacs 시뮬레이션 시작"
echo "PDB 파일: $PDB_FILE"
echo "펩타이드 잔기: $PEPTIDE_RES"
echo "단백질 잔기: $PROTEIN_RES"
echo "출력 디렉토리: $OUTPUT_DIR"
echo "거리 임계값: $DISTANCE_THRESHOLD"
echo "시스테인 처리: $PROCESS_CYS"

# 시스테인 잔기 처리
PROCESSED_PDB="$PDB_FILE"
if [ "$PROCESS_CYS" = "true" ]; then
    echo "시스테인 잔기 처리 중..."
    PROCESSED_PDB="$OUTPUT_DIR/processed_cys.pdb"
    python3 /app/scripts/process_cystein.py --input "$PDB_FILE" --output "$PROCESSED_PDB"
fi

# SuMD Gromacs 스크립트 실행
python3 /app/scripts/sumd_gromacs.py \
    --pdb "$PROCESSED_PDB" \
    --peptide_res $PEPTIDE_RES \
    --protein_res $PROTEIN_RES \
    --output_dir "$OUTPUT_DIR" \
    --distance_threshold "$DISTANCE_THRESHOLD"

# MD 시뮬레이션 결과가 있는지 확인하고 에너지 분석 실행
if [ -f "$OUTPUT_DIR/md.xtc" ]; then
    echo "MD 시뮬레이션 결과 발견. 에너지 분석 시작..."
    
    python3 /app/scripts/energy_analysis.py \
        --output_dir "$OUTPUT_DIR" \
        --peptide_res $PEPTIDE_RES \
        --protein_res $PROTEIN_RES
        
    echo "에너지 분석 완료"
    echo "결과 파일: $OUTPUT_DIR/interaction_energy.png"
    echo "요약 파일: $OUTPUT_DIR/energy_summary.txt"
else
    echo "MD 시뮬레이션 결과가 없습니다. 에너지 분석을 건너뜁니다."
    echo "무게중심 거리가 임계값 $DISTANCE_THRESHOLD nm보다 큽니다."
    echo "에너지 최소화 결과는 $OUTPUT_DIR/em.gro 파일에서 확인할 수 있습니다."
fi

echo "모든 작업이 완료되었습니다."
