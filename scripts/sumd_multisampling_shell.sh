#!/bin/bash
# run_sumd_gromacs_multisampling_process.sh - 프로세스 기반 다중 샘플 SuMD Gromacs 시뮬레이션 실행 스크립트

set -e  # 오류 발생 시 스크립트 중단

# 기본 설정
PDB_FILE=$1
CHAIN1=$2
CHAIN2=$3
MAX_ITERATIONS=${4:-10}        # 최대 반복 횟수 (기본값: 10)
DISTANCE_THRESHOLD=${5:-0.5}   # 거리 임계값 (nm, 기본값: 0.5)
CLEAN_PDB=${6:-"true"}         # PDB 파일 정리 여부 (기본값: true)
PROCESS_CYS=${7:-"true"}       # 시스테인 처리 여부 (기본값: true)
NUM_SAMPLES=${8:-5}            # 각 반복당 샘플 수 (기본값: 5)
SIMULATION_TIME=${9:-2}        # 각 샘플의 시뮬레이션 시간 (ns, 기본값: 2)
SIMULATION_THRESHOLD=${10:-5.0} # MD 실행 거리 임계값 (Å, 기본값: 5.0)
MAX_WORKERS=${11:-0}           # 최대 워커 프로세스 수 (0: 자동)
OUTPUT_BASE_DIR=${12:-"./sumd_multisampling_output"}  # 출력 기본 디렉토리

# 인자 확인
if [ -z "$PDB_FILE" ] || [ -z "$CHAIN1" ] || [ -z "$CHAIN2" ]; then
    echo "사용법: $0 <PDB_파일> <체인1_ID> <체인2_ID> [최대반복횟수] [거리_임계값] [PDB_정리] [시스테인_처리] [샘플수] [시뮬레이션_시간] [MD_거리_임계값] [최대워커수] [출력_디렉토리]"
    echo "예시: $0 complex.pdb A B 10 0.5 true true 5 2 5.0 4 ./results"
    echo ""
    echo "프로세스 기반 다중 샘플 SuMD 설정:"
    echo "  - 각 반복마다 ${NUM_SAMPLES}개의 샘플을 병렬 프로세스로 실행"
    echo "  - 초기 거리 ≤ ${SIMULATION_THRESHOLD}Å: EM + MD (${SIMULATION_TIME}ns) 수행"
    echo "  - 초기 거리 > ${SIMULATION_THRESHOLD}Å: EM만 수행"
    echo "  - 최대 워커 프로세스: ${MAX_WORKERS} (0=자동)"
    echo "  - 가장 좋은 결과(거리가 가장 짧은)를 다음 반복의 시작점으로 사용"
    echo "  - 수렴 조건: 거리 ≤ ${DISTANCE_THRESHOLD} nm"
    exit 1
fi

# 시스템 정보 확인
CPU_COUNT=$(nproc)
if [ "$MAX_WORKERS" -eq 0 ]; then
    EFFECTIVE_WORKERS=$CPU_COUNT
else
    EFFECTIVE_WORKERS=$(($MAX_WORKERS < $CPU_COUNT ? $MAX_WORKERS : $CPU_COUNT))
fi

# 타임스탬프
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
OUTPUT_DIR="${OUTPUT_BASE_DIR}_${TIMESTAMP}"
mkdir -p "$OUTPUT_DIR"

# 로그 파일 설정
LOG_FILE="${OUTPUT_DIR}/sumd_multisampling_process_run.log"
exec > >(tee -a "$LOG_FILE") 2>&1

echo "==================================================="
echo "프로세스 기반 다중 샘플 SuMD Gromacs 시뮬레이션 시작"
echo "==================================================="
echo "PDB 파일: $PDB_FILE"
echo "체인1 ID: $CHAIN1"
echo "체인2 ID: $CHAIN2"
echo "최대 반복 횟수: $MAX_ITERATIONS"
echo "거리 임계값: $DISTANCE_THRESHOLD nm"
echo "PDB 파일 정리: $CLEAN_PDB"
echo "시스테인 처리: $PROCESS_CYS"
echo "각 반복당 샘플 수: $NUM_SAMPLES"
echo "각 샘플 시뮬레이션 시간: $SIMULATION_TIME ns"
echo "MD 실행 거리 임계값: $SIMULATION_THRESHOLD Å"
echo "수렴 거리 임계값: $DISTANCE_THRESHOLD nm"
echo "시스템 CPU 개수: $CPU_COUNT"
echo "효과적 워커 프로세스: $EFFECTIVE_WORKERS"
echo "출력 디렉토리: $OUTPUT_DIR"
echo "==================================================="

# 1. PDB 파일 정리 (표준 아미노산만 추출)
CURRENT_PDB="$PDB_FILE"
if [ "$CLEAN_PDB" = "true" ]; then
    echo "PDB 파일에서 표준 아미노산만 추출 중..."
    CLEANED_PDB="${OUTPUT_DIR}/cleaned_$(basename "$PDB_FILE")"
    
    if [ -f "/app/scripts/clean_pdb.py" ]; then
        python3 /app/scripts/clean_pdb.py --input "$PDB_FILE" --output "$CLEANED_PDB"
    elif [ -f "./scripts/clean_pdb.py" ]; then
        python3 ./scripts/clean_pdb.py --input "$PDB_FILE" --output "$CLEANED_PDB"
    elif [ -f "./clean_pdb.py" ]; then
        python3 ./clean_pdb.py --input "$PDB_FILE" --output "$CLEANED_PDB"
    else
        echo "경고: clean_pdb.py 스크립트를 찾을 수 없습니다. 원본 PDB 파일을 계속 사용합니다."
    fi
    
    # 파일이 성공적으로 생성되었는지 확인
    if [ -f "$CLEANED_PDB" ]; then
        echo "PDB 파일 정리 완료: $CLEANED_PDB"
        CURRENT_PDB="$CLEANED_PDB"
    else
        echo "경고: PDB 파일 정리 실패, 원본 파일을 계속 사용합니다."
    fi
fi

# 2. 시스테인 잔기 처리
if [ "$PROCESS_CYS" = "true" ]; then
    echo "시스테인 잔기 처리 중..."
    PROCESSED_PDB="${OUTPUT_DIR}/processed_cys_$(basename "$CURRENT_PDB")"
    
    if [ -f "/app/scripts/process_cystein.py" ]; then
        python3 /app/scripts/process_cystein.py --input "$CURRENT_PDB" --output "$PROCESSED_PDB"
    elif [ -f "./scripts/process_cystein.py" ]; then
        python3 ./scripts/process_cystein.py --input "$CURRENT_PDB" --output "$PROCESSED_PDB"
    elif [ -f "./process_cystein.py" ]; then
        python3 ./process_cystein.py --input "$CURRENT_PDB" --output "$PROCESSED_PDB"
    else
        echo "경고: process_cystein.py 스크립트를 찾을 수 없습니다. 이전 파일을 계속 사용합니다."
    fi
    
    # 파일이 성공적으로 생성되었는지 확인
    if [ -f "$PROCESSED_PDB" ]; then
        echo "시스테인 잔기 처리 완료: $PROCESSED_PDB"
        CURRENT_PDB="$PROCESSED_PDB"
    else
        echo "경고: 시스테인 잔기 처리 실패, 이전 파일을 계속 사용합니다."
    fi
fi

# 3. 프로세스 기반 다중 샘플 시뮬레이션 스크립트 경로 확인
MULTISAMPLING_SCRIPT=""
if [ -f "/app/scripts/sumd_multisampling.py" ]; then
    MULTISAMPLING_SCRIPT="/app/scripts/sumd_multisampling.py"
elif [ -f "./scripts/sumd_multisampling.py" ]; then
    MULTISAMPLING_SCRIPT="./scripts/sumd_multisampling.py"
elif [ -f "./sumd_multisampling.py" ]; then
    MULTISAMPLING_SCRIPT="./sumd_multisampling.py"
else
    echo "오류: 다중 샘플링 스크립트를 찾을 수 없습니다."
    exit 1
fi

# 4. 프로세스 기반 다중 샘플 SuMD 시뮬레이션 실행
echo "프로세스 기반 다중 샘플 SuMD 시뮬레이션 시작..."
echo "각 반복마다 ${NUM_SAMPLES}개의 샘플을 ${EFFECTIVE_WORKERS}개 프로세스로 병렬 실행합니다."

# 시뮬레이션 결과 요약을 저장할 파일
SUMMARY_FILE="${OUTPUT_DIR}/sumd_multisampling_process_summary.txt"
echo "=== 프로세스 기반 다중 샘플 SuMD 시뮬레이션 요약 ===" > "$SUMMARY_FILE"
echo "시작 구조: $CURRENT_PDB" >> "$SUMMARY_FILE"
echo "체인1: $CHAIN1, 체인2: $CHAIN2" >> "$SUMMARY_FILE"
echo "거리 임계값: $DISTANCE_THRESHOLD nm" >> "$SUMMARY_FILE"
echo "각 반복당 샘플 수: $NUM_SAMPLES" >> "$SUMMARY_FILE"
echo "각 샘플 시뮬레이션 시간: $SIMULATION_TIME ns" >> "$SUMMARY_FILE"
echo "최대 반복 횟수: $MAX_ITERATIONS" >> "$SUMMARY_FILE"
echo "효과적 워커 프로세스: $EFFECTIVE_WORKERS" >> "$SUMMARY_FILE"
echo "----------------------------------------" >> "$SUMMARY_FILE"

# Python 스크립트 실행 (워커 수 옵션 추가)
echo "프로세스 기반 다중 샘플 SuMD 파이썬 스크립트 실행 중..."

# 워커 수 인자 구성
WORKER_ARG=""
if [ "$MAX_WORKERS" -ne 0 ]; then
    WORKER_ARG="--max_workers $MAX_WORKERS"
fi

RESULT=$(python3 "$MULTISAMPLING_SCRIPT" \
    --input "$CURRENT_PDB" \
    --chain1 "$CHAIN1" \
    --chain2 "$CHAIN2" \
    --output_dir "$OUTPUT_DIR" \
    --distance_threshold "$DISTANCE_THRESHOLD" \
    --max_iterations "$MAX_ITERATIONS" \
    --num_samples "$NUM_SAMPLES" \
    --simulation_time "$SIMULATION_TIME" \
    --simulation_threshold "$SIMULATION_THRESHOLD" \
    $WORKER_ARG 2>&1 | grep "SUMD_RESULT:")

# 결과 파싱
if [[ $RESULT == SUMD_RESULT:* ]]; then
    # SUMD_RESULT:[파일 경로]:[거리] 형식에서 추출
    IFS=':' read -ra RESULT_PARTS <<< "$RESULT"
    
    if [ ${#RESULT_PARTS[@]} -ge 3 ]; then
        FINAL_FILE="${RESULT_PARTS[1]}"
        for i in $(seq 2 $(( ${#RESULT_PARTS[@]} - 2 ))); do
            FINAL_FILE="$FINAL_FILE:${RESULT_PARTS[$i]}"
        done
        FINAL_DISTANCE="${RESULT_PARTS[${#RESULT_PARTS[@]}-1]}"
        
        echo "프로세스 기반 다중 샘플 SuMD 시뮬레이션 완료"
        echo "최종 구조 파일: $FINAL_FILE"
        echo "최종 거리: $FINAL_DISTANCE nm"
        
        # 요약 파일에 최종 결과 추가
        echo "" >> "$SUMMARY_FILE"
        echo "=== 최종 결과 ===" >> "$SUMMARY_FILE"
        echo "최종 구조: $FINAL_FILE" >> "$SUMMARY_FILE"
        echo "최종 거리: $FINAL_DISTANCE nm" >> "$SUMMARY_FILE"
        
        # 수렴 여부 판단
        if (( $(echo "$FINAL_DISTANCE <= $DISTANCE_THRESHOLD" | bc -l) )); then
            echo "상태: 수렴 완료 (거리: $FINAL_DISTANCE nm <= 임계값: $DISTANCE_THRESHOLD nm)"
            echo "상태: 수렴 완료" >> "$SUMMARY_FILE"
        else
            echo "상태: 최대 반복 횟수 도달 (거리: $FINAL_DISTANCE nm > 임계값: $DISTANCE_THRESHOLD nm)"
            echo "상태: 최대 반복 횟수 도달 (수렴하지 않음)" >> "$SUMMARY_FILE"
        fi
        
        # 에너지 분석 실행 (선택적)
        if [ -f "/app/scripts/energy_analysis.py" ] || [ -f "./scripts/energy_analysis.py" ] || [ -f "./energy_analysis.py" ]; then
            echo ""
            echo "에너지 분석 실행 중..."
            
            # 최종 구조가 있는 디렉토리에서 에너지 분석 수행
            FINAL_DIR=$(dirname "$FINAL_FILE")
            
            if [ -f "/app/scripts/energy_analysis.py" ]; then
                ENERGY_SCRIPT="/app/scripts/energy_analysis.py"
            elif [ -f "./scripts/energy_analysis.py" ]; then
                ENERGY_SCRIPT="./scripts/energy_analysis.py"
            else
                ENERGY_SCRIPT="./energy_analysis.py"
            fi
            
            # 에너지 분석 실행 (체인 정보를 잔기 번호로 변환 필요)
            # 여기서는 간단히 체인 전체를 분석한다고 가정
            python3 "$ENERGY_SCRIPT" \
                --output_dir "$FINAL_DIR" \
                --peptide_res "1 2 3 4 5" \
                --protein_res "6 7 8 9 10" 2>/dev/null || echo "에너지 분석 건너뜀 (선택적 기능)"
        fi
        
    else
        echo "오류: 시뮬레이션 결과 형식이 잘못되었습니다: $RESULT"
        echo "  오류: 부적절한 결과 형식" >> "$SUMMARY_FILE"
    fi
else
    echo "오류: 시뮬레이션 결과를 찾을 수 없습니다."
    echo "스크립트 출력:"
    python3 "$MULTISAMPLING_SCRIPT" \
        --input "$CURRENT_PDB" \
        --chain1 "$CHAIN1" \
        --chain2 "$CHAIN2" \
        --output_dir "$OUTPUT_DIR" \
        --distance_threshold "$DISTANCE_THRESHOLD" \
        --max_iterations "$MAX_ITERATIONS" \
        --num_samples "$NUM_SAMPLES" \
        --simulation_time "$SIMULATION_TIME" \
        --simulation_threshold "$SIMULATION_THRESHOLD" \
        $WORKER_ARG
    echo "  오류: 결과 없음" >> "$SUMMARY_FILE"
fi

# 최종 결과 요약
echo ""
echo "==================================================="
echo "프로세스 기반 다중 샘플 SuMD 시뮬레이션 완료"
echo "==================================================="
echo "시뮬레이션 방식: 프로세스 기반 거리 적응형 다중 샘플"
echo "  - 초기 거리 ≤ ${SIMULATION_THRESHOLD}Å: EM + MD (${SIMULATION_TIME}ns) 수행"
echo "  - 초기 거리 > ${SIMULATION_THRESHOLD}Å: EM만 수행"
echo "각 반복당 샘플 수: ${NUM_SAMPLES}개"
echo "사용된 워커 프로세스: ${EFFECTIVE_WORKERS}개"
echo "요약 파일: $SUMMARY_FILE"
echo "로그 파일: $LOG_FILE"
echo "결과 디렉토리: $OUTPUT_DIR"
echo ""
echo "프로세스 기반 병렬 처리의 장점:"
echo "  - 각 샘플이 독립적인 프로세스에서 실행"
echo "  - 메모리 충돌 없는 안정적인 병렬 처리"
echo "  - GROMACS의 OpenMP 충돌 방지"
echo "  - 시스템 리소스 효율적 활용"
echo "  - 프로세스별 격리된 작업 환경"
echo "==================================================="

echo "프로세스 기반 다중 샘플 SuMD 시뮬레이션이 완료되었습니다."