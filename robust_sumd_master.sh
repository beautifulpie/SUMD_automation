#!/bin/bash
# golden_standard_sumd_master.sh - Golden Standard 기반 SuMD 시뮬레이션 마스터 실행 스크립트

set -e  # 오류 발생 시 스크립트 중단

# 기본 설정
INPUT_PDB=$1
GOLDEN_STANDARD_PDB=$2           # 새로 추가된 파라미터
RECEPTOR_CHAIN=$3
LIGAND_CHAIN=$4
SIMULATION_TIME=${5:-1.0}          # 시뮬레이션 시간 (ns, 기본값: 1.0)
DISTANCE_THRESHOLD=${6:-5.0}       # 동작 임계거리 (Å, 기본값: 5.0)
RMSD_THRESHOLD=${7:-1.5}           # Golden Standard와의 RMSD 임계값 (Å, 기본값: 1.5)
MAX_ITERATIONS=${8:-10}            # 최대 반복 횟수 (기본값: 10)
NUM_SAMPLES=${9:-5}                # 각 반복당 샘플 수 (기본값: 5)
OUTPUT_DIR=${10:-"/app/output"}    # 출력 디렉토리
JOB_ID=${11:-""}                   # 작업 ID (기본값: 자동생성)
SKIP_PREPROCESSING=${12:-"false"}  # 전처리 건너뛰기 (기본값: false)

# 인자 확인
if [ -z "$INPUT_PDB" ] || [ -z "$GOLDEN_STANDARD_PDB" ] || [ -z "$RECEPTOR_CHAIN" ] || [ -z "$LIGAND_CHAIN" ]; then
    echo "사용법: $0 <입력_PDB> <Golden_Standard_PDB> <수용체_체인> <리간드_체인> [시뮬레이션_시간] [거리_임계값] [RMSD_임계값] [최대_반복] [샘플수] [출력_디렉토리] [작업_ID] [전처리_건너뛰기]"
    echo "예시: $0 displaced_complex.pdb native_complex.pdb A B 2.0 5.0 1.5 10 5 /app/output \"\" false"
    echo ""
    echo "Golden Standard 기반 다중 샘플링 SuMD 시뮬레이션 설정:"
    echo "  - 입력 PDB: 변형된 구조 (chain이 이동/회전된 상태)"
    echo "  - Golden Standard PDB: 원래 네이티브 구조"
    echo "  - 각 반복마다 ${NUM_SAMPLES}개의 샘플을 순차 실행"
    echo "  - 거리 ≤ ${DISTANCE_THRESHOLD}Å: EM + MD (${SIMULATION_TIME}ns) 수행"
    echo "  - 거리 > ${DISTANCE_THRESHOLD}Å: EM만 수행"
    echo "  - 수렴 조건: Golden Standard와의 RMSD ≤ ${RMSD_THRESHOLD}Å (마지막 5회 평균)"
    echo "  - 최대 반복: ${MAX_ITERATIONS}회"
    echo ""
    echo "새로운 수렴 기준:"
    echo "  - Golden Standard PDB와 비교하여 절대적 수렴 판정"
    echo "  - 마지막 5회 iteration의 평균 RMSD가 ${RMSD_THRESHOLD}Å 이하면 수렴"
    echo "  - 네이티브 구조로의 복귀 정도를 정량적으로 측정"
    exit 1
fi

# 파일 존재 확인
if [ ! -f "$INPUT_PDB" ]; then
    echo "오류: 입력 PDB 파일 '$INPUT_PDB'를 찾을 수 없습니다."
    exit 1
fi

if [ ! -f "$GOLDEN_STANDARD_PDB" ]; then
    echo "오류: Golden Standard PDB 파일 '$GOLDEN_STANDARD_PDB'를 찾을 수 없습니다."
    exit 1
fi

# PDB 파일의 원자 수 확인 (시스템 크기 추정)
INPUT_ATOM_COUNT=$(grep "^ATOM" "$INPUT_PDB" | wc -l)
GOLDEN_ATOM_COUNT=$(grep "^ATOM" "$GOLDEN_STANDARD_PDB" | wc -l)

if [ "$INPUT_ATOM_COUNT" -gt 10000 ]; then
    SYSTEM_SIZE="large"
    EXPECTED_TIME="30-60분"
elif [ "$INPUT_ATOM_COUNT" -lt 3000 ]; then
    SYSTEM_SIZE="small"
    EXPECTED_TIME="5-15분"
else
    SYSTEM_SIZE="medium"
    EXPECTED_TIME="10-30분"
fi

echo "PDB 파일 분석:"
echo "  입력 PDB 원자 수: $INPUT_ATOM_COUNT"
echo "  Golden Standard PDB 원자 수: $GOLDEN_ATOM_COUNT"
echo "  시스템 크기: $SYSTEM_SIZE"
echo "  예상 실행 시간: $EXPECTED_TIME"

# 원자 수 차이 확인
ATOM_DIFF=$((INPUT_ATOM_COUNT - GOLDEN_ATOM_COUNT))
ATOM_DIFF_ABS=${ATOM_DIFF#-}  # 절댓값

if [ "$ATOM_DIFF_ABS" -gt 100 ]; then
    echo "⚠️  경고: 입력 PDB와 Golden Standard PDB의 원자 수 차이가 큽니다 ($ATOM_DIFF)"
    echo "  동일한 시스템인지 확인해주세요."
    echo ""
    read -p "계속 진행하시겠습니까? (y/N): " CONTINUE_DIFF
    if [ "$CONTINUE_DIFF" != "y" ] && [ "$CONTINUE_DIFF" != "Y" ]; then
        echo "시뮬레이션이 취소되었습니다."
        exit 1
    fi
fi

# PDB 파일의 체인 정보 미리 확인
echo "PDB 파일 체인 정보 확인 중..."
INPUT_CHAINS=$(grep "^ATOM" "$INPUT_PDB" | awk '{print $5}' | sort | uniq | tr '\n' ' ')
GOLDEN_CHAINS=$(grep "^ATOM" "$GOLDEN_STANDARD_PDB" | awk '{print $5}' | sort | uniq | tr '\n' ' ')

echo "입력 PDB 체인: $INPUT_CHAINS"
echo "Golden Standard PDB 체인: $GOLDEN_CHAINS"

# 체인 존재 여부 간단 확인
INPUT_CHAIN1_EXISTS=$(grep "^ATOM.*[[:space:]]${RECEPTOR_CHAIN}[[:space:]]" "$INPUT_PDB" | head -1)
INPUT_CHAIN2_EXISTS=$(grep "^ATOM.*[[:space:]]${LIGAND_CHAIN}[[:space:]]" "$INPUT_PDB" | head -1)
GOLDEN_CHAIN1_EXISTS=$(grep "^ATOM.*[[:space:]]${RECEPTOR_CHAIN}[[:space:]]" "$GOLDEN_STANDARD_PDB" | head -1)
GOLDEN_CHAIN2_EXISTS=$(grep "^ATOM.*[[:space:]]${LIGAND_CHAIN}[[:space:]]" "$GOLDEN_STANDARD_PDB" | head -1)

if [ -z "$INPUT_CHAIN1_EXISTS" ] || [ -z "$INPUT_CHAIN2_EXISTS" ]; then
    echo "경고: 입력 PDB에서 지정된 체인을 찾을 수 없습니다."
fi

if [ -z "$GOLDEN_CHAIN1_EXISTS" ] || [ -z "$GOLDEN_CHAIN2_EXISTS" ]; then
    echo "경고: Golden Standard PDB에서 지정된 체인을 찾을 수 없습니다."
fi

# 큰 시스템의 경우 사용자 확인
if [ "$SYSTEM_SIZE" = "large" ]; then
    echo ""
    echo "⚠️  큰 시스템 감지됨!"
    echo "  원자 수: $INPUT_ATOM_COUNT (>10,000)"
    echo "  예상 실행 시간: $EXPECTED_TIME"
    echo "  메모리 사용량이 높을 수 있습니다."
    echo ""
    read -p "계속 진행하시겠습니까? (y/N): " CONTINUE_LARGE
    if [ "$CONTINUE_LARGE" != "y" ] && [ "$CONTINUE_LARGE" != "Y" ]; then
        echo "시뮬레이션이 취소되었습니다."
        exit 1
    fi
fi

# 출력 디렉토리 생성
mkdir -p "$OUTPUT_DIR"

# 작업 ID가 없으면 자동 생성
if [ -z "$JOB_ID" ]; then
    TIMESTAMP=$(date +%Y%m%d_%H%M%S)
    JOB_ID="golden_${TIMESTAMP}_$$"
fi

# 최종 출력 디렉토리
FINAL_OUTPUT_DIR="${OUTPUT_DIR}/sumd_${JOB_ID}"
mkdir -p "$FINAL_OUTPUT_DIR"

# 로그 파일 설정
LOG_FILE="${FINAL_OUTPUT_DIR}/sumd_golden_master.log"

# Process substitution 대신 간단한 로깅 방식 사용
{

echo "=========================================="
echo "Golden Standard 기반 SuMD 다중 샘플링 시뮬레이션 시작"
echo "=========================================="
echo "작업 ID: $JOB_ID"
echo "입력 PDB: $INPUT_PDB"
echo "Golden Standard PDB: $GOLDEN_STANDARD_PDB"
echo "수용체 체인: $RECEPTOR_CHAIN"
echo "리간드 체인: $LIGAND_CHAIN"
echo "시뮬레이션 시간: $SIMULATION_TIME ns"
echo "거리 임계값: $DISTANCE_THRESHOLD Å"
echo "RMSD 임계값 (vs Golden): $RMSD_THRESHOLD Å"
echo "최대 반복 횟수: $MAX_ITERATIONS"
echo "각 반복당 샘플 수: $NUM_SAMPLES"
echo "출력 디렉토리: $FINAL_OUTPUT_DIR"
echo "시스템 크기: $SYSTEM_SIZE (원자수: $INPUT_ATOM_COUNT)"
echo "전처리 건너뛰기: $SKIP_PREPROCESSING"
echo "예상 실행 시간: $EXPECTED_TIME"
echo "=========================================="

# Python 스크립트 경로 확인
PYTHON_SCRIPT=""
if [ -f "/app/scripts/robust_sumd_master.py" ]; then
    PYTHON_SCRIPT="/app/scripts/robust_sumd_master.py"
else
    echo "오류: sumd_master 스크립트를 찾을 수 없습니다."
    exit 1
fi

echo "사용할 Python 스크립트: $PYTHON_SCRIPT"

# 시스템 정보 출력
echo ""
echo "시스템 정보:"
echo "  CPU 코어 수: $(nproc)"
echo "  메모리: $(free -h | grep '^Mem:' | awk '{print $2}')"
echo "  디스크 여유공간: $(df -h . | tail -1 | awk '{print $4}')"
echo ""

# Golden Standard SuMD 시뮬레이션 실행
echo "Golden Standard 기반 SuMD 마스터 스크립트 실행 중..."
START_TIME=$(date +%s)

# 전처리 건너뛰기 옵션 처리
SKIP_PREPROCESSING_FLAG=""
if [ "$SKIP_PREPROCESSING" = "true" ] || [ "$SKIP_PREPROCESSING" = "True" ] || [ "$SKIP_PREPROCESSING" = "1" ]; then
    SKIP_PREPROCESSING_FLAG="--skip_preprocessing"
fi

# Golden Standard 버전인지 확인하고 적절한 파라미터 사용
RESULT=$(python3 "$PYTHON_SCRIPT" \
    --input_pdb "$INPUT_PDB" \
    --golden_standard_pdb "$GOLDEN_STANDARD_PDB" \
    --output_dir "$OUTPUT_DIR" \
    --simulation_time "$SIMULATION_TIME" \
    --receptor_chain "$RECEPTOR_CHAIN" \
    --ligand_chain "$LIGAND_CHAIN" \
    --distance_threshold "$DISTANCE_THRESHOLD" \
    --rmsd_threshold "$RMSD_THRESHOLD" \
    --max_iterations "$MAX_ITERATIONS" \
    --num_samples "$NUM_SAMPLES" \
    --job_id "$JOB_ID" \
    $SKIP_PREPROCESSING_FLAG 2>&1 | grep "SUMD_RESULT:")


END_TIME=$(date +%s)
DURATION=$((END_TIME - START_TIME))

# 결과 파싱
if [[ $RESULT == SUMD_RESULT:* ]]; then
    # SUMD_RESULT:[파일 경로]:[수렴 여부] 형식에서 추출
    IFS=':' read -ra RESULT_PARTS <<< "$RESULT"
    
    if [ ${#RESULT_PARTS[@]} -ge 3 ]; then
        FINAL_PDB="${RESULT_PARTS[1]}"
        for i in $(seq 2 $(( ${#RESULT_PARTS[@]} - 2 ))); do
            FINAL_PDB="$FINAL_PDB:${RESULT_PARTS[$i]}"
        done
        CONVERGED="${RESULT_PARTS[${#RESULT_PARTS[@]}-1]}"
        
        echo ""
        echo "=========================================="
        echo "Golden Standard 기반 SuMD 시뮬레이션 완료"
        echo "=========================================="
        echo "최종 구조 파일: $FINAL_PDB"
        echo "수렴 여부: $CONVERGED"
        echo "실행 시간: ${DURATION}초 ($(date -d@$DURATION -u +%H:%M:%S))"
        echo "시스템 크기: $SYSTEM_SIZE"
        
        # 결과 요약 파일 생성
        SUMMARY_FILE="${FINAL_OUTPUT_DIR}/simulation_summary.txt"
        cat > "$SUMMARY_FILE" << EOF
=== Golden Standard 기반 SuMD 시뮬레이션 요약 ===
작업 ID: $JOB_ID
입력 PDB: $INPUT_PDB
Golden Standard PDB: $GOLDEN_STANDARD_PDB
수용체 체인: $RECEPTOR_CHAIN
리간드 체인: $LIGAND_CHAIN
시뮬레이션 시간: $SIMULATION_TIME ns
거리 임계값: $DISTANCE_THRESHOLD Å
RMSD 임계값 (vs Golden): $RMSD_THRESHOLD Å
최대 반복 횟수: $MAX_ITERATIONS
시스템 크기: $SYSTEM_SIZE (원자수: $INPUT_ATOM_COUNT)

=== 실행 결과 ===
최종 구조: $FINAL_PDB
수렴 여부: $CONVERGED
실행 시간: ${DURATION}초 ($(date -d@$DURATION -u +%H:%M:%S))

=== 파일 위치 ===
결과 디렉토리: $FINAL_OUTPUT_DIR
로그 파일: $LOG_FILE
결과 JSON: ${FINAL_OUTPUT_DIR}/results.json
EOF
        
        # 수렴 상태에 따른 메시지
        if [ "$CONVERGED" = "True" ] || [ "$CONVERGED" = "true" ]; then
            echo "상태: Golden Standard로 성공적으로 수렴됨"
            echo "상태: Golden Standard로 수렴 완료" >> "$SUMMARY_FILE"
            EXIT_CODE=0
        else
            echo "상태: 최대 반복 횟수 도달 (부분적 성공)"
            echo "상태: 최대 반복 도달" >> "$SUMMARY_FILE"
            EXIT_CODE=0  # 부분적 성공도 정상 종료로 처리
        fi
        
        # 추가 분석 실행 (선택적)
        if [ -f "${FINAL_OUTPUT_DIR}/results.json" ]; then
            echo ""
            echo "추가 분석 수행 중..."
            
            # 반복별 결과 분석 (jq 오류 방지)
            if command -v jq >/dev/null 2>&1; then
                # JSON 파일이 유효한지 먼저 확인
                if jq empty "${FINAL_OUTPUT_DIR}/results.json" 2>/dev/null; then
                    TOTAL_ITERATIONS=$(jq '.iterations | length' "${FINAL_OUTPUT_DIR}/results.json" 2>/dev/null || echo "0")
                    GOLDEN_STANDARD_FILE=$(jq -r '.golden_standard_pdb // "N/A"' "${FINAL_OUTPUT_DIR}/results.json" 2>/dev/null || echo "N/A")
                    echo "총 반복 횟수: $TOTAL_ITERATIONS"
                    echo "Golden Standard 파일: $GOLDEN_STANDARD_FILE"
                    echo "총 반복 횟수: $TOTAL_ITERATIONS" >> "$SUMMARY_FILE"
                    echo "Golden Standard PDB: $GOLDEN_STANDARD_FILE" >> "$SUMMARY_FILE"
                    
                    # 각 반복의 Golden Standard RMSD 변화 추출
                    echo "" >> "$SUMMARY_FILE"
                    echo "=== 반복별 Golden Standard RMSD 변화 ===" >> "$SUMMARY_FILE"
                    
                    jq -r '.iterations[]? | select(. != null) | 
                        "반복 \(.iteration // "?"):  거리 \(.initial_distance // "N/A")Å → \(.final_distance // "N/A")Å, Golden RMSD: \(.rmsd_vs_golden // "N/A")Å (\(.simulation_type // "Unknown"))"' \
                        "${FINAL_OUTPUT_DIR}/results.json" 2>/dev/null >> "$SUMMARY_FILE" || \
                        echo "반복별 상세 정보 추출 실패" >> "$SUMMARY_FILE"
                    
                else
                    echo "JSON 파일이 손상되었거나 유효하지 않습니다." >> "$SUMMARY_FILE"
                fi
            else
                echo "jq가 설치되지 않아 상세 분석을 건너뜁니다." >> "$SUMMARY_FILE"
            fi
        fi
        
    else
        echo "오류: 시뮬레이션 결과 형식이 잘못되었습니다: $RESULT"
        EXIT_CODE=1
    fi
else
    echo "오류: 시뮬레이션 결과를 찾을 수 없습니다."
    echo "Python 스크립트 직접 실행 결과:"
    
    python3 "$PYTHON_SCRIPT" \
        --input_pdb "$INPUT_PDB" \
        --golden_standard_pdb "$GOLDEN_STANDARD_PDB" \
        --output_dir "$OUTPUT_DIR" \
        --simulation_time "$SIMULATION_TIME" \
        --receptor_chain "$RECEPTOR_CHAIN" \
        --ligand_chain "$LIGAND_CHAIN" \
        --distance_threshold "$DISTANCE_THRESHOLD" \
        --rmsd_threshold "$RMSD_THRESHOLD" \
        --max_iterations "$MAX_ITERATIONS" \
        --num_samples "$NUM_SAMPLES" \
        --job_id "$JOB_ID" \
        $SKIP_PREPROCESSING_FLAG
    EXIT_CODE=1
fi

# 최종 결과 출력
echo ""
echo "=========================================="
echo "Golden Standard 기반 SuMD 다중 샘플링 시뮬레이션 완료"
echo "=========================================="
echo "작업 ID: $JOB_ID"
echo "실행 시간: ${DURATION}초 ($(date -d@$DURATION -u +%H:%M:%S))"
echo "시스템 크기: $SYSTEM_SIZE (원자수: $INPUT_ATOM_COUNT)"
echo "결과 디렉토리: $FINAL_OUTPUT_DIR"
echo "로그 파일: $LOG_FILE"
echo ""
echo "주요 출력 파일:"
echo "  - 결과 JSON: ${FINAL_OUTPUT_DIR}/results.json"
echo "  - 요약 파일: ${FINAL_OUTPUT_DIR}/simulation_summary.txt"
echo "  - 각 반복 결과: ${FINAL_OUTPUT_DIR}/iteration_*/sample_*/result_*.pdb"
echo ""
echo "Golden Standard 기반 수렴 알고리즘 특징:"
echo "  - Golden Standard PDB와의 절대적 RMSD 비교"
echo "  - 마지막 5회 iteration 평균 RMSD ≤ ${RMSD_THRESHOLD}Å로 수렴 판정"
echo "  - 네이티브 구조로의 복귀 정도를 정량적 측정"
echo "  - 체인 이동/회전된 구조의 복원 능력 평가"
echo "  - 각 반복마다 ${NUM_SAMPLES}개 샘플 순차 실행"
echo "  - 거리 기반 적응형 시뮬레이션 (EM/EM+MD)"
echo "=========================================="

exit $EXIT_CODE
} 2>&1 | tee -a "$LOG_FILE"