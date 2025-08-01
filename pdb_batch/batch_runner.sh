#!/bin/bash
# run_pdb_pipeline.sh - CSV 기반 PDB 처리 파이프라인 실행 스크립트
# 파일 위치: /app/SUMD_automation/pdb_preprocessing/

set -e  # 오류 발생 시 스크립트 중단

# 기본 설정
CSV_FILE=$1
OUTPUT_DIR=${2:-"pdb_processing_results"}
MAX_ENTRIES=${3:-""}
CONFIG_FILE=${4:-""}
VALIDATE_ONLY=${5:-"false"}

# 스크립트 디렉토리 설정
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PARENT_DIR="$(dirname "$SCRIPT_DIR")"  # /app/SUMD_automation
ROOT_DIR="$(dirname "$PARENT_DIR")"    # /app

# GROMACS 설정 파일 기본값
if [ -z "$CONFIG_FILE" ]; then
    if [ -f "$PARENT_DIR/gromacs_commands_config.json" ]; then
        CONFIG_FILE="$PARENT_DIR/gromacs_commands_config.json"
    fi
fi

# 도움말 출력
show_help() {
    echo "CSV 기반 PDB 처리 파이프라인"
    echo ""
    echo "사용법: $0 <CSV파일> [출력디렉토리] [최대항목수] [설정파일] [검증만]"
    echo ""
    echo "파라미터:"
    echo "  CSV파일      : 처리할 CSV 파일 (필수)"
    echo "  출력디렉토리 : 결과 저장 디렉토리 (기본: pdb_processing_results)"
    echo "  최대항목수   : 처리할 최대 항목 수 (기본: 전체)"
    echo "  설정파일     : GROMACS 설정 JSON 파일 (기본: 자동 탐지)"
    echo "  검증만       : true시 검증만 수행하고 종료"
    echo ""
    echo "예시:"
    echo "  $0 data.csv results 5"
    echo "  $0 data.csv results \"\" ../gromacs_commands_config.json"
    echo "  $0 data.csv \"\" \"\" \"\" true  # 검증만"
    echo ""
    echo "CSV 형식 요구사항:"
    echo "  - PDB: PDB ID (예: 1a07)"
    echo "  - Peptide Chain: 리간드 체인 ID (예: C)"
    echo "  - Receptor Chain: 수용체 체인 ID (예: A)"
    echo ""
    echo "파일 구조:"
    echo "  $ROOT_DIR/"
    echo "  └── SUMD_automation/ (기존 SUMD 파일들)"
    echo "      └── pdb_preprocessing/ (PDB 전처리 파일들)"
    echo ""
    echo "처리 과정:"
    echo "  1. CSV 검증 → 2. PDB 다운로드 → 3. Chain 추출"
    echo "  4. PDB 전처리 → 5. GROMACS pdb2gmx 실행"
}

# 인자 확인
if [ -z "$CSV_FILE" ] || [ "$CSV_FILE" = "-h" ] || [ "$CSV_FILE" = "--help" ]; then
    show_help
    exit 1
fi

if [ ! -f "$CSV_FILE" ]; then
    echo "❌ CSV 파일을 찾을 수 없습니다: $CSV_FILE"
    exit 1
fi

# 출력 디렉토리 생성
mkdir -p "$OUTPUT_DIR"

# 로그 파일 설정
LOG_FILE="${OUTPUT_DIR}/pipeline.log"
TIMESTAMP=$(date +%Y%m%d_%H%M%S)

echo "=========================================="
echo "CSV 기반 PDB 처리 파이프라인"
echo "=========================================="
echo "스크립트 위치: $SCRIPT_DIR"
echo "상위 디렉토리: $PARENT_DIR"
echo "CSV 파일: $CSV_FILE"
echo "출력 디렉토리: $OUTPUT_DIR"
echo "최대 항목 수: ${MAX_ENTRIES:-'전체'}"
echo "설정 파일: ${CONFIG_FILE:-'기본값'}"
echo "검증만: $VALIDATE_ONLY"
echo "로그 파일: $LOG_FILE"
echo "시작 시간: $(date)"
echo "=========================================="

# Process substitution 대신 간단한 로깅 방식 사용
{

# 환경 확인
echo ""
echo "=== 환경 확인 ==="
echo "Python 버전: $(python3 --version 2>/dev/null || echo '❌ Python3 없음')"
echo "현재 작업 디렉토리: $(pwd)"

# GROMACS 확인
if command -v gmx >/dev/null 2>&1; then
    echo "GROMACS: ✅ $(gmx -version 2>/dev/null | head -1 || echo 'version unknown')"
elif command -v gmx_mpi >/dev/null 2>&1; then
    echo "GROMACS: ✅ MPI 버전"
else
    echo "GROMACS: ⚠️  설치되지 않음 (PDB 전처리만 가능)"
fi

# 필수 파일들 확인
echo ""
echo "=== 필수 파일 확인 ==="

REQUIRED_SCRIPTS=(
    "pdb_csv_pipeline.py"
    "csv_validator.py"
    "simple_chain_extractor.py"
)

MISSING_SCRIPTS=()

for script in "${REQUIRED_SCRIPTS[@]}"; do
    if [ -f "$SCRIPT_DIR/$script" ]; then
        echo "✅ $script"
    else
        echo "❌ $script (누락)"
        MISSING_SCRIPTS+=("$script")
    fi
done

if [ ${#MISSING_SCRIPTS[@]} -gt 0 ]; then
    echo ""
    echo "❌ 필수 스크립트가 누락되었습니다:"
    for script in "${MISSING_SCRIPTS[@]}"; do
        echo "  - $script"
    done
    echo ""
    echo "setup_pdb_pipeline.sh를 먼저 실행하여 설치를 완료하세요."
    exit 1
fi

# 상위 디렉토리의 기존 파일들 확인
SUMD_FILES=(
    "gromacs_runner_class.py"
    "mda_pdb_processor.py"
    "interface_distance_calculator.py"
)

echo ""
echo "=== 기존 SUMD 파일 확인 ==="
MISSING_SUMD=()

for file in "${SUMD_FILES[@]}"; do
    if [ -f "$PARENT_DIR/$file" ]; then
        echo "✅ $file"
    else
        echo "⚠️  $file (누락 - 일부 기능 제한됨)"
        MISSING_SUMD+=("$file")
    fi
done

# 1단계: CSV 검증
echo ""
echo "=== 1단계: CSV 파일 검증 ==="

CSV_VALIDATOR="$SCRIPT_DIR/csv_validator.py"
if [ -f "$CSV_VALIDATOR" ]; then
    echo "CSV 검증 중..."
    
    VALIDATION_REPORT="${OUTPUT_DIR}/validation_report.json"
    CLEANED_CSV="${OUTPUT_DIR}/cleaned_data.csv"
    
    # CSV 검증 실행
    cd "$SCRIPT_DIR"  # 작업 디렉토리를 스크립트 디렉토리로 변경
    
    if python3 "$CSV_VALIDATOR" \
        --input "$CSV_FILE" \
        --clean "$CLEANED_CSV" \
        --report "$VALIDATION_REPORT"; then
        
        echo "✅ CSV 검증 완료"
        
        # 정제된 파일이 있고 유효하면 사용
        if [ -f "$CLEANED_CSV" ] && [ -s "$CLEANED_CSV" ]; then
            echo "정제된 CSV 파일 사용: $CLEANED_CSV"
            CSV_TO_USE="$CLEANED_CSV"
        else
            echo "정제된 파일이 없어서 원본 사용: $CSV_FILE"
            CSV_TO_USE="$CSV_FILE"
        fi
    else
        echo "⚠️  CSV 검증 실패, 원본 파일로 계속 진행"
        CSV_TO_USE="$CSV_FILE"
    fi
else
    echo "⚠️  CSV 검증기 없음, 원본 파일 사용"
    CSV_TO_USE="$CSV_FILE"
fi

# 검증만 수행하는 경우 여기서 종료
if [ "$VALIDATE_ONLY" = "true" ]; then
    echo ""
    echo "검증만 수행하고 종료합니다."
    echo "검증 결과: $VALIDATION_REPORT"
    exit 0
fi

# CSV 파일 통계 출력
if command -v python3 >/dev/null 2>&1; then
    TOTAL_LINES=$(python3 -c "
import pandas as pd
try:
    df = pd.read_csv('$CSV_TO_USE')
    print(len(df))
except:
    print('Unknown')
" 2>/dev/null || echo "Unknown")
    echo "처리할 항목 수: $TOTAL_LINES"
fi

# 2단계: 파이프라인 실행
echo ""
echo "=== 2단계: PDB 처리 파이프라인 실행 ==="

PIPELINE_SCRIPT="$SCRIPT_DIR/pdb_csv_pipeline.py"

if [ ! -f "$PIPELINE_SCRIPT" ]; then
    echo "❌ 파이프라인 스크립트를 찾을 수 없습니다: $PIPELINE_SCRIPT"
    exit 1
fi

echo "파이프라인 스크립트: $PIPELINE_SCRIPT"

# 시스템 정보 출력
echo ""
echo "시스템 정보:"
echo "  CPU 코어 수: $(nproc)"
echo "  메모리: $(free -h | grep '^Mem:' | awk '{print $2}')"
echo "  디스크 여유공간: $(df -h . | tail -1 | awk '{print $4}')"
echo ""

# 파이프라인 실행 파라미터 준비
echo "PDB 처리 파이프라인 시작..."
START_TIME=$(date +%s)

# 최대 항목 수 파라미터
MAX_ENTRIES_PARAM=""
if [ -n "$MAX_ENTRIES" ]; then
    MAX_ENTRIES_PARAM="--max-entries $MAX_ENTRIES"
fi

# 설정 파일 파라미터
CONFIG_PARAM=""
if [ -n "$CONFIG_FILE" ] && [ -f "$CONFIG_FILE" ]; then
    CONFIG_PARAM="--config $CONFIG_FILE"
    echo "GROMACS 설정 파일 사용: $CONFIG_FILE"
else
    echo "⚠️  GROMACS 설정 파일 없음, 기본 설정 사용"
fi

# 출력 디렉토리를 절대 경로로 변환
OUTPUT_DIR_ABS=$(cd "$OUTPUT_DIR" && pwd)

# 파이프라인 실행 (스크립트 디렉토리에서)
cd "$SCRIPT_DIR"

if python3 "$PIPELINE_SCRIPT" \
    --csv "$CSV_TO_USE" \
    --output "$OUTPUT_DIR_ABS" \
    $MAX_ENTRIES_PARAM \
    $CONFIG_PARAM \
    --verbose; then
    
    PIPELINE_SUCCESS=true
else
    PIPELINE_SUCCESS=false
fi

END_TIME=$(date +%s)
DURATION=$((END_TIME - START_TIME))

# 3단계: 결과 요약
echo ""
echo "=== 3단계: 결과 요약 ==="

if [ "$PIPELINE_SUCCESS" = true ]; then
    echo "✅ 파이프라인 실행 완료"
    
    # 결과 파일들 확인
    SUMMARY_FILE="${OUTPUT_DIR_ABS}/processing_summary.json"
    
    if [ -f "$SUMMARY_FILE" ]; then
        echo "📋 처리 요약:"
        
        # jq가 있으면 JSON 파싱
        if command -v jq >/dev/null 2>&1; then
            TOTAL=$(jq -r '.total_entries // 0' "$SUMMARY_FILE")
            SUCCESSFUL=$(jq -r '.successful // 0' "$SUMMARY_FILE")
            FAILED=$(jq -r '.failed // 0' "$SUMMARY_FILE")
            
            echo "  총 항목: $TOTAL"
            echo "  성공: $SUCCESSFUL"
            echo "  실패: $FAILED"
            
            if [ "$TOTAL" -gt 0 ]; then
                SUCCESS_RATE=$(echo "scale=1; $SUCCESSFUL * 100 / $TOTAL" | bc -l 2>/dev/null || echo "0")
                echo "  성공률: ${SUCCESS_RATE}%"
            fi
            
            # 실패한 항목들 표시
            if [ "$FAILED" -gt 0 ]; then
                echo ""
                echo "실패한 항목들:"
                jq -r '.results[] | select(.success == false) | "  - \(.pdb_id): \(.error // "Unknown error")"' "$SUMMARY_FILE" 2>/dev/null | head -5
                
                if [ "$FAILED" -gt 5 ]; then
                    echo "  ... 및 $((FAILED - 5))개 더"
                fi
            fi
        else
            echo "  상세 정보는 $SUMMARY_FILE 참조"
        fi
    fi
    
    # 주요 출력 디렉토리들
    echo ""
    echo "📁 결과 디렉토리 구조:"
    echo "  $OUTPUT_DIR_ABS/"
    echo "  ├── downloaded_pdbs/      # 다운로드된 PDB 파일들"
    echo "  ├── extracted_chains/     # 추출된 receptor/ligand chains"
    echo "  ├── processed_pdbs/       # GROMACS 호환성 처리"
    echo "  ├── gromacs_results/      # pdb2gmx 실행 결과"
    echo "  │   ├── [pdb_id]/"
    echo "  │   │   ├── input.pdb"
    echo "  │   │   ├── processed.gro"
    echo "  │   │   └── topol.top"
    echo "  └── processing_summary.json  # 전체 처리 요약"
    
    # 성공한 결과가 있으면 예시 표시
    if [ -d "$OUTPUT_DIR_ABS/gromacs_results" ]; then
        FIRST_RESULT=$(ls "$OUTPUT_DIR_ABS/gromacs_results" 2>/dev/null | head -1)
        if [ -n "$FIRST_RESULT" ]; then
            echo ""
            echo "📄 GROMACS 결과 예시 ($FIRST_RESULT):"
            ls -la "$OUTPUT_DIR_ABS/gromacs_results/$FIRST_RESULT/" 2>/dev/null | head -5
        fi
    fi
    
else
    echo "❌ 파이프라인 실행 실패"
    echo ""
    echo "문제 해결 방법:"
    echo "  1. 로그 파일 확인: $LOG_FILE"
    echo "  2. --verbose 옵션으로 재실행"
    echo "  3. --max-entries 1로 단일 항목 테스트"
    echo "  4. CSV 파일 형식 재확인"
fi

echo ""
echo "실행 시간: ${DURATION}초 ($(date -d@$DURATION -u +%H:%M:%S))"
echo "완료 시간: $(date)"

# 최종 상태 코드
if [ "$PIPELINE_SUCCESS" = true ]; then
    echo ""
    echo "🎉 모든 작업이 완료되었습니다!"
    echo ""
    echo "다음 단계:"
    echo "  - 결과 확인: $OUTPUT_DIR_ABS"
    echo "  - GROMACS 시뮬레이션: processed.gro, topol.top 사용"
    echo "  - SuMD 실행: ../robust_sumd_master.py 활용"
    EXIT_CODE=0
else
    echo ""
    echo "⚠️  일부 작업이 실패했습니다."
    echo ""
    echo "문제 해결:"
    echo "  - 로그 확인: $LOG_FILE"
    echo "  - 인터넷 연결 확인 (PDB 다운로드용)"
    echo "  - GROMACS 설치 확인"
    echo "  - CSV 형식 확인"
    EXIT_CODE=1
fi

echo "=========================================="

exit $EXIT_CODE

} 2>&1 | tee -a "$LOG_FILE"
