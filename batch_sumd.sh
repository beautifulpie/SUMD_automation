#!/bin/bash

# Batch SuMD 실행 스크립트
# 사용법: ./batch_sumd.sh <입력_폴더> <출력_폴더> [simple_sumd.py_경로]

# 색상 정의
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# 로그 함수
log() {
    echo -e "[$(date '+%Y-%m-%d %H:%M:%S')] $1"
}

log_info() {
    log "${BLUE}[INFO]${NC} $1"
}

log_success() {
    log "${GREEN}[SUCCESS]${NC} $1"
}

log_warning() {
    log "${YELLOW}[WARNING]${NC} $1"
}

log_error() {
    log "${RED}[ERROR]${NC} $1"
}

# 사용법 출력
usage() {
    echo "사용법: $0 <입력_폴더> <출력_폴더> [simple_sumd.py_경로]"
    echo ""
    echo "매개변수:"
    echo "  입력_폴더        PDB 파일들이 있는 폴더"
    echo "  출력_폴더        결과를 저장할 폴더"
    echo "  simple_sumd.py_경로  simple_sumd.py 파일의 경로 (기본값: ./simple_sumd.py)"
    echo ""
    echo "PDB 파일 명명 규칙:"
    echo "  (pdb_코드)_(receptor_chain)_(ligand_chain)_processed.pdb"
    echo "  예: 1abi_H_L_processed.pdb"
    echo ""
    exit 1
}

# 매개변수 확인
if [ $# -lt 2 ]; then
    log_error "매개변수가 부족합니다."
    usage
fi

INPUT_DIR="$1"
OUTPUT_DIR="$2"
SUMD_SCRIPT="${3:-./simple_sumd.py}"

# 파일 및 폴더 존재 확인
if [ ! -d "$INPUT_DIR" ]; then
    log_error "입력 폴더가 존재하지 않습니다: $INPUT_DIR"
    exit 1
fi

if [ ! -f "$SUMD_SCRIPT" ]; then
    log_error "simple_sumd.py 파일을 찾을 수 없습니다: $SUMD_SCRIPT"
    exit 1
fi

# 출력 폴더 생성
if [ ! -d "$OUTPUT_DIR" ]; then
    log_info "출력 폴더 생성: $OUTPUT_DIR"
    mkdir -p "$OUTPUT_DIR"
fi

# PDB 파일명에서 정보 추출 함수
parse_pdb_filename() {
    local filename="$1"
    local basename=$(basename "$filename" .pdb)
    
    # _processed 제거
    basename=${basename%_processed}
    
    # 언더스코어로 분할
    IFS='_' read -ra PARTS <<< "$basename"
    
    if [ ${#PARTS[@]} -ne 3 ]; then
        log_error "잘못된 파일명 형식: $filename"
        log_error "올바른 형식: (pdb_코드)_(receptor_chain)_(ligand_chain)_processed.pdb"
        return 1
    fi
    
    PDB_CODE="${PARTS[0]}"
    RECEPTOR_CHAIN="${PARTS[1]}"
    LIGAND_CHAIN="${PARTS[2]}"
    
    return 0
}

# 실행 시작
log_info "=== Batch SuMD 실행 시작 ==="
log_info "입력 폴더: $INPUT_DIR"
log_info "출력 폴더: $OUTPUT_DIR"
log_info "SuMD 스크립트: $SUMD_SCRIPT"

# PDB 파일 찾기
PDB_FILES=($(find "$INPUT_DIR" -name "*_processed.pdb" -type f | sort))

if [ ${#PDB_FILES[@]} -eq 0 ]; then
    log_error "처리할 PDB 파일을 찾을 수 없습니다."
    log_error "파일명 형식: *_processed.pdb"
    exit 1
fi

log_info "발견된 PDB 파일: ${#PDB_FILES[@]}개"

# 전체 결과 요약을 위한 변수
TOTAL_FILES=${#PDB_FILES[@]}
SUCCESS_COUNT=0
FAILED_COUNT=0
FAILED_FILES=()

# 전체 실행 시작 시간
BATCH_START_TIME=$(date '+%Y-%m-%d %H:%M:%S')
BATCH_START_TIMESTAMP=$(date +%s)

# 전체 로그 파일
BATCH_LOG="$OUTPUT_DIR/batch_execution.log"
echo "=== Batch SuMD 실행 로그 ===" > "$BATCH_LOG"
echo "시작 시간: $BATCH_START_TIME" >> "$BATCH_LOG"
echo "총 파일 수: $TOTAL_FILES" >> "$BATCH_LOG"
echo "" >> "$BATCH_LOG"

# 각 PDB 파일 처리
for i in "${!PDB_FILES[@]}"; do
    PDB_FILE="${PDB_FILES[$i]}"
    FILE_NUM=$((i + 1))
    
    log_info "[$FILE_NUM/$TOTAL_FILES] 처리 중: $(basename "$PDB_FILE")"
    
    # 파일명에서 정보 추출
    if ! parse_pdb_filename "$PDB_FILE"; then
        log_error "파일명 파싱 실패: $PDB_FILE"
        FAILED_COUNT=$((FAILED_COUNT + 1))
        FAILED_FILES+=("$(basename "$PDB_FILE") - 파일명 파싱 실패")
        continue
    fi
    
    log_info "  PDB 코드: $PDB_CODE"
    log_info "  Receptor Chain: $RECEPTOR_CHAIN"
    log_info "  Ligand Chain: $LIGAND_CHAIN"
    
    # 개별 출력 폴더 생성
    JOB_OUTPUT_DIR="$OUTPUT_DIR/${PDB_CODE}_${RECEPTOR_CHAIN}_${LIGAND_CHAIN}"
    if [ -d "$JOB_OUTPUT_DIR" ]; then
        log_warning "기존 출력 폴더 발견, 삭제 후 재생성: $JOB_OUTPUT_DIR"
        rm -rf "$JOB_OUTPUT_DIR"
    fi
    mkdir -p "$JOB_OUTPUT_DIR"
    
    # 실행 시작 시간 기록
    JOB_START_TIME=$(date '+%Y-%m-%d %H:%M:%S')
    JOB_START_TIMESTAMP=$(date +%s)
    
    log_info "  SuMD 실행 시작: $JOB_START_TIME"
    
    # simple_sumd.py 실행
    log_info "  명령어: python3 $SUMD_SCRIPT \"$PDB_FILE\" $RECEPTOR_CHAIN $LIGAND_CHAIN \"$JOB_OUTPUT_DIR\""
    
    # 개별 작업 로그 파일
    JOB_LOG="$JOB_OUTPUT_DIR/execution.log"
    
    if python3 "$SUMD_SCRIPT" "$PDB_FILE" "$RECEPTOR_CHAIN" "$LIGAND_CHAIN" "$JOB_OUTPUT_DIR" > "$JOB_LOG" 2>&1; then
        # 실행 완료 시간 계산
        JOB_END_TIMESTAMP=$(date +%s)
        JOB_DURATION=$((JOB_END_TIMESTAMP - JOB_START_TIMESTAMP))
        JOB_END_TIME=$(date '+%Y-%m-%d %H:%M:%S')
        
        log_success "  SuMD 실행 완료: $JOB_END_TIME (소요시간: ${JOB_DURATION}초)"
        SUCCESS_COUNT=$((SUCCESS_COUNT + 1))
        
        # 결과 파일 확인
        if [ -f "$JOB_OUTPUT_DIR/final_results.json" ]; then
            log_success "  결과 파일 확인됨: final_results.json"
        else
            log_warning "  결과 파일이 생성되지 않았습니다"
        fi
        
        # 배치 로그에 기록
        echo "[$FILE_NUM/$TOTAL_FILES] SUCCESS: $(basename "$PDB_FILE") (${JOB_DURATION}초)" >> "$BATCH_LOG"
        
    else
        JOB_END_TIME=$(date '+%Y-%m-%d %H:%M:%S')
        log_error "  SuMD 실행 실패: $JOB_END_TIME"
        FAILED_COUNT=$((FAILED_COUNT + 1))
        FAILED_FILES+=("$(basename "$PDB_FILE") - 실행 오류")
        
        # 배치 로그에 기록
        echo "[$FILE_NUM/$TOTAL_FILES] FAILED: $(basename "$PDB_FILE")" >> "$BATCH_LOG"
    fi
    
    echo "" >> "$BATCH_LOG"
    
    # 진행률 표시
    PROGRESS=$((FILE_NUM * 100 / TOTAL_FILES))
    log_info "  진행률: $PROGRESS% ($FILE_NUM/$TOTAL_FILES)"
    echo ""
done

# 전체 실행 완료 시간 계산
BATCH_END_TIMESTAMP=$(date +%s)
BATCH_DURATION=$((BATCH_END_TIMESTAMP - BATCH_START_TIMESTAMP))
BATCH_END_TIME=$(date '+%Y-%m-%d %H:%M:%S')

# 최종 결과 요약
log_info "=== Batch SuMD 실행 완료 ==="
log_info "종료 시간: $BATCH_END_TIME"
log_info "총 소요시간: ${BATCH_DURATION}초 ($(($BATCH_DURATION / 60))분 $(($BATCH_DURATION % 60))초)"
log_info "총 파일 수: $TOTAL_FILES"
log_success "성공: $SUCCESS_COUNT"
if [ $FAILED_COUNT -gt 0 ]; then
    log_error "실패: $FAILED_COUNT"
    log_error "실패한 파일들:"
    for failed_file in "${FAILED_FILES[@]}"; do
        log_error "  - $failed_file"
    done
fi

# 성공률 계산
if [ $TOTAL_FILES -gt 0 ]; then
    SUCCESS_RATE=$((SUCCESS_COUNT * 100 / TOTAL_FILES))
    log_info "성공률: $SUCCESS_RATE%"
fi

# 배치 로그 파일에 최종 요약 추가
echo "" >> "$BATCH_LOG"
echo "=== 실행 완료 요약 ===" >> "$BATCH_LOG"
echo "종료 시간: $BATCH_END_TIME" >> "$BATCH_LOG"
echo "총 소요시간: ${BATCH_DURATION}초" >> "$BATCH_LOG"
echo "성공: $SUCCESS_COUNT / $TOTAL_FILES" >> "$BATCH_LOG"
echo "실패: $FAILED_COUNT / $TOTAL_FILES" >> "$BATCH_LOG"
if [ $TOTAL_FILES -gt 0 ]; then
    echo "성공률: $SUCCESS_RATE%" >> "$BATCH_LOG"
fi

if [ $FAILED_COUNT -gt 0 ]; then
    echo "" >> "$BATCH_LOG"
    echo "실패한 파일들:" >> "$BATCH_LOG"
    for failed_file in "${FAILED_FILES[@]}"; do
        echo "  - $failed_file" >> "$BATCH_LOG"
    done
fi

log_info "배치 실행 로그: $BATCH_LOG"

# 결과 요약 JSON 파일 생성
SUMMARY_JSON="$OUTPUT_DIR/batch_summary.json"
cat > "$SUMMARY_JSON" << EOF
{
  "batch_execution_summary": {
    "start_time": "$BATCH_START_TIME",
    "end_time": "$BATCH_END_TIME",
    "duration_seconds": $BATCH_DURATION,
    "total_files": $TOTAL_FILES,
    "successful_runs": $SUCCESS_COUNT,
    "failed_runs": $FAILED_COUNT,
    "success_rate_percent": $SUCCESS_RATE,
    "input_directory": "$INPUT_DIR",
    "output_directory": "$OUTPUT_DIR",
    "sumd_script": "$SUMD_SCRIPT"
  }
}
EOF

log_info "배치 요약 JSON: $SUMMARY_JSON"

if [ $FAILED_COUNT -eq 0 ]; then
    log_success "모든 작업이 성공적으로 완료되었습니다!"
    exit 0
else
    log_warning "일부 작업이 실패했습니다. 로그를 확인해주세요."
    exit 1
fi