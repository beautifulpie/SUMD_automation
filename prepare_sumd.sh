#!/bin/bash
# SuMD 사전 처리 스크립트
# 이 스크립트는 SuMD 시뮬레이션을 위해 PDB 파일을 준비합니다.

# 기본 값
PDB_FILE=""
OUTPUT_DIR="./prepared_pdb"
PROTEIN_CHAIN="A"
PEPTIDE_CHAIN="B"
CLEAN_ONLY=false
VERBOSE=false

# 도움말 표시
show_help() {
    echo "사용법: $0 [옵션] -p PDB_FILE"
    echo "옵션:"
    echo "  -p, --pdb PDB_FILE       처리할 PDB 파일 (필수)"
    echo "  -o, --output OUTPUT_DIR  출력 디렉토리 (기본값: ./prepared_pdb)"
    echo "  -pc, --protein CHAIN     단백질 체인 ID (기본값: A)"
    echo "  -lc, --peptide CHAIN     펩티드 체인 ID (기본값: B)"
    echo "  -c, --clean-only         단백질만 추출 (기본값: false)"
    echo "  -v, --verbose            상세 출력 모드"
    echo "  -h, --help               도움말 표시"
    exit 1
}

# 파라미터 파싱
while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        -p|--pdb)
            PDB_FILE="$2"
            shift 2
            ;;
        -o|--output)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        -pc|--protein)
            PROTEIN_CHAIN="$2"
            shift 2
            ;;
        -lc|--peptide)
            PEPTIDE_CHAIN="$2"
            shift 2
            ;;
        -c|--clean-only)
            CLEAN_ONLY=true
            shift
            ;;
        -v|--verbose)
            VERBOSE=true
            shift
            ;;
        -h|--help)
            show_help
            ;;
        *)
            echo "알 수 없는 옵션: $1"
            show_help
            ;;
    esac
done

# PDB 파일 확인
if [ -z "$PDB_FILE" ]; then
    echo "오류: PDB 파일이 지정되지 않았습니다."
    show_help
fi

if [ ! -f "$PDB_FILE" ]; then
    echo "오류: PDB 파일 '$PDB_FILE'을 찾을 수 없습니다."
    exit 1
fi

# 출력 디렉토리 생성
mkdir -p "$OUTPUT_DIR"

# 파일명에서 확장자 제거
PDB_NAME=$(basename "$PDB_FILE" .pdb)

# 로그 함수
log() {
    if [ "$VERBOSE" = true ]; then
        echo "$1"
    fi
}

log "처리 중인 PDB 파일: $PDB_FILE"
log "출력 디렉토리: $OUTPUT_DIR"

# 단백질만 추출
if [ "$CLEAN_ONLY" = true ]; then
    OUTPUT_PDB="${OUTPUT_DIR}/${PDB_NAME}_protein_only.pdb"
    grep "^ATOM\|^TER\|^END\|^CRYST1\|^HEADER\|^TITLE" "$PDB_FILE" > "$OUTPUT_PDB"
    echo "단백질만 포함된 PDB 파일이 생성되었습니다: $OUTPUT_PDB"
    exit 0
fi

# 1. 체인 정보 추출
OUTPUT_CHAINS="${OUTPUT_DIR}/${PDB_NAME}_chains.txt"
log "체인 정보 추출 중..."

CHAINS=$(grep "^ATOM\|^HETATM" "$PDB_FILE" | awk '{print substr($0, 22, 1)}' | sort -u)
echo "발견된 체인: $CHAINS" > "$OUTPUT_CHAINS"

# 각 체인의 잔기 수 계산
for CHAIN in $CHAINS; do
    COUNT=$(grep "^ATOM\|^HETATM" "$PDB_FILE" | awk -v chain="$CHAIN" '{if (substr($0, 22, 1) == chain) print substr($0, 23, 4)}' | sort -u | wc -l)
    echo "체인 $CHAIN: $COUNT 잔기" >> "$OUTPUT_CHAINS"
done

log "체인 정보가 저장되었습니다: $OUTPUT_CHAINS"

# 2. 단백질과 펩티드 분리
PROTEIN_PDB="${OUTPUT_DIR}/${PDB_NAME}_protein.pdb"
PEPTIDE_PDB="${OUTPUT_DIR}/${PDB_NAME}_peptide.pdb"

log "단백질 체인($PROTEIN_CHAIN) 추출 중..."
grep "^ATOM\|^HETATM" "$PDB_FILE" | awk -v chain="$PROTEIN_CHAIN" '{if (substr($0, 22, 1) == chain) print $0}' > "$PROTEIN_PDB"
echo "TER" >> "$PROTEIN_PDB"
echo "END" >> "$PROTEIN_PDB"

log "펩티드 체인($PEPTIDE_CHAIN) 추출 중..."
grep "^ATOM\|^HETATM" "$PDB_FILE" | awk -v chain="$PEPTIDE_CHAIN" '{if (substr($0, 22, 1) == chain) print $0}' > "$PEPTIDE_PDB"
echo "TER" >> "$PEPTIDE_PDB"
echo "END" >> "$PEPTIDE_PDB"

# 3. 표준 잔기만 포함하는 버전 생성
CLEAN_PDB="${OUTPUT_DIR}/${PDB_NAME}_clean.pdb"
STANDARD_RESIDUES="ALA|ARG|ASN|ASP|CYS|GLN|GLU|GLY|HIS|ILE|LEU|LYS|MET|PHE|PRO|SER|THR|TRP|TYR|VAL|HOH|WAT|TIP"

log "표준 잔기만 포함하는 PDB 생성 중..."
grep "^ATOM\|^TER\|^END\|^CRYST1\|^HEADER\|^TITLE" "$PDB_FILE" > "$CLEAN_PDB"
grep "^HETATM" "$PDB_FILE" | grep -E "$STANDARD_RESIDUES" >> "$CLEAN_PDB"

# 4. pdb2gmx 테스트 실행
TEST_GRO="${OUTPUT_DIR}/${PDB_NAME}_test.gro"
TEST_TOP="${OUTPUT_DIR}/${PDB_NAME}_test.top"
TEST_LOG="${OUTPUT_DIR}/${PDB_NAME}_pdb2gmx_test.log"

log "GROMACS pdb2gmx 테스트 실행 중..."
gmx pdb2gmx -f "$CLEAN_PDB" -o "$TEST_GRO" -p "$TEST_TOP" -ff charmm36 -water tip3p -ignh -merge all > "$TEST_LOG" 2>&1
TEST_RESULT=$?

if [ $TEST_RESULT -eq 0 ]; then
    echo "pdb2gmx 테스트 성공!"
    echo "생성된 파일들:"
    echo "  단백질 체인: $PROTEIN_PDB"
    echo "  펩티드 체인: $PEPTIDE_PDB"
    echo "  정리된 PDB: $CLEAN_PDB"
    echo "  테스트 GRO: $TEST_GRO"
    echo "  테스트 TOP: $TEST_TOP"
else
    echo "pdb2gmx 테스트 실패."
    echo "로그 파일을 확인하세요: $TEST_LOG"
    
    # 설정 파일 생성
    CONFIG_FILE="${OUTPUT_DIR}/${PDB_NAME}_sumd_config.txt"
    echo "# SuMD 설정 파일" > "$CONFIG_FILE"
    echo "PDB_FILE=$CLEAN_PDB" >> "$CONFIG_FILE"
    echo "PROTEIN_CHAIN=$PROTEIN_CHAIN" >> "$CONFIG_FILE"
    echo "PEPTIDE_CHAIN=$PEPTIDE_CHAIN" >> "$CONFIG_FILE"
    echo "FORCE_FIELD=charmm36" >> "$CONFIG_FILE"
    echo "WATER_MODEL=tip3p" >> "$CONFIG_FILE"
    
    echo "SuMD 설정 파일이 생성되었습니다: $CONFIG_FILE"
fi

echo "전처리가 완료되었습니다."
