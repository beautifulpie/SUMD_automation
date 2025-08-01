#!/bin/bash
# setup_pdb_pipeline.sh - CSV 기반 PDB 처리 파이프라인 설치 스크립트

set -e

echo "=========================================="
echo "CSV 기반 PDB 처리 파이프라인 설치"
echo "=========================================="

# 현재 디렉토리 확인 (/app/SUMD_automation/pdb_preprocessing)
CURRENT_DIR=$(pwd)
echo "현재 디렉토리: $CURRENT_DIR"

# 상위 디렉토리 (/app/SUMD_automation)
PARENT_DIR=$(dirname "$CURRENT_DIR")
echo "상위 디렉토리: $PARENT_DIR"

# 프로젝트 루트 디렉토리 (/app)
ROOT_DIR=$(dirname "$PARENT_DIR")
echo "프로젝트 루트: $ROOT_DIR"

echo ""
echo "파일 구조 확인:"
echo "  📁 $ROOT_DIR/"
echo "  📁 $PARENT_DIR/ (기존 파일들)"
echo "  📁 $CURRENT_DIR/ (새 파일들)"

echo ""
echo "필요한 파일들 확인 중..."

# 상위 디렉토리의 기존 파일들 확인
EXISTING_FILES=(
    "interface_distance_calculator.py"
    "gromacs_runner_class.py"
    "mda_pdb_processor.py"
    "robust_sumd_master.py"
    "MDPGenerator.py"
    "utils.py"
)

FOUND_EXISTING=()
MISSING_EXISTING=()

echo ""
echo "=== 기존 파일들 (상위 디렉토리) ==="
for file in "${EXISTING_FILES[@]}"; do
    if [ -f "$PARENT_DIR/$file" ]; then
        FOUND_EXISTING+=("$file")
        echo "✅ $file"
    else
        MISSING_EXISTING+=("$file")
        echo "❌ $file (누락)"
    fi
done

# 현재 디렉토리의 새 파일들 확인
NEW_FILES=(
    "pdb_csv_pipeline.py"
    "csv_validator.py"
    "simple_chain_extractor.py"
)

echo ""
echo "=== 새 파일들 (현재 디렉토리) ==="
FOUND_NEW=()
MISSING_NEW=()

for file in "${NEW_FILES[@]}"; do
    if [ -f "$CURRENT_DIR/$file" ]; then
        FOUND_NEW+=("$file")
        echo "✅ $file"
    else
        MISSING_NEW+=("$file")
        echo "❌ $file (생성 필요)"
    fi
done

# GROMACS 설정 파일 확인
GROMACS_CONFIG="$PARENT_DIR/gromacs_commands_config.json"
if [ -f "$GROMACS_CONFIG" ]; then
    echo "✅ GROMACS 설정 파일: gromacs_commands_config.json"
else
    echo "⚠️  GROMACS 설정 파일 없음: gromacs_commands_config.json"
fi

echo ""
echo "파일 상태 요약:"
echo "  기존 파일: ${#FOUND_EXISTING[@]}/${#EXISTING_FILES[@]} 발견"
echo "  새 파일: ${#FOUND_NEW[@]}/${#NEW_FILES[@]} 발견"

if [ ${#MISSING_EXISTING[@]} -gt 0 ]; then
    echo ""
    echo "⚠️  누락된 기존 파일들:"
    for file in "${MISSING_EXISTING[@]}"; do
        echo "  - $file"
    done
    echo ""
    echo "일부 기능이 제한될 수 있습니다."
fi

if [ ${#MISSING_NEW[@]} -gt 0 ]; then
    echo ""
    echo "❌ 누락된 새 파일들:"
    for file in "${MISSING_NEW[@]}"; do
        echo "  - $file"
    done
    echo ""
    echo "이 파일들을 먼저 생성해야 합니다."
    echo "계속 진행하려면 Enter를 누르세요..."
    read -r
fi

# Python 의존성 확인
echo ""
echo "=== Python 의존성 확인 ==="

PYTHON_MODULES=(
    "pandas"
    "requests" 
    "Bio"
    "numpy"
    "pathlib"
    "json"
    "logging"
)

MISSING_MODULES=()

for module in "${PYTHON_MODULES[@]}"; do
    if python3 -c "import $module" 2>/dev/null; then
        echo "✅ $module"
    else
        echo "❌ $module (누락)"
        MISSING_MODULES+=("$module")
    fi
done

if [ ${#MISSING_MODULES[@]} -gt 0 ]; then
    echo ""
    echo "누락된 모듈 설치 명령:"
    echo "pip3 install ${MISSING_MODULES[*]}"
    echo ""
fi

# 선택적 의존성 확인
echo ""
echo "=== 선택적 의존성 확인 ==="

OPTIONAL_MODULES=(
    "MDAnalysis"
    "scipy"
    "matplotlib"
)

for module in "${OPTIONAL_MODULES[@]}"; do
    if python3 -c "import $module" 2>/dev/null; then
        echo "✅ $module (선택적)"
    else
        echo "⚠️  $module (선택적, 없어도 동작)"
    fi
done

# GROMACS 확인
echo ""
echo "=== GROMACS 확인 ==="

if command -v gmx >/dev/null 2>&1; then
    echo "✅ GROMACS 설치됨"
    GMX_VERSION=$(gmx -version 2>/dev/null | head -1 || echo "버전 확인 실패")
    echo "  $GMX_VERSION"
else
    echo "❌ GROMACS가 설치되지 않음"
    echo "  PDB 전처리만 가능하고 pdb2gmx는 실행되지 않습니다."
fi

if command -v gmx_mpi >/dev/null 2>&1; then
    echo "✅ GROMACS MPI 버전 사용 가능"
else
    echo "⚠️  GROMACS MPI 버전 없음 (단일 프로세서로 실행)"
fi

# 테스트 CSV 파일 생성
echo ""
echo "=== 테스트 데이터 생성 ==="

TEST_CSV="$CURRENT_DIR/test_data.csv"
cat > "$TEST_CSV" << 'EOF'
PDB,Protein Name,Resolution,Classification,Peptide Chain,Peptide Size,Peptide Sequence,Peptide Description,Peptide Organism,Peptide Interface Area,Peptide Molecular Weight,Peptide Aromaticity,Peptide Instability,Peptide Isoelectric Point,Receptor Chain,Receptor Size,Receptor Sequence,Receptor Description,Receptor Organism
148l,A COVALENT ENZYME-SUBSTRATE INTERMEDIATE WITH SACCHARIDE DISTORTION IN A MUTANT T4 LYSOZYME,1.9,HYDROLASE/HYDROLASE SUBSTRATE,S,5,AXXXX,SUBSTRATE CLEAVED FROM CELL WALL OF ESCHERICHIA COLI,Escherichia coli,151.7,-,0,-,5.57,E,163,MNIFEMLRIDEGLRLKIYKDTEGYYEIGIGHLLTKSPSLNAAKSELDKAIGRNTNGVITKDEAEKLFNQDVDAAVRGILRNAKLKPVYDSLDAVRRAALINMVFQMGETGVAGFTNSLRMLQQKRWDEAAVNLAKSRWYNQTPNRAKRVITTFRTGTWDAYKE,T4 LYSOZYME,Escherichia virus T4
1a07,C-SRC COMPLEXED WITH PEPTIDE,2.2,COMPLEX,C,3,XEX,PEPTIDE,-,104.76,-,0,-,4,A,105,SIQAEEWYFGKITRRESERLLLNAENPRGTFLVRESETTKGAYCLSVSDFDNAKGLNVKHYKIRKLDSGGFYITSRTQFNSLQQLVAYYSKHADGLCHRLTTVCP,C-SRC TYROSINE KINASE,Homo sapiens
EOF

echo "✅ 테스트 CSV 파일 생성: test_data.csv"

# 실행 권한 부여
echo ""
echo "=== 실행 권한 설정 ==="

SCRIPT_FILES=(
    "run_pdb_pipeline.sh"
    "setup_pdb_pipeline.sh"
)

for script in "${SCRIPT_FILES[@]}"; do
    if [ -f "$CURRENT_DIR/$script" ]; then
        chmod +x "$CURRENT_DIR/$script"
        echo "✅ $script 실행 권한 부여"
    fi
done

# 출력 디렉토리 생성
echo ""
echo "=== 출력 디렉토리 생성 ==="

OUTPUT_DIR="$CURRENT_DIR/test_output"
mkdir -p "$OUTPUT_DIR"
echo "✅ 테스트 출력 디렉토리: $OUTPUT_DIR"

# 설치 완료 안내
echo ""
echo "=========================================="
echo "설치 완료!"
echo "=========================================="
echo ""
echo "📂 파일 구조:"
echo "  $PARENT_DIR/ (기존 SUMD 파일들)"
echo "  │"
echo "  └── pdb_preprocessing/ (새 PDB 전처리 파일들)"
echo "      ├── pdb_csv_pipeline.py      # 메인 파이프라인"
echo "      ├── csv_validator.py         # CSV 검증 도구"
echo "      ├── simple_chain_extractor.py # 체인 추출기"
echo "      ├── run_pdb_pipeline.sh      # 배치 실행"
echo "      ├── test_data.csv           # 테스트 데이터"
echo "      └── test_output/            # 테스트 결과"
echo ""
echo "🚀 빠른 시작 가이드:"
echo ""
echo "1. CSV 파일 검증:"
echo "   python3 csv_validator.py --input test_data.csv --validate"
echo ""
echo "2. 체인 추출 테스트 (PDB 다운로드 후):"
echo "   python3 simple_chain_extractor.py --input sample.pdb --output extracted.pdb --chains A,B --info"
echo ""
echo "3. 간단한 파이프라인 테스트 (1개 항목만):"
echo "   python3 pdb_csv_pipeline.py --csv test_data.csv --output test_output --max-entries 1 --verbose"
echo ""
echo "4. 배치 스크립트로 실행:"
echo "   ./run_pdb_pipeline.sh test_data.csv test_output 1"
echo ""
echo "📊 처리 과정:"
echo "  1. CSV 파싱 → PDB ID, Receptor/Ligand Chain 추출"
echo "  2. PDB 다운로드 → RCSB에서 자동 다운로드"
echo "  3. Chain 추출 → 지정된 체인만 분리"
echo "  4. PDB 전처리 → GROMACS 호환성 수정"
echo "  5. pdb2gmx 실행 → 토폴로지 생성"
echo ""
echo "📁 결과 구조:"
echo "  test_output/"
echo "  ├── downloaded_pdbs/     # 다운로드된 원본 PDB"
echo "  ├── extracted_chains/    # 추출된 체인들"
echo "  ├── processed_pdbs/      # GROMACS 호환 처리"
echo "  ├── gromacs_results/     # pdb2gmx 결과"
echo "  └── processing_summary.json # 전체 요약"
echo ""
echo "⚠️  주의사항:"
echo "  - 인터넷 연결 필요 (PDB 다운로드)"
echo "  - 처리 시간이 오래 걸릴 수 있음"
echo "  - 대량 처리시 단계별로 테스트 권장"
echo ""
echo "🔧 문제 해결:"
echo "  - 로그 파일: test_output/pipeline.log"
echo "  - 상세 로그: --verbose 옵션 사용"
echo "  - 단계별 테스트: --max-entries 1 옵션"
echo ""
echo "💡 고급 사용법:"
echo "  - 설정 파일: --config $PARENT_DIR/gromacs_commands_config.json"
echo "  - CSV 정제: csv_validator.py --clean 옵션"
echo "  - 체인 정보 확인: simple_chain_extractor.py --info 옵션"
echo "=========================================="

# 시스템 체크 요약
echo ""
echo "🔍 시스템 준비도 체크:"
echo "  Python 모듈: ${#PYTHON_MODULES[@]}/${#PYTHON_MODULES[@]} ✅"
echo "  기존 SUMD 파일: ${#FOUND_EXISTING[@]}/${#EXISTING_FILES[@]} $([ ${#FOUND_EXISTING[@]} -eq ${#EXISTING_FILES[@]} ] && echo '✅' || echo '⚠️')"
echo "  새 파이프라인 파일: ${#FOUND_NEW[@]}/${#NEW_FILES[@]} $([ ${#FOUND_NEW[@]} -eq ${#NEW_FILES[@]} ] && echo '✅' || echo '❌')"
echo "  GROMACS: $(command -v gmx >/dev/null 2>&1 && echo '✅' || echo '❌')"
echo ""

if [ ${#FOUND_NEW[@]} -eq ${#NEW_FILES[@]} ] && [ ${#MISSING_MODULES[@]} -eq 0 ]; then
    echo "🎉 모든 준비가 완료되었습니다! 파이프라인을 실행할 수 있습니다."
else
    echo "⚠️  일부 파일이나 모듈이 누락되었습니다. 위의 안내를 참조하여 설치를 완료하세요."
fi

echo ""
echo "설치 스크립트 완료!"
