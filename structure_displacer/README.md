# Structure Displacer - 50Å 거리 기반 단백질 복합체 변형 시스템

정확한 50Å 체인 간 거리를 가진 변형된 단백질 복합체 구조를 생성하는 시스템입니다. 필요에 따라 SUMD(Supervised Molecular Dynamics) 시뮬레이션과 연동할 수 있습니다.

## 🎯 개요

Structure Displacer 시스템은 지정된 두 체인 간의 최소 거리가 정확히 50Å가 되도록 수정된 단백질 복합체 구조를 생성합니다. 이 시스템의 주요 목적:

- 정확한 50Å 거리로 체인 간 분리된 구조 생성
- 변형 중 구조적 무결성 유지
- 충돌 없고 경로가 막히지 않은 변형 보장
- 포괄적인 품질 분석 및 검증 제공
- 선택적으로 SUMD 시뮬레이션을 위한 시작 구조 준비

## 📁 시스템 구조

### 핵심 모듈

1. **`structure_displacer.py`** - 메인 변형 생성 엔진
   - 50Å 거리 기반 변형 알고리즘 구현
   - 충돌 검사 및 경로 차단 검사 처리
   - 적응형 변형 거리 계산 지원
   - 다중 변형 변형체 생성

2. **`displacement_utils.py`** - 품질 분석 및 검증 유틸리티
   - 다중 메트릭을 사용한 변형 품질 분석
   - 거리, 충돌, 경로 기준에 대한 구조 검증
   - 포괄적인 보고서 및 시각화 생성
   - 품질 점수에 따른 구조 순위 매기기

3. **`example_workflow.py`** - 완전 자동화된 워크플로우 관리자
   - 변형 생성, 분석 통합
   - 선택적 SUMD 시뮬레이션 실행
   - 구성 관리 및 결과 추적 처리

## 🚀 빠른 시작

### 기본 사용법 (구조 생성 및 분석만)

변형된 구조 생성 및 품질 분석 (SUMD 없이):

```bash
python example_workflow.py \
    --golden_standard native_complex.pdb \
    --target_chains A,B \
    --num_variants 5 \
    --output_dir my_50A_workflow \
    --skip_sumd
```

### SUMD 포함 전체 워크플로우

변형 생성 + 분석 + SUMD 시뮬레이션:

```bash
python example_workflow.py \
    --golden_standard native_complex.pdb \
    --target_chains A,B \
    --num_variants 5 \
    --output_dir my_50A_workflow
```

### 단계별 사용법

1. **50Å 변형 구조 생성:**
```bash
python structure_displacer.py \
    --input complex.pdb \
    --target_chains A,B \
    --num_variants 5 \
    --target_distance 50.0 \
    --output_dir displaced_structures_50A
```

2. **변형 품질 분석:**
```bash
python displacement_utils.py \
    --golden_standard complex.pdb \
    --displaced_dir displaced_structures_50A \
    --target_chains A,B \
    --plot \
    --output_dir analysis_results
```

3. **선택적: SUMD 시뮬레이션 실행**
```bash
# 최적 구조로 SUMD 실행 (수동으로 선택된 구조 사용)
python /app/SUMD_automation/robust_sumd_master.py \
    --input_pdb analysis_results/best_structure.pdb \
    --golden_standard_pdb complex.pdb \
    --receptor_chain A \
    --ligand_chain B \
    --output_dir sumd_results
```

## 🔧 핵심 클래스 및 함수

### DisplacementConfig
변형 매개변수를 위한 구성 클래스:
- `target_distance_for_validation`: 목표 표면 거리 (50.0Å)
- `clash_threshold`: 충돌 검사 임계값 (2.0Å)
- `max_attempts`: 최대 생성 시도 횟수 (200)
- `use_adaptive_displacement`: 적응형 거리 계산 (True)
- `path_check_threshold`: 경로 차단 검사 임계값 (5.0Å)

### AdvancedStructureDisplacer
메인 변형 생성 클래스:

**주요 메서드:**
- `displace_structure(input_pdb, target_chains, output_pdb)`: 단일 변형 구조 생성
- `calculate_displacement_vector(chain1, chain2)`: 정확한 변형 벡터 계산
- `validate_displaced_structure()`: 생성된 구조의 품질 검증

### MultipleDisplacer
다중 변형체 생성 클래스:
- `generate_multiple_variants()`: 여러 변형 구조 동시 생성

### Advanced50AAnalyzer
포괄적인 품질 분석:
- `analyze_displaced_structure()`: 단일 구조 분석
- `analyze_batch_structures()`: 배치 분석 및 순위 매기기
- `create_analysis_plots()`: 시각화 플롯 생성

## 📊 품질 메트릭

시스템은 다음 세 가지 주요 메트릭을 사용하여 변형된 구조를 평가합니다:

### 1. 거리 달성도 (50% 가중치)
- **목표**: 정확히 50Å 체인 간 최소 거리
- **허용 오차**: ±5Å (15Å까지는 관대하게 허용)
- **점수**: 거리 정확도 기반 (1.0 = 완벽, 0.0 = >10Å 오차)

### 2. 충돌 검사 (30% 가중치)
- **기준**: 원자 충돌 없음 (거리 < 2.0Å)
- **점수**: 충돌 없으면 1.0, 그렇지 않으면 최소 충돌 거리로 스케일링

### 3. 경로 차단 (20% 가중치)
- **기준**: 변형 경로를 막는 단백질 체인 없음
- **해상도**: 1.0Å 경로 검사 해상도
- **점수**: 경로가 깨끗하면 1.0, 차단되면 0.0

### 전체 품질 점수
가중 결합 점수: `0.5 × 거리달성도 + 0.3 × 충돌점수 + 0.2 × 경로점수`

## 🛠️ 주요 사용 옵션

### 기본 변형 생성 옵션
```bash
# 기본 50Å 변형 (추천)
python structure_displacer.py \
    --input complex.pdb \
    --target_chains A,B \
    --num_variants 5

# 더 많은 변형체와 엄격한 조건
python structure_displacer.py \
    --input complex.pdb \
    --target_chains A,B \
    --num_variants 10 \
    --target_distance 50.0 \
    --clash_threshold 2.5 \
    --max_attempts 500
```

### 품질 분석 옵션
```bash
# 기본 분석 (그래프 포함)
python displacement_utils.py \
    --golden_standard native.pdb \
    --displaced_dir displaced_structures_50A \
    --target_chains A,B \
    --plot

# 사용자 정의 분석
python displacement_utils.py \
    --golden_standard native.pdb \
    --displaced_dir displaced_structures_50A \
    --target_chains A,B \
    --target_distance 50.0 \
    --output_dir custom_analysis \
    --plot
```

### 워크플로우 옵션
```bash
# 기본: 구조 생성 + 분석만 (SUMD 없음)
python example_workflow.py \
    --golden_standard native.pdb \
    --target_chains A,B \
    --num_variants 5 \
    --skip_sumd

# SUMD 포함 전체 워크플로우
python example_workflow.py \
    --golden_standard native.pdb \
    --target_chains A,B \
    --num_variants 5 \
    --simulation_time 2.0 \
    --rmsd_threshold 1.5

# 관대한 조건으로 실행
python example_workflow.py \
    --golden_standard native.pdb \
    --target_chains A,B \
    --num_variants 10 \
    --allow_invalid \
    --min_quality_score 0.5 \
    --skip_sumd
```

## 📈 출력 구조

### 기본 출력 (SUMD 없이)
```
output_directory/
├── displaced_structures_50A/           # 생성된 변형 구조들
│   ├── complex_displaced_001.pdb
│   ├── complex_displaced_002.pdb
│   └── complex_displaced_003.pdb
├── quality_analysis_50A/               # 품질 분석 결과
│   ├── displacement_analysis_50A.csv   # 상세 분석 데이터
│   ├── displacement_analysis_50A.json  # JSON 형식 결과
│   ├── analysis_summary_50A.json       # 요약 통계
│   └── displacement_analysis_50A_plots.png  # 시각화 플롯
└── workflow_result_50A.json            # 워크플로우 요약
```

### SUMD 포함 출력
```
output_directory/
├── displaced_structures_50A/           # 생성된 변형 구조들
├── quality_analysis_50A/               # 품질 분석 결과  
├── sumd_simulation_50A/                # SUMD 시뮬레이션 결과
│   ├── sumd_[job_id]/                 # SUMD 작업 디렉토리
│   │   ├── iteration_1/
│   │   ├── iteration_2/
│   │   └── results.json
│   └── final_structure.pdb
└── workflow_result_50A.json            # 완전한 워크플로우 요약
```

### 결과 해석

**구조 생성 성공:**
```
=== 변형 구조 생성 결과 ===
목표: 5개
성공: 5개
출력 디렉토리: displaced_structures_50A

생성된 구조들:
  1. complex_displaced_001.pdb
     최종 거리: 50.12Å
     시도 횟수: 23회
     거리 오차: 0.12Å
```

**품질 분석 결과:**
```
=== 50Å 거리 기반 구조 변형 분석 결과 ===
분석된 파일: 5개
유효한 구조: 4개 (80.0%)
무효한 구조: 1개 (20.0%)

거리 달성 통계:
  평균 거리: 50.08±0.15Å
  목표 거리 달성률: 100.0%

품질 통계:
  Clash 없음: 80.0%
  경로 차단 없음: 100.0%
  평균 품질 점수: 0.823±0.145

최고 품질 유효 구조 TOP 3:
  1위: complex_displaced_003.pdb
       품질점수: 0.952, 거리: 50.03Å (오차: 0.03Å)
  2위: complex_displaced_001.pdb
       품질점수: 0.887, 거리: 49.87Å (오차: 0.13Å)
```

## 🔍 알고리즘 세부사항

### 변형 방법
1. **표면 중심 계산**: 목표 체인들의 마주보는 표면 식별
2. **적응형 거리 계산**: 50Å 목표를 위한 정확한 이동 거리 계산
3. **방향 벡터 생성**: 체인 중심으로부터 변형 방향 계산
4. **랜덤 변화**: 다양성을 위한 제어된 각도 변화 (±15°) 적용
5. **검증**: 충돌, 경로 차단, 거리 정확도 검사

### 품질 평가
1. **거리 검증**: 50Å ± 5Å 달성 확인 (최대 15Å까지 허용)
2. **충돌 검사**: < 2.0Å 원자 겹침 식별
3. **경로 분석**: 변형 궤적의 차단 검사
4. **점수 매기기**: 가중 공식으로 메트릭 결합

## 🚨 문제 해결

### 일반적인 문제들

**생성 실패:**
```bash
# 시도 횟수 증가 및 충돌 임계값 완화
python structure_displacer.py \
    --input complex.pdb \
    --target_chains A,B \
    --clash_threshold 1.8 \
    --max_attempts 1000 \
    --verbose
```

**유효한 구조 없음:**
```bash
# 더 낮은 품질 임계값으로 무효한 구조도 허용
python example_workflow.py \
    --golden_standard complex.pdb \
    --target_chains A,B \
    --allow_invalid \
    --min_quality_score 0.5 \
    --skip_sumd
```

**디버깅 모드:**
```bash
# 진단을 위한 상세 로깅
python example_workflow.py \
    --golden_standard complex.pdb \
    --target_chains A,B \
    --verbose \
    --skip_sumd
```

**체인 ID 확인:**
```bash
# PDB 파일의 체인 정보 확인
grep "^ATOM" complex.pdb | awk '{print $5}' | sort | uniq
```

## 📋 의존성

### 필수 Python 패키지
- `numpy` - 수치 계산
- `scipy` - 과학 계산  
- `pandas` - 데이터 분석
- `matplotlib` - 플롯팅 (선택적: --plot 사용시)
- `seaborn` - 통계 시각화 (선택적: --plot 사용시)
- `Biopython` - PDB 구조 처리
- `MDAnalysis` - 분자 구조 분석 (선택적)

### 시스템 요구사항
- Python 3.7+
- GROMACS (SUMD 시뮬레이션용, 선택적)
- Unix 계열 환경 (Linux/macOS)

## 🎯 주요 사용 사례

### 1. 구조 분석 목적 (기본)
50Å 분리된 구조 생성 및 품질 분석:
```bash
python example_workflow.py \
    --golden_standard native.pdb \
    --target_chains A,B \
    --num_variants 10 \
    --skip_sumd \
    --plot
```

### 2. SUMD 시뮬레이션 준비
구조 생성 후 최적 구조로 SUMD 실행:
```bash
# 1단계: 구조 생성 및 분석
python example_workflow.py \
    --golden_standard native.pdb \
    --target_chains A,B \
    --num_variants 5 \
    --skip_sumd

# 2단계: 수동으로 최적 구조 선택 후 SUMD 실행
python /app/SUMD_automation/robust_sumd_master.py \
    --input_pdb workflow_output_50A/displaced_structures_50A/best_structure.pdb \
    --golden_standard_pdb native.pdb \
    --receptor_chain A \
    --ligand_chain B
```

### 3. 전체 자동화 (고급 사용자)
구조 생성부터 SUMD까지 전체 자동화:
```bash
python example_workflow.py \
    --golden_standard native.pdb \
    --target_chains A,B \
    --num_variants 5 \
    --simulation_time 2.0 \
    --rmsd_threshold 1.5
```

### 4. 배치 처리
여러 PDB 파일 처리:
```bash
#!/bin/bash
for pdb_file in *.pdb; do
    echo "Processing $pdb_file..."
    python structure_displacer.py \
        --input "$pdb_file" \
        --target_chains A,B \
        --num_variants 3 \
        --output_dir "results_$(basename $pdb_file .pdb)"
done
```

## 🔬 기술 사양

### 변형 알고리즘
- **방법**: 표면 기반 적응형 기하학적 변형
- **정밀도**: 50Å 목표에 대한 서브-앙스트롬 정확도
- **검증**: 다중 기준 품질 평가 (거리 + 충돌 + 경로)
- **다양성**: 제어된 랜덤 방향 변화로 다중 변형체 생성

### 성능 특성
- **속도**: 구조당 1-30초 (복잡도에 따라)
- **성공률**: 표준 단백질 복합체 80-95%
- **메모리**: 구조 크기에 선형 비례
- **확장성**: ~100k 원자까지 지원

### 검증 기준
- **거리 정확도**: 50Å ± 5Å (최대 15Å까지 허용)
- **충돌 임계값**: 최소 2.0Å 원자간 거리
- **경로 해상도**: 1.0Å 간격 차단 검사
- **품질 임계값**: 기본 0.7 (조정 가능)

---

**최종 업데이트**: 2024년  
**버전**: 2.0 (50Å 거리 기반)  
**호환성**: SUMD 자동화 시스템 v2.x  
**기본 모드**: 구조 생성 + 분석 (SUMD 선택적)