# SUMD 자동화 시스템

GROMACS를 사용한 단백질-리간드 복합체 시뮬레이션을 위한 포괄적인 SUMD (Supervised Molecular Dynamics) 자동화 시스템입니다. 거리 기반 수렴 기준을 통해 단백질-리간드 결합을 연구하는 반복적인 분자동역학 시뮬레이션을 수행합니다.

## 개요

SUMD 자동화 시스템은 두 가지 주요 구성 요소로 이루어져 있습니다:
1. **robust_sumd_master** - 핵심 SUMD 시뮬레이션 오케스트레이터
2. **structure_displacer** - 고급 구조 변형 및 분석 도구

## 핵심 구성 요소

### 주요 진입점
- `robust_sumd_master.py` - SUMD 시뮬레이션용 메인 Python 오케스트레이터
- `robust_sumd_master.sh` - 포괄적인 매개변수 처리 및 로깅이 포함된 셸 래퍼 스크립트
- `scripts/` - 메인 자동화 스크립트의 심볼릭 링크 버전들

### 핵심 모듈
- `gromacs_runner_class.py` - GROMACS 명령 실행 및 관리
- `interface_distance_calculator.py` - 단백질 체인 간 거리 계산
- `mda_pdb_processor.py` - MDAnalysis를 사용한 PDB 파일 전처리
- `MDPGenerator.py` - GROMACS .mdp 매개변수 파일 생성
- `utils.py` - NumPy 배열용 JSON 인코딩을 포함한 유틸리티 함수

### 설정
- `gromacs_commands_config.json` - 포괄적인 GROMACS 명령 템플릿 및 시스템별 매개변수
- `templates/` - 다양한 시뮬레이션 유형을 위한 .mdp 템플릿 파일들

## robust_sumd_master 사용법

### 기본 SUMD 시뮬레이션

```bash
# 기본 매개변수로 기본 SUMD 실행
./robust_sumd_master.sh input.pdb A B

# 사용자 정의 매개변수로 실행
./robust_sumd_master.sh input.pdb A B 2.0 5.0 1.5 10 5 /app/output
```

### Golden Standard SUMD (네이티브 구조로의 수렴)

```bash
# 네이티브 구조와 비교하는 골든 스탠다드 시뮬레이션 실행
./robust_sumd_master.sh displaced.pdb native.pdb A B 2.0 5.0 1.5 10 5 /app/output
```

### Python 인터페이스

```bash
# 전체 매개변수 제어로 직접 Python 실행
python3 robust_sumd_master.py \
  --input_pdb complex.pdb \
  --receptor_chain A \
  --ligand_chain B \
  --simulation_time 2.0 \
  --distance_threshold 5.0 \
  --rmsd_threshold 1.5 \
  --max_iterations 10 \
  --num_samples 5 \
  --output_dir /path/to/output
```

### 매개변수

| 매개변수 | 설명 | 기본값 | 셸 위치 |
|----------|------|--------|---------|
| `input_pdb` | 복합체 구조 파일 | 필수 | 1 |
| `golden_standard_pdb` | 네이티브 참조 구조 (선택사항) | None | 2 (제공된 경우) |
| `receptor_chain` | 수용체 체인 식별자 | 필수 | 2/3 |
| `ligand_chain` | 리간드 체인 식별자 | 필수 | 3/4 |
| `simulation_time` | 지속 시간 (ns) | 1.0 | 4/5 |
| `distance_threshold` | 활성화 거리 (Å) | 5.0 | 5/6 |
| `rmsd_threshold` | 수렴 기준 (Å) | 1.5 | 6/7 |
| `max_iterations` | 최대 시뮬레이션 사이클 | 10 | 7/8 |
| `num_samples` | 반복당 병렬 샘플 수 | 5 | 8/9 |
| `output_dir` | 출력 디렉토리 | 자동 생성 | 9/10 |

### 사용 예시

#### 표준 SUMD 시뮬레이션
```bash
# 기본 단백질-리간드 해리 연구
./robust_sumd_master.sh complex.pdb A B 2.0 5.0 1.5 10 5

# 사용자 정의 출력 디렉토리와 함께
./robust_sumd_master.sh complex.pdb A B 2.0 5.0 1.5 10 5 /app/results
```

#### Golden Standard SUMD
```bash
# 변위된 구조를 네이티브 결합 포즈와 비교
./robust_sumd_master.sh displaced_complex.pdb native_complex.pdb A B 2.0 5.0 1.5 10 5
```

#### Python API 예시
```python
# Python API 가져오고 사용하기
from robust_sumd_master import SUMDRunner

runner = SUMDRunner(
    input_pdb="complex.pdb",
    receptor_chain="A",
    ligand_chain="B",
    simulation_time=2.0,
    distance_threshold=5.0,
    rmsd_threshold=1.5,
    max_iterations=10,
    num_samples=5
)

result = runner.run_sumd_workflow()
```

## 시스템 테스트

```bash
# 환경 및 의존성 테스트
./test_sumd.sh

# 테스트 PDB 구조 다운로드
./download_test_pdb.sh
```

## 시스템 설정

### GROMACS 설정

시스템은 GPU 가속화를 지원하는 표준 및 MPI 버전의 GROMACS를 모두 지원합니다. 설정은 다음을 포함하는 `gromacs_commands_config.json`을 통해 처리됩니다:

- **Force Field**: CHARMM36, AMBER, GROMOS, OPLS-AA
- **GPU 설정**: 기본 GPU ID: 3
- **MPI 병렬화**: 멀티코어 최적화
- **시스템 크기 최적화**: 자동 감지 및 최적화
- **타임아웃 관리**: 시뮬레이션 단계별 다른 타임아웃 설정

### 디렉토리 구조

```
SUMD_automation/
├── robust_sumd_master.py          # 메인 Python 오케스트레이터
├── robust_sumd_master.sh          # 셸 래퍼
├── structure_displacer/           # 고급 변위 도구
├── gromacs_commands_config.json   # GROMACS 설정
├── templates/                     # .mdp 매개변수 템플릿
├── scripts/                       # 심볼릭 링크된 스크립트
└── output/                        # 결과 (자동 생성)
    └── sumd_[job_id]/
        └── iteration_N/
            └── sample_N/           # 개별 시뮬레이션 결과
```

## 주요 기능

### 적응형 시뮬레이션 프로토콜
- 거리 기반 시뮬레이션 선택 (EM-only vs EM+MD)
- 네이티브 구조 비교를 위한 골든 스탠다드 수렴
- 반복당 다중 샘플 병렬 실행
- 포괄적인 로깅 및 진행 상황 추적

### 출력 파일
- `results.json` - 완전한 시뮬레이션 메타데이터 및 결과
- `progress_tracker.txt` - 실시간 진행 상황 업데이트
- `simulation_summary.txt` - 사람이 읽을 수 있는 요약
- 각 반복/샘플별 개별 PDB 구조

### 오류 처리
- 대체 옵션이 있는 강력한 전처리
- 시스템 크기 감지 및 최적화
- 포괄적인 타임아웃 관리
- 디버깅을 위한 상세한 로깅

## 의존성

시스템은 다음을 필요로 합니다:
- **GROMACS** (MPI 및 GPU 지원)
- **Python 3** 패키지:
  - numpy, scipy, matplotlib, pandas
  - Bio (Biopython), MDAnalysis
- **표준 Unix 도구**: bash, awk, grep 등

모든 분자동역학 시뮬레이션은 기본적으로 GPU ID 3에서 GPU 가속화를 사용하는 GROMACS로 수행됩니다.

## Structure Displacer

고급 구조 변위 및 분석 기능에 대해서는 [structure_displacer/README.md](structure_displacer/README.md)를 참조하세요.

## 지원

문제, 질문 또는 기여에 대해서는 개별 모듈 문서를 참조하거나 개발팀에 문의하세요.