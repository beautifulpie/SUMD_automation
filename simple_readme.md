# Simple SuMD

핵심 기능만 포함한 간단한 SuMD (Supervised Molecular Dynamics) 시뮬레이션 도구

## 🎯 특징

- **단순함**: 복잡한 클래스 구조 없이 함수 기반의 간단한 구조
- **핵심 기능**: GROMACS 실행, attempt/iteration 제어, 기본 로깅만 포함
- **JSON 출력**: 각 단계별 상세 결과를 JSON으로 기록
- **긴 MD 지원**: 체인 간 거리가 가까워지면 자동으로 긴 MD 시뮬레이션 실행
- **체인 복원**: topology 파일을 참조하여 최종 PDB의 체인 정보 자동 복원
- **유지보수 용이**: 코드 수정이 간단하고 직관적

## 📁 파일 구조

```
simple_sumd/
├── simple_sumd.py          # 메인 Python 스크립트
├── simple_sumd.sh          # Shell 실행 스크립트
├── simple_config.py        # 설정 파일
├── analyze_results.py      # 결과 분석 도구
└── README.md              # 이 파일
```

## 🚀 사용법

### 기본 실행

```bash
# 실행 권한 부여
chmod +x simple_sumd.sh

# 기본 실행
./simple_sumd.sh complex.pdb A B

# 출력 디렉토리 지정
./simple_sumd.sh complex.pdb A B my_output
```

### Python 직접 실행

```bash
python3 simple_sumd.py complex.pdb A B output_dir
```

## ⚙️ 설정

`simple_config.py` 파일에서 다음 설정들을 수정할 수 있습니다:

### 기본 시뮬레이션 설정
```python
MAX_ITERATIONS = 10          # 최대 iteration 수
MAX_ATTEMPTS = 100           # iteration당 최대 attempt 수  
SIMULATION_TIME_NS = 0.3     # 기본 MD 시뮬레이션 시간 (300ps)
SLOPE_THRESHOLD = -0.001     # 기울기 임계값
```

### 긴 MD 설정 (NEW!)
```python
CLOSE_DISTANCE_THRESHOLD = 10.0    # 긴 MD 실행 거리 임계값 (Å)
LONG_MD_TIME_NS = 10.0            # 긴 MD 시뮬레이션 시간 (ns)
ENABLE_LONG_MD = True             # 긴 MD 기능 활성화
TIMEOUT_LONG_MD = 7200            # 긴 MD 타임아웃 (초)
```

### GROMACS 설정
```python
GPU_ID = "3"                 # 사용할 GPU ID
MPI_RANKS = "4"              # MPI 랭크 수
FORCE_FIELD = "charmm36-jul2022"  # Force field
```

## 🧬 긴 MD 기능

### 동작 원리
1. 기본 300ps MD 시뮬레이션 실행
2. 궤적에서 체인 간 최소 거리 계산
3. **거리가 10Å 이내로 가까워지면 자동으로 10ns 긴 MD 실행**
4. 긴 MD 결과를 최종 기울기 판정에 사용

### 장점
- 근접 접촉 시 더 정확한 결합/해리 예측
- 초기 접촉 후 안정화 과정 관찰
- 더 신뢰할 수 있는 기울기 계산

### 주의사항
- 긴 MD 실행 시 시뮬레이션 시간 대폭 증가 (30분-2시간)
- 테스트 시에는 `ENABLE_LONG_MD = False` 권장
- GPU 메모리 및 디스크 공간 충분히 확보

## 📊 출력 결과

시뮬레이션 완료 후 다음 파일들이 생성됩니다:

```
output_dir/
├── final_results.json                    # 전체 결과 요약
├── iteration_1_summary.json              # Iteration 1 요약
├── iteration_1_attempt_1.json            # Iteration 1, Attempt 1 상세
├── iteration_1_attempt_2.json            # Iteration 1, Attempt 2 상세
├── ...
├── accepted_result.pdb                   # 최종 채택된 구조
└── attempt_*/                            # 각 attempt의 GROMACS 파일들
    ├── md.xtc                           # 기본 MD 궤적
    ├── long_md.xtc                      # 긴 MD 궤적 (실행된 경우)
    └── ...
```

### JSON 결과 형식

#### attempt_detail.json (업데이트됨)
```json
{
  "attempt": 1,
  "success": true,
  "slope": -0.0025,
  "distances": [25.3, 24.8, 18.7, 12.1, 8.5, ...],
  "min_distance": 8.5,
  "final_distance": 8.2,
  "close_contact_detected": true,
  "long_md_executed": true,
  "long_md_distances": [8.5, 8.3, 8.1, 8.2, ...],
  "stages": [...]
}
```

## 🔧 알고리즘

1. **Iteration 루프** (최대 10회)
   - 각 iteration마다 최대 100번 attempt 시도
   
2. **Attempt 처리**
   - GROMACS 파이프라인 실행: pdb2gmx → editconf → solvate → genion → EM → NVT → NPT → MD
   - 300ps 기본 MD 시뮬레이션 수행
   - 궤적에서 체인 간 거리 추출
   - **거리 < 10Å이면 10ns 긴 MD 추가 실행**
   - 최종 거리 기울기 계산
   
3. **채택 기준**
   - 기울기 < -0.001이면 채택
   - 긴 MD가 실행된 경우 긴 MD 기울기 사용
   - 채택되면 다음 iteration의 입력으로 사용

## 🛠️ 요구사항

- Python 3.6+
- GROMACS 2021+
- BioPython
- NumPy
- Matplotlib (분석용)
- CUDA 지원 GPU (권장)

## 📈 성능 팁

1. **빠른 테스트**: 
   ```python
   MAX_ATTEMPTS = 10
   ENABLE_LONG_MD = False
   ```

2. **GPU 최적화**: 
   ```python
   MPI_RANKS = "1"  # 단일 프로세스
   ```

3. **대형 시스템**: 
   ```python
   LONG_MD_TIME_NS = 5.0  # 긴 MD 시간 단축
   ```

4. **디스크 공간**: 
   ```python
   KEEP_FAILED_ATTEMPTS = False
   ```

## 📈 결과 분석

```bash
# 상세 분석 실행
python3 analyze_results.py output_directory

# 생성되는 파일들
output_directory/
├── analysis_plots.png      # 거리 변화 그래프
├── detailed_report.txt     # 상세 분석 보고서
```

## 🔍 문제 해결

### 일반적인 오류

1. **GROMACS 경로 문제**
   ```bash
   which gmx
   which gmx_mpi
   ```

2. **긴 MD 타임아웃**
   - `TIMEOUT_LONG_MD` 값 증가
   - 또는 `LONG_MD_TIME_NS` 값 감소

3. **GPU 메모리 부족**
   - `GPU_ID` 변경
   - `MPI_RANKS` 줄이기
   - `ENABLE_LONG_MD = False`로 설정

### 로그 확인

- 실시간 로그: 터미널 출력
- 상세 결과: JSON 파일들
- GROMACS 로그: 각 attempt 디렉토리 내
- 긴 MD 로그: `attempt_*/long_md.log`

## 🆕 업데이트 내역

### v2.0 (현재 버전)
- **긴 MD 기능 추가**: 10Å 이내 근접 시 자동 10ns MD 실행
- 거리 기반 적응적 시뮬레이션
- 결과 JSON에 긴 MD 정보 추가
- 타임아웃 설정 개선

### v1.0
- 기본 SuMD 기능
- JSON 결과 출력
- 기본 분석 도구

## 🤝 기여

버그 신고나 기능 요청은 이슈로 등록해주세요.

---

**Note**: 이 도구는 연구 목적으로 개발되었습니다. 상업적 사용 시 해당 기관의 라이선스 정책을 확인하세요.