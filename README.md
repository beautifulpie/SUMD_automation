# GROMACS 기반 SuMD (Supervised Molecular Dynamics) 시스템

GROMACS를 SuMD(Supervised Molecular Dynamics)의 엔진으로 활용하여 단백질-펩타이드 결합 시뮬레이션을 수행하는 자동화 시스템입니다.

## 개요

SuMD(Supervised Molecular Dynamics)는 일반적인 MD 시뮬레이션보다 효율적으로 분자 결합 과정을 시뮬레이션할 수 있는 방법입니다. 이 구현에서는 GROMACS를 사용하여 SuMD를 구현하였습니다.

이 시스템의 주요 특징:
- GROMACS를 사용한 SuMD 구현
- 복수의 PDB 구조 일괄 처리를 위한 자동화 시스템
- 에너지 최소화와 짧은 MD 시뮬레이션을 조합하여 결합 과정 탐색
- 최적의 결합 포즈 자동 저장

## 요구 사항

- GROMACS (2020 이상 버전 권장)
- Python 3.6 이상
- 필요 Python 패키지: 
  - tqdm
  - numpy

## 파일 구성

- `SuMD_script.py`: 단일 PDB에 대해 SuMD 시뮬레이션을 수행하는 메인 스크립트
- `Running_SuMD_Automation.py`: 여러 PDB를 일괄 처리하기 위한 자동화 스크립트
- `README.md`: 사용 설명서 (현재 파일)

## ff

https://mackerell.umaryland.edu/download.php?filename=CHARMM_ff_params_files/charmm36-jul2022.ff.tgz


## 사용 방법

### 단일 PDB 시뮬레이션 (SuMD_script.py)

```bash
python SuMD_script.py --pdb <PDB 파일 경로> \
                      --peptide_res <펩타이드 잔기 번호들> \
                      --protein_res <단백질 잔기 번호들> \
                      [--force_field <힘장>] \
                      [--water_model <물 모델>] \
                      [--max_emin_iter <에너지 최소화 반복 횟수>] \
                      [--max_segments <최대 MD 세그먼트 수>] \
                      [--distance_threshold <결합 기준 거리(nm)>] \
                      [--output_dir <결과 저장 디렉토리>]
```

예시:
```bash
python SuMD_script.py --pdb complex.pdb --peptide_res 45 46 47 --protein_res 123 124 125 --output_dir results
```

### 여러 PDB 일괄 처리 (Running_SuMD_Automation.py)

```bash
python Running_SuMD_Automation.py --dataset_path <PDB 파일 디렉토리> \
                                  --save_path <결과 저장 디렉토리> \
                                  --sumd_script_path <SuMD_script.py 경로> \
                                  --pdb_list_path <처리할 PDB 목록 파일> \
                                  --protein_chain <단백질 체인 ID> \
                                  --peptide_chain <펩타이드 체인 ID> \
                                  [--box_size <시뮬레이션 박스 크기(nm)>]
```

예시:
```bash
python Running_SuMD_Automation.py --dataset_path ./pdbs/ \
                                  --save_path ./results/ \
                                  --sumd_script_path ./SuMD_script.py \
                                  --pdb_list_path ./pdb_list.txt \
                                  --protein_chain A \
                                  --peptide_chain B
```

## PDB 목록 파일 형식

PDB 목록 파일은 한 줄에 하나의 PDB ID를 포함해야 합니다:
```
1a2b
3c4d
5e6f
```

확장자(.pdb)가 있어도 처리 가능합니다:
```
1a2b.pdb
3c4d.pdb
5e6f.pdb
```

## 결과 구성

각 PDB 시뮬레이션 결과는 지정된 저장 디렉토리 아래 PDB ID 이름의 폴더에 저장되며, 다음과 같은 파일들이 포함됩니다:

- `best_trajectories/`: 시뮬레이션 중 거리가 가장 가까운 상위 5개 구조
  - `best_structure_1.gro`, `best_structure_1.tpr`: 가장 좋은 구조 파일
  - `best_structure_1_info.txt`: 구조 정보 (거리 등)
- `emin_*.gro`, `emin_*.tpr`: 에너지 최소화 결과 파일
- `segment_*.gro`, `segment_*.xtc`, `segment_*.tpr`: MD 세그먼트 결과 파일
- `distance.xvg`: 거리 측정 결과

## 주의 사항

1. 실행 전 GROMACS가 올바르게 설치되어 있는지 확인하세요.
2. 시뮬레이션은 계산 자원을 많이 사용할 수 있으므로, 충분한 CPU와 메모리가 있는 환경에서 실행하는 것이 좋습니다.
3. 큰 단백질-펩타이드 복합체의 경우 시뮬레이션 시간이 길어질 수 있습니다.

## 알고리즘 설명

SuMD 알고리즘은 다음과 같은 과정으로 동작합니다:

1. **초기 시스템 준비**:
   - PDB 구조 로드 및 처리
   - 시뮬레이션 박스 생성 및 물 분자 추가
   - 이온 추가로 시스템 중성화
   - 펩타이드와 단백질 결합 부위 그룹 정의

2. **에너지 최소화 단계**:
   - 여러 개의 에너지 최소화 샘플 생성
   - 각 샘플에서 펩타이드-단백질 간 거리 측정
   - 거리가 가장 가까운 샘플 선택하여 다음 단계 진행

3. **MD 시뮬레이션 단계**:
   - 짧은 MD 시뮬레이션 세그먼트 실행
   - 각 세그먼트 후 펩타이드-단백질 간 거리 측정
   - 거리가 감소하는 경우 다음 세그먼트 진행
   - 지정된 거리(cutoff) 이하로 도달하면 시뮬레이션 종료

4. **결과 분석 및 저장**:
   - 시뮬레이션 중 가장 좋은 구조 선택 및 저장
   - 거리 데이터 및 궤적 저장

## 예제 실행

### 1. 단일 PDB 시뮬레이션 예제

```bash
# 폴더 생성 및 이동
mkdir -p ~/sumd_test
cd ~/sumd_test

# 스크립트 복사
cp /path/to/SuMD_script.py .

# PDB 파일 준비 (1a2b.pdb 파일이 있다고 가정)
# SuMD 실행
python SuMD_script.py --pdb 1a2b.pdb --peptide_res 10 11 12 --protein_res 100 101 102 --output_dir results_1a2b
```

### 2. 여러 PDB 일괄 처리 예제

```bash
# PDB 목록 파일 생성
echo "1a2b.pdb" > pdb_list.txt
echo "3c4d.pdb" >> pdb_list.txt

# 자동화 스크립트 실행
python Running_SuMD_Automation.py --dataset_path ./pdbs/ \
                                 --save_path ./results/ \
                                 --sumd_script_path ./SuMD_script.py \
                                 --pdb_list_path ./pdb_list.txt \
                                 --protein_chain A \
                                 --peptide_chain B
```

## 문제 해결

- **GROMACS 오류 발생 시**: GROMACS 설치가 올바른지 확인하고, 필요한 환경 변수가 설정되어 있는지 확인하세요.
- **메모리 부족 오류**: 큰 시스템의 경우 더 많은 메모리가 필요할 수 있습니다. 시스템 크기를 줄이거나 더 많은 메모리가 있는 시스템에서 실행하세요.
- **시뮬레이션이 너무 오래 걸림**: `--max_emin_iter`와 `--max_segments` 값을 줄여 총 시뮬레이션 시간을 단축할 수 있습니다.