# GROMACS-GPU 기반 SuMD 테스트 환경 설정 가이드

이 가이드는 기존 GROMACS-GPU 도커 이미지를 기반으로 SuMD(Supervised Molecular Dynamics) 시뮬레이션을 위한 테스트 환경을 구축하는 방법을 설명합니다.

## 1. 필요한 파일

먼저 다음 파일들을 준비하세요:

- `Dockerfile`: GROMACS-GPU 이미지 기반 SuMD 환경 정의
- `docker-compose.yml`: 컨테이너 설정 및 GPU 지원 정의
- `SuMD_script.py`: 개별 SuMD 시뮬레이션 스크립트
- `Running_SuMD_Automation.py`: 일괄 처리 자동화 스크립트
- `prepare_sumd.sh`: PDB 파일 전처리 스크립트

## 2. 디렉토리 구조 설정

```
sumd-test/
├── Dockerfile
├── docker-compose.yml
├── SuMD_script.py
├── Running_SuMD_Automation.py
├── prepare_sumd.sh
├── pdbs/
│   ├── protein1.pdb
│   ├── protein2.pdb
│   └── pdb_list.txt
└── results/
```

다음 명령으로 필요한 디렉토리를 생성하세요:

```bash
mkdir -p sumd-test/pdbs sumd-test/results
```

## 3. 도커 이미지 빌드

디렉토리로 이동하여 도커 이미지를 빌드합니다:

```bash
cd sumd-test
docker build -t sumd-gromacs-gpu .
```

## 4. PDB 파일 준비

테스트할 PDB 파일을 `pdbs` 디렉토리에 복사하고, 처리할 PDB 파일 목록을 작성합니다:

```bash
# PDB 파일 복사
cp path/to/your/pdb/files/*.pdb pdbs/

# PDB 목록 생성
cd pdbs
ls *.pdb > pdb_list.txt
cd ..
```

## 5. PDB 전처리 (선택 사항)

특수 잔기(HEM, GOL, SO4 등)가 포함된 PDB 파일은 전처리가 필요할 수 있습니다. 컨테이너 내에서 전처리 스크립트를 실행하세요:

```bash
# 대화형 모드로 컨테이너 실행
docker run -it --rm -v $(pwd)/pdbs:/sumd/data -v $(pwd)/results:/sumd/results sumd-gromacs-gpu

# 컨테이너 내부에서 실행
cd /sumd
./prepare_sumd.sh --pdb /sumd/data/your_protein.pdb --output /sumd/data/prepared
```

또는 호스트에서 직접 실행할 수도 있습니다:

```bash
chmod +x prepare_sumd.sh
./prepare_sumd.sh --pdb pdbs/your_protein.pdb --output pdbs/prepared
```

## 6. Docker Compose로 실행

### 6.1 자동 모드

docker-compose.yml 파일을 사용하여 자동으로 모든 PDB 파일을 처리합니다:

```bash
docker-compose up
```

### 6.2 대화형 모드

docker-compose.yml 파일에서 `command` 섹션을 주석 처리하고 `stdin_open`과 `tty` 옵션을 활성화한 후:

```bash
docker-compose run --rm sumd-gpu
```

## 7. 단일 PDB 시뮬레이션

컨테이너 내부에서 특정 PDB 파일에 대해 SuMD 시뮬레이션을 실행합니다:

```bash
python SuMD_script.py --pdb /sumd/data/your_protein.pdb \
                       --peptide_res 45 46 47 \
                       --protein_res 123 124 125 \
                       --output_dir /sumd/results/your_protein
```

## 8. 고급 사용법: 특수 잔기 처리

### 8.1 단백질만 사용하는 방법

특수 잔기가 포함된 PDB 파일에서 단백질 부분만 추출:

```bash
grep "^ATOM" /sumd/data/your_protein.pdb > /sumd/data/your_protein_protein_only.pdb
```

### 8.2 -merge 옵션 사용

알려진 잔기만 처리하고 나머지는 무시:

```bash
python SuMD_script.py --pdb /sumd/data/your_protein.pdb \
                       --peptide_res 45 46 47 \
                       --protein_res 123 124 125 \
                       --output_dir /sumd/results/your_protein \
                       --merge
```

## 9. GPU 사용 확인

시뮬레이션이 GPU를 사용하고 있는지 확인하려면:

```bash
# 컨테이너 내부에서
nvidia-smi

# 또는 호스트에서
docker exec -it sumd-test_sumd-gpu_1 nvidia-smi
```

## 10. 문제 해결

### 10.1 GROMACS 버전 확인

```bash
gmx --version
```

### 10.2 GPU 지원 확인

```bash
gmx mdrun -version
```
출력에 GPU 지원 관련 정보가 포함되어 있는지 확인하세요.

### 10.3 일반적인 오류

- **"Residue 'XXX' not found in residue topology database"**: 특수 잔기가 포스필드에 정의되어 있지 않음. 해결책: 단백질만 추출하거나 `-merge` 옵션 사용
- **"CUDA error: out of memory"**: GPU 메모리 부족. 해결책: 박스 크기를 줄이거나 작은 단백질로 테스트
- **"Fatal error: Invalid box vector"**: 박스 정의 오류. 해결책: 박스 크기 매개변수 조정

## 11. SuMD 스크립트 수정 (선택 사항)

필요한 경우 SuMD 스크립트를 수정하여 성능을 최적화하거나 GPU 특화 기능을 활용할 수 있습니다:

```bash
vi /sumd/SuMD_script.py
```

주요 수정 포인트:
- `-nb gpu` 옵션 추가하여 비결합 계산을 GPU에서 수행
- 스레드 및 MPI 프로세스 수 최적화
- GPU 특화 GROMACS 매개변수 조정

이제 GROMACS-GPU 이미지 기반 SuMD 테스트 환경이 준비되었습니다!
