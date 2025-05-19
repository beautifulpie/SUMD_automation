# SuMD를 활용한 Gromacs 분자 동역학 자동화 시스템

이 프로젝트는 웹서버 없이 직접 실행 방식으로 SuMD(Supervised Molecular Dynamics) 접근법과 Gromacs 엔진을 활용한 분자 동역학 시뮬레이션 자동화 시스템입니다.

## 주요 특징

1. **직접 실행 방식**
   - 웹서버 없이 Docker 컨테이너 내에서 직접 스크립트 실행
   - 간단한 명령어로 복잡한 MD 시뮬레이션 수행

2. **자동화된 워크플로우**
   - PDB 파일 및 관심 잔기 지정만으로 전체 시뮬레이션 자동 수행
   - 무게중심 거리에 따른 시뮬레이션 분기 처리

3. **시간 최적화**
   - dt=0.004, nsteps=2500000 설정으로 10ns 시뮬레이션 빠르게 수행
   - 무게중심 거리가 임계값보다 큰 경우 MD 스킵 (EM만 수행)

4. **특수 잔기 처리**
   - 시스테인(CYS) 잔기에 자동으로 HG 원자 추가
   - 시스테인 이황화 결합 처리 지원

5. **결과 분석**
   - Coulomb 및 LJ 상호작용 에너지 자동 분석
   - 결합 강도 평가 및 결과 시각화

## 필요 환경

- Docker 및 Docker Compose
- NVIDIA GPU 및 NVIDIA Container Toolkit (권장)

## 디렉토리 구조

```
.
├── Dockerfile                # 도커 이미지 빌드 파일
├── docker-compose.yml        # 도커 컴포즈 설정 파일
├── mdp_templates/            # Gromacs MDP 템플릿 파일
│   ├── ions.mdp              # 이온화 설정
│   ├── em.mdp                # 에너지 최소화 설정
│   └── md.mdp                # MD 시뮬레이션 설정
├── scripts/                  # 스크립트 파일
│   ├── sumd_gromacs.py       # 주요 SuMD 자동화 스크립트
│   ├── energy_analysis.py    # 에너지 분석 스크립트
│   ├── process_cystein.py    # 시스테인 잔기 처리 스크립트
│   └── run_sumd_gromacs.sh   # 실행 스크립트
├── test_data/                # 테스트 데이터 디렉토리
└── output/                   # 결과 출력 디렉토리
```

## 빠른 시작

1. **환경 구축**
   ```bash
   # 저장소 클론 또는 디렉토리 생성 후 파일 복사
   git clone <repository-url>
   cd sumd-gromacs
   
   # Docker 이미지 빌드 및 컨테이너 실행
   docker-compose up -d --build
   ```

2. **시뮬레이션 실행**

   ```
   # 기본 체인 사용 (첫 번째 체인)
   docker exec -it sumd-gromacs /app/scripts/run_sumd_gromacs.sh \
    /app/test_data/example.pdb \
    "45 46 47" \
    "123 124 125" \
    /app/output \
    0.5 \
    true
   ```
```
# 체인 A를 사용하여 시뮬레이션 실행
docker exec -it sumd-gromacs /app/scripts/run_sumd_gromacs.sh \
    /app/test_data/example.pdb \
    "45 46 47" \
    "123 124 125" \
    /app/output \
    0.5 \
    true \
    "A"
```

```
# 체인 B를 사용하여 시뮬레이션 실행
docker exec -it sumd-gromacs /app/scripts/run_sumd_gromacs.sh \
    /app/test_data/example.pdb \
    "45 46 47" \
    "123 124 125" \
    /app/output \
    0.5 \
    true \
    "B"
```
3. **결과 확인**
   ```bash
   # 결과 파일 확인
   ls -la output/
   ```

## 매개변수 설명

```
/app/scripts/run_sumd_gromacs.sh <PDB_파일> <펩타이드_잔기_번호> <단백질_잔기_번호> [출력_디렉토리] [거리_임계값] [시스테인_처리]
```

- `<PDB_파일>`: 단백질-펩타이드 복합체 PDB 파일 경로
- `<펩타이드_잔기_번호>`: 펩타이드 잔기 번호 (공백으로 구분)
- `<단백질_잔기_번호>`: 단백질 잔기 번호 (공백으로 구분)
- `[출력_디렉토리]`: 결과 저장 디렉토리 (기본값: ./output)
- `[거리_임계값]`: 무게중심 거리 임계값 (nm, 기본값: 0.5)
- `[시스테인_처리]`: 시스테인 잔기에 HG 원자 추가 여부 (기본값: true)

## 성능 최적화

1. **시뮬레이션 매개변수 최적화**
   - dt=0.004, nsteps=2500000 설정으로 10ns 시뮬레이션 빠르게 수행
   - LINCS 알고리즘 매개변수 조정 (lincs_iter=2, lincs_order=6)
   - 출력 빈도 최적화 (nstenergy=5000, nstlog=5000, nstxout-compressed=5000)

2. **하드웨어 가속**
   - NVIDIA GPU 활용
   - 멀티 GPU 지원 (NVIDIA_VISIBLE_DEVICES 환경 변수 설정)

3. **메모리 사용량 최적화**
   - nstlist, rlist 매개변수 조정으로 메모리 사용량 조절 가능

## 주요 프로세스

1. **시스테인 잔기 처리**
   - PDB 파일에서 시스테인 잔기 식별
   - 필요한 경우 HG 원자 추가 (SG-HG 결합)
   - 처리된 PDB 파일 생성

2. **시스템 준비**
   - PDB 파일에서 토폴로지 생성
   - 시뮬레이션 박스 생성 및 솔베이션
   - 이온 추가로 전하 중화

3. **거리 기반 분기 처리**
   - 펩타이드와 단백질 잔기 간 무게중심 거리 계산
   - 에너지 최소화(EM) 수행
   - 거리 임계값 기준 MD 시뮬레이션 실행 여부 결정

4. **에너지 분석**
   - Coulomb 및 LJ 상호작용 에너지 계산
   - 시간에 따른 에너지 변화 그래프 생성
   - 결합 강도 평가 (강, 중, 약)

## 시스테인 잔기 처리

시스테인 잔기 처리는 다음과 같은 방식으로 수행됩니다:

1. **HG 원자 추가**: SG 원자에 연결된 HG 원자가 없을 경우 자동으로 추가합니다.
2. **위치 계산**: CA와 SG 원자를 이용하여 적절한 HG 원자 위치를 계산합니다.
3. **결합 길이**: S-H 결합 길이는 1.33 Å을 사용합니다.

시스테인 잔기 처리는 다음 명령어로 수동으로 실행할 수도 있습니다:

```bash
python3 /app/scripts/process_cystein.py --input input.pdb --output output.pdb
```

## 결과 파일

- `processed_cys.pdb`: 시스테인 처리된 PDB 파일
- `em.gro`: 에너지 최소화 후 구조
- `md.gro`: MD 시뮬레이션 후 최종 구조 (MD 실행 시)
- `md.xtc`: MD 트라젝토리 (MD 실행 시)
- `interaction_energy.png`: 상호작용 에너지 그래프 (MD 실행 시)
- `interaction_energy.csv`: 상호작용 에너지 데이터 (MD 실행 시)
- `energy_summary.txt`: 에너지 분석 요약 (MD 실행 시)

## 시뮬레이션 시간 조정

MD 시뮬레이션 시간은 mdp_templates/md.mdp 파일에서 조정 가능합니다:
- `dt * nsteps = 총 시뮬레이션 시간 (ps)`
- 예: 0.004 ps * 2,500,000 = 10,000 ps = 10 ns

더 짧은 시뮬레이션을 위해 nsteps 값을 줄일 수 있습니다:
- 5 ns: nsteps = 1,250,000
- 2 ns: nsteps = 500,000
- 1 ns: nsteps = 250,000

## 트러블슈팅

1. **메모리 부족 오류**
   - mdp_templates/md.mdp 파일에서 nstlist 값 증가
   - 더 작은 시뮬레이션 박스 사용 (editconf에서 -d 값 감소)

2. **GPU 관련 오류**
   - NVIDIA 드라이버 및 CUDA 설치 확인
   - docker-compose.yml의 GPU 관련 설정 점검

3. **시뮬레이션 충돌 문제**
   - dt 값을 줄여 시간 간격 감소 (0.003 또는 0.002)
   - constraints 설정 변경 (h-bonds에서 all-bonds로)

4. **시스테인 잔기 처리 오류**
   - PDB 파일의 시스테인 잔기 형식 확인
   - 처리 오류 시 시스테인 처리 단계 건너뛰기 (`false` 옵션 사용)

## 추가 자원

- [Gromacs 공식 문서](http://manual.gromacs.org/)
- [SuMD 관련 논문](https://pubs.acs.org/doi/10.1021/ci400552t)
- CHARMM36m 포스필드: inhyeoksong/gromacs-new Docker 이미지에 포함됨
