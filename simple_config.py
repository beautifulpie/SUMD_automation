# Simple SuMD 설정 파일 (체인 복원 기능 포함)

# ===== 시뮬레이션 기본 설정 =====
MAX_ITERATIONS = 2          # 최대 iteration 수
MAX_ATTEMPTS = 10           # iteration당 최대 attempt 수
SIMULATION_TIME_NS = 0.3     # MD 시뮬레이션 시간 (ns) - 300ps
ENABLE_CHAIN_RESTORATION = True

# 기울기 임계값 (음수여야 채택)
SLOPE_THRESHOLD = -0.001

# ===== 복합체 분석 설정 =====
BINDING_SITE_CUTOFF = 4.0        # Binding site 정의 거리 임계값 (Å)
SEPARATION_DISTANCE = 30.0       # 구조 이격 거리 (Å)
AUTO_DETECT_LIGAND_RECEPTOR = True  # 자동 ligand/receptor 판별

# ===== Binding site 설정 =====
BINDING_SITE_RESIDUES = []       # 수동으로 지정된 binding site (빈 리스트면 자동 탐지)

# ===== 구조 이격 설정 =====
ENABLE_MULTI_DIRECTION_SEPARATION = True  # 다방향 이격 활성화
SEPARATION_DISTANCE = 30.0                # 이격 거리 (Angstrom)
CONE_ANGLES = [60,120]                           # 기준 벡터에서 기울기 [] 배열 입력.
ENABLE_RANDOM_DIRECTION=True # 다방향 이격 랜덤 활성화
ROTATION_STEP = 1                  # 이격 평면 회전 각도 (degrees) ENABLE_RANDOM_DIRECTION False 시 활용
MAX_SEPARATION_VARIANTS=1           # 구조 이격 최대 수 (con_angle당) + 1(기준이격구조) = 최종 구조 수. ENABLE_RANDOM_DIRECTION True 시 활용

# ===== 회전 변형 설정 =====
ENABLE_ROTATIONAL_VARIANTS = True  # 회전 변형 활성화
MAX_ROTATION_VARIANTS = 2           # 최대 회전 변형 수
STRUCTURE_ROTATION = 60

# ===== GROMACS 설정 =====
GPU_ID = "1"                 #  사용할 GPU ID
MPI_RANKS = "4"              # MPI 랭크 수
NTOMP = "2"                  # OpenMP 스레드 수
FORCE_FIELD = "charmm36-jul2022"  # Force field
WATER_MODEL = "tip3p"        # 물 모델
TIMEOUT_GROMACS = 3600       # 기본 GROMACS 명령어 타임아웃
NPME = "1"

# ===== 시스템 설정 =====
BOX_DISTANCE = 1.5           # 박스 거리 (nm)
MAX_WARNINGS = 1             # GROMACS 최대 경고 수

# ===== 긴 MD 설정 =====
CLOSE_DISTANCE_THRESHOLD = 10.0  # 긴 MD 실행 거리 임계값 (Å)
LONG_MD_TIME_NS = 10.0           # 긴 MD 시뮬레이션 시간 (ns)
ENABLE_LONG_MD = True            # 긴 MD 기능 활성화
TIMEOUT_LONG_MD = 72000           # 긴 MD 타임아웃 (초)

# ===== 로그 설정 =====
VERBOSE = True               # 상세 로그 출력
KEEP_FAILED_ATTEMPTS = False # 실패한 attempt 디렉토리 보존

# ===== MDP 템플릿 설정 =====
MDP_SETTINGS = {
    "em": {
        "integrator": "steep",
        "nsteps": 50000,
        "emtol": 100.0,
        "emstep": 0.01
    },
    "nvt": {
        "integrator": "md",
        "dt": 0.002,
        "nsteps": 25000,  # 50ps
        "temperature": 300
    },
    "npt": {
        "integrator": "md", 
        "dt": 0.002,
        "nsteps": 25000,  # 50ps
        "temperature": 300,
        "pressure": 1.0
    },
    "npt2": {
        "integrator": "md", 
        "dt": 0.002,
        "nsteps": 25000,  # 50ps
        "temperature": 300,
        "pressure": 1.0
    },
    "md": {
        "integrator": "sd",
        "dt": 0.002,
        "temperature": 300,
        "pressure": 1.0
    },
    "long_md": {
        "integrator": "sd",
        "dt": 0.002,
        "temperature": 300,
        "pressure": 1.0,
        "output_frequency_multiplier": 10  # 긴 MD에서는 출력 빈도 증가
    }
}

# ===== 출력 빈도 설정 =====
OUTPUT_FREQUENCY = {
    "energy": 5000,
    "log": 5000, 
    "trajectory": 5000
}

# ===== 시스템 크기별 권장 설정 =====
"""
소형 시스템 (<3000 원자):
- ENABLE_CHAIN_RESTORATION = True
- LONG_MD_TIME_NS = 5.0
- MPI_RANKS = "2"

중형 시스템 (3000-10000 원자):
- ENABLE_CHAIN_RESTORATION = True  
- LONG_MD_TIME_NS = 10.0
- MPI_RANKS = "4"

대형 시스템 (>10000 원자):
- ENABLE_CHAIN_RESTORATION = True
- LONG_MD_TIME_NS = 5.0  # 시간 단축
- MPI_RANKS = "4"
"""