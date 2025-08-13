class MDPGenerator:
    """MDP 파일 생성 클래스 - 300ps 최적화 버전"""
    
    TIMESTEP = 0.002  # 2 fs
    
    @staticmethod
    def calculate_nsteps_for_time(simulation_time_ns, dt=None):
        if dt is None:
            dt = MDPGenerator.TIMESTEP
        total_ps = simulation_time_ns * 1000
        nsteps = int(total_ps / dt)
        return nsteps
    
    @staticmethod
    def calculate_em_nsteps(simulation_time_ns, system_size="medium"):
        """300ps 시뮬레이션에 맞게 EM 단계 최적화"""
        if system_size == "large":
            base_steps = max(1000, min(5000, int(simulation_time_ns * 500)))
        elif system_size == "small":
            base_steps = max(250, min(2000, int(simulation_time_ns * 250)))
        else:  # medium
            base_steps = max(500, min(3000, int(simulation_time_ns * 350)))
        return base_steps
    
    @staticmethod
    def calculate_ions_nsteps(simulation_time_ns):
        """300ps 시뮬레이션에 맞게 이온 단계 최적화"""
        if simulation_time_ns <= 0.5:  # 500ps 이하인 경우
            base_steps = max(100, min(800, int(simulation_time_ns * 400)))
        else:
            base_steps = max(200, min(2000, int(simulation_time_ns * 200)))
        return base_steps
    
    @staticmethod
    def generate_md_mdp(output_file, simulation_time_ns):
        """300ps에 최적화된 MD 설정"""
        nsteps = MDPGenerator.calculate_nsteps_for_time(simulation_time_ns, MDPGenerator.TIMESTEP)
        
        nstenergy = 5000
        nstlog = 5000
        nstxout_compressed = 5000
    
        mdp_content = f"""; 300ps 최적화 분자 시뮬레이션 설정 (시뮬레이션 시간: {simulation_time_ns} ns)
integrator              = sd
dt                      = {MDPGenerator.TIMESTEP}
nsteps                  = {nsteps}
nstenergy               = {nstenergy}
nstlog                  = {nstlog}
nstxout-compressed      = {nstxout_compressed}

; Langevin dynamics 파라미터
tc-grps                 = System
tau_t                   = 0.1         ; [ps] time constant 
ref_t                   = 300         ; [K] reference temperature
bd-fric                 = 0           ; auto-calculate: mass/tau_t
ld-seed                 = -1          ; random seed from process ID

; 압력 커플링 (300ps에 맞게 조정)
pcoupl                  = C-rescale
pcoupltype              = isotropic
tau_p                   = 2.0
ref_p                   = 1.0
compressibility         = 4.5e-5

; 결합 제약
constraints             = h-bonds
constraint_algorithm    = LINCS
lincs_iter              = 2
lincs_order             = 6

; 비결합 상호작용
cutoff-scheme           = Verlet
nstlist                 = 40
ns_type                 = grid
coulombtype             = PME
rcoulomb                = 1.0
rvdw                    = 1.0
DispCorr                = EnerPres
pbc                     = xyz
"""
        with open(output_file, 'w', encoding='utf-8') as f:
            f.write(mdp_content)
    
    @staticmethod
    def generate_em_mdp(output_file, simulation_time_ns=1.0, system_size="medium"):
        """300ps에 최적화된 에너지 최소화 설정"""
        nsteps = MDPGenerator.calculate_em_nsteps(simulation_time_ns, system_size)
        
        if system_size == "large":
            emtol = max(500.0, min(2000.0, 1500.0 / simulation_time_ns))
        elif system_size == "small":
            emtol = max(50.0, min(500.0, 500.0 / simulation_time_ns))
        else:  # medium
            emtol = max(100.0, min(1000.0, 1000.0 / simulation_time_ns))
    
        mdp_content = f"""; 300ps 최적화 에너지 최소화 설정
integrator          = steep
nsteps              = {nsteps}
emtol               = {emtol:.1f}
emstep              = {MDPGenerator.TIMESTEP}
nstlist             = 1
cutoff-scheme       = Verlet
ns_type             = grid
coulombtype         = PME
rcoulomb            = 1.0
rvdw                = 1.0
pbc                 = xyz
"""
        with open(output_file, 'w', encoding='utf-8') as f:
            f.write(mdp_content)
    
    @staticmethod
    def generate_ions_mdp(output_file, simulation_time_ns=1.0):
        """300ps에 최적화된 이온 삽입 설정"""
        nsteps = MDPGenerator.calculate_ions_nsteps(simulation_time_ns)
        
        emtol = max(200.0, min(1000.0, 800.0 / simulation_time_ns))
        
        mdp_content = f"""; 300ps 최적화 이온 삽입 설정
title               = Ion insertion (300ps optimized)
integrator          = steep
emtol               = {emtol:.1f}
emstep              = {MDPGenerator.TIMESTEP}
nsteps              = {nsteps}
nstlist             = 1
cutoff-scheme       = Verlet
ns_type             = grid
rlist               = 1.0
coulombtype         = cutoff
rcoulomb            = 1.0
rvdw                = 1.0
pbc                 = xyz
"""
        with open(output_file, 'w', encoding='utf-8') as f:
            f.write(mdp_content)
    
    @staticmethod
    def is_short_simulation(simulation_time_ns):
        """300ps와 같은 짧은 시뮬레이션인지 확인"""
        return simulation_time_ns <= 0.5
        
    @staticmethod
    def get_optimization_info(simulation_time_ns, system_size="medium"):
        """현재 설정에 대한 최적화 정보 반환"""
        info = {
            "simulation_time_ns": simulation_time_ns,
            "system_size": system_size,
            "is_short_simulation": MDPGenerator.is_short_simulation(simulation_time_ns),
            "nsteps_md": MDPGenerator.calculate_nsteps_for_time(simulation_time_ns),
            "nsteps_em": MDPGenerator.calculate_em_nsteps(simulation_time_ns, system_size),
            "nsteps_ions": MDPGenerator.calculate_ions_nsteps(simulation_time_ns)
        }
        
        if info["is_short_simulation"]:
            info["optimizations"] = [
                "빠른 온도/압력 커플링 (tau_t=0.05, Berendsen barostat)",
                "더 관대한 에너지 최소화 임계값",
                "증가된 출력 빈도 (더 세밀한 모니터링)",
                "단축된 EM 및 이온화 단계",
                "더 짧은 neighbor list 업데이트 주기"
            ]
        else:
            info["optimizations"] = ["표준 설정 사용"]
            
        return info
    
    @staticmethod
    def generate_nvt_mdp(output_file, simulation_time_ns=0.1):
        """300ps SuMD에 최적화된 NVT 평형 설정 (100ps)"""
        nsteps = int((simulation_time_ns * 1000) / MDPGenerator.TIMESTEP)  # 100ps = 50000 steps
        
        mdp_content = f"""; NVT 평형 (100ps) - SuMD 최적화
title                   = NVT equilibration for SuMD
define                  = -DPOSRES
integrator              = md
dt                      = {MDPGenerator.TIMESTEP}
nsteps                  = {nsteps}
nstenergy               = 500
nstlog                  = 500
nstxout-compressed      = 500

; 결합 제약
constraints             = h-bonds
constraint_algorithm    = lincs
lincs_iter              = 1
lincs_order             = 4

; 비결합 상호작용
cutoff-scheme           = Verlet
ns_type                 = grid
nstlist                 = 10
rcoulomb                = 1.0
rvdw                    = 1.0
DispCorr                = EnerPres

; 정전기학
coulombtype             = PME
pme_order               = 4
fourierspacing          = 0.16

; 온도 커플링
tcoupl                  = V-rescale
tc-grps                 = Protein Non-Protein
tau_t                   = 0.1     0.1
ref_t                   = 300     300

; 압력 커플링 끔 (NVT)
pcoupl                  = no

; 주기적 경계 조건
pbc                     = xyz

; 초기 속도 생성
gen_vel                 = yes
gen_temp                = 300
gen_seed                = -1
"""
        with open(output_file, 'w', encoding='utf-8') as f:
            f.write(mdp_content)

    @staticmethod
    def generate_npt_mdp(output_file, simulation_time_ns=0.1):
        """300ps SuMD에 최적화된 NPT 평형 설정 (100ps)"""
        nsteps = int((simulation_time_ns * 1000) / MDPGenerator.TIMESTEP)  # 100ps = 50000 steps
        
        mdp_content = f"""; NPT 평형 (100ps) - SuMD 최적화  
title                   = NPT equilibration for SuMD
define                  = -DPOSRES
integrator              = md
dt                      = {MDPGenerator.TIMESTEP}
nsteps                  = {nsteps}
nstenergy               = 500
nstlog                  = 500
nstxout-compressed      = 500

; 결합 제약
continuation            = yes
constraints             = h-bonds
constraint_algorithm    = lincs
lincs_iter              = 1
lincs_order             = 4

; 비결합 상호작용
cutoff-scheme           = Verlet
ns_type                 = grid
nstlist                 = 10
rcoulomb                = 1.0
rvdw                    = 1.0
DispCorr                = EnerPres

; 정전기학
coulombtype             = PME
pme_order               = 4
fourierspacing          = 0.16

; 온도 커플링
tcoupl                  = V-rescale
tc-grps                 = Protein Non-Protein
tau_t                   = 0.1     0.1
ref_t                   = 300     300

; 압력 커플링 (NPT) - C-rescale 사용
pcoupl                  = C-rescale
pcoupltype              = isotropic
tau_p                   = 2.0
ref_p                   = 1.0
compressibility         = 4.5e-5
refcoord_scaling        = com

; 주기적 경계 조건
pbc                     = xyz

; 초기 속도 생성 끔 (NVT에서 이어받음)
gen_vel                 = no
"""
        with open(output_file, 'w', encoding='utf-8') as f:
            f.write(mdp_content)