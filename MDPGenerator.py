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
        if simulation_time_ns <= 0.5:  # 500ps 이하인 경우 (300ps 포함)
            if system_size == "large":
                base_steps = max(500, min(2000, int(simulation_time_ns * 1000)))  # 더 빠른 EM
            elif system_size == "small":
                base_steps = max(200, min(1000, int(simulation_time_ns * 500)))
            else:  # medium
                base_steps = max(300, min(1500, int(simulation_time_ns * 800)))
        else:
            # 기존 로직 유지 (긴 시뮬레이션용)
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
        
        # 300ps의 경우 출력 빈도를 더 자주 설정
        if simulation_time_ns <= 0.5:  # 500ps 이하
            nstenergy = min(1000, max(100, nsteps // 150))    # 약 150회 출력
            nstlog = min(1000, max(100, nsteps // 150))
            nstxout_compressed = min(1000, max(100, nsteps // 100))  # 약 100회 출력
        else:
            # 기존 설정 유지
            nstenergy = 5000
            nstlog = 5000
            nstxout_compressed = 5000
        
        mdp_content = f"""; 300ps 최적화 분자 시뮬레이션 설정 (시뮬레이션 시간: {simulation_time_ns} ns)
integrator              = md
dt                      = {MDPGenerator.TIMESTEP}
nsteps                  = {nsteps}
nstenergy               = {nstenergy}
nstlog                  = {nstlog}
nstxout-compressed      = {nstxout_compressed}

; 온도 커플링 (300ps에 맞게 더 빠른 반응)
tcoupl                  = V-rescale
tc-grps                 = System
tau_t                   = {'0.05' if simulation_time_ns <= 0.5 else '0.1'}
ref_t                   = 300

; 압력 커플링 (300ps에 맞게 조정)
pcoupl                  = {'Berendsen' if simulation_time_ns <= 0.5 else 'Parrinello-Rahman'}
pcoupltype              = isotropic
tau_p                   = {'1.0' if simulation_time_ns <= 0.5 else '2.0'}
ref_p                   = 1.0
compressibility         = 4.5e-5

; 결합 제약
constraints             = h-bonds
constraint_algorithm    = LINCS
lincs_iter              = 2
lincs_order             = 6

; 비결합 상호작용
cutoff-scheme           = Verlet
nstlist                 = {'20' if simulation_time_ns <= 0.5 else '40'}
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
        
        # 300ps 시뮬레이션의 경우 더 관대한 emtol 사용 (빠른 수렴)
        if simulation_time_ns <= 0.5:  # 500ps 이하
            if system_size == "large":
                emtol = max(1000.0, min(3000.0, 2000.0 / simulation_time_ns))
            elif system_size == "small":
                emtol = max(200.0, min(800.0, 600.0 / simulation_time_ns))
            else:  # medium
                emtol = max(500.0, min(1500.0, 1200.0 / simulation_time_ns))
        else:
            # 기존 로직 유지
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
        
        # 300ps의 경우 더 관대한 설정
        if simulation_time_ns <= 0.5:  # 500ps 이하
            emtol = max(500.0, min(1500.0, 1000.0 / simulation_time_ns))
        else:
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