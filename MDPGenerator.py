

class MDPGenerator:
    """MDP 파일 생성 클래스"""
    
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
        if system_size == "large":
            base_steps = max(1000, min(5000, int(simulation_time_ns * 500)))
        elif system_size == "small":
            base_steps = max(250, min(2000, int(simulation_time_ns * 250)))
        else:  # medium
            base_steps = max(500, min(3000, int(simulation_time_ns * 350)))
        return base_steps
    
    @staticmethod
    def calculate_ions_nsteps(simulation_time_ns):
        base_steps = max(200, min(2000, int(simulation_time_ns * 200)))
        return base_steps
    
    @staticmethod
    def generate_md_mdp(output_file, simulation_time_ns):
        nsteps = MDPGenerator.calculate_nsteps_for_time(simulation_time_ns, MDPGenerator.TIMESTEP)
        
        mdp_content = f"""; 분자 시뮬레이션 설정 (시뮬레이션 시간: {simulation_time_ns} ns)
integrator              = md
dt                      = {MDPGenerator.TIMESTEP}
nsteps                  = {nsteps}
nstenergy               = 5000
nstlog                  = 5000
nstxout-compressed      = 5000

; 온도 커플링
tcoupl                  = V-rescale
tc-grps                 = System
tau_t                   = 0.1
ref_t                   = 300

; 압력 커플링
pcoupl                  = Parrinello-Rahman
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
        nsteps = MDPGenerator.calculate_em_nsteps(simulation_time_ns, system_size)
        
        if system_size == "large":
            emtol = max(500.0, min(2000.0, 1500.0 / simulation_time_ns))
        elif system_size == "small":
            emtol = max(50.0, min(500.0, 500.0 / simulation_time_ns))
        else:  # medium
            emtol = max(100.0, min(1000.0, 1000.0 / simulation_time_ns))
        
        mdp_content = f"""; 에너지 최소화 설정
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
        nsteps = MDPGenerator.calculate_ions_nsteps(simulation_time_ns)
        emtol = max(200.0, min(1000.0, 800.0 / simulation_time_ns))
        
        mdp_content = f"""; 이온 삽입 설정
title               = Ion insertion
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
