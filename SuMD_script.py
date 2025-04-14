import subprocess
import argparse
import os

def run_command(cmd, input_str=None):
    print(f"Running command: {' '.join(cmd)}")
    subprocess.run(cmd, input=input_str.encode() if input_str else None, check=True)

def create_short_mdp(filename="short.mdp"):
    if not os.path.exists(filename):
        mdp_content = """
integrator  = md
nsteps      = 2500      ; 5 ps with 2 fs timestep
dt          = 0.002
nstxout     = 500
nstvout     = 500
nstenergy   = 500
nstlog      = 500
continuation    = no
constraint_algorithm = lincs
constraints = all-bonds
lincs_iter  = 1
lincs_order = 4
; Temperature coupling
tcoupl      = V-rescale
tc-grps     = System
tau_t       = 0.1
ref_t       = 300
; Pressure coupling
pcoupl      = Parrinello-Rahman
pcoupltype  = isotropic
tau_p       = 2.0
ref_p       = 1.0
compressibility = 4.5e-5
; Periodic boundary conditions
pbc         = xyz
; Dispersion correction
DispCorr    = EnerPres
"""
        with open(filename, "w") as f:
            f.write(mdp_content)
        print(f"Created default {filename}")

def create_emin_mdp(filename="emin.mdp"):
    if not os.path.exists(filename):
        mdp_content = """
integrator  = steep
nsteps      = 5000
emtol       = 1000.0
energygrps  = System
"""
        with open(filename, "w") as f:
            f.write(mdp_content)
        print(f"Created default {filename}")

def prepare_system(pdb_file, peptide_res_list, protein_res_list):
    # (1) pdb2gmx: AMBER14SB force field와 TIP3P water 모델 사용
    run_command(["gmx", "pdb2gmx", "-f", pdb_file, "-o", "processed.gro", "-p", "topol.top", "-ff", "charmm27", "-water", "tip3p", "-ignh"])
    
    # (2) editconf: 시뮬레이션 박스 생성 (중심 배치, 1.0 nm 여유, cubic)
    run_command(["gmx", "editconf", "-f", "processed.gro", "-o", "newbox.gro", "-c", "-d", "1.0", "-bt", "cubic"])

    # (3) solvate: 물 추가
    run_command(["gmx", "solvate", "-cp", "newbox.gro", "-cs", "spc216.gro", "-o", "solvated.gro", "-p", "topol.top"])

    # (4) make_ndx: peptide와 protein의 결합하는 residue 선택 → union 그룹 생성
    # 주의: "r 45" 등으로 선택하면 생성된 그룹 이름은 보통 "45"가 됩니다.
    index_commands_lines = []
    for res in peptide_res_list:
        index_commands_lines.append(f"r {res}")
    # union: peptide_grp = 45 | 46 | 47 (입력 번호에 따라)
    index_commands_lines.append("peptide_grp = " + " | ".join(f"{res}" for res in peptide_res_list))
    
    for res in protein_res_list:
        index_commands_lines.append(f"r {res}")
    index_commands_lines.append("protein_grp = " + " | ".join(f"{res}" for res in protein_res_list))
    
    index_commands_lines.append("q")
    index_commands = "\n".join(index_commands_lines) + "\n"
    
    run_command(["gmx", "make_ndx", "-f", "solvated.gro", "-o", "index.ndx"], input_str=index_commands)

def run_energy_minimization(iteration, start_coord):
    emin_prefix = f"emin_{iteration}"
    # emin.mdp를 이용하여 에너지 최소화 시뮬레이션 준비
    run_command(["gmx", "grompp", "-f", "emin.mdp", "-c", start_coord, "-p", "topol.top", "-n", "index.ndx", "-o", f"{emin_prefix}.tpr"])
    # 에너지 최소화 실행
    run_command(["gmx", "mdrun", "-deffnm", emin_prefix])
    # 최종 구조는 emin_<iteration>.gro 파일로 출력됨
    return f"{emin_prefix}.tpr", f"{emin_prefix}.gro"

def run_segment(segment_number, start_coord):
    seg_prefix = f"segment_{segment_number}"
    run_command(["gmx", "grompp", "-f", "short.mdp", "-c", start_coord, "-p", "topol.top", "-n", "index.ndx", "-o", f"{seg_prefix}.tpr"])
    run_command(["gmx", "mdrun", "-deffnm", seg_prefix])
    # 반환: tpr, 최종 구조 gro, trajectory xtc
    return f"{seg_prefix}.tpr", f"{seg_prefix}.gro", f"{seg_prefix}.xtc"

def check_distance(sim_file, tpr_file, peptide_res_list, protein_res_list):
    """
    –select 옵션에서 inline으로 여러 residue를 union하여 계산하도록 합니다.
    예를 들어, peptide_res_list = ['45','46','47']인 경우,
      peptide_sel = "r 45 | r 46 | r 47"
    최종 선택식:
      "com of ( r 45 | r 46 | r 47 ) plus com of ( r 123 | r 124 | r 125 )"
    """
    peptide_sel = " | ".join(f"r {res}" for res in peptide_res_list)
    protein_sel = " | ".join(f"r {res}" for res in protein_res_list)
    select_string = f"com of ( {peptide_sel} ) plus com of ( {protein_sel} )"
    
    distance_output = "distance.xvg"
    run_command([
        "gmx", "distance",
        "-s", tpr_file,
        "-f", sim_file,
        "-n", "index.ndx",
        "-select", select_string,
        "-oall", distance_output
    ])
    
    final_distance = None
    with open(distance_output, "r") as f:
        lines = f.readlines()
        for line in reversed(lines):
            if not line.startswith(("#", "@")):
                parts = line.strip().split()
                if len(parts) >= 2:
                    final_distance = float(parts[1])
                    break
    if final_distance is None:
        raise RuntimeError("Failed to extract distance from distance.xvg")
    
    print(f"Final distance between peptide and protein groups: {final_distance} nm")
    return final_distance

def main():
    parser = argparse.ArgumentParser(
        description="Gromacs와 AMBER14SB force field를 이용한 peptide-protein SuMD 시뮬레이션 자동화\n"
                    "여러 residue 번호를 공백으로 구분해서 입력할 수 있습니다."
    )
    parser.add_argument("--pdb", required=True, help="Complex PDB file")
    parser.add_argument("--peptide_res", nargs="+", required=True,
                        help="Peptide의 결합하는 residue 번호들 (여러 개 가능)")
    parser.add_argument("--protein_res", nargs="+", required=True,
                        help="Protein의 결합하는 residue 번호들 (여러 개 가능)")
    parser.add_argument("--max_emin_iter", type=int, default=5,
                        help="최대 Energy Minimization 반복 횟수")
    parser.add_argument("--max_segments", type=int, default=10,
                        help="최대 MD 세그먼트 수")
    parser.add_argument("--distance_threshold", type=float, default=0.5,
                        help="Cutoff 거리 (nm); 기본 0.5 nm (5 Å)")
    args = parser.parse_args()
    
    # MD 및 에너지 최소화용 mdp 파일 생성
    create_short_mdp("short.mdp")
    create_emin_mdp("emin.mdp")
    
    # 시스템 준비: Topology 생성, 박스, 솔베이션, 기본 index 파일 생성
    prepare_system(args.pdb)
    
    # 초기 구조: solvated.gro 파일을 사용
    current_coord = "solvated.gro"
    
    # Energy Minimization 단계: peptide와 protein의 COM 사이 거리가 cutoff (0.5 nm) 이하가 될 때까지 반복
    emin_iter = 0
    current_distance = None
    while emin_iter < args.max_emin_iter:
        emin_iter += 1
        print(f"\n--- Energy Minimization Iteration {emin_iter} (4 samples) ---")

        sample_results = []
        for sample_id in range(1, 5):  # 4개의 sample 실행
            tag = f"{emin_iter}_{sample_id}"
            tpr_file, coord_file = run_energy_minimization(tag, current_coord)
            try:
                dist = check_distance(coord_file, tpr_file, args.peptide_res, args.protein_res)
                sample_results.append((dist, coord_file, tpr_file))
            except Exception as e:
                print(f"Sample {tag} 실패: {e}")

        if not sample_results:
            print("모든 에너지 최소화 sample이 실패했습니다. 시뮬레이션을 종료합니다.")
            break

        # 가장 짧은 거리의 sample 선택
        sample_results.sort(key=lambda x: x[0])  # distance 기준 정렬
        best_distance, best_coord, best_tpr = sample_results[0]

        print(f"가장 짧은 거리: {best_distance:.3f} nm (선택된 sample)")

        if best_distance <= args.distance_threshold:
            print("Cutoff 도달! Energy Minimization 단계 완료.")
            current_coord = best_coord
            break
        else:
            print("Cutoff에 도달하지 않았으므로, 에너지 최소화를 반복합니다.")
            current_coord = best_coord  # 다음 iteration의 시작 구조
    else:
        print("최대 Energy Minimization 반복 횟수에 도달하였습니다.")

    
    # Energy Minimization 후 cutoff 상태에서 MD 세그먼트 실행 (최종 pose 형성 탐색)
    previous_distance = None
    for segment in range(1, args.max_segments + 1):
        print(f"\n--- Running MD Segment {segment} ---")
        md_tpr, md_coord, md_xtc = run_segment(segment, current_coord)
        current_distance = check_distance(md_xtc, md_tpr, args.peptide_res, args.protein_res)
        current_coord = md_coord  # 다음 세그먼트 시작 좌표로 업데이트
        
        if current_distance < args.distance_threshold:
            print("최종 pose가 cutoff 내에 도달하였습니다. MD 시뮬레이션 완료.")
            break
        
        if previous_distance is not None and current_distance >= previous_distance:
            print("거리 감소가 없으므로 MD 시뮬레이션을 종료합니다.")
            break
        previous_distance = current_distance

if __name__ == "__main__":
    main()