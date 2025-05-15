import subprocess
import argparse
import os
import numpy as np
import shutil
import time
from pathlib import Path

def run_command(cmd, input_str=None, verbose=True):
    """
    명령어를 실행하고 결과를 반환합니다.
    verbose가 True이면 실행 명령어를 출력합니다.
    """
    if verbose:
        print(f"실행 명령어: {' '.join(cmd)}")
    
    try:
        result = subprocess.run(
            cmd, 
            input=input_str.encode() if input_str else None, 
            check=True,
            capture_output=True,
            text=True
        )
        if verbose and result.stdout:
            print(result.stdout)
        return result
    except subprocess.CalledProcessError as e:
        print(f"명령어 실행 오류: {e}")
        if e.stdout:
            print(f"표준 출력: {e.stdout}")
        if e.stderr:
            print(f"에러 출력: {e.stderr}")
        raise

def create_short_mdp(filename="short.mdp"):
    """짧은 MD 시뮬레이션을 위한 mdp 파일 생성"""
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
        print(f"{filename} 기본 파일 생성됨")

def create_emin_mdp(filename="emin.mdp"):
    """에너지 최소화를 위한 mdp 파일 생성"""
    if not os.path.exists(filename):
        mdp_content = """
integrator  = steep
nsteps      = 5000
emtol       = 1000.0
emstep      = 0.01
nstlist     = 1
nstenergy   = 100
nstxout     = 100
pbc         = xyz
rlist       = 1.0
coulombtype = PME
rcoulomb    = 1.0
vdwtype     = Cut-off
rvdw        = 1.0
"""
        with open(filename, "w") as f:
            f.write(mdp_content)
        print(f"{filename} 기본 파일 생성됨")

def prepare_system(pdb_file, peptide_res_list, protein_res_list, force_field="amber99sb-ildn", water_model="tip3p"):
    """
    GROMACS를 사용하여 시스템 준비: 
    1. PDB 처리
    2. 박스 생성
    3. 솔베이션
    4. 이온 추가
    5. Index 그룹 생성
    """
    print(f"\n=== 시스템 준비 시작: {pdb_file} ===")
    
    # (1) pdb2gmx: 힘장(force field)과 물 모델 사용
    run_command(["gmx", "pdb2gmx", 
                 "-f", pdb_file, 
                 "-o", "processed.gro", 
                 "-p", "topol.top", 
                 "-ff", force_field, 
                 "-water", water_model, 
                 "-ignh"])
    
    # (2) editconf: 시뮬레이션 박스 생성 (중심 배치, 1.0 nm 여유, cubic)
    run_command(["gmx", "editconf", 
                 "-f", "processed.gro", 
                 "-o", "newbox.gro", 
                 "-c", "-d", "1.0", 
                 "-bt", "cubic"])

    # (3) solvate: 물 추가
    run_command(["gmx", "solvate", 
                 "-cp", "newbox.gro", 
                 "-cs", "spc216.gro", 
                 "-o", "solvated.gro", 
                 "-p", "topol.top"])
    
    # (4) 중성화를 위한 이온 추가
    # 먼저 ions.mdp 생성
    with open("ions.mdp", "w") as f:
        f.write("""
integrator  = steep
nsteps      = 5000
emtol       = 1000.0
emstep      = 0.01
nstlist     = 1
nstenergy   = 100
pbc         = xyz
rlist       = 1.0
coulombtype = PME
rcoulomb    = 1.0
vdwtype     = Cut-off
rvdw        = 1.0
""")
    
    # grompp 실행 (topology 처리)
    run_command(["gmx", "grompp", 
                 "-f", "ions.mdp", 
                 "-c", "solvated.gro", 
                 "-p", "topol.top", 
                 "-o", "ions.tpr"])
    
    # genion 실행 (중성화)
    ion_input = "SOL\n"  # 물 그룹 선택
    run_command(["gmx", "genion", 
                 "-s", "ions.tpr", 
                 "-o", "neutralized.gro", 
                 "-p", "topol.top", 
                 "-pname", "NA", 
                 "-nname", "CL", 
                 "-neutral"], 
                input_str=ion_input)

    # (5) make_ndx: peptide와 protein의 결합하는 residue 선택 → 그룹 생성
    index_commands_lines = []
    
    # Peptide residue 그룹 생성
    for res in peptide_res_list:
        index_commands_lines.append(f"r {res}")
    
    # 모든 peptide residue를 결합
    if len(peptide_res_list) > 1:
        peptide_sel = " | ".join(f"r {res}" for res in peptide_res_list)
        index_commands_lines.append(f"name 0 Peptide_group")
        index_commands_lines.append(peptide_sel)
        index_commands_lines.append("name 1 Peptide_group")
    else:
        index_commands_lines.append(f"r {peptide_res_list[0]}")
        index_commands_lines.append("name 0 Peptide_group")
    
    # Protein residue 그룹 생성
    for res in protein_res_list:
        index_commands_lines.append(f"r {res}")
    
    # 모든 protein residue를 결합
    if len(protein_res_list) > 1:
        protein_sel = " | ".join(f"r {res}" for res in protein_res_list)
        index_commands_lines.append(protein_sel)
        index_commands_lines.append("name 1 Protein_group")
    else:
        index_commands_lines.append(f"r {protein_res_list[0]}")
        index_commands_lines.append("name 1 Protein_group")
    
    # 작업 완료
    index_commands_lines.append("q")
    index_commands = "\n".join(index_commands_lines) + "\n"
    
    run_command(["gmx", "make_ndx", 
                 "-f", "neutralized.gro", 
                 "-o", "index.ndx"], 
                input_str=index_commands)
    
    print("=== 시스템 준비 완료 ===\n")
    return "neutralized.gro"

def run_energy_minimization(iteration, start_coord):
    """에너지 최소화 실행"""
    emin_prefix = f"emin_{iteration}"
    
    # emin.mdp를 이용하여 에너지 최소화 시뮬레이션 준비
    run_command(["gmx", "grompp", 
                 "-f", "emin.mdp", 
                 "-c", start_coord, 
                 "-p", "topol.top", 
                 "-n", "index.ndx", 
                 "-o", f"{emin_prefix}.tpr"])
    
    # 에너지 최소화 실행
    run_command(["gmx", "mdrun", "-deffnm", emin_prefix])
    
    # 최종 구조는 emin_<iteration>.gro 파일로 출력됨
    return f"{emin_prefix}.tpr", f"{emin_prefix}.gro"

def run_segment(segment_number, start_coord):
    """MD 시뮬레이션 세그먼트 실행"""
    seg_prefix = f"segment_{segment_number}"
    
    run_command(["gmx", "grompp", 
                 "-f", "short.mdp", 
                 "-c", start_coord, 
                 "-p", "topol.top", 
                 "-n", "index.ndx", 
                 "-o", f"{seg_prefix}.tpr"])
    
    run_command(["gmx", "mdrun", "-deffnm", seg_prefix])
    
    # 반환: tpr, 최종 구조 gro, trajectory xtc
    return f"{seg_prefix}.tpr", f"{seg_prefix}.gro", f"{seg_prefix}.xtc"

def check_distance(sim_file, tpr_file, output_prefix="distance"):
    """
    Peptide_group와 Protein_group 사이의 중심 간 거리 계산
    index.ndx에서 이미 생성한 그룹을 사용함
    """
    distance_output = f"{output_prefix}.xvg"
    
    # 두 그룹 선택 (이미 index에 정의되어 있음)
    select_input = "Peptide_group\nProtein_group\n"
    
    run_command([
        "gmx", "pairdist",
        "-s", tpr_file,
        "-f", sim_file,
        "-n", "index.ndx",
        "-type", "com",  # center of mass
        "-o", distance_output
    ], input_str=select_input, verbose=False)
    
    # 결과 파일에서 마지막 거리 값 추출
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
        raise RuntimeError(f"거리 파일({distance_output})에서 값을 추출할 수 없습니다.")
    
    print(f"Peptide와 Protein 그룹 간 최종 거리: {final_distance:.3f} nm")
    return final_distance

def save_best_trajectory(best_structures, output_dir="best_trajectories"):
    """최고의 구조들을 저장"""
    os.makedirs(output_dir, exist_ok=True)
    
    for i, (dist, coord_file, tpr_file) in enumerate(best_structures, 1):
        dest_coord = os.path.join(output_dir, f"best_structure_{i}.gro")
        dest_tpr = os.path.join(output_dir, f"best_structure_{i}.tpr")
        
        shutil.copy2(coord_file, dest_coord)
        shutil.copy2(tpr_file, dest_tpr)
        
        with open(os.path.join(output_dir, f"best_structure_{i}_info.txt"), "w") as f:
            f.write(f"Distance: {dist:.3f} nm\nSource files: {coord_file}, {tpr_file}\n")
    
    print(f"\n최상의 구조들이 {output_dir} 디렉토리에 저장되었습니다.")

def run_sumd(args):
    """SuMD 메인 함수"""
    # 작업 디렉토리 생성 (지정된 경우)
    if args.output_dir:
        os.makedirs(args.output_dir, exist_ok=True)
        os.chdir(args.output_dir)
    
    # MD 및 에너지 최소화용 mdp 파일 생성
    create_short_mdp("short.mdp")
    create_emin_mdp("emin.mdp")
    
    # PDB 파일 복사 (필요한 경우)
    pdb_basename = os.path.basename(args.pdb)
    local_pdb = pdb_basename
    if args.pdb != local_pdb:
        shutil.copy2(args.pdb, local_pdb)
    
    # 시스템 준비: Topology 생성, 박스, 솔베이션, 기본 index 파일 생성
    current_coord = prepare_system(local_pdb, args.peptide_res, args.protein_res, 
                                   force_field=args.force_field, water_model=args.water_model)
    
    # 모든 과정에서의 최상 결과를 저장
    all_best_structures = []
    
    # Energy Minimization 단계
    emin_iter = 0
    while emin_iter < args.max_emin_iter:
        emin_iter += 1
        print(f"\n--- Energy Minimization Iteration {emin_iter} (4 samples) ---")
        
        sample_results = []
        for sample_id in range(1, 5):  # 4개의 sample 실행
            tag = f"{emin_iter}_{sample_id}"
            try:
                print(f"  Sample {tag} 진행 중...")
                tpr_file, coord_file = run_energy_minimization(tag, current_coord)
                dist = check_distance(coord_file, tpr_file)
                sample_results.append((dist, coord_file, tpr_file))
                print(f"  Sample {tag} 완료, 거리: {dist:.3f} nm")
            except Exception as e:
                print(f"  Sample {tag} 실패: {e}")
        
        if not sample_results:
            print("모든 에너지 최소화 sample이 실패했습니다. 시뮬레이션을 종료합니다.")
            break
        
        # 가장 짧은 거리의 sample 선택
        sample_results.sort(key=lambda x: x[0])  # distance 기준 정렬
        best_distance, best_coord, best_tpr = sample_results[0]
        all_best_structures.append((best_distance, best_coord, best_tpr))
        
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
    
    # Energy Minimization 후 MD 세그먼트 실행 (최종 pose 형성 탐색)
    print("\n--- MD 세그먼트 시뮬레이션 시작 ---")
    
    previous_distance = None
    for segment in range(1, args.max_segments + 1):
        print(f"\n--- MD 세그먼트 {segment} 실행 중 ---")
        md_tpr, md_coord, md_xtc = run_segment(segment, current_coord)
        current_distance = check_distance(md_xtc, md_tpr)
        all_best_structures.append((current_distance, md_coord, md_tpr))
        
        current_coord = md_coord  # 다음 세그먼트 시작 좌표로 업데이트
        
        if current_distance <= args.distance_threshold:
            print(f"최종 pose가 cutoff({args.distance_threshold} nm) 내에 도달하였습니다. MD 시뮬레이션 완료.")
            break
        
        if previous_distance is not None and current_distance >= previous_distance:
            print("거리 감소가 없으므로 MD 시뮬레이션을 종료합니다.")
            break
        
        previous_distance = current_distance
    
    # 결과 저장
    all_best_structures.sort(key=lambda x: x[0])  # 거리 기준 정렬
    best_structures = all_best_structures[:5]  # 상위 5개
    
    print("\n=== SuMD 시뮬레이션 결과 ===")
    print(f"최종 Peptide-Protein 거리: {all_best_structures[0][0]:.3f} nm")
    print(f"총 {len(all_best_structures)} 개의 샘플 중 상위 5개를 저장합니다.")
    
    save_best_trajectory(best_structures)
    
    return all_best_structures[0][0]  # 최종 거리 반환

def main():
    parser = argparse.ArgumentParser(
        description="GROMACS를 이용한 peptide-protein SuMD(Supervised Molecular Dynamics) 시뮬레이션 자동화"
    )
    parser.add_argument("--pdb", required=True, help="Complex PDB 파일 경로")
    parser.add_argument("--peptide_res", nargs="+", required=True,
                        help="Peptide의 결합하는 residue 번호들 (여러 개 가능)")
    parser.add_argument("--protein_res", nargs="+", required=True,
                        help="Protein의 결합하는 residue 번호들 (여러 개 가능)")
    parser.add_argument("--force_field", default="amber99sb-ildn",
                        help="사용할 힘장(Force field) (기본값: amber99sb-ildn)")
    parser.add_argument("--water_model", default="tip3p",
                        help="사용할 물 모델 (기본값: tip3p)")
    parser.add_argument("--max_emin_iter", type=int, default=5,
                        help="최대 Energy Minimization 반복 횟수 (기본값: 5)")
    parser.add_argument("--max_segments", type=int, default=10,
                        help="최대 MD 세그먼트 수 (기본값: 10)")
    parser.add_argument("--distance_threshold", type=float, default=0.5,
                        help="Cutoff 거리 (nm); 기본 0.5 nm (5 Å)")
    parser.add_argument("--output_dir", help="결과 저장 디렉토리 (기본값: 현재 디렉토리)")
    
    args = parser.parse_args()
    
    # 시작 시간
    start_time = time.time()
    
    try:
        final_distance = run_sumd(args)
        elapsed_time = time.time() - start_time
        print(f"\n총 실행 시간: {elapsed_time:.1f}초")
        print(f"최종 거리: {final_distance:.3f} nm")
        if final_distance <= args.distance_threshold:
            print("SuMD 시뮬레이션 성공!")
        else:
            print("SuMD 목표 거리에 도달하지 못했습니다.")
    except Exception as e:
        print(f"SuMD 시뮬레이션 오류: {e}")
        return 1
    
    return 0

if __name__ == "__main__":
    exit(main())
