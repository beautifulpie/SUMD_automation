import os
import subprocess
from tqdm import tqdm
from os.path import join
import json
import re
import select
import sys
import time

debug_mode = True

configuration = {
    "dataset_path" : "/root/SMC_PPI_MD/Dataset/test/",
    "save_path" : "/root/SMC_PPI_MD/Result/Test/",
    "working_path" : f"/root/SMC_PPI_MD/SuMD/",
    # "csv_path" : "",
    "SuMD_path" : "/root/SMC_PPI_MD/SuMD/",
    "pdb_list_path" : "/root/SMC_PPI_MD/Dataset/test_pdb_list",
    "pdb2amber_path" : "/root/SMC_PPI_MD/pdb2amber/",
    "box_size" : 35,   #35 angstrom
    "ions_mdp_path" : "/root/SMC_PPI_MD/SuMD/ions.mdp"
}

def run_script(script, path):
    try :
        subprocess.run(script, cwd = path,  shell =True)
    except Exception as e :
        print(f"Running Failed : {e}  || {script} ||")

class MD_automation():
    def __init__(self, config ):
        self.dataset_path = config["dataset_path"]
        self.save_path = config["save_path"]
        self.pdb_ID_list = []
        with open(config["pdb_list_path"], 'r') as f:
            for line in f:  # iterate over each line in the file
                self.pdb_ID_list.append(line.strip()[:-4])  # .strip() removes any trailing newline characters
        self.SuMD_path = config["SuMD_path"]
        self.pdb2amber = config["pdb2amber_path"]
        self.box_len = config["box_size"]
        self.ions_mdp = config["ions_mdp_path"]
        
    def inputfile_generation(self, complex_name):
        protein_resid =[]
        ligand_resid = []
        protein = ""
        peptide = ""
        protein_ID = complex_name[-1]
        peptide_chain_ID = complex_name[-3]

        with open(f"{self.dataset_path}{complex_name}.pdb", 'r') as complex:
            for line in complex.readlines():
                line = str(line).split()
                if line[4] == protein_ID:
                    protein_resid.append(line[5])
                elif line[4] == peptide_chain_ID:
                    ligand_resid.append(line[5])

        protein_resid = list(set(protein_resid))
        ligand_resid = list(set(ligand_resid))

        for resid in protein_resid:
            protein += f"{resid} "
        protein = protein[:-1]
        
        for resid in ligand_resid:
            peptide += f"{resid} "
        peptide = peptide[:-1]

        script = f"""### SYSTEM SETTING
structre={complex_name}_sumd.pdb
parameters={complex_name}_sumd.prmtop
ForceField=AMBER

#RECEPTOR
main_chain=protein
resid={protein}

#LIGAND
ligand=peptide
ligand_chain=peptide
ligand_cm={peptide}

randomize=no
constrain=no
slope=0.00

#### SIMULATION SETTINGS
n_device=0
n_steps=150000
timestep=2
MaxFailed=17
opt_dist=5
bound_dist=10
meta_dist=15
"""

        with open(f"{self.save_path}{complex_name}/{complex_name}.dat", 'w') as f:
            f.write(script)

        with open(f"{self.save_path}{complex_name}/{complex_name}_sumd.xsc", 'w') as f:
            f.write(f"1000000 3.500000e+01 0.000000e+00 0.000000e+00 0.000000e+00 3.500000e+01 0.000000e+00 0.000000e+00 0.000000e+00 3.500000e+01 0.000000e+00 0.000000e+00 0.000000e+00")

    def generate_prmtop(self, complex_name):
        pdb_path = self.save_path + complex_name + "/"
        box_len = self.box_len

        with open(f"{pdb_path}{complex_name}.pdb", 'r') as infile, open(f"{pdb_path}H_fixed_{complex_name}.pdb", 'w') as outfile:   ## HETATM chain 어떻게 처리할지 
            for line in infile:
                if not line.startswith("HETATM"):
                    outfile.write(line)

        subprocess.run(f"cp {self.SuMD_path}ions.mdp {pdb_path}", cwd = self.SuMD_path, shell=True)
        subprocess.run(f"pdbfixer H_fixed_{complex_name}.pdb --replace-nonstandard --add-atoms=all --output=fixed_{complex_name}.pdb", cwd = pdb_path, shell=True)

        #editconf
        script = f"gmx editconf -f fixed_{complex_name}.pdb -o processed_{complex_name}.pdb -box {box_len} {box_len} {box_len}" 
        run_script(script=script, path = pdb_path)
        time.sleep(0.5)
        script = f"gmx pdb2gmx -f processed_{complex_name}.pdb -o processed_{complex_name}.gro -water tip3p -ignh"
        subprocess.run(script, input="1\n", cwd= pdb_path , text=True, shell =True)
        time.sleep(0.5)
        #solvate
        script = f"gmx solvate -cp processed_{complex_name}.gro -cs spc216.gro -o {complex_name}_solvate.gro -p topol.top"
        run_script(script=script, path = pdb_path)

        script = f"gmx grompp -f ions.mdp -c {complex_name}_solvate.gro -p topol.top -o ions.tpr"
        run_script(script=script, path = pdb_path)

        #ionize
        master_fd, slave_fd = os.openpty()
        script = f"gmx genion -s {pdb_path}ions.tpr -o {pdb_path}{complex_name}_ionized.pdb -p {pdb_path}topol.top -pname NA -nname CL -conc 0.15"

        proc = subprocess.Popen(
            script,
            stdin=slave_fd,
            stdout=slave_fd,
            stderr=slave_fd,
            text=True,
            bufsize=0,
            shell = True
        )
        os.close(slave_fd)

        sol_group = None
        buffer = ""
        prompt_sent = False

        while True:
            # master_fd에서 1초 동안 데이터가 도착하는지 확인
            r, _, _ = select.select([master_fd], [], [], 4.0)
            if master_fd in r:
                try:
                    data = os.read(master_fd, 1024).decode()
                except OSError:
                    break
                if not data:
                    break
                # 받은 데이터를 즉시 출력
                sys.stdout.write(data)
                sys.stdout.flush()
                buffer += data

                # 버퍼를 줄 단위로 나눔
                lines = buffer.splitlines(keepends=True)
                for line in lines:
                    # "SOL"이 포함된 줄에서 그룹 번호 추출 (최초 한 번만)
                    if "SOL" in line and sol_group is None:
                        match = re.search(r"Group\s+(\d+)\s+\(.*SOL.*\)", line)
                        if match:
                            sol_group = match.group(1)
                        else :
                            match_2 = re.search(r"Group\s+(\d+)\s+\(.*Water.*\)", line)
                            if match_2:
                                sol_group = match_2.group(1)
                    # "Select a group:" 프롬프트가 나오면 입력 전송
                    if "Select a group:" in line and not prompt_sent:
                        if sol_group is None:
                            # SOL 그룹 번호를 찾지 못했다면 빈 문자열을 입력하거나 원하는 기본값 사용
                            sol_group = ""
                        os.write(master_fd, (sol_group + "\n").encode())
                        prompt_sent = True
                        # 프롬프트 이후 버퍼는 초기화 (또는 필요에 따라 유지)
                        buffer = ""
            else:
                # 지정한 시간 동안 더 이상 출력이 없으면 루프 종료
                break

        # 프로세스 종료 대기
        proc.wait()
        os.close(master_fd)

        script = "rm '#topol.top.1#'; rm '#topol.top.2#' ; rm *.gro ; rm *.itp ; rm *.tpr; rm *.top; rm mdout.mdp; rm fixed_*.pdb; rm processed_*.pdb; rm H_fixed_*.pdb"
        
        if not debug_mode:
            subprocess.run(script, cwd = pdb_path, shell = True)

        script = f"cp {complex_name}_ionized.pdb {self.pdb2amber}"
        subprocess.run(script, cwd = pdb_path, shell = True)

        script = f"pdbfixer {complex_name}_ionized.pdb --replace-nonstandard --add-atoms=all --output={complex_name}_sumd.pdb"
        run_script(script=script, path = self.pdb2amber)

        #inputfile 
        data = {
            "fname_pdb": f"{complex_name}_sumd.pdb",
            "fname_prmtop": f"{complex_name}_sumd.prmtop",
            "fname_ff": [
                "/root/SMC_PPI_MD/SuMD/forcefield/protein.ff14SB.xml",
                "/root/SMC_PPI_MD/SuMD/forcefield/tip3p.xml"
            ]
        }

        with open(f"{self.pdb2amber}{complex_name}_input.json", 'w') as f:
            json.dump(data, f)

        script = f"python pdb2amber.py -i {complex_name}_input.json"
        run_script(script=script, path = self.pdb2amber)

        time.sleep(0.5)

        if not debug_mode:
            script = f"rm {complex_name}_final.pdb; rm {complex_name}_input.json; cd {complex_name}_sumd.prmtop {pdb_path}"
        else : 
            script = f"cd {complex_name}_sumd.prmtop {pdb_path}"
        run_script(script=script, path = self.pdb2amber)

    def MD_run(self, inputfile, pdb_ID):
        subprocess.run(f"python3 {self.SuMD_path}suMD {inputfile}", cwd = f"{self.save_path}{pdb_ID}", shell=True)
        if not debug_mode:
            subprocess.run(f'rm {pdb_ID}.pdb ; rm {pdb_ID}.dat', cwd = f"{self.save_path}{pdb_ID}", shell = True)

    def MD_automation_start(self):
        for pdb_ID in tqdm(self.pdb_ID_list):
            if os.path.exists(f"{self.save_path}{pdb_ID}"):
                continue
            else :
                os.mkdir(f"{self.save_path}{pdb_ID}")
                
            subprocess.run(f"cp {pdb_ID}.pdb {self.save_path}{pdb_ID}/", cwd = self.dataset_path, shell = True)   
            try :
                self.generate_prmtop(pdb_ID)     
                self.inputfile_generation(pdb_ID)
                self.MD_run(f"{self.save_path}{pdb_ID}.dat", pdb_ID)
            except Exception as e :
                with open(f"{self.save_path}{pdb_ID}/error.txt", 'w') as f:
                    f.write(f"Error : {e}")
                continue

if __name__ == '__main__':
    MD = MD_automation(configuration)
    MD.MD_automation_start()

    