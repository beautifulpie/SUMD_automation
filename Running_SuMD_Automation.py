import os
import subprocess
import argparse
import json
import time
import shutil
from pathlib import Path
from tqdm import tqdm

class Running_SuMD_Automation:
    """
    GROMACS 기반 SuMD 자동화 시스템
    여러 PDB에 대한 SuMD 시뮬레이션을 일괄 처리
    """
    
    def __init__(self, config):
        """설정 초기화"""
        self.dataset_path = config["dataset_path"]
        self.save_path = config["save_path"]
        self.sumd_script_path = config["sumd_script_path"]
        self.pdb_list_path = config["pdb_list_path"]
        self.box_size = config["box_size"]
        self.pdb_list = []
        
        # PDB 목록 로드
        with open(self.pdb_list_path, 'r') as f:
            for line in f:
                pdb_id = line.strip()
                if pdb_id.endswith('.pdb'):
                    pdb_id = pdb_id[:-4]  # .pdb 확장자 제거
                self.pdb_list.append(pdb_id)
                
        # 기본 디렉토리 생성
        os.makedirs(self.save_path, exist_ok=True)
        
        # SuMD_script.py 복사 (없는 경우)
        if not os.path.exists(os.path.join(self.save_path, "SuMD_script.py")):
            shutil.copy2(self.sumd_script_path, os.path.join(self.save_path, "SuMD_script.py"))
    
    def identify_chain_residues(self, pdb_file):
        """
        PDB 파일에서 chain과 residue 정보를 분석
        """
        print(f"PDB 파일 분석 중: {pdb_file}")
        
        chains = {}  # chain ID별 residue 목록
        
        with open(pdb_file, 'r') as f:
            for line in f:
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    try:
                        chain_id = line[21].strip()
                        res_id = line[22:26].strip()
                        
                        if chain_id not in chains:
                            chains[chain_id] = set()
                        
                        chains[chain_id].add(res_id)
                    except:
                        continue
        
        # Chain별 통계
        print("Chain 정보:")
        for chain_id, residues in chains.items():
            print(f"Chain {chain_id}: {len(residues)}개 residue")
        
        return chains
    
    def prepare_input_file(self, pdb_id, protein_chain, peptide_chain):
        """
        SuMD 입력 파일 준비
        """
        pdb_file = f"{self.dataset_path}/{pdb_id}.pdb"
        chains = self.identify_chain_residues(pdb_file)
        
        if protein_chain not in chains or peptide_chain not in chains:
            print(f"ERROR: 지정한 chain ID가 PDB에 존재하지 않습니다.")
            return False
        
        # 각 chain의 residue 목록
        protein_residues = list(chains[protein_chain])
        peptide_residues = list(chains[peptide_chain])
        
        # 출력 디렉토리 생성
        output_dir = os.path.join(self.save_path, pdb_id)
        os.makedirs(output_dir, exist_ok=True)
        
        # 입력 파일 생성
        with open(os.path.join(output_dir, f"{pdb_id}_input.json"), 'w') as f:
            json.dump({
                "pdb_file": pdb_file,
                "output_dir": output_dir,
                "protein_chain": protein_chain,
                "peptide_chain": peptide_chain,
                "protein_residues": protein_residues,
                "peptide_residues": peptide_residues,
                "box_size": self.box_size
            }, f, indent=2)
        
        return True
    
    def run_sumd(self, pdb_id):
        """
        단일 PDB에 대해 SuMD 실행
        """
        input_file = os.path.join(self.save_path, pdb_id, f"{pdb_id}_input.json")
        if not os.path.exists(input_file):
            print(f"ERROR: {input_file} 입력 파일이 존재하지 않습니다.")
            return False
        
        # 입력 파일 로드
        with open(input_file, 'r') as f:
            config = json.load(f)
        
        # SuMD 명령어 조합
        sumd_script = os.path.join(self.save_path, "SuMD_script.py")
        cmd = [
            "python", sumd_script,
            "--pdb", config["pdb_file"],
            "--protein_res"
        ]
        
        # protein_residues 추가
        cmd.extend(config["protein_residues"])
        cmd.extend(["--peptide_res"])
        
        # peptide_residues 추가
        cmd.extend(config["peptide_residues"])
        
        # 출력 디렉토리 추가
        cmd.extend(["--output_dir", config["output_dir"]])
        
        print(f"SuMD 실행 중: {' '.join(cmd)}")
        
        # 실행
        try:
            subprocess.run(cmd, check=True)
            print(f"SuMD 완료: {pdb_id}")
            return True
        except subprocess.CalledProcessError as e:
            print(f"SuMD 실행 오류: {e}")
            return False
    
    def process_all(self, protein_chain, peptide_chain):
        """
        모든 PDB 파일에 대해 SuMD 시뮬레이션 실행
        """
        success_count = 0
        fail_count = 0
        
        for pdb_id in tqdm(self.pdb_list, desc="PDB 처리"):
            print(f"\n===== PDB 처리 중: {pdb_id} =====")
            
            # 이미 처리된 경우 확인
            results_dir = os.path.join(self.save_path, pdb_id, "best_trajectories")
            if os.path.exists(results_dir) and len(os.listdir(results_dir)) > 0:
                print(f"{pdb_id}는 이미 처리되었습니다. 건너뜁니다.")
                success_count += 1
                continue
            
            # 입력 파일 준비
            if not self.prepare_input_file(pdb_id, protein_chain, peptide_chain):
                print(f"{pdb_id} 입력 파일 준비 실패")
                fail_count += 1
                continue
            
            # SuMD 실행
            if self.run_sumd(pdb_id):
                success_count += 1
            else:
                fail_count += 1
        
        print(f"\n===== 처리 완료 =====")
        print(f"총 PDB: {len(self.pdb_list)}")
        print(f"성공: {success_count}")
        print(f"실패: {fail_count}")

def main():
    """명령행 인터페이스"""
    parser = argparse.ArgumentParser(
        description="GROMACS 기반 SuMD 자동화 시스템"
    )
    parser.add_argument("--dataset_path", required=True, 
                        help="PDB 파일이 저장된 디렉토리 경로")
    parser.add_argument("--save_path", required=True, 
                        help="결과 저장 디렉토리 경로")
    parser.add_argument("--sumd_script_path", required=True, 
                        help="SuMD_script.py 파일 경로")
    parser.add_argument("--pdb_list_path", required=True, 
                        help="처리할 PDB 목록 파일 경로")
    parser.add_argument("--box_size", type=float, default=10.0,
                        help="시뮬레이션 박스 크기 (nm)")
    parser.add_argument("--protein_chain", required=True,
                        help="단백질 체인 ID")
    parser.add_argument("--peptide_chain", required=True,
                        help="펩타이드 체인 ID")
    
    args = parser.parse_args()
    
    # 설정 구성
    config = {
        "dataset_path": args.dataset_path,
        "save_path": args.save_path,
        "sumd_script_path": args.sumd_script_path,
        "pdb_list_path": args.pdb_list_path,
        "box_size": args.box_size
    }
    
    # 자동화 시스템 실행
    automation = Running_SuMD_Automation(config)
    automation.process_all(args.protein_chain, args.peptide_chain)

if __name__ == "__main__":
    main()