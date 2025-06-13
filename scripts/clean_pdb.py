from Bio import PDB
from Bio.PDB.PDBIO import PDBIO
import argparse
import os
import sys
import numpy as np

def extract_amino_acids_only(input_pdb, output_pdb):
    """
    PDB 파일에서 표준 아미노산만 추출하여 새 PDB 파일로 저장합니다.
    
    Parameters:
    input_pdb (str): 입력 PDB 파일 경로
    output_pdb (str): 출력 PDB 파일 경로
    """
    # 표준 아미노산 목록
    standard_aa = [
        "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", 
        "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", 
        "THR", "TRP", "TYR", "VAL"
    ]
    
    # N 말단과 C 말단 변형 아미노산도 포함
    terminal_aa = ["ACE", "NME", "NHE", "NH2"]
    
    # 히스티딘 변형
    his_variants = ["HID", "HIE", "HIP", "HSE", "HSD", "HSP"]
    
    # 시스테인 변형
    cys_variants = ["CYX", "CYM"]
    
    # 모든 허용 아미노산
    allowed_residues = standard_aa + terminal_aa + his_variants + cys_variants
    
    # PDB 파서 초기화
    parser = PDB.PDBParser(QUIET=True)
    
    # PDB 파일 읽기
    structure_id = os.path.basename(input_pdb).split('.')[0]
    structure = parser.get_structure(structure_id, input_pdb)
    
    # 비표준 잔기 제거를 위한 클래스 정의
    class StandardAminoAcidSelect(PDB.Select):
        def accept_residue(self, residue):
            # 표준 아미노산만 선택
            resname = residue.get_resname()
            return resname in allowed_residues
    
    # PDB 파일 작성을 위한 IO 객체
    io = PDBIO()
    io.set_structure(structure)
    
    # 표준 아미노산만 선택하여 파일 저장
    io.save(output_pdb, StandardAminoAcidSelect())
    
    print(f"표준 아미노산만 포함된 PDB 파일이 '{output_pdb}'에 저장되었습니다.")
    return True

def main():
    parser = argparse.ArgumentParser(description="PDB 파일에서 표준 아미노산만 추출하여 새 PDB 파일로 저장합니다.")
    parser.add_argument("--input", "-i", required=True, help="입력 PDB 파일 경로")
    parser.add_argument("--output", "-o", required=True, help="출력 PDB 파일 경로")
    
    args = parser.parse_args()
    
    # 입력 파일 확인
    if not os.path.exists(args.input):
        print(f"오류: 입력 파일 {args.input}을 찾을 수 없습니다.")
        return 1
    
    # 출력 디렉토리 확인
    output_dir = os.path.dirname(args.output)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)
    
    # 표준 아미노산 추출 및 pdb 저장.
    stdout = extract_amino_acids_only(args.input, args.output)
    
    return 0

# 사용 예시
if __name__ == "__main__":
    # 명령줄 인자를 사용하도록 수정
    sys.exit(main())