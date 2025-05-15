#!/usr/bin/env python3
"""
Cysteine 잔기에 수소 원자(HG)를 자동으로 추가하는 스크립트

이 스크립트는 PDB 파일에서 Cysteine(CYS) 잔기를 찾아 HG 원자가 없는 경우
SG-CA 벡터 방향으로 1.33 Å 거리에 HG 원자를 추가합니다.
"""

from Bio.PDB import PDBParser, PDBIO
import numpy as np
from Bio.PDB.Atom import Atom
import argparse
import os
import sys

def add_hydrogens_to_cysteines(input_file, output_file):
    """
    PDB 파일의 CYS 잔기에 HG 원자를 추가합니다.
    
    Parameters:
    ----------
    input_file : str
        입력 PDB 파일 경로
    output_file : str
        출력 PDB 파일 경로
        
    Returns:
    -------
    int
        추가된 HG 원자 수
    """
    if not os.path.exists(input_file):
        print(f"오류: 입력 파일 '{input_file}'이 존재하지 않습니다.")
        return 0
    
    # 1. PDB 파일 파싱
    parser = PDBParser(QUIET=True)
    try:
        structure = parser.get_structure("structure", input_file)
    except Exception as e:
        print(f"PDB 파일 파싱 오류: {e}")
        return 0
    
    added_count = 0
    missing_atoms = []
    
    # 2. 모델, 체인, 잔기를 순회하면서 CYS 잔기 찾기
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.get_resname() == "CYS":
                    # 이미 HG가 존재하면 건너뜁니다.
                    if "HG" in residue:
                        continue
                        
                    # SG와 CA 원자가 있는지 확인
                    if "SG" in residue and "CA" in residue:
                        sg_atom = residue["SG"]
                        ca_atom = residue["CA"]
                        
                        # 3. SG와 CA 사이의 벡터 계산 및 단위 벡터 산출
                        vector = sg_atom.coord - ca_atom.coord
                        norm = np.linalg.norm(vector)
                        if norm < 0.001:  # 영벡터 방지 (오류 발생 가능)
                            print(f"경고: 잔기 {residue.get_id()}의 SG-CA 원자가 같은 위치에 있습니다.")
                            missing_atoms.append(residue.get_id())
                            continue
                            
                        unit_vector = vector / norm
                        
                        # S-H 결합 길이 (Å): 일반적으로 약 1.33 Å
                        bond_length = 1.33
                        hg_coord = sg_atom.coord + unit_vector * bond_length
                        
                        # 4. 새로운 HG 원자 생성
                        try:
                            # 시리얼 번호를 구조에서 가장 큰 값 + 1로 설정
                            max_serial = max([a.serial_number for a in structure.get_atoms()])
                            
                            hg_atom = Atom(name="HG",
                                        coord=hg_coord,
                                        bfactor=sg_atom.get_bfactor(),
                                        occupancy=sg_atom.get_occupancy(),
                                        altloc=sg_atom.get_altloc(),
                                        fullname=" HG ",
                                        serial_number=max_serial + 1,
                                        element="H")
                            
                            # HG 원자를 잔기에 추가합니다.
                            residue.add(hg_atom)
                            added_count += 1
                        except Exception as e:
                            print(f"원자 추가 오류 (잔기 {residue.get_id()}): {e}")
                            missing_atoms.append(residue.get_id())
                    else:
                        missing_msg = f"Residue {residue.get_id()}는 SG 또는 CA 원자가 없어 HG를 추가하지 않습니다."
                        print(missing_msg)
                        missing_atoms.append(residue.get_id())
    
    # 5. 변경된 구조를 새로운 PDB 파일로 저장
    if added_count > 0:
        try:
            io = PDBIO()
            io.set_structure(structure)
            io.save(output_file)
            print(f"총 {added_count}개의 HG 원자가 추가되었습니다.")
            print(f"결과가 {output_file}에 저장되었습니다.")
        except Exception as e:
            print(f"파일 저장 오류: {e}")
            return 0
    else:
        print("추가된 HG 원자가 없습니다.")
        
    # 누락된 원자 목록 출력
    if missing_atoms:
        print(f"주의: {len(missing_atoms)}개 CYS 잔기에 HG 원자를 추가할 수 없었습니다.")
        
    return added_count

def main():
    parser = argparse.ArgumentParser(description="PDB 파일의 CYS 잔기에 HG 원자를 추가합니다.")
    parser.add_argument("-i", "--input", required=True, help="입력 PDB 파일 경로")
    parser.add_argument("-o", "--output", help="출력 PDB 파일 경로 (기본값: input_cys.pdb)")
    parser.add_argument("-v", "--verbose", action="store_true", help="상세 출력 모드 활성화")
    
    args = parser.parse_args()
    
    # 출력 파일이 지정되지 않은 경우 기본값 설정
    if not args.output:
        input_base = os.path.splitext(args.input)[0]
        args.output = f"{input_base}_cys.pdb"
    
    # HG 원자 추가 실행
    count = add_hydrogens_to_cysteines(args.input, args.output)
    
    return 0 if count > 0 else 1

if __name__ == "__main__":
    sys.exit(main())