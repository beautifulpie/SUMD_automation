#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import os
import sys
import numpy as np
from Bio.PDB import PDBParser, PDBIO, Atom

def process_cystein(input_pdb, output_pdb):
    """
    PDB 파일의 시스테인(CYS) 잔기에 수소 원자(HG)를 추가합니다.
    
    Args:
        input_pdb: 입력 PDB 파일 경로
        output_pdb: 출력 PDB 파일 경로
    
    Returns:
        처리된 CYS 잔기 수
    """
    print(f"시스테인 잔기 처리 중: {input_pdb} -> {output_pdb}")
    
    # PDB 파일 파싱
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("structure", input_pdb)
    
    # 처리된 CYS 잔기 수 카운트
    processed_count = 0
    
    # 모델, 체인, 잔기를 순회하면서 CYS 잔기 찾기
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.get_resname() == "CYS":
                    # 이미 HG가 존재하면, 건너뜀
                    if "HG" not in residue:
                        # SG와 CA 원자가 있는지 확인
                        if "SG" in residue and "CA" in residue:
                            sg_atom = residue["SG"]
                            ca_atom = residue["CA"]
                            
                            # SG와 CA 사이의 벡터 계산 및 단위 벡터 산출
                            vector = sg_atom.coord - ca_atom.coord
                            norm = np.linalg.norm(vector)
                            if norm == 0:
                                continue  # 0 나누기 방지
                            unit_vector = vector / norm
                            
                            # S-H 결합 길이 (Å): 일반적으로 약 1.33 Å
                            bond_length = 1.33
                            hg_coord = sg_atom.coord + unit_vector * bond_length
                            
                            # 새로운 HG 원자 생성
                            hg_atom = Atom(name="HG",
                                       coord=hg_coord,
                                       bfactor=sg_atom.get_bfactor(),
                                       occupancy=sg_atom.get_occupancy(),
                                       altloc=sg_atom.get_altloc(),
                                       fullname=" HG ",
                                       serial_number=0,
                                       element="H")
                            
                            # HG 원자를 잔기에 추가
                            residue.add(hg_atom)
                            processed_count += 1
                        else:
                            print(f"경고: 잔기 {residue.get_id()}에 SG 또는 CA 원자가 없어 HG를 추가하지 않습니다.")
    
    # 변경된 구조를 새로운 PDB 파일로 저장
    io = PDBIO()
    io.set_structure(structure)
    io.save(output_pdb)
    
    print(f"처리 완료: {processed_count}개의 시스테인 잔기에 HG 원자 추가")
    return processed_count

def main():
    parser = argparse.ArgumentParser(description="PDB 파일의 시스테인 잔기에 HG 원자 추가")
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
    
    # 시스테인 잔기 처리
    processed_count = process_cystein(args.input, args.output)
    
    print(f"총 {processed_count}개의 시스테인 잔기가 처리되었습니다.")
    return 0

if __name__ == "__main__":
    sys.exit(main())
