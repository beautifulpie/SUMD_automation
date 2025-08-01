#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
즉시 사용 가능한 PDB 수정 스크립트

GLN CD 원자 누락 문제 및 기타 GROMACS 호환성 문제 해결
"""

import os
import sys
import numpy as np
from Bio.PDB import PDBParser, PDBIO, Select
from Bio.PDB.Atom import Atom
import logging

# 로깅 설정
logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
logger = logging.getLogger(__name__)

def fix_problematic_pdb(input_pdb, output_pdb, target_chains=['A', 'B']):
    """
    GROMACS pdb2gmx 오류를 해결하는 PDB 수정
    
    Args:
        input_pdb: 문제가 있는 PDB 파일
        output_pdb: 수정된 PDB 파일
        target_chains: 처리할 체인 리스트
    """
    
    logger.info(f"PDB 수정 시작: {input_pdb} -> {output_pdb}")
    
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("structure", input_pdb)
    
    # 필수 원자 정의
    backbone_atoms = {'N', 'CA', 'C', 'O'}
    critical_sidechain = {
        'GLN': {'CB', 'CG', 'CD'},  # GLN에는 CD가 필수!
        'GLU': {'CB', 'CG', 'CD'},
        'ARG': {'CB', 'CG', 'CD'},
        'LYS': {'CB', 'CG', 'CD'},
        'ASN': {'CB', 'CG'},
        'ASP': {'CB', 'CG'},
        'CYS': {'CB', 'SG'},
        'HIS': {'CB', 'CG'},
        'PRO': {'CB', 'CG', 'CD'},
    }
    
    # 표준 아미노산 리스트
    standard_aa = {
        'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY',
        'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER',
        'THR', 'TRP', 'TYR', 'VAL'
    }
    
    removed_residues = []
    substituted_residues = []
    
    for model in structure:
        for chain in model:
            if chain.id not in target_chains:
                continue
                
            residues_to_remove = []
            
            for residue in chain:
                resname = residue.get_resname().strip()
                res_id = f"{chain.id}:{residue.get_id()[1]}{resname}"
                
                # 비표준 아미노산 제거
                if resname not in standard_aa:
                    logger.warning(f"비표준 잔기 제거: {res_id}")
                    residues_to_remove.append(residue)
                    removed_residues.append(res_id)
                    continue
                
                # 현재 원자들 확인
                present_atoms = {atom.get_name().strip() for atom in residue}
                
                # 백본 원자 확인
                missing_backbone = backbone_atoms - present_atoms
                if missing_backbone:
                    logger.error(f"백본 원자 누락으로 제거: {res_id} (누락: {missing_backbone})")
                    residues_to_remove.append(residue)
                    removed_residues.append(res_id)
                    continue
                
                # 중요한 사이드체인 원자 확인
                if resname in critical_sidechain:
                    required_sidechain = critical_sidechain[resname]
                    missing_sidechain = required_sidechain - present_atoms
                    
                    if missing_sidechain:
                        # GLN이나 GLU는 ALA로 대체 시도
                        if resname in ['GLN', 'GLU']:
                            logger.warning(f"중요 원자 누락으로 ALA 대체: {res_id} (누락: {missing_sidechain})")
                            try:
                                substitute_to_alanine(residue)
                                substituted_residues.append(f"{res_id} -> ALA")
                            except Exception as e:
                                logger.error(f"ALA 대체 실패, 제거: {res_id} - {e}")
                                residues_to_remove.append(residue)
                                removed_residues.append(res_id)
                        else:
                            logger.error(f"중요 원자 누락으로 제거: {res_id} (누락: {missing_sidechain})")
                            residues_to_remove.append(residue)
                            removed_residues.append(res_id)
            
            # 문제 있는 잔기들 제거
            for residue in residues_to_remove:
                try:
                    chain.detach_child(residue.get_id())
                except Exception as e:
                    logger.error(f"잔기 제거 실패: {residue.get_id()} - {e}")
    
    # 수정된 구조 저장
    io = PDBIO()
    io.set_structure(structure)
    io.save(output_pdb)
    
    logger.info(f"PDB 수정 완료!")
    logger.info(f"제거된 잔기: {len(removed_residues)}개")
    logger.info(f"대체된 잔기: {len(substituted_residues)}개")
    
    if removed_residues:
        logger.warning(f"제거된 잔기들: {removed_residues[:10]}")  # 처음 10개만 표시
    
    if substituted_residues:
        logger.info(f"대체된 잔기들: {substituted_residues}")
    
    return output_pdb


def substitute_to_alanine(residue):
    """잔기를 ALA(알라닌)으로 대체"""
    
    # 잔기 이름 변경
    residue.resname = 'ALA'
    
    # ALA에 필요한 원자: N, CA, C, O, CB
    required_atoms = {'N', 'CA', 'C', 'O', 'CB'}
    
    # 불필요한 원자들 제거
    atoms_to_remove = []
    ca_atom = None
    
    for atom in residue:
        atom_name = atom.get_name().strip()
        if atom_name == 'CA':
            ca_atom = atom
        if atom_name not in required_atoms:
            atoms_to_remove.append(atom.get_id())
    
    for atom_id in atoms_to_remove:
        try:
            residue.detach_child(atom_id)
        except:
            pass
    
    # CB 원자가 없으면 추가
    present_atoms = {atom.get_name().strip() for atom in residue}
    if 'CB' not in present_atoms and ca_atom is not None:
        # CB 위치 추정 (CA에서 약간 떨어진 곳)
        ca_coord = ca_atom.get_coord()
        cb_coord = ca_coord + np.array([1.54, 0.0, 0.0])  # 간단한 추정
        
        cb_atom = Atom(
            name='CB',
            coord=cb_coord,
            bfactor=ca_atom.get_bfactor(),
            occupancy=1.0,
            altloc=' ',
            fullname=' CB ',
            serial_number=0,
            element='C'
        )
        
        residue.add(cb_atom)


def main():
    """메인 함수 - 명령행에서 사용"""
    
    if len(sys.argv) < 3:
        print("사용법: python quick_pdb_fix.py <입력PDB> <출력PDB> [체인1,체인2,...]")
        print("예시: python quick_pdb_fix.py input.pdb fixed.pdb A,B")
        return
    
    input_pdb = sys.argv[1]
    output_pdb = sys.argv[2]
    
    # 체인 지정 (선택적)
    if len(sys.argv) > 3:
        target_chains = [c.strip() for c in sys.argv[3].split(',')]
    else:
        target_chains = ['A', 'B']  # 기본값
    
    if not os.path.exists(input_pdb):
        logger.error(f"입력 파일을 찾을 수 없습니다: {input_pdb}")
        return
    
    try:
        fixed_pdb = fix_problematic_pdb(input_pdb, output_pdb, target_chains)
        print(f"\n✅ 성공! 수정된 PDB: {fixed_pdb}")
        print(f"이제 GROMACS pdb2gmx를 다시 시도해보세요:")
        print(f"gmx pdb2gmx -f {fixed_pdb} -o processed.gro -p topol.top -water tip3p -ff charmm36-jul2022 -ignh")
        
    except Exception as e:
        logger.error(f"PDB 수정 실패: {e}")
        return

if __name__ == "__main__":
    main()
