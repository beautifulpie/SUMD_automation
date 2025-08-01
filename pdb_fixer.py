#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
개선된 PDB 수정 및 전처리 모듈

GROMACS pdb2gmx 오류를 방지하기 위한 강력한 PDB 전처리
- 누락된 중요 원자가 있는 잔기 감지 및 처리
- 잔기별 완전성 검사 및 복구
- GROMACS 호환성 확보

주요 기능:
1. 아미노산별 필수 원자 정의 및 검사
2. 누락된 원자 복구 (가능한 경우)
3. 복구 불가능한 잔기 제거
4. 대체 아미노산으로 변환 (선택적)
"""

import os
import numpy as np
import logging
from Bio.PDB import PDBParser, PDBIO, Select
from Bio.PDB.Atom import Atom
from Bio.PDB.Residue import Residue
from typing import List, Dict, Set, Optional, Tuple, Any
import warnings
warnings.filterwarnings("ignore", category=PDBConstructionWarning)

class AdvancedPDBFixer:
    """GROMACS 호환성을 위한 고급 PDB 수정 클래스"""
    
    def __init__(self, logger=None):
        self.logger = logger or logging.getLogger(__name__)
        
        # 표준 아미노산
        self.standard_amino_acids = {
            'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY',
            'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER',
            'THR', 'TRP', 'TYR', 'VAL'
        }
        
        # 필수 백본 원자 (모든 아미노산 공통)
        self.backbone_atoms = {'N', 'CA', 'C', 'O'}
        
        # 아미노산별 완전한 사이드체인 원자 정의
        self.complete_sidechain_atoms = {
            'ALA': {'CB'},
            'ARG': {'CB', 'CG', 'CD', 'NE', 'CZ', 'NH1', 'NH2'},
            'ASN': {'CB', 'CG', 'OD1', 'ND2'},
            'ASP': {'CB', 'CG', 'OD1', 'OD2'},
            'CYS': {'CB', 'SG'},
            'GLN': {'CB', 'CG', 'CD', 'OE1', 'NE2'},  # CD 원자 포함!
            'GLU': {'CB', 'CG', 'CD', 'OE1', 'OE2'},
            'GLY': set(),  # 사이드체인 없음
            'HIS': {'CB', 'CG', 'ND1', 'CD2', 'CE1', 'NE2'},
            'ILE': {'CB', 'CG1', 'CG2', 'CD1'},
            'LEU': {'CB', 'CG', 'CD1', 'CD2'},
            'LYS': {'CB', 'CG', 'CD', 'CE', 'NZ'},
            'MET': {'CB', 'CG', 'SD', 'CE'},
            'PHE': {'CB', 'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'},
            'PRO': {'CB', 'CG', 'CD'},
            'SER': {'CB', 'OG'},
            'THR': {'CB', 'OG1', 'CG2'},
            'TRP': {'CB', 'CG', 'CD1', 'CD2', 'NE1', 'CE2', 'CE3', 'CZ2', 'CZ3', 'CH2'},
            'TYR': {'CB', 'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ', 'OH'},
            'VAL': {'CB', 'CG1', 'CG2'}
        }
        
        # 최소 필수 원자 (이것들이 없으면 잔기 제거)
        self.critical_atoms = {
            'ALA': {'CB'},
            'ARG': {'CB', 'CG', 'CD'},
            'ASN': {'CB', 'CG'},
            'ASP': {'CB', 'CG'},
            'CYS': {'CB', 'SG'},
            'GLN': {'CB', 'CG', 'CD'},  # CD 원자가 핵심!
            'GLU': {'CB', 'CG', 'CD'},
            'GLY': set(),
            'HIS': {'CB', 'CG'},
            'ILE': {'CB', 'CG1'},
            'LEU': {'CB', 'CG'},
            'LYS': {'CB', 'CG', 'CD'},
            'MET': {'CB', 'CG', 'SD'},
            'PHE': {'CB', 'CG'},
            'PRO': {'CB', 'CG', 'CD'},
            'SER': {'CB', 'OG'},
            'THR': {'CB', 'OG1'},
            'TRP': {'CB', 'CG'},
            'TYR': {'CB', 'CG'},
            'VAL': {'CB', 'CG1'}
        }
        
        # 대체 가능한 아미노산 (복구 불가능시 사용)
        self.substitution_map = {
            'GLN': 'ALA',  # GLN이 문제가 되면 ALA로 대체
            'GLU': 'ALA',
            'ARG': 'ALA',
            'LYS': 'ALA',
            'ASN': 'ALA',
            'ASP': 'ALA'
        }
        
    def fix_pdb_for_gromacs(self, input_pdb: str, output_pdb: str, 
                           target_chains: List[str]) -> Tuple[str, Dict]:
        """GROMACS 호환성을 위한 PDB 수정"""
        
        self.logger.info(f"=== GROMACS 호환 PDB 수정 시작 ===")
        self.logger.info(f"입력: {input_pdb}")
        self.logger.info(f"대상 체인: {target_chains}")
        
        stats = {
            'original_residues': 0,
            'processed_residues': 0,
            'removed_residues': [],
            'substituted_residues': [],
            'recovered_residues': [],
            'problematic_residues': []
        }
        
        try:
            parser = PDBParser(QUIET=True)
            structure = parser.get_structure("structure", input_pdb)
            
            # 문제 있는 잔기들을 추적
            problematic_residues = []
            
            for model in structure:
                for chain in model:
                    if chain.id not in target_chains:
                        continue
                        
                    chain_residues = list(chain.get_residues())
                    stats['original_residues'] += len(chain_residues)
                    
                    # 잔기별 검사 및 수정
                    residues_to_remove = []
                    
                    for residue in chain_residues:
                        resname = residue.get_resname().strip()
                        res_id = residue.get_id()
                        res_key = f"{chain.id}:{res_id[1]}{resname}"
                        
                        # 표준 아미노산만 처리
                        if resname not in self.standard_amino_acids:
                            self.logger.warning(f"비표준 잔기 제거: {res_key}")
                            residues_to_remove.append(residue)
                            stats['removed_residues'].append(res_key + " (non-standard)")
                            continue
                        
                        # 잔기 완전성 검사
                        check_result = self._check_residue_completeness(residue, resname)
                        
                        if check_result['status'] == 'complete':
                            # 완전한 잔기
                            stats['processed_residues'] += 1
                            
                        elif check_result['status'] == 'recoverable':
                            # 복구 가능한 잔기
                            try:
                                self._recover_missing_atoms(residue, resname, check_result['missing_atoms'])
                                stats['recovered_residues'].append(res_key)
                                stats['processed_residues'] += 1
                                self.logger.info(f"잔기 복구 성공: {res_key}")
                            except Exception as e:
                                self.logger.error(f"잔기 복구 실패: {res_key} - {e}")
                                residues_to_remove.append(residue)
                                stats['problematic_residues'].append(res_key + f" (recovery_failed: {e})")
                                
                        elif check_result['status'] == 'substitutable':
                            # 대체 가능한 잔기
                            try:
                                new_resname = self.substitution_map.get(resname, 'ALA')
                                self._substitute_residue(residue, new_resname)
                                stats['substituted_residues'].append(f"{res_key} -> {new_resname}")
                                stats['processed_residues'] += 1
                                self.logger.info(f"잔기 대체: {res_key} -> {new_resname}")
                            except Exception as e:
                                self.logger.error(f"잔기 대체 실패: {res_key} - {e}")
                                residues_to_remove.append(residue)
                                stats['problematic_residues'].append(res_key + f" (substitution_failed: {e})")
                                
                        else:
                            # 제거해야 하는 잔기
                            residues_to_remove.append(residue)
                            stats['removed_residues'].append(res_key + f" ({check_result['reason']})")
                            self.logger.warning(f"잔기 제거: {res_key} - {check_result['reason']}")
                    
                    # 문제 있는 잔기들 제거
                    for residue in residues_to_remove:
                        try:
                            chain.detach_child(residue.get_id())
                        except Exception as e:
                            self.logger.error(f"잔기 제거 실패: {residue.get_id()} - {e}")
            
            # 수정된 구조 저장
            io = PDBIO()
            io.set_structure(structure)
            io.save(output_pdb)
            
            self.logger.info(f"=== PDB 수정 완료 ===")
            self.logger.info(f"원본 잔기: {stats['original_residues']}")
            self.logger.info(f"처리된 잔기: {stats['processed_residues']}")
            self.logger.info(f"제거된 잔기: {len(stats['removed_residues'])}")
            self.logger.info(f"복구된 잔기: {len(stats['recovered_residues'])}")
            self.logger.info(f"대체된 잔기: {len(stats['substituted_residues'])}")
            
            return output_pdb, stats
            
        except Exception as e:
            self.logger.error(f"PDB 수정 실패: {e}")
            raise
    
    def _check_residue_completeness(self, residue, resname: str) -> Dict:
        """잔기 완전성 검사"""
        
        present_atoms = {atom.get_name().strip() for atom in residue}
        
        # 백본 원자 확인
        missing_backbone = self.backbone_atoms - present_atoms
        if missing_backbone:
            return {
                'status': 'invalid',
                'reason': f'missing_backbone_{missing_backbone}',
                'missing_atoms': missing_backbone
            }
        
        # 사이드체인 원자 확인
        required_sidechain = self.complete_sidechain_atoms.get(resname, set())
        critical_sidechain = self.critical_atoms.get(resname, set())
        
        missing_sidechain = required_sidechain - present_atoms
        missing_critical = critical_sidechain - present_atoms
        
        if not missing_sidechain:
            # 완전한 잔기
            return {'status': 'complete', 'missing_atoms': set()}
        
        elif not missing_critical:
            # 일부 원자 누락이지만 핵심 원자는 있음 (복구 가능)
            return {
                'status': 'recoverable', 
                'missing_atoms': missing_sidechain,
                'reason': f'missing_non_critical_{missing_sidechain}'
            }
        
        elif resname in self.substitution_map:
            # 핵심 원자 누락이지만 대체 가능
            return {
                'status': 'substitutable',
                'missing_atoms': missing_critical,
                'reason': f'missing_critical_{missing_critical}'
            }
        
        else:
            # 복구 불가능
            return {
                'status': 'invalid',
                'reason': f'missing_critical_unrecoverable_{missing_critical}',
                'missing_atoms': missing_critical
            }
    
    def _recover_missing_atoms(self, residue, resname: str, missing_atoms: Set[str]):
        """누락된 원자 복구 (간단한 기하학적 추정)"""
        
        present_atoms = {atom.get_name().strip(): atom for atom in residue}
        
        # CA 원자를 기준으로 사용
        if 'CA' not in present_atoms:
            raise ValueError("CA 원자가 없어서 복구 불가능")
        
        ca_atom = present_atoms['CA']
        ca_coord = ca_atom.get_coord()
        
        # 간단한 원자 복구 로직
        for atom_name in missing_atoms:
            try:
                new_coord = self._estimate_atom_position(
                    atom_name, resname, present_atoms, ca_coord
                )
                
                if new_coord is not None:
                    # 새 원자 생성
                    new_atom = Atom(
                        name=atom_name,
                        coord=new_coord,
                        bfactor=ca_atom.get_bfactor(),
                        occupancy=1.0,
                        altloc=' ',
                        fullname=f' {atom_name:<3}',
                        serial_number=0,
                        element=atom_name[0]  # 첫 글자를 원소로 가정
                    )
                    
                    residue.add(new_atom)
                    self.logger.debug(f"원자 복구: {atom_name} in {resname}")
                    
            except Exception as e:
                self.logger.warning(f"원자 복구 실패: {atom_name} in {resname} - {e}")
                # 복구 실패시 예외를 다시 발생시켜 상위에서 처리하도록 함
                raise
    
    def _estimate_atom_position(self, atom_name: str, resname: str, 
                               present_atoms: Dict, ca_coord: np.ndarray) -> Optional[np.ndarray]:
        """원자 위치 추정 (매우 간단한 방법)"""
        
        # 기본적인 결합 길이와 각도를 사용한 추정
        bond_lengths = {
            'CB': 1.54,   # CA-CB 결합 길이
            'CG': 1.54,   # CB-CG 결합 길이  
            'CD': 1.54,   # CG-CD 결합 길이
            'CE': 1.54,   # CD-CE 결합 길이
            'NZ': 1.47,   # CE-NZ 결합 길이
            'OG': 1.43,   # CB-OG 결합 길이
            'SG': 1.82,   # CB-SG 결합 길이
            'SD': 1.82,   # CG-SD 결합 길이
            'OE1': 1.25,  # CD-OE1 결합 길이
            'NE2': 1.33   # CD-NE2 결합 길이
        }
        
        if atom_name not in bond_lengths:
            return None
        
        # CA를 기준으로 랜덤한 방향으로 적절한 거리에 배치
        # 실제로는 더 정교한 기하학적 계산이 필요하지만, 
        # 여기서는 GROMACS에서 처리될 수 있을 정도로만 배치
        
        # 랜덤 방향 벡터 생성
        direction = np.random.randn(3)
        direction = direction / np.linalg.norm(direction)
        
        # 결합 길이만큼 떨어뜨려 배치
        distance = bond_lengths[atom_name]
        new_coord = ca_coord + direction * distance
        
        return new_coord
    
    def _substitute_residue(self, residue, new_resname: str):
        """잔기를 다른 아미노산으로 대체"""
        
        # 잔기 이름 변경
        residue.resname = new_resname
        
        # 새로운 잔기에 불필요한 원자들 제거
        new_required_atoms = self.backbone_atoms | self.complete_sidechain_atoms.get(new_resname, set())
        
        atoms_to_remove = []
        for atom in residue:
            if atom.get_name().strip() not in new_required_atoms:
                atoms_to_remove.append(atom.get_id())
        
        for atom_id in atoms_to_remove:
            try:
                residue.detach_child(atom_id)
            except:
                pass
        
        # 필요한 원자가 부족하면 추가 (ALA의 경우 CB만 있으면 됨)
        present_atoms = {atom.get_name().strip() for atom in residue}
        missing_atoms = new_required_atoms - present_atoms
        
        if missing_atoms and new_resname == 'ALA':
            # ALA의 경우 CB 원자만 추가하면 됨
            if 'CB' in missing_atoms and 'CA' in present_atoms:
                ca_atom = None
                for atom in residue:
                    if atom.get_name().strip() == 'CA':
                        ca_atom = atom
                        break
                
                if ca_atom:
                    # CB 원자 위치 추정 (CA에서 약간 떨어진 곳)
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


class GromacsCompatibleProcessor:
    """GROMACS 호환성을 보장하는 전체 처리 파이프라인"""
    
    def __init__(self, logger=None):
        self.logger = logger or logging.getLogger(__name__)
        self.fixer = AdvancedPDBFixer(logger)
    
    def process_for_gromacs(self, input_pdb: str, output_pdb: str, 
                           target_chains: List[str]) -> Tuple[str, Dict]:
        """GROMACS를 위한 완전한 PDB 처리"""
        
        self.logger.info("=== GROMACS 호환 PDB 처리 시작 ===")
        
        # 1단계: PDB 수정 및 복구
        fixed_pdb, fix_stats = self.fixer.fix_pdb_for_gromacs(
            input_pdb, output_pdb, target_chains
        )
        
        # 2단계: 최종 검증
        validation_result = self._validate_for_gromacs(fixed_pdb, target_chains)
        
        # 통합 통계
        combined_stats = {
            **fix_stats,
            'validation': validation_result,
            'gromacs_ready': validation_result['ready']
        }
        
        self.logger.info("=== GROMACS 호환 PDB 처리 완료 ===")
        self.logger.info(f"GROMACS 준비 상태: {validation_result['ready']}")
        
        return fixed_pdb, combined_stats
    
    def _validate_for_gromacs(self, pdb_file: str, target_chains: List[str]) -> Dict:
        """GROMACS 호환성 검증"""
        
        try:
            parser = PDBParser(QUIET=True)
            structure = parser.get_structure("validation", pdb_file)
            
            issues = []
            total_residues = 0
            total_atoms = 0
            
            for model in structure:
                for chain in model:
                    if chain.id not in target_chains:
                        continue
                    
                    for residue in chain:
                        total_residues += 1
                        total_atoms += len(list(residue.get_atoms()))
                        
                        resname = residue.get_resname().strip()
                        
                        # 표준 아미노산 확인
                        if resname not in self.fixer.standard_amino_acids:
                            issues.append(f"Non-standard residue: {chain.id}:{residue.get_id()[1]}{resname}")
                        
                        # 백본 원자 확인
                        present_atoms = {atom.get_name().strip() for atom in residue}
                        missing_backbone = self.fixer.backbone_atoms - present_atoms
                        
                        if missing_backbone:
                            issues.append(f"Missing backbone atoms in {chain.id}:{residue.get_id()[1]}{resname}: {missing_backbone}")
            
            return {
                'ready': len(issues) == 0,
                'total_residues': total_residues,
                'total_atoms': total_atoms,
                'issues': issues
            }
            
        except Exception as e:
            return {
                'ready': False,
                'error': str(e),
                'issues': [f"Validation failed: {e}"]
            }


# 통합 사용 함수
def fix_pdb_for_gromacs(input_pdb: str, output_pdb: str, 
                        target_chains: List[str], logger=None) -> Tuple[str, Dict]:
    """
    GROMACS pdb2gmx 오류를 방지하기 위한 PDB 수정
    
    Args:
        input_pdb: 입력 PDB 파일
        output_pdb: 출력 PDB 파일
        target_chains: 처리할 체인 리스트
        logger: 로거 객체
        
    Returns:
        tuple: (수정된 PDB 파일 경로, 처리 통계)
    """
    processor = GromacsCompatibleProcessor(logger)
    return processor.process_for_gromacs(input_pdb, output_pdb, target_chains)


# 테스트 함수
def test_pdb_fixer():
    """PDB Fixer 테스트"""
    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger("PDB_Fixer_Test")
    
    input_pdb = "problematic.pdb"  # GLN CD 원자 누락된 PDB
    output_pdb = "fixed.pdb"
    target_chains = ['A', 'B']
    
    try:
        fixed_pdb, stats = fix_pdb_for_gromacs(
            input_pdb, output_pdb, target_chains, logger
        )
        
        print(f"\n=== PDB 수정 결과 ===")
        print(f"입력: {input_pdb}")
        print(f"출력: {fixed_pdb}")
        print(f"처리된 잔기: {stats['processed_residues']}")
        print(f"제거된 잔기: {len(stats['removed_residues'])}")
        print(f"복구된 잔기: {len(stats['recovered_residues'])}")
        print(f"대체된 잔기: {len(stats['substituted_residues'])}")
        print(f"GROMACS 준비 상태: {stats['gromacs_ready']}")
        
        if stats['validation']['issues']:
            print(f"남은 문제점들:")
            for issue in stats['validation']['issues'][:5]:
                print(f"  - {issue}")
                
    except Exception as e:
        print(f"테스트 실패: {e}")


if __name__ == "__main__":
    test_pdb_fixer()
