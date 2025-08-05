
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
개선된 PDB 수정 및 전처리 모듈

GROMACS pdb2gmx 오류를 방지하기 위한 강력한 PDB 전처리
- 누락된 중요 원자가 있는 잔기 감지 및 처리
- 잔기별 완전성 검사 및 복구
- GROMACS 호환성 확보
"""

import os
import numpy as np
import logging
from Bio.PDB import PDBParser, PDBIO, Select
from Bio.PDB.Atom import Atom
from Bio.PDB.Residue import Residue
from typing import List, Dict, Set, Optional, Tuple, Any
import warnings

# PDBConstructionWarning import 추가
try:
    from Bio.PDB.PDBExceptions import PDBConstructionWarning
    warnings.filterwarnings("ignore", category=PDBConstructionWarning)
except ImportError:
    # Bio.PDB 버전에 따른 호환성
    try:
        from Bio.PDB import PDBConstructionWarning
        warnings.filterwarnings("ignore", category=PDBConstructionWarning)
    except ImportError:
        pass

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
        
        # 아미노산별 완전한 사이드체인 원자 정의 (개선됨)
        self.complete_sidechain_atoms = {
            'ALA': {'CB'},
            'ARG': {'CB', 'CG', 'CD', 'NE', 'CZ', 'NH1', 'NH2'},
            'ASN': {'CB', 'CG', 'OD1', 'ND2'},
            'ASP': {'CB', 'CG', 'OD1', 'OD2'},
            'CYS': {'CB', 'SG'},
            'GLN': {'CB', 'CG', 'CD', 'OE1', 'NE2'},  # CD 원자가 핵심!
            'GLU': {'CB', 'CG', 'CD', 'OE1', 'OE2'},  # CD 원자가 핵심!
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
            'ARG': {'CB', 'CG'},  # CD까지는 요구하지 않음
            'ASN': {'CB', 'CG'},
            'ASP': {'CB', 'CG'},
            'CYS': {'CB', 'SG'},
            'GLN': {'CB', 'CG', 'CD'},  # CD 원자는 반드시 필요!
            'GLU': {'CB', 'CG', 'CD'},  # CD 원자는 반드시 필요!
            'GLY': set(),
            'HIS': {'CB', 'CG'},
            'ILE': {'CB', 'CG1'},
            'LEU': {'CB', 'CG'},
            'LYS': {'CB', 'CG'},  # CD까지는 요구하지 않음
            'MET': {'CB', 'CG', 'SD'},
            'PHE': {'CB', 'CG'},
            'PRO': {'CB', 'CG', 'CD'},
            'SER': {'CB', 'OG'},
            'THR': {'CB', 'OG1'},
            'TRP': {'CB', 'CG'},
            'TYR': {'CB', 'CG'},
            'VAL': {'CB', 'CG1'}
        }
        
        # 표준 결합 길이 (Å)
        self.bond_lengths = {
            'CA-CB': 1.54,
            'CB-CG': 1.54,
            'CG-CD': 1.54,
            'CD-CE': 1.54,
            'CD-NE': 1.47,
            'CD-OE1': 1.25,
            'CD-OE2': 1.25,
            'CG-OD1': 1.25,
            'CG-OD2': 1.25,
            'CB-SG': 1.82,
            'CG-SD': 1.82,
            'CB-OG': 1.43,
            'CB-OG1': 1.43
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
            'problematic_residues': [],
            'gln_glu_issues': 0  # GLN/GLU CD 원자 문제 추적
        }
        
        try:
            parser = PDBParser(QUIET=True)
            structure = parser.get_structure("structure", input_pdb)
            
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
                            # 복구 시도
                            try:
                                success = self._smart_atom_recovery(residue, resname, check_result['missing_atoms'])
                                if success:
                                    stats['recovered_residues'].append(res_key)
                                    stats['processed_residues'] += 1
                                    self.logger.info(f"잔기 복구 성공: {res_key}")
                                else:
                                    # 복구 실패시 대체 시도
                                    if resname in ['GLN', 'GLU']:
                                        stats['gln_glu_issues'] += 1
                                        self._substitute_to_alanine(residue)
                                        stats['substituted_residues'].append(f"{res_key} -> ALA (recovery_failed)")
                                        stats['processed_residues'] += 1
                                    else:
                                        residues_to_remove.append(residue)
                                        stats['removed_residues'].append(res_key + " (recovery_failed)")
                            except Exception as e:
                                self.logger.error(f"잔기 처리 실패: {res_key} - {e}")
                                residues_to_remove.append(residue)
                                stats['problematic_residues'].append(res_key + f" (exception: {e})")
                                
                        elif check_result['status'] == 'substitutable':
                            # GLN/GLU CD 원자 문제 등으로 ALA 대체
                            if resname in ['GLN', 'GLU']:
                                stats['gln_glu_issues'] += 1
                            try:
                                self._substitute_to_alanine(residue)
                                stats['substituted_residues'].append(f"{res_key} -> ALA")
                                stats['processed_residues'] += 1
                                self.logger.info(f"잔기 대체: {res_key} -> ALA")
                            except Exception as e:
                                self.logger.error(f"잔기 대체 실패: {res_key} - {e}")
                                residues_to_remove.append(residue)
                                stats['problematic_residues'].append(res_key + f" (substitution_failed: {e})")
                                
                        else:
                            # 제거
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
            self.logger.info(f"GLN/GLU CD 문제: {stats['gln_glu_issues']}개")
            
            return output_pdb, stats
            
        except Exception as e:
            self.logger.error(f"PDB 수정 실패: {e}")
            raise
    
    def _smart_atom_recovery(self, residue, resname: str, missing_atoms: Set[str]) -> bool:
        """개선된 원자 복구 로직"""
        
        present_atoms = {atom.get_name().strip(): atom for atom in residue}
        
        # CA 원자가 없으면 복구 불가능
        if 'CA' not in present_atoms:
            return False
        
        ca_atom = present_atoms['CA']
        ca_coord = ca_atom.get_coord()
        
        recovery_success = True
        
        # 원자별 복구 시도
        for atom_name in missing_atoms:
            try:
                new_coord = self._estimate_atom_position_improved(
                    atom_name, resname, present_atoms, ca_coord
                )
                
                if new_coord is not None:
                    # 새 원자 생성
                    element = atom_name[0] if atom_name[0] in ['C', 'N', 'O', 'S'] else 'C'
                    
                    new_atom = Atom(
                        name=atom_name,
                        coord=new_coord,
                        bfactor=ca_atom.get_bfactor(),
                        occupancy=1.0,
                        altloc=' ',
                        fullname=f' {atom_name:<3}',
                        serial_number=0,
                        element=element
                    )
                    
                    residue.add(new_atom)
                    self.logger.debug(f"원자 복구: {atom_name} in {resname}")
                else:
                    recovery_success = False
                    break
                    
            except Exception as e:
                self.logger.warning(f"원자 복구 실패: {atom_name} in {resname} - {e}")
                recovery_success = False
                break
        
        return recovery_success
    
    def _estimate_atom_position_improved(self, atom_name: str, resname: str, 
                                       present_atoms: Dict, ca_coord: np.ndarray) -> Optional[np.ndarray]:
        """개선된 원자 위치 추정"""
        
        # GLN/GLU의 CD 원자는 특별히 처리
        if resname in ['GLN', 'GLU'] and atom_name == 'CD':
            return self._estimate_cd_position(resname, present_atoms, ca_coord)
        
        # 기본 결합 길이 사용
        bond_key = self._get_parent_bond_key(atom_name, resname)
        if bond_key not in self.bond_lengths:
            return None
        
        bond_length = self.bond_lengths[bond_key]
        
        # 부모 원자 찾기
        parent_atom_name = self._get_parent_atom(atom_name, resname)
        if parent_atom_name not in present_atoms:
            return None
        
        parent_coord = present_atoms[parent_atom_name].get_coord()
        
        # 간단한 방향 추정 (개선된 버전)
        if parent_atom_name == 'CA':
            # CA를 기준으로 하는 경우, 백본과의 각도 고려
            direction = self._get_sidechain_direction(present_atoms, ca_coord)
        else:
            # 다른 원자를 기준으로 하는 경우
            direction = self._get_chain_direction(parent_atom_name, present_atoms)
        
        if direction is None:
            # fallback: 랜덤 방향
            direction = np.random.randn(3)
            direction = direction / np.linalg.norm(direction)
        
        new_coord = parent_coord + direction * bond_length
        return new_coord
    
    def _estimate_cd_position(self, resname: str, present_atoms: Dict, ca_coord: np.ndarray) -> Optional[np.ndarray]:
        """GLN/GLU의 CD 원자 위치 특별 추정"""
        
        # CB와 CG가 모두 있어야 함
        if 'CB' not in present_atoms or 'CG' not in present_atoms:
            return None
        
        cb_coord = present_atoms['CB'].get_coord()
        cg_coord = present_atoms['CG'].get_coord()
        
        # CB -> CG 방향으로 연장
        cb_cg_vector = cg_coord - cb_coord
        if np.linalg.norm(cb_cg_vector) == 0:
            return None
        
        cb_cg_direction = cb_cg_vector / np.linalg.norm(cb_cg_vector)
        
        # CD는 CG에서 CB-CG 방향으로 약 1.54Å 떨어진 곳
        cd_coord = cg_coord + cb_cg_direction * 1.54
        
        return cd_coord
    
    def _substitute_to_alanine(self, residue):
        """잔기를 ALA로 안전하게 대체"""
        
        # 잔기 이름 변경
        residue.resname = 'ALA'
        
        # ALA에 필요한 원자만 남기기: N, CA, C, O, CB
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
            cb_coord = self._estimate_cb_for_alanine(ca_atom, residue)
            if cb_coord is not None:
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

    def _estimate_cb_for_alanine(self, ca_atom, residue) -> Optional[np.ndarray]:
        """ALA의 CB 원자 위치 추정 (개선된 버전)"""
        
        ca_coord = ca_atom.get_coord()
        
        # N과 C 원자 찾기
        n_atom = None
        c_atom = None
        
        for atom in residue:
            if atom.get_name().strip() == 'N':
                n_atom = atom
            elif atom.get_name().strip() == 'C':
                c_atom = atom
        
        if n_atom is None or c_atom is None:
            # N이나 C가 없으면 간단한 추정
            return ca_coord + np.array([1.54, 0.0, 0.0])
        
        # N-CA-C 평면에 수직인 방향으로 CB 배치
        n_coord = n_atom.get_coord()
        c_coord = c_atom.get_coord()
        
        ca_n = n_coord - ca_coord
        ca_c = c_coord - ca_coord
        
        # 외적으로 수직 벡터 구하기
        perpendicular = np.cross(ca_n, ca_c)
        if np.linalg.norm(perpendicular) == 0:
            return ca_coord + np.array([1.54, 0.0, 0.0])
        
        perpendicular = perpendicular / np.linalg.norm(perpendicular)
        cb_coord = ca_coord + perpendicular * 1.54
        
        return cb_coord

    # 나머지 헬퍼 메서드들...
    def _get_parent_bond_key(self, atom_name: str, resname: str) -> str:
        """원자와 부모 원자 간의 결합 키 반환"""
        if atom_name == 'CB':
            return 'CA-CB'
        elif atom_name == 'CG':
            return 'CB-CG'
        elif atom_name == 'CD':
            return 'CG-CD'
        elif atom_name in ['OE1', 'OE2', 'NE2']:
            return 'CD-OE1'  # 기본값
        elif atom_name in ['OD1', 'OD2', 'ND2']:
            return 'CG-OD1'  # 기본값
        elif atom_name == 'SG':
            return 'CB-SG'
        elif atom_name == 'OG':
            return 'CB-OG'
        else:
            return 'CA-CB'  # fallback
    
    def _get_parent_atom(self, atom_name: str, resname: str) -> str:
        """원자의 부모 원자 이름 반환"""
        if atom_name == 'CB':
            return 'CA'
        elif atom_name == 'CG':
            return 'CB'
        elif atom_name == 'CD':
            return 'CG'
        elif atom_name in ['OE1', 'OE2', 'NE2']:
            return 'CD'
        elif atom_name in ['OD1', 'OD2', 'ND2']:
            return 'CG'
        elif atom_name in ['SG', 'OG', 'OG1']:
            return 'CB'
        else:
            return 'CA'  # fallback

    def _get_sidechain_direction(self, present_atoms: Dict, ca_coord: np.ndarray) -> Optional[np.ndarray]:
        """사이드체인 방향 추정"""
        # 백본 원자들을 이용해 사이드체인 방향 추정
        if 'N' in present_atoms and 'C' in present_atoms:
            n_coord = present_atoms['N'].get_coord()
            c_coord = present_atoms['C'].get_coord()
            
            # N-CA와 CA-C의 이등분선에 수직인 방향
            n_ca = ca_coord - n_coord
            ca_c = c_coord - ca_coord
            
            if np.linalg.norm(n_ca) > 0 and np.linalg.norm(ca_c) > 0:
                n_ca = n_ca / np.linalg.norm(n_ca)
                ca_c = ca_c / np.linalg.norm(ca_c)
                
                bisector = n_ca + ca_c
                if np.linalg.norm(bisector) > 0:
                    bisector = bisector / np.linalg.norm(bisector)
                    # 백본 평면에 수직인 방향
                    perpendicular = np.cross(n_ca, ca_c)
                    if np.linalg.norm(perpendicular) > 0:
                        perpendicular = perpendicular / np.linalg.norm(perpendicular)
                        return perpendicular
        
        return None

    def _get_chain_direction(self, parent_atom_name: str, present_atoms: Dict) -> Optional[np.ndarray]:
        """체인 방향 추정"""
        # 간단한 체인 방향 추정
        if parent_atom_name == 'CB' and 'CA' in present_atoms:
            ca_coord = present_atoms['CA'].get_coord()
            cb_coord = present_atoms[parent_atom_name].get_coord()
            direction = cb_coord - ca_coord
            if np.linalg.norm(direction) > 0:
                return direction / np.linalg.norm(direction)
        
        return None


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
