#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Interface 기반 거리 계산 모듈

SuMD 시뮬레이션을 위한 개선된 거리 측정 시스템
- Interface residue 기반 정밀 거리 계산
- 기존 center of mass 방식과 호환성 유지
- 적응형 측정 방법 (Interface/Center of mass)

사용법:
    from interface_distance_calculator import DistanceCalculator
    
    distance = DistanceCalculator.calculate_chain_distance(
        "complex.pdb", "A", "B"
    )
"""

import numpy as np
from Bio.PDB import PDBParser
from typing import List, Tuple, Dict, Optional
import logging

class InterfaceDistanceCalculator:
    """Interface 기반 체인 간 거리 계산 클래스"""
    
    def __init__(self, interface_cutoff: float = 5.0, min_interface_residues: int = 2):
        """
        Interface 기반 거리 계산기 초기화
        
        Args:
            interface_cutoff: Interface 정의를 위한 거리 임계값 (Å)
            min_interface_residues: Interface로 인정하기 위한 최소 residue 수
        """
        self.interface_cutoff = interface_cutoff
        self.min_interface_residues = min_interface_residues
        
    def calculate_chain_distance(self, pdb_file: str, chain1_id: str, chain2_id: str) -> float:
        """
        두 체인 간의 interface 기반 거리 계산 (Å 단위)
        
        Interface가 충분하지 않으면 자동으로 center of mass 방식으로 전환
        
        Args:
            pdb_file: PDB 파일 경로
            chain1_id: 첫 번째 체인 ID (수용체)
            chain2_id: 두 번째 체인 ID (리간드)
            
        Returns:
            float: Interface 간 거리 (Å)
            
        Raises:
            Exception: PDB 파일 읽기 실패 또는 체인 없음
        """
        try:
            parser = PDBParser(QUIET=True)
            structure = parser.get_structure("structure", pdb_file)
            
            # 체인별 원자 및 residue 정보 수집
            chain1_atoms, chain1_residues = self._get_chain_atoms_and_residues(structure, chain1_id)
            chain2_atoms, chain2_residues = self._get_chain_atoms_and_residues(structure, chain2_id)
            
            if not chain1_atoms:
                raise Exception(f"체인 {chain1_id}를 찾을 수 없거나 원자가 없습니다")
            if not chain2_atoms:
                raise Exception(f"체인 {chain2_id}를 찾을 수 없거나 원자가 없습니다")
            
            # Interface residue 찾기
            interface1_residues, interface2_residues = self._find_interface_residues(
                chain1_residues, chain2_residues
            )
            
            # Interface가 충분하지 않으면 center of mass 방식으로 fallback
            if (len(interface1_residues) < self.min_interface_residues or 
                len(interface2_residues) < self.min_interface_residues):
                
                self._log_warning(f"Interface residue가 부족합니다 (체인 {chain1_id}: {len(interface1_residues)}, "
                                 f"체인 {chain2_id}: {len(interface2_residues)}). Center of mass 방식으로 fallback.")
                
                return self._calculate_center_of_mass_distance(chain1_atoms, chain2_atoms)
            
            # Interface 기반 거리 계산
            interface_distance = self._calculate_interface_distance(
                interface1_residues, interface2_residues
            )
            
            self._log_info(f"Interface 기반 거리: {interface_distance:.2f}Å "
                          f"(Interface residue: {len(interface1_residues)}+{len(interface2_residues)})")
            
            return interface_distance
            
        except Exception as e:
            raise Exception(f"거리 계산 실패: {e}")
    
    def get_interface_analysis(self, pdb_file: str, chain1_id: str, chain2_id: str) -> Dict:
        """
        Interface 상세 분석 정보 반환
        
        Args:
            pdb_file: PDB 파일 경로
            chain1_id: 첫 번째 체인 ID
            chain2_id: 두 번째 체인 ID
            
        Returns:
            dict: Interface 분석 결과
                - total_residues: 각 체인의 총 residue 수
                - interface_residues: 각 체인의 interface residue 수
                - interface_percentage: Interface 비율
                - interface_residue_list: Interface residue 목록
                - interface_distance: Interface 간 거리
                - min_interface_distance: 최소 interface 거리
                - avg_interface_distance: 평균 interface 거리
        """
        try:
            parser = PDBParser(QUIET=True)
            structure = parser.get_structure("structure", pdb_file)
            
            chain1_atoms, chain1_residues = self._get_chain_atoms_and_residues(structure, chain1_id)
            chain2_atoms, chain2_residues = self._get_chain_atoms_and_residues(structure, chain2_id)
            
            interface1, interface2 = self._find_interface_residues(chain1_residues, chain2_residues)
            
            # 상세 분석
            analysis = {
                'chain1_id': chain1_id,
                'chain2_id': chain2_id,
                'total_residues': {
                    chain1_id: len(chain1_residues),
                    chain2_id: len(chain2_residues)
                },
                'interface_residues': {
                    chain1_id: len(interface1),
                    chain2_id: len(interface2)
                },
                'interface_percentage': {
                    chain1_id: len(interface1) / len(chain1_residues) * 100 if chain1_residues else 0,
                    chain2_id: len(interface2) / len(chain2_residues) * 100 if chain2_residues else 0
                },
                'interface_residue_list': {
                    chain1_id: [f"{key[1]}{key[2]}" for key in interface1.keys()],
                    chain2_id: [f"{key[1]}{key[2]}" for key in interface2.keys()]
                }
            }
            
            if interface1 and interface2:
                analysis['interface_distance'] = self._calculate_interface_distance(interface1, interface2)
                analysis['min_interface_distance'] = self._calculate_min_interface_distance(interface1, interface2)
                analysis['avg_interface_distance'] = self._calculate_average_interface_distance(interface1, interface2)
            else:
                analysis['interface_distance'] = self._calculate_center_of_mass_distance(chain1_atoms, chain2_atoms)
                analysis['fallback_method'] = 'center_of_mass'
            
            return analysis
            
        except Exception as e:
            return {'error': str(e)}
    
    def _get_chain_atoms_and_residues(self, structure, chain_id: str) -> Tuple[List, Dict]:
        """체인의 원자와 residue 정보 수집"""
        atoms = []
        residues = {}
        
        for model in structure:
            for chain in model:
                if chain.id == chain_id:
                    for residue in chain:
                        res_key = (chain_id, residue.id[1], residue.get_resname())
                        residue_atoms = []
                        
                        for atom in residue:
                            atoms.append(atom.coord)
                            residue_atoms.append({
                                'name': atom.get_name(),
                                'coord': atom.coord,
                                'element': atom.element
                            })
                        
                        if residue_atoms:  # 원자가 있는 residue만 저장
                            residues[res_key] = {
                                'atoms': residue_atoms,
                                'center': np.mean([atom['coord'] for atom in residue_atoms], axis=0)
                            }
        
        return atoms, residues
    
    def _find_interface_residues(self, chain1_residues: Dict, chain2_residues: Dict) -> Tuple[Dict, Dict]:
        """Interface에 속하는 residue들 찾기"""
        interface1 = {}
        interface2 = {}
        
        # 각 체인1의 residue에 대해 체인2와의 최소 거리 확인
        for res1_key, res1_data in chain1_residues.items():
            min_distance = float('inf')
            
            for res2_key, res2_data in chain2_residues.items():
                min_atom_distance = self._calculate_min_atom_distance(
                    res1_data['atoms'], res2_data['atoms']
                )
                
                if min_atom_distance < min_distance:
                    min_distance = min_atom_distance
            
            # Interface 임계값 이내이면 interface residue로 분류
            if min_distance <= self.interface_cutoff:
                interface1[res1_key] = res1_data
                interface1[res1_key]['min_distance_to_other_chain'] = min_distance
        
        # 각 체인2의 residue에 대해 체인1과의 최소 거리 확인
        for res2_key, res2_data in chain2_residues.items():
            min_distance = float('inf')
            
            for res1_key, res1_data in chain1_residues.items():
                min_atom_distance = self._calculate_min_atom_distance(
                    res2_data['atoms'], res1_data['atoms']
                )
                
                if min_atom_distance < min_distance:
                    min_distance = min_atom_distance
            
            if min_distance <= self.interface_cutoff:
                interface2[res2_key] = res2_data
                interface2[res2_key]['min_distance_to_other_chain'] = min_distance
        
        return interface1, interface2
    
    def _calculate_min_atom_distance(self, atoms1: List[Dict], atoms2: List[Dict]) -> float:
        """두 residue 간 최소 원자 거리 계산"""
        min_distance = float('inf')
        
        for atom1 in atoms1:
            for atom2 in atoms2:
                distance = np.linalg.norm(atom1['coord'] - atom2['coord'])
                if distance < min_distance:
                    min_distance = distance
        
        return min_distance
    
    def _calculate_interface_distance(self, interface1: Dict, interface2: Dict) -> float:
        """Interface 간 거리 계산 - 적응형 방법 선택"""
        
        # Interface residue들의 center of mass 간 거리
        interface1_center = self._calculate_interface_center(interface1)
        interface2_center = self._calculate_interface_center(interface2)
        center_distance = np.linalg.norm(interface1_center - interface2_center)
        
        # Interface 원자들 간 최소 거리
        min_interface_distance = self._calculate_min_interface_distance(interface1, interface2)
        
        # Interface residue들 간 평균 거리
        avg_interface_distance = self._calculate_average_interface_distance(interface1, interface2)
        
        # 적응형 방법 선택
        # 최소 거리가 매우 작으면 (접촉 상태) 평균 거리 사용
        # 그렇지 않으면 center 간 거리 사용 (더 안정적)
        if min_interface_distance < 3.0:
            return avg_interface_distance
        else:
            return center_distance
    
    def _calculate_interface_center(self, interface_residues: Dict) -> np.ndarray:
        """Interface residue들의 center of mass 계산"""
        all_coords = []
        for res_data in interface_residues.values():
            for atom in res_data['atoms']:
                all_coords.append(atom['coord'])
        
        return np.mean(all_coords, axis=0)
    
    def _calculate_min_interface_distance(self, interface1: Dict, interface2: Dict) -> float:
        """Interface 원자들 간 최소 거리"""
        min_distance = float('inf')
        
        for res1_data in interface1.values():
            for atom1 in res1_data['atoms']:
                for res2_data in interface2.values():
                    for atom2 in res2_data['atoms']:
                        distance = np.linalg.norm(atom1['coord'] - atom2['coord'])
                        if distance < min_distance:
                            min_distance = distance
        
        return min_distance
    
    def _calculate_average_interface_distance(self, interface1: Dict, interface2: Dict) -> float:
        """Interface residue들 간 평균 거리"""
        distances = []
        
        for res1_key, res1_data in interface1.items():
            for res2_key, res2_data in interface2.items():
                center_distance = np.linalg.norm(res1_data['center'] - res2_data['center'])
                distances.append(center_distance)
        
        return np.mean(distances) if distances else float('inf')
    
    def _calculate_center_of_mass_distance(self, chain1_atoms: List, chain2_atoms: List) -> float:
        """Fallback: 기존 center of mass 방식"""
        chain1_center = np.mean(chain1_atoms, axis=0)
        chain2_center = np.mean(chain2_atoms, axis=0)
        return float(np.linalg.norm(chain1_center - chain2_center))
    
    def _log_info(self, message: str):
        """안전한 정보 로깅"""
        try:
            logging.info(message)
        except:
            pass  # 로깅 설정이 없어도 오류 없이 계속 진행
    
    def _log_warning(self, message: str):
        """안전한 경고 로깅"""
        try:
            logging.warning(message)
        except:
            pass  # 로깅 설정이 없어도 오류 없이 계속 진행


class DistanceCalculator:
    """
    체인 간 거리 계산 클래스 - Interface 기반으로 개선됨
    
    기존 SuMD 코드와의 호환성을 위한 래퍼 클래스
    """
    
    @staticmethod
    def calculate_chain_distance(pdb_file: str, chain1_id: str, chain2_id: str) -> float:
        """
        두 체인 간의 interface 기반 거리 계산 (Å 단위)
        
        기존 API 호환성을 위한 정적 메서드
        Interface가 충분하지 않으면 자동으로 center of mass 방식으로 전환
        
        Args:
            pdb_file: PDB 파일 경로
            chain1_id: 첫 번째 체인 ID (수용체)
            chain2_id: 두 번째 체인 ID (리간드)
            
        Returns:
            float: 체인 간 거리 (Å)
            
        Raises:
            Exception: PDB 파일 읽기 실패 또는 체인 없음
        """
        calculator = InterfaceDistanceCalculator(
            interface_cutoff=5.0,
            min_interface_residues=2
        )
        return calculator.calculate_chain_distance(pdb_file, chain1_id, chain2_id)


def compare_distance_methods(pdb_file: str, chain1_id: str, chain2_id: str) -> Optional[Dict]:
    """
    기존 center of mass 방식과 새로운 interface 방식 비교
    
    디버깅 및 검증 목적으로 사용
    
    Args:
        pdb_file: PDB 파일 경로
        chain1_id: 첫 번째 체인 ID
        chain2_id: 두 번째 체인 ID
        
    Returns:
        dict: 비교 결과 또는 None (실패 시)
    """
    
    def old_center_of_mass_method():
        """기존 center of mass 방식"""
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("structure", pdb_file)
        
        chain1_atoms = []
        chain2_atoms = []
        
        for model in structure:
            for chain in model:
                if chain.id == chain1_id:
                    for residue in chain:
                        for atom in residue:
                            chain1_atoms.append(atom.coord)
                elif chain.id == chain2_id:
                    for residue in chain:
                        for atom in residue:
                            chain2_atoms.append(atom.coord)
        
        if chain1_atoms and chain2_atoms:
            chain1_center = np.mean(chain1_atoms, axis=0)
            chain2_center = np.mean(chain2_atoms, axis=0)
            return float(np.linalg.norm(chain1_center - chain2_center))
        return float('inf')
    
    try:
        # 기존 방식
        old_distance = old_center_of_mass_method()
        
        # 새로운 방식
        calculator = InterfaceDistanceCalculator()
        new_distance = calculator.calculate_chain_distance(pdb_file, chain1_id, chain2_id)
        analysis = calculator.get_interface_analysis(pdb_file, chain1_id, chain2_id)
        
        print(f"=== 거리 측정 방식 비교 ===")
        print(f"기존 (Center of Mass): {old_distance:.2f}Å")
        print(f"신규 (Interface 기반): {new_distance:.2f}Å")
        print(f"차이: {abs(old_distance - new_distance):.2f}Å")
        
        if 'interface_residues' in analysis:
            print(f"\n=== Interface 분석 ===")
            print(f"Interface residue 수: {analysis['interface_residues']}")
            print(f"Interface 비율: {analysis['interface_percentage']}")
            print(f"Interface residue 목록:")
            print(f"  {chain1_id}: {analysis['interface_residue_list'][chain1_id]}")
            print(f"  {chain2_id}: {analysis['interface_residue_list'][chain2_id]}")
        
        return {
            'old_distance': old_distance,
            'new_distance': new_distance,
            'difference': abs(old_distance - new_distance),
            'analysis': analysis
        }
        
    except Exception as e:
        print(f"비교 실패: {e}")
        return None


if __name__ == "__main__":
    # 테스트 코드
    print("Interface Distance Calculator 모듈")
    print("사용법:")
    print("  from interface_distance_calculator import DistanceCalculator")
    print("  distance = DistanceCalculator.calculate_chain_distance('complex.pdb', 'A', 'B')")
    print("  # 또는")
    print("  from interface_distance_calculator import compare_distance_methods")
    print("  compare_distance_methods('complex.pdb', 'A', 'B')")
