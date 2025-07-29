#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
구조 변형 생성기 V2 - 50Å 거리 기반 변형

Golden Standard PDB 파일에서 특정 체인을 정확히 50Å 거리로 이동시켜
변형된 구조들을 생성하는 도구입니다.

주요 특징:
1. 두 체인 간 최소 거리가 정확히 50Å이 되도록 이동
2. 이동 방향만 랜덤 (벡터 랜덤)
3. 회전 없음 (이동만)
4. Clash 검사 및 폐기
5. 이동 경로상 다른 단백질 존재 시 폐기

사용법:
    python structure_displacer.py --input complex.pdb --target_chains A,B --num_variants 5
"""

import os
import sys
import argparse
import numpy as np
import logging
import multiprocessing as mp
import random
import time
from typing import List, Tuple, Dict, Optional, Union
from Bio.PDB import PDBParser, PDBIO, Atom
from Bio.PDB.Structure import Structure
from Bio.PDB.Model import Model
from Bio.PDB.Chain import Chain
import shutil
from datetime import datetime


class DisplacementConfig:
    """구조 변형 설정 클래스"""
    
    def __init__(self):
        # 실제 이동할 거리 (Angstrom) - 60Å로 이동하여 50Å 근처에 도달
        self.actual_displacement_distance = 60.0
        
        # 적응형 이동 거리 사용 여부
        self.use_adaptive_displacement = True   # True: 기하학적 계산 기반, False: 고정 거리
        
        # 이동 방향 다양화 각도 (라디안) - 15도 단위 다양화
        self.rotation_angle_range = np.pi / 12  # ±15도 방향 다양화
        
        # Clash 검사 파라미터
        self.clash_threshold = 2.0              # Clash 판정 거리 (Å)
        self.max_attempts = 200                 # 최대 시도 횟수
        
        # 이동 체인 선택 (첫 번째 체인을 이동)
        self.move_first_chain = True
        
        # 경로 차단 검사 파라미터
        self.path_check_threshold = 5.0         # 경로상 장애물 판정 거리 (Å)
        self.path_check_resolution = 1.0        # 경로 체크 해상도 (Å)
        
        # 유효성 검사 파라미터 - 최종 목표는 50Å
        self.target_distance_for_validation = 50.0  # 유효성 검사용 목표 거리


class MinDistanceCalculator:
    """최소 거리 계산기"""
    
    @staticmethod
    def calculate_min_distance_between_chains(chain1: Chain, chain2: Chain) -> Tuple[float, Tuple[np.ndarray, np.ndarray]]:
        """
        두 체인 간의 최소 거리와 해당 원자 좌표들 반환 (단백질 원자만, HETATM 제외)
        
        Returns:
            Tuple[float, Tuple[np.ndarray, np.ndarray]]: (최소거리, (chain1_atom_coord, chain2_atom_coord))
        """
        min_distance = float('inf')
        closest_coords = (None, None)
        
        # 표준 아미노산 잔기만 선택 (물 분자 및 HETATM 제외)
        atoms1 = []
        atoms2 = []
        
        for residue in chain1:
            if residue.id[0] == ' ':  # 표준 아미노산만 (HETATM 제외)
                atoms1.extend(list(residue.get_atoms()))
        
        for residue in chain2:
            if residue.id[0] == ' ':  # 표준 아미노산만 (HETATM 제외)
                atoms2.extend(list(residue.get_atoms()))
        
        for atom1 in atoms1:
            for atom2 in atoms2:
                distance = np.linalg.norm(atom1.coord - atom2.coord)
                if distance < min_distance:
                    min_distance = distance
                    closest_coords = (atom1.coord.copy(), atom2.coord.copy())
        
        return min_distance, closest_coords
    
    @staticmethod
    def calculate_chain_center(chain: Chain) -> np.ndarray:
        """
        체인의 중심 좌표 계산 (단백질 원자만)
        
        Returns:
            np.ndarray: 중심 좌표
        """
        atoms = []
        for residue in chain:
            if residue.id[0] == ' ':  # 표준 아미노산만
                atoms.extend([atom.coord for atom in residue.get_atoms()])
        
        if not atoms:
            return np.array([0.0, 0.0, 0.0])
        
        return np.mean(atoms, axis=0)
    
    @staticmethod
    def calculate_directional_surface_distance(chain: Chain, center: np.ndarray, direction: np.ndarray) -> float:
        """
        특정 방향에서 중심으로부터 표면까지의 최대 거리 계산
        
        Args:
            chain: 체인
            center: 체인 중심 좌표
            direction: 방향 벡터
            
        Returns:
            float: 방향별 표면 거리
        """
        max_distance = 0.0
        
        # 표준 아미노산 원자만 선택
        atoms = []
        for residue in chain:
            if residue.id[0] == ' ':
                atoms.extend(list(residue.get_atoms()))
        
        if not atoms:
            return 0.0
        
        # 방향 벡터 정규화
        direction_norm = direction / np.linalg.norm(direction)
        
        for atom in atoms:
            # 중심에서 원자로의 벡터
            atom_vector = atom.coord - center
            # 방향 벡터와의 내적 (투영 길이)
            projection = np.dot(atom_vector, direction_norm)
            
            # 해당 방향으로의 거리만 고려 (양수만)
            if projection > max_distance:
                max_distance = projection
        
        return float(max_distance)
    
    @staticmethod
    def calculate_required_displacement_distance(chain1: Chain, chain2: Chain, target_surface_distance: float) -> float:
        """
        목표 표면간 거리를 달성하기 위해 필요한 이동 거리 계산
        
        Args:
            chain1: 이동할 체인
            chain2: 고정된 체인  
            target_surface_distance: 목표 표면간 거리
            
        Returns:
            float: 필요한 이동 거리
        """
        # 각 체인의 중심 계산
        center1 = MinDistanceCalculator.calculate_chain_center(chain1)
        center2 = MinDistanceCalculator.calculate_chain_center(chain2)
        
        # 방향 벡터 계산
        direction_1_to_2 = center2 - center1
        current_center_distance = np.linalg.norm(direction_1_to_2)
        
        if current_center_distance == 0:
            return target_surface_distance
        
        # 각 체인의 방향별 표면 거리 계산
        surface_dist_1 = MinDistanceCalculator.calculate_directional_surface_distance(
            chain1, center1, direction_1_to_2
        )
        surface_dist_2 = MinDistanceCalculator.calculate_directional_surface_distance(
            chain2, center2, -direction_1_to_2
        )
        
        # 목표 표면간 거리를 위한 필요한 중심간 거리
        required_center_distance = target_surface_distance + surface_dist_1 + surface_dist_2
        
        # 추가 이동 필요 거리
        additional_displacement = required_center_distance - current_center_distance
        
        return max(0.0, additional_displacement)
    
    @staticmethod
    def calculate_surface_centers(chain1: Chain, chain2: Chain) -> Tuple[np.ndarray, np.ndarray]:
        """
        두 체인의 서로를 향한 표면 중심점 계산
        """
        atoms1 = np.array([atom.coord for atom in chain1.get_atoms()])
        atoms2 = np.array([atom.coord for atom in chain2.get_atoms()])
        
        center1 = np.mean(atoms1, axis=0)
        center2 = np.mean(atoms2, axis=0)
        
        # 각 체인에서 상대방에 가장 가까운 표면 영역 찾기
        direction_12 = center2 - center1
        direction_12 = direction_12 / np.linalg.norm(direction_12)
        
        direction_21 = center1 - center2
        direction_21 = direction_21 / np.linalg.norm(direction_21)
        
        # 각 체인에서 상대방 방향의 원자들만 선택
        surface_atoms1 = []
        surface_atoms2 = []
        
        for atom_coord in atoms1:
            if np.dot(atom_coord - center1, direction_12) > 0:  # 상대방 방향
                surface_atoms1.append(atom_coord)
        
        for atom_coord in atoms2:
            if np.dot(atom_coord - center2, direction_21) > 0:  # 상대방 방향
                surface_atoms2.append(atom_coord)
        
        if surface_atoms1 and surface_atoms2:
            surface_center1 = np.mean(surface_atoms1, axis=0)
            surface_center2 = np.mean(surface_atoms2, axis=0)
        else:
            surface_center1 = center1
            surface_center2 = center2
        
        return surface_center1, surface_center2


class ClashChecker:
    """충돌 검사 클래스"""
    
    def __init__(self, threshold: float = 2.0):
        self.threshold = threshold
    
    def check_chain_clash(self, moved_chain: Chain, other_chains: List[Chain]) -> Tuple[bool, float, str]:
        """
        이동된 체인과 다른 체인들 간의 충돌 검사
        
        Returns:
            Tuple[bool, float, str]: (충돌 여부, 최소 거리, 세부 정보)
        """
        min_distance = float('inf')
        clash_info = ""
        
        moved_atoms = list(moved_chain.get_atoms())
        
        for other_chain in other_chains:
            other_atoms = list(other_chain.get_atoms())
            
            for moved_atom in moved_atoms:
                for other_atom in other_atoms:
                    distance = np.linalg.norm(moved_atom.coord - other_atom.coord)
                    
                    if distance < min_distance:
                        min_distance = distance
                        clash_info = f"Chain {moved_chain.id} atom {moved_atom.get_fullname()} vs Chain {other_chain.id} atom {other_atom.get_fullname()}"
                    
                    if distance < self.threshold:
                        return True, distance, clash_info
        
        return False, min_distance, clash_info


class PathChecker:
    """이동 경로 차단 검사 클래스"""
    
    def __init__(self, threshold: float = 5.0, resolution: float = 1.0):
        self.threshold = threshold
        self.resolution = resolution
    
    def check_path_obstruction(self, start_center: np.ndarray, end_center: np.ndarray, 
                             other_chains: List[Chain]) -> Tuple[bool, str]:
        """
        이동 경로상에 다른 체인이 있는지 검사
        
        Args:
            start_center: 시작 위치 (이동 전 체인 중심)
            end_center: 끝 위치 (이동 후 체인 중심)
            other_chains: 검사할 다른 체인들
            
        Returns:
            Tuple[bool, str]: (경로 차단 여부, 차단 정보)
        """
        if not other_chains:
            return False, "No other chains to check"
        
        # 이동 벡터와 거리
        movement_vector = end_center - start_center
        movement_distance = np.linalg.norm(movement_vector)
        
        if movement_distance == 0:
            return False, "No movement"
        
        movement_direction = movement_vector / movement_distance
        
        # 경로상의 점들을 체크
        num_steps = int(movement_distance / self.resolution)
        
        for step in range(1, num_steps):  # 시작점과 끝점 제외
            check_point = start_center + (step * self.resolution) * movement_direction
            
            # 각 다른 체인과 이 점 사이의 최소 거리 확인
            for other_chain in other_chains:
                other_atoms = list(other_chain.get_atoms())
                
                for atom in other_atoms:
                    distance = np.linalg.norm(atom.coord - check_point)
                    
                    if distance < self.threshold:
                        obstruction_info = f"Chain {other_chain.id} blocks path at distance {distance:.2f}Å"
                        return True, obstruction_info
        
        return False, "Path clear"


class AdvancedStructureDisplacer:
    """고급 구조 변형 클래스 - 50Å 거리 기반"""
    
    def __init__(self, config: DisplacementConfig = None, logger=None):
        self.config = config or DisplacementConfig()
        self.logger = logger or logging.getLogger(__name__)
        self.parser = PDBParser(QUIET=True)
        self.clash_checker = ClashChecker(self.config.clash_threshold)
        self.path_checker = PathChecker(self.config.path_check_threshold, self.config.path_check_resolution)
        self.min_dist_calc = MinDistanceCalculator()
    
    def generate_random_direction_vector(self) -> np.ndarray:
        """랜덤 방향 벡터 생성 (단위벡터)"""
        # 구면 좌표계에서 균등 분포 생성
        theta = np.random.uniform(0, 2 * np.pi)  # 방위각
        phi = np.random.uniform(0, np.pi)        # 극각
        
        # 구면 좌표를 카르테시안 좌표로 변환
        x = np.sin(phi) * np.cos(theta)
        y = np.sin(phi) * np.sin(theta)
        z = np.cos(phi)
        
        return np.array([x, y, z])
    
    def calculate_displacement_vector(self, chain1: Chain, chain2: Chain) -> Tuple[np.ndarray, Dict]:
        """
        표면에서부터 고정된 거리로 이동하는 변위 벡터 계산
        
        Args:
            chain1: 이동할 체인
            chain2: 고정된 체인
            
        Returns:
            Tuple[np.ndarray, Dict]: (변위 벡터, 계산 정보)
        """
        # 현재 최소 거리와 가장 가까운 원자들의 좌표
        current_min_dist, (closest_coord1, closest_coord2) = self.min_dist_calc.calculate_min_distance_between_chains(chain1, chain2)
        
        # 표면 중심점들 계산 
        surface_center1, surface_center2 = self.min_dist_calc.calculate_surface_centers(chain1, chain2)
        
        # 기본 방향: chain1에서 chain2로부터 멀어지는 방향
        basic_direction = surface_center1 - surface_center2
        if np.linalg.norm(basic_direction) > 0:
            basic_direction = basic_direction / np.linalg.norm(basic_direction)
        else:
            # fallback: 랜덤 방향
            basic_direction = self.generate_random_direction_vector()
        
        # 회전 각도 적용 (config에서 설정, 기본값 0)
        if self.config.rotation_angle_range > 0:
            random_angle = np.random.uniform(-self.config.rotation_angle_range, self.config.rotation_angle_range)
        else:
            random_angle = 0.0  # 회전 없음
        
        # 회전 적용 (random_angle이 0이면 회전하지 않음)
        if random_angle != 0.0:
            # 랜덤 회전축 생성 (기본 방향에 수직)
            perpendicular = np.cross(basic_direction, np.array([0, 0, 1]))
            if np.linalg.norm(perpendicular) < 1e-6:  # 기본 방향이 z축과 평행한 경우
                perpendicular = np.cross(basic_direction, np.array([1, 0, 0]))
            perpendicular = perpendicular / np.linalg.norm(perpendicular)
            
            # Rodriguez 회전 공식으로 기본 방향을 회전
            cos_angle = np.cos(random_angle)
            sin_angle = np.sin(random_angle)
            
            final_direction = (basic_direction * cos_angle + 
                             np.cross(perpendicular, basic_direction) * sin_angle +
                             perpendicular * np.dot(perpendicular, basic_direction) * (1 - cos_angle))
        else:
            final_direction = basic_direction
        
        # 정규화
        final_direction = final_direction / np.linalg.norm(final_direction)
        
        # 적응형 이동 거리 계산 (기하학적 계산 기반)
        if self.config.use_adaptive_displacement:
            # 목표 표면간 거리를 달성하기 위한 정확한 이동 거리 계산
            required_displacement = self.min_dist_calc.calculate_required_displacement_distance(
                chain1, chain2, self.config.target_distance_for_validation
            )
            actual_target_distance = current_min_dist + required_displacement
            
            self.logger.info(f"적응형 이동: 현재 거리 {current_min_dist:.2f}Å → 목표 {self.config.target_distance_for_validation}Å, 이동 필요 {required_displacement:.2f}Å")
        else:
            # 기존 고정 거리 방식
            actual_target_distance = self.config.actual_displacement_distance
            required_displacement = actual_target_distance - current_min_dist
        
        # 변위 벡터 계산
        displacement_vector = final_direction * required_displacement
        
        calc_info = {
            'current_min_distance': current_min_dist,
            'actual_displacement_distance': self.config.actual_displacement_distance,
            'validation_target_distance': self.config.target_distance_for_validation,
            'required_displacement': required_displacement,
            'displacement_magnitude': np.linalg.norm(displacement_vector),
            'displacement_direction': final_direction.tolist(),
            'basic_direction': basic_direction.tolist(),
            'rotation_angle_deg': np.degrees(random_angle),
            'rotation_applied': random_angle != 0.0,
            'surface_center1': surface_center1.tolist(),
            'surface_center2': surface_center2.tolist()
        }
        
        return displacement_vector, calc_info
    
    def apply_displacement(self, chain: Chain, displacement_vector: np.ndarray) -> Chain:
        """체인에 변위 적용"""
        new_chain = Chain(chain.id)
        
        for residue in chain:
            new_residue = residue.copy()
            new_residue.detach_parent()
            
            for atom in new_residue:
                atom.coord = atom.coord + displacement_vector
            
            new_chain.add(new_residue)
        
        return new_chain
    
    def validate_displaced_structure(self, displaced_chain: Chain, fixed_chain: Chain, 
                                   other_chains: List[Chain], original_chain_center: np.ndarray) -> Tuple[bool, Dict]:
        """
        변형된 구조의 유효성 검사
        
        Returns:
            Tuple[bool, Dict]: (유효 여부, 검사 결과)
        """
        validation_result = {
            'valid': True,
            'clash_check': {},
            'path_check': {},
            'distance_check': {}
        }
        
        # 1. Clash 검사
        has_clash, min_clash_dist, clash_info = self.clash_checker.check_chain_clash(
            displaced_chain, [fixed_chain] + other_chains
        )
        
        validation_result['clash_check'] = {
            'has_clash': has_clash,
            'min_distance': min_clash_dist,
            'clash_info': clash_info
        }
        
        if has_clash:
            validation_result['valid'] = False
            validation_result['failure_reason'] = f"Clash detected: {clash_info}"
            return False, validation_result
        
        # 2. 목표 거리 달성 확인
        final_min_dist, _ = self.min_dist_calc.calculate_min_distance_between_chains(displaced_chain, fixed_chain)
        distance_error = abs(final_min_dist - self.config.target_distance_for_validation)
        
        validation_result['distance_check'] = {
            'final_distance': final_min_dist,
            'target_distance': self.config.target_distance_for_validation,
            'distance_error': distance_error,
            'acceptable_error': 5.0  # 5Å 오차 허용
        }
        
        if distance_error > 15.0:  # 15Å 이상 오차는 실패로 간주 (매우 관대하게)
            validation_result['valid'] = False
            validation_result['failure_reason'] = f"Distance error too large: {distance_error:.2f}Å"
            return False, validation_result
        
        # 3. 경로 차단 검사
        displaced_center = np.mean([atom.coord for atom in displaced_chain.get_atoms()], axis=0)
        
        is_obstructed, obstruction_info = self.path_checker.check_path_obstruction(
            original_chain_center, displaced_center, other_chains
        )
        
        validation_result['path_check'] = {
            'is_obstructed': is_obstructed,
            'obstruction_info': obstruction_info,
            'original_center': original_chain_center.tolist(),
            'displaced_center': displaced_center.tolist()
        }
        
        if is_obstructed:
            validation_result['valid'] = False
            validation_result['failure_reason'] = f"Path obstructed: {obstruction_info}"
            return False, validation_result
        
        return True, validation_result
    
    def displace_structure(self, input_pdb: str, target_chains: List[str], 
                         output_pdb: str = None) -> Tuple[str, Dict]:
        """
        구조 변형 실행
        
        Args:
            input_pdb: 입력 PDB 파일
            target_chains: [이동할 체인, 고정된 체인] 
            output_pdb: 출력 PDB 파일
            
        Returns:
            Tuple[str, Dict]: (출력 파일 경로, 변형 정보)
        """
        if len(target_chains) != 2:
            raise ValueError("정확히 2개의 체인이 필요합니다 (이동할 체인, 고정된 체인)")
        
        move_chain_id, fixed_chain_id = target_chains
        
        try:
            # PDB 파일 로드
            structure = self.parser.get_structure("structure", input_pdb)
            
            # 출력 파일명 자동 생성
            if output_pdb is None:
                base_name = os.path.splitext(os.path.basename(input_pdb))[0]
                timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
                output_pdb = f"{base_name}_displaced_{timestamp}.pdb"
            
            # 체인 수집
            move_chain = None
            fixed_chain = None
            other_chains = []
            
            for model in structure:
                for chain in model:
                    if chain.id == move_chain_id:
                        move_chain = chain
                    elif chain.id == fixed_chain_id:
                        fixed_chain = chain
                    else:
                        other_chains.append(chain)
            
            if move_chain is None:
                raise ValueError(f"이동할 체인 {move_chain_id}를 찾을 수 없습니다")
            if fixed_chain is None:
                raise ValueError(f"고정된 체인 {fixed_chain_id}를 찾을 수 없습니다")
            
            # 원래 이동 체인의 중심점 저장
            original_center = np.mean([atom.coord for atom in move_chain.get_atoms()], axis=0)
            
            # 초기 거리 확인
            initial_distance, _ = self.min_dist_calc.calculate_min_distance_between_chains(move_chain, fixed_chain)
            
            displacement_info = {
                'input_pdb': input_pdb,
                'target_chains': target_chains,
                'move_chain': move_chain_id,
                'fixed_chain': fixed_chain_id,
                'other_chains': [chain.id for chain in other_chains],
                'initial_distance': initial_distance,
                'actual_displacement_distance': self.config.actual_displacement_distance,
                'target_distance_for_validation': self.config.target_distance_for_validation,
                'attempts': [],
                'final_result': None
            }
            
            self.logger.info(f"구조 변형 시작: {move_chain_id} 체인을 {self.config.actual_displacement_distance}Å 거리로 이동 (목표: {self.config.target_distance_for_validation}Å)")
            self.logger.info(f"초기 거리: {initial_distance:.2f}Å")
            
            # 여러 번 시도
            successful = False
            
            for attempt in range(self.config.max_attempts):
                try:
                    # 변위 벡터 계산
                    displacement_vector, calc_info = self.calculate_displacement_vector(
                        move_chain, fixed_chain
                    )
                    
                    # 변위 적용
                    displaced_chain = self.apply_displacement(move_chain, displacement_vector)
                    
                    # 유효성 검사
                    is_valid, validation_result = self.validate_displaced_structure(
                        displaced_chain, fixed_chain, other_chains, original_center
                    )
                    
                    attempt_info = {
                        'attempt_number': attempt + 1,
                        'calculation_info': calc_info,
                        'validation_result': validation_result,
                        'success': is_valid
                    }
                    
                    displacement_info['attempts'].append(attempt_info)
                    
                    if is_valid:
                        # 성공: 새로운 구조 생성 및 저장
                        new_structure = Structure("displaced")
                        new_model = Model(0)
                        new_structure.add(new_model)
                        
                        # 변형된 체인 추가
                        displaced_chain.detach_parent()
                        new_model.add(displaced_chain)
                        
                        # 고정된 체인 추가
                        fixed_chain_copy = fixed_chain.copy()
                        fixed_chain_copy.detach_parent()
                        new_model.add(fixed_chain_copy)
                        
                        # 다른 체인들 추가
                        for other_chain in other_chains:
                            other_chain_copy = other_chain.copy()
                            other_chain_copy.detach_parent()
                            new_model.add(other_chain_copy)
                        
                        # PDB 저장
                        io = PDBIO()
                        io.set_structure(new_structure)
                        io.save(output_pdb)
                        
                        # 최종 결과 정보
                        final_distance = validation_result['distance_check']['final_distance']
                        
                        displacement_info['final_result'] = {
                            'success': True,
                            'output_file': output_pdb,
                            'final_distance': final_distance,
                            'distance_error': validation_result['distance_check']['distance_error'],
                            'total_attempts': attempt + 1,
                            'displacement_vector': displacement_vector.tolist(),
                            'displacement_magnitude': np.linalg.norm(displacement_vector)
                        }
                        
                        successful = True
                        self.logger.info(f"구조 변형 성공 (시도 {attempt + 1}회): 최종 거리 {final_distance:.2f}Å")
                        break
                    
                    else:
                        failure_reason = validation_result.get('failure_reason', 'Unknown')
                        if (attempt + 1) % 5 == 0 or attempt == 0:  # 첫 번째와 5회마다 로그 출력
                            # 자세한 정보 출력
                            distance_info = validation_result.get('distance_check', {})
                            clash_info = validation_result.get('clash_check', {})
                            path_info = validation_result.get('path_check', {})
                            
                            self.logger.info(f"시도 {attempt + 1} 실패: {failure_reason}")
                            self.logger.info(f"  거리: {distance_info.get('final_distance', 'N/A'):.2f}Å (목표: {distance_info.get('target_distance', 'N/A')}Å, 오차: {distance_info.get('distance_error', 'N/A'):.2f}Å)")
                            self.logger.info(f"  Clash: {clash_info.get('has_clash', 'N/A')} (최소거리: {clash_info.get('min_distance', 'N/A'):.2f}Å)")
                            self.logger.info(f"  경로차단: {path_info.get('is_obstructed', 'N/A')}")
                
                except Exception as e:
                    self.logger.info(f"시도 {attempt + 1} 중 오류: {e}")
                    continue
            
            if not successful:
                displacement_info['final_result'] = {
                    'success': False,
                    'error': f"최대 시도 횟수 {self.config.max_attempts}회 도달",
                    'total_attempts': self.config.max_attempts
                }
                raise Exception(f"구조 변형 실패: 최대 시도 횟수 {self.config.max_attempts}회 도달")
            
            return output_pdb, displacement_info
            
        except Exception as e:
            self.logger.error(f"구조 변형 실패: {e}")
            raise


class MultipleDisplacer:
    """다중 변형 구조 생성 클래스"""
    
    def __init__(self, config: DisplacementConfig = None, logger=None):
        self.config = config or DisplacementConfig()
        self.logger = logger or logging.getLogger(__name__)
        self.single_displacer = AdvancedStructureDisplacer(config, logger)
    
    def generate_multiple_variants(self, input_pdb: str, target_chains: List[str], 
                                 num_variants: int = 5, output_dir: str = None,
                                 output_prefix: str = None) -> List[Tuple[str, Dict]]:
        """여러 변형된 구조 생성"""
        if output_dir is None:
            output_dir = os.path.dirname(input_pdb) or "."
        
        if output_prefix is None:
            base_name = os.path.splitext(os.path.basename(input_pdb))[0]
            output_prefix = f"{base_name}_displaced"
        
        os.makedirs(output_dir, exist_ok=True)
        
        results = []
        successful_count = 0
        
        self.logger.info(f"다중 변형 구조 생성 시작: {num_variants}개 목표")
        
        for i in range(1, num_variants + 1):
            try:
                output_pdb = os.path.join(output_dir, f"{output_prefix}_{i:03d}.pdb")
                
                # 각 변형마다 서로 다른 랜덤 시드
                seed_value = (int(time.time() * 1000) + i) % (2**32 - 1)
                np.random.seed(seed_value)
                random.seed(seed_value)
                
                self.logger.info(f"변형 구조 {i}/{num_variants} 생성 중...")
                
                result = self.single_displacer.displace_structure(
                    input_pdb, target_chains, output_pdb
                )
                
                results.append(result)
                successful_count += 1
                
                # 결과 요약
                _, info = result
                final_result = info['final_result']
                
                self.logger.info(f"변형 구조 {i} 완료: "
                               f"최종 거리 {final_result['final_distance']:.2f}Å, "
                               f"시도 횟수 {final_result['total_attempts']}회")
                
            except Exception as e:
                self.logger.error(f"변형 구조 {i} 생성 실패: {e}")
                continue
        
        self.logger.info(f"다중 변형 구조 생성 완료: 성공 {successful_count}/{num_variants}개")
        
        return results


def main():
    """메인 함수"""
    parser = argparse.ArgumentParser(
        description="고급 구조 변형 생성기 - 50Å 거리 기반 변형",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
특징:
- 두 체인 간 최소 거리가 정확히 50Å이 되도록 이동
- 이동 방향만 랜덤 (벡터 랜덤)
- 회전 없음 (이동만)
- Clash 검사 및 폐기
- 이동 경로상 다른 단백질 존재 시 폐기

사용 예시:
  python structure_displacer.py --input complex.pdb --target_chains A,B --num_variants 5
  python structure_displacer.py --input complex.pdb --target_chains A,B --target_distance 60 --clash_threshold 2.5
        """
    )
    
    # 필수 인수
    parser.add_argument("--input", required=True, help="입력 PDB 파일")
    parser.add_argument("--target_chains", required=True, 
                       help="타겟 체인 ID (쉼표로 구분, 첫 번째가 이동할 체인, 예: A,B)")
    
    # 기본 파라미터
    parser.add_argument("--num_variants", type=int, default=5,
                       help="생성할 변형 구조 수 (기본값: 5)")
    parser.add_argument("--output_dir", default="displaced_structures_50A",
                       help="출력 디렉토리 (기본값: displaced_structures_50A)")
    
    # 변형 파라미터
    parser.add_argument("--target_distance", type=float, default=50.0,
                       help="목표 거리 (Å, 기본값: 50.0)")
    parser.add_argument("--clash_threshold", type=float, default=2.0,
                       help="Clash 판정 거리 (Å, 기본값: 2.0)")
    parser.add_argument("--max_attempts", type=int, default=200,
                       help="최대 시도 횟수 (기본값: 200)")
    parser.add_argument("--path_threshold", type=float, default=5.0,
                       help="경로 차단 판정 거리 (Å, 기본값: 5.0)")
    
    # 기타
    parser.add_argument("--verbose", "-v", action="store_true", help="상세 로그 출력")
    parser.add_argument("--seed", type=int, help="랜덤 시드")
    
    args = parser.parse_args()
    
    # 로깅 설정
    level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )
    logger = logging.getLogger("AdvancedDisplacer")
    
    # 입력 파일 확인
    if not os.path.exists(args.input):
        logger.error(f"입력 파일을 찾을 수 없습니다: {args.input}")
        return 1
    
    # 타겟 체인 파싱
    target_chains = [chain.strip().upper() for chain in args.target_chains.split(',')]
    if len(target_chains) != 2:
        logger.error("정확히 2개의 체인이 필요합니다 (이동할 체인, 고정된 체인)")
        return 1
    
    logger.info(f"이동할 체인: {target_chains[0]}, 고정된 체인: {target_chains[1]}")
    
    # 랜덤 시드 설정
    if args.seed is not None:
        np.random.seed(args.seed)
        random.seed(args.seed)
        logger.info(f"랜덤 시드 설정: {args.seed}")
    
    # 설정 생성
    config = DisplacementConfig()
    config.actual_displacement_distance = 60.0  # 60Å로 이동하여 여유를 둠
    config.target_distance_for_validation = 50.0  # 유효성 검사는 50Å 기준
    config.clash_threshold = args.clash_threshold
    config.max_attempts = args.max_attempts
    config.path_check_threshold = args.path_threshold
    
    logger.info(f"설정: 실제이동거리={config.actual_displacement_distance}Å, 유효성검사목표={config.target_distance_for_validation}Å, Clash임계값={config.clash_threshold}Å, "
               f"최대시도={config.max_attempts}회, 경로차단임계값={config.path_check_threshold}Å")
    
    try:
        # 다중 변형 생성기
        displacer = MultipleDisplacer(config, logger)
        
        # 변형 구조 생성
        results = displacer.generate_multiple_variants(
            args.input, target_chains, args.num_variants, args.output_dir
        )
        
        # 결과 요약
        print(f"\n=== 변형 구조 생성 결과 ===")
        print(f"목표: {args.num_variants}개")
        print(f"성공: {len(results)}개")
        print(f"출력 디렉토리: {args.output_dir}")
        
        if results:
            print(f"\n생성된 구조들:")
            for i, (output_file, info) in enumerate(results, 1):
                final_result = info['final_result']
                print(f"  {i}. {os.path.basename(output_file)}")
                print(f"     최종 거리: {final_result['final_distance']:.2f}Å")
                print(f"     시도 횟수: {final_result['total_attempts']}회")
                print(f"     거리 오차: {final_result['distance_error']:.2f}Å")
        
        return 0
        
    except Exception as e:
        logger.error(f"실행 중 오류: {e}")
        return 1


if __name__ == "__main__":
    exit(main())
