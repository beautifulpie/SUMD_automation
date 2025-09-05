#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import json
import glob
import shutil
import subprocess, math
import numpy as np
import multiprocessing as mp
from multiprocessing import Process, Queue, Manager
import time
from datetime import datetime
from Bio.PDB import Structure, Model, Chain as PDBChain, Residue as PDBResidue, Atom
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.PDBIO import PDBIO,Select

# ===== 설정 불러오기 =====
try:
    from simple_config import *
except ImportError:
    # 기본 설정
    MAX_ITERATIONS = 10
    MAX_ATTEMPTS = 100
    SIMULATION_TIME_NS = 0.3
    SLOPE_THRESHOLD = -0.001
    GPU_ID = "1"
    MPI_RANKS = "4"
    NTOMP = "2"
    ENABLE_CHAIN_RESTORATION = True
    FORCE_FIELD = "charmm36-jul2022"
    WATER_MODEL = "tip3p"
    CLOSE_DISTANCE_THRESHOLD = 10.0
    LONG_MD_TIME_NS = 10.0
    ENABLE_LONG_MD = True
    TIMEOUT_LONG_MD = 7200
    BOX_DISTANCE = 1.5
    MAX_WARNINGS = 1
    ENABLE_MULTI_DIRECTION_SEPARATION = True
    SEPARATION_DISTANCE = 30.0
    ENABLE_ROTATIONAL_VARIANTS = True
    ROTATION_STEP = 60
    MAX_ROTATION_VARIANTS = 3
    print("설정 파일을 찾을 수 없어 기본 설정을 사용합니다.")

def log(message):
    """간단한 로깅"""
    timestamp = datetime.now().strftime("%H:%M:%S")
    print(f"[{timestamp}] {message}")

def parse_gpu_ids(gpu_id_string):
    """GPU ID 문자열을 개별 GPU 리스트로 파싱"""
    try:
        # "0123" -> ["0", "1", "2", "3"]
        gpu_list = [gpu_id_string[i] for i in range(len(gpu_id_string))]
        log(f"사용 가능한 GPU: {gpu_list} (총 {len(gpu_list)}개)")
        return gpu_list
    except Exception as e:
        log(f"GPU ID 파싱 실패: {e}")
        return ["0"]  # 기본값

# def save_original_chain_info(input_pdb, chain1, chain2):
#     """원본 체인 정보 저장 - 딕셔너리 반환"""
#     try:
#         original_chains = {chain1: chain1, chain2: chain2}
#         log(f"원본 체인 정보 저장: {original_chains}")
#         return original_chains
#     except Exception as e:
#         log(f"원본 체인 정보 저장 실패: {e}")
#         return None

def extract_target_chains_pdb(input_pdb, output_pdb, chain1, chain2):
    """타겟 체인만 추출하여 새로운 PDB 생성"""
    try:
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("complex", input_pdb)
        
        if structure is None:
            log(f"PDB 구조가 None: {input_pdb}")
            return False
        
        # 새로운 구조 생성
        new_structure = Structure.Structure("target")
        new_model = Model.Model(0)
        new_structure.add(new_model)
        
        # 타겟 체인들만 복사
        for model in structure:
            for chain in model:
                if chain.id in [chain1, chain2]:
                    new_chain = PDBChain.Chain(chain.id)
                    for residue in chain:
                        new_residue = PDBResidue.Residue(residue.id, residue.resname, residue.segid)
                        for atom in residue:
                            new_atom = Atom.Atom(atom.name, atom.coord, atom.bfactor, 
                                           atom.occupancy, atom.altloc, atom.fullname, 
                                           atom.serial_number, atom.element)
                            new_residue.add(new_atom)
                        new_chain.add(new_residue)
                    new_model.add(new_chain)
        
        # 저장
        io = PDBIO()
        io.set_structure(new_structure)
        io.save(output_pdb)
        
        log(f"타겟 체인 추출 완료: {output_pdb} (체인: {chain1}, {chain2})")
        return True
        
    except Exception as e:
        log(f"타겟 체인 추출 실패: {e}")
        return False

def calculate_separation_axis(pdb_file, chain1, chain2):
    """두 체인 간의 분리 축 계산"""
    try:
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("structure", pdb_file)
        
        chain_centers = {}
        
        for model in structure:
            for chain in model:
                if chain.id in [chain1, chain2]:
                    atoms = []
                    elements = []
                    for residue in chain:
                        for atom in residue:
                            atoms.append(atom.coord)
                            elements.append(atom.element)
                    
                    if atoms:
                        center = calculate_mass_weighted_center(atoms, elements)
                        if center is not None:
                            chain_centers[chain.id] = center
        
        if len(chain_centers) == 2:
            centers = list(chain_centers.values())
            # 체인 간 벡터 (chain2에서 chain1 방향)
            separation_vector = centers[1] - centers[0]
            # 정규화
            separation_vector = separation_vector / np.linalg.norm(separation_vector)
            
            log(f"분리 축 계산 완료: {separation_vector}")
            return separation_vector, chain_centers
        
        log("체인 중심 계산 실패")
        return None, None
        
    except Exception as e:
        log(f"분리 축 계산 오류: {e}")
        return None, None

def generate_multi_direction_vectors(base_vector):
    """원뿔형 벡터 생성 - ROTATION_STEP 설정에 따라 분할 각도 조절"""
    vectors = []
    
    # 1. 기준 벡터 (역방향 30Å)
    vectors.append(-base_vector)  # 반대 방향
    log("기준 역방향 벡터 추가")
    
    # 2. 기준 벡터에 수직인 두 벡터 찾기 (원뿔의 기저 평면용)
    if abs(base_vector[2]) < 0.9:  # Z 성분이 작으면
        perp1 = np.cross(base_vector, np.array([0, 0, 1]))
    else:  # Z 성분이 크면
        perp1 = np.cross(base_vector, np.array([1, 0, 0]))
    
    perp1 = perp1 / np.linalg.norm(perp1)
    perp2 = np.cross(base_vector, perp1)
    perp2 = perp2 / np.linalg.norm(perp2)
    
    # 3. 원뿔 표면에 점들 생성 (CONE_ANGLES 원뿔 밑넓이 결정하는 각도)
    for cone_angle in CONE_ANGLES:
        cone_angle_deg = cone_angle  # 원뿔 밑넓이 결정하는 각도
        cone_angle_rad = math.radians(cone_angle_deg)
        
        # 360도를 ROTATION_STEP으로 나눠서 균등 분할
        
        if  ENABLE_RANDOM_DIRECTION:
            num_directions=MAX_SEPARATION_VARIANTS
            directions=np.random.choice(range(360),size=num_directions, replace=False) 
        else:
            num_directions=360 // ROTATION_STEP

        for i in range(num_directions):

            plane_angle_deg = directions[i] if ENABLE_RANDOM_DIRECTION else i * ROTATION_STEP
            plane_angle_rad = math.radians(plane_angle_deg)
            
            # 원뿔 기저 평면에서의 방향 벡터
            direction_in_plane = perp1 * math.cos(plane_angle_rad) + perp2 * math.sin(plane_angle_rad)
            direction_in_plane = direction_in_plane / np.linalg.norm(direction_in_plane)
            
            # 원뿔 표면의 점 계산: 기준벡터에서 60도 기울어진 방향
            # 원뿔 축(기준벡터) 성분과 반지름 성분 조합
            cone_vector = (-base_vector * math.cos(cone_angle_rad) + 
                        direction_in_plane * math.sin(cone_angle_rad))
            cone_vector = cone_vector / np.linalg.norm(cone_vector)
            
            vectors.append(cone_vector)
        
        log(f"원뿔형 배치 완료: 총 {len(vectors)}개 방향 벡터 생성")
        log(f"  - 기준 역방향: 1개")
        log(f"  - 원뿔 표면 ({cone_angle_deg}도 반각): {len(vectors)-1}개")
        log(f"  - 분할 각도: {ROTATION_STEP}도 간격 ({num_directions}개 방향)")
    
    return vectors

def apply_structure_separation(input_pdb, output_pdb, separation_vector, distance, chain_to_move):
    """구조에 이격 적용"""
    try:
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("structure", input_pdb)
        
        # 이동할 거리 벡터 계산 (Angstrom -> Angstrom)
        move_vector = separation_vector * distance
        
        log(f"체인 {chain_to_move}를 {distance:.1f}Å 이동: {move_vector}")
        
        for model in structure:
            for chain in model:
                if chain.id == chain_to_move:
                    for residue in chain:
                        for atom in residue:
                            atom.coord = atom.coord + move_vector
        
        # 수정된 구조 저장
        io = PDBIO()
        io.set_structure(structure)
        io.save(output_pdb)
        
        log(f"이격된 구조 저장: {output_pdb}")
        return True
        
    except Exception as e:
        log(f"구조 이격 적용 실패: {e}")
        return False

def create_rotation_matrix(axis, angle_degrees):
    """축 기준 회전 행렬 생성 (Rodrigues' rotation formula)"""
    angle = math.radians(angle_degrees)
    axis = axis / np.linalg.norm(axis)  # 정규화
    
    cos_angle = math.cos(angle)
    sin_angle = math.sin(angle)
    
    # 외적 행렬
    cross_matrix = np.array([[0, -axis[2], axis[1]],
                            [axis[2], 0, -axis[0]], 
                            [-axis[1], axis[0], 0]])
    
    # Rodrigues' formula
    rotation_matrix = (np.eye(3) + sin_angle * cross_matrix + 
                      (1 - cos_angle) * np.dot(cross_matrix, cross_matrix))
    
    return rotation_matrix

def apply_rotational_transform(input_pdb, output_pdb, rotation_matrix, center_point, chain_to_move):
    """구조에 회전 변형 적용"""
    try:
        log(f"회전 변형 적용 시작: {os.path.basename(input_pdb)} -> {os.path.basename(output_pdb)}")
        
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("structure", input_pdb)
        
        if structure is None:
            log(f"구조 로드 실패: {input_pdb}")
            return False
        
        atom_count = 0
        for model in structure:
            for chain in model:
                if chain.id == chain_to_move:
                    for residue in chain:
                        for atom in residue:
                            # 중심점 기준으로 이동
                            relative_coord = atom.coord - center_point
                            # 회전 적용
                            rotated_coord = np.dot(rotation_matrix, relative_coord)
                            # 다시 원래 위치로 이동
                            atom.coord = rotated_coord + center_point
                            atom_count += 1
        
        log(f"회전 적용된 원자 수: {atom_count}")
        
        # 수정된 구조 저장
        io = PDBIO()
        io.set_structure(structure)
        io.save(output_pdb)
        
        # 파일 생성 확인
        if os.path.exists(output_pdb) and os.path.getsize(output_pdb) > 0:
            log(f"회전된 구조 저장 완료: {os.path.basename(output_pdb)} ({os.path.getsize(output_pdb)} bytes)")
            return True
        else:
            log(f"회전된 구조 저장 실패: 파일이 생성되지 않음")
            return False
        
    except Exception as e:
        log(f"회전 변형 적용 실패: {e}")
        import traceback
        log(f"상세 오류: {traceback.format_exc()}")
        return False

def generate_structure_variants(base_pdb, output_dir, base_name, chain_to_move):
    """기본 구조에서 회전 변형들 생성"""
    variants = [base_pdb]  # 원본 포함
    
    if not ENABLE_ROTATIONAL_VARIANTS:
        log(f"회전 변형 비활성화됨 (ENABLE_ROTATIONAL_VARIANTS = False)")
        return variants
    
    log(f"회전 변형 생성 시작: {base_name} (최대 {MAX_ROTATION_VARIANTS}개)")
    
    try:
        # 구조의 중심점 계산
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("structure", base_pdb)
        
        all_coords = []
        for model in structure:
            for chain in model:
                if chain.id==chain_to_move:
                    for residue in chain:
                        for atom in residue:
                            all_coords.append(atom.coord)
        
        if not all_coords:
            log(f"회전 변형 실패: 원자 좌표를 찾을 수 없음")
            return variants
        
        center_point = np.mean(all_coords, axis=0)
        log(f"구조 중심점: {center_point}")
        
        # 회전 축들 (X, Y, Z)
        rotation_axes = [
            np.array([1, 0, 0]),  # X축
            np.array([0, 1, 0]),  # Y축  
            np.array([0, 0, 1])   # Z축
        ]
        
        variant_count = 1
        for axis in rotation_axes:
            if variant_count >= MAX_ROTATION_VARIANTS:
                log(f"최대 회전 변형 수 도달: {MAX_ROTATION_VARIANTS}")
                break
                
            angle = STRUCTURE_ROTATION
            log(f"회전 변형 {variant_count}: 축{axis}, 각도{angle}도")
            
            rotation_matrix = create_rotation_matrix(axis, angle)
            
            variant_pdb = os.path.join(output_dir, f"{base_name}_rot{variant_count}.pdb")
            
            if apply_rotational_transform(base_pdb, variant_pdb, rotation_matrix, center_point, chain_to_move):
                variants.append(variant_pdb)
                variant_count += 1
                log(f"✓ 회전 변형 생성 성공: {os.path.basename(variant_pdb)}")
            else:
                log(f"✗ 회전 변형 생성 실패: {os.path.basename(variant_pdb)}")
        
        log(f"회전 변형 완료: 원본 1개 + 회전 {variant_count-1}개 = 총 {len(variants)}개")
        return variants
        
    except Exception as e:
        log(f"회전 변형 생성 중 오류: {e}")
        import traceback
        log(f"상세 오류: {traceback.format_exc()}")
        return variants

def create_initial_structure_pool(input_pdb, chain1, chain2, output_dir):
    """초기 구조 풀 생성 - 구조 리스트 반환"""
    structure_pool = []
    
    if not ENABLE_MULTI_DIRECTION_SEPARATION:
        # 기본 모드: 원본 구조만 사용
        structure_pool = [input_pdb]
        log("기본 모드: 원본 구조만 사용")
        return structure_pool
    
    try:
        # 구조 풀 디렉토리 생성
        pool_dir = os.path.join(output_dir, "structure_pool")
        if os.path.exists(pool_dir):
            shutil.rmtree(pool_dir)
        os.makedirs(pool_dir)
        
        # 1단계: 분리 축 계산
        separation_vector, chain_centers = calculate_separation_axis(input_pdb, chain1, chain2)
        if separation_vector is None:
            log("분리 축 계산 실패 - 원본 구조 사용")
            structure_pool = [input_pdb]
            return structure_pool
        
        # 2단계: 다방향 벡터 생성 (원뿔형)
        direction_vectors = generate_multi_direction_vectors(separation_vector)
        
        # 3단계: 각 방향으로 이격된 구조 생성
        base_structures = []
        
        # 더 작은 체인을 이동시키기 위해 체인 크기 비교
        chain_sizes = {}
        try:
            parser = PDBParser(QUIET=True)
            structure = parser.get_structure("structure", input_pdb)
            for model in structure:
                for chain in model:
                    if chain.id in [chain1, chain2]:
                        atom_count = sum(1 for residue in chain for atom in residue)
                        chain_sizes[chain.id] = atom_count
        except:
            chain_sizes = {chain1: 1000, chain2: 1000}  # 기본값
        
        # 더 작은 체인 결정
        chain_to_move = chain1 if chain_sizes.get(chain1, 0) <= chain_sizes.get(chain2, 0) else chain2
        log(f"이동할 체인: {chain_to_move} (크기: {chain_sizes.get(chain_to_move, 0)} 원자)")
        
        for i, direction in enumerate(direction_vectors):
            separated_pdb = os.path.join(pool_dir, f"separated_{i}.pdb")
            
            if apply_structure_separation(input_pdb, separated_pdb, direction, 
                                        SEPARATION_DISTANCE, chain_to_move):
                base_structures.append(separated_pdb)
                log(f"이격 구조 {i} 생성: {separated_pdb}")
        
        # 4단계: 각 이격 구조에 대해 회전 변형 생성
        for i, base_struct in enumerate(base_structures):
            variants = generate_structure_variants(base_struct, pool_dir, f"sep{i}", chain_to_move)
            structure_pool.extend(variants)
        
        log(f"구조 풀 생성 완료: 이 {len(structure_pool)}개 구조")
        
        # 구조 풀 정보 저장
        pool_info = {
            "total_structures": len(structure_pool),
            "separation_directions": len(direction_vectors),
            "rotation_variants_per_structure": len(generate_structure_variants(base_structures[0], pool_dir, "test", chain_to_move)) if base_structures else 0,
            "structures": [os.path.basename(s) for s in structure_pool]
        }
        
        with open(os.path.join(pool_dir, "pool_info.json"), "w") as f:
            json.dump(pool_info, f, indent=2)
        
        return structure_pool
        
    except Exception as e:
        log(f"구조 풀 생성 실패: {e}")
        structure_pool = [input_pdb]
        return structure_pool

def identify_ligand_receptor(input_pdb):
    """복합체에서 ligand와 receptor 판별 (크기 기준)"""
    try:
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("complex", input_pdb)

        if structure is None:
            log(f"PDB 구조가 None: {input_pdb}")
            raise

        chain_sizes = {}
        for model in structure:
            for chain in model:
                atom_count = len(list(chain.get_atoms()))
                chain_sizes[chain.id] = atom_count
        
        # 크기순 정렬
        sorted_chains = sorted(chain_sizes.items(), key=lambda x: x[1])
        ligand_chain = sorted_chains[0][0]  # 가장 작은 chain
        receptor_chain = sorted_chains[-1][0]  # 가장 큰 chain
        
        log(f"Ligand (작은 분자): Chain {ligand_chain} ({chain_sizes[ligand_chain]} atoms)")
        log(f"Receptor (큰 분자): Chain {receptor_chain} ({chain_sizes[receptor_chain]} atoms)")
        
        return ligand_chain, receptor_chain
    except Exception as e:
        log(f"체인 판별 실패: {e}")
        return None, None

def define_binding_site(input_pdb, ligand_chain, receptor_chain, distance_cutoff=4.0):
    """Binding site 정의 (<4Å 거리 기준)"""
    try:
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("complex", input_pdb)

        if structure is None:
            log(f"PDB 구조가 None: {input_pdb}")
            raise

        ligand_atoms = []
        receptor_residues = set()
        
        for model in structure:
            for chain in model:
                if chain.id == ligand_chain:
                    ligand_atoms = [atom for atom in chain.get_atoms()]
        
        for model in structure:
            for chain in model:
                if chain.id == receptor_chain:
                    for residue in chain:
                        for atom in residue:
                            # ligand의 모든 원자와 거리 계산
                            for lig_atom in ligand_atoms:
                                distance = np.linalg.norm(atom.coord - lig_atom.coord)
                                # log(f"{lig_atom.coord}, {atom.coord}, {distance}")
                                if distance < distance_cutoff:
                                    residue_id = residue.get_id()
                                    full_id = f"{residue_id[1]}{residue_id[2].strip()}"
                                    receptor_residues.add(full_id)
                                    break
        
        binding_site_residues = list(receptor_residues)
        log(f"Binding site residues: {binding_site_residues}")
        
        return binding_site_residues
    except Exception as e:
        log(f"Binding site 정의 실패: {e}")
        return None

def calculate_mass_weighted_center(atoms_coords, atom_types):
    """질량 가중 중심 계산"""
    # 원자 타입별 대략적인 질량 (원자량)
    mass_dict = {'C': 12.01, 'N': 14.01, 'O': 16.00, 'S': 32.06, 'H': 1.008}
    
    total_mass = 0
    weighted_coords = np.zeros(3)
    try:
        for coord, atom_type in zip(atoms_coords, atom_types):
            mass = mass_dict.get(atom_type, 12.01)  # 첫 글자로 원소 판단
            weighted_coords += coord * mass
            total_mass += mass
        weighted_coords/=total_mass

        return weighted_coords
    except Exception as e:
        log(f"질량 가중 중심 계산 실패 {e}")
        return weighted_coords

def calculate_distance_binding_site(pdb_file, binding_site_residues):
    """체인 간 거리 계산"""
    try:
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("s", pdb_file)
        if structure is None:
            log(f"PDB 구조가 None: {pdb_file}")
            raise

        for model in structure:
            chains = list(model)
            if len(list(chains[0])) < len(list(chains[1])):
                ligand_chain = chains[0]
                receptor_chain = chains[1]
            else:
                ligand_chain = chains[1]
                receptor_chain = chains[0]
            
            # Ligand center mass
            ligand_atoms = [atom.coord for atom in ligand_chain.get_atoms()]
            ligand_elements = [atom.element for atom in ligand_chain.get_atoms()]
            ligand_center = calculate_mass_weighted_center(ligand_atoms, ligand_elements)
            
            # Binding site center mass
            binding_site_atoms = []
            binding_site_elements=[]
            for residue in receptor_chain:
                residue_id = residue.get_id()
                full_id = f"{residue_id[1]}{residue_id[2].strip()}"
                if full_id in binding_site_residues:
                    binding_site_atoms.extend([atom.coord for atom in residue.get_atoms()])
                    binding_site_elements.extend([atom.element for atom in residue.get_atoms()])
            
            if len(binding_site_atoms)>0:
                binding_site_center = calculate_mass_weighted_center(binding_site_atoms, binding_site_elements)
                return float(np.linalg.norm(ligand_center - binding_site_center))

        return float('inf')
    except:
        return float('inf')

# ===== GROMACS 관련 함수들 =====

def run_command_with_output_check(cmd, cwd=None, input_text=None, expected_output=None, timeout=3600):
    """명령어 실행 및 목적 파일 생성 확인"""
    try:
        log(f"실행: {cmd}")
        result = subprocess.run(
            cmd, shell=True, cwd=cwd, 
            input=input_text.encode() if input_text else None,
            capture_output=True, text=True, timeout=timeout
        )
        
        if result.returncode != 0:
            log(f"Return code 오류 ({result.returncode}): {result.stderr}")
            return False
        
        if expected_output:
            if isinstance(expected_output, str):
                expected_output = [expected_output]
            
            for output_file in expected_output:
                full_path = os.path.join(cwd, output_file) if cwd else output_file
                if not os.path.exists(full_path):
                    log(f"목적 파일 생성 실패: {output_file}")
                    return False
                elif os.path.getsize(full_path) == 0:
                    log(f"빈 파일 생성: {output_file}")
                    return False
        
        return True
    except subprocess.TimeoutExpired:
        log(f"명령어 실행 시간 초과: {cmd}")
        return False
    except Exception as e:
        log(f"명령어 실행 실패: {e}")
        return False
    
def run_mdrun_with_checkpoint_recovery(cmd, cwd="", input_text="", expected_output=[], timeout=3600, max_retries=2, is_long_md=False):
    """GROMACS mdrun을 checkpoint 복구 기능과 함께 실행"""
    attempt = 0
    
    while attempt <= max_retries:
        try:
            if attempt == 0:
                # 첫 번째 시도 - 일반 실행
                log(f"MD 실행 시도 {attempt + 1}/{max_retries + 1}")
                result = subprocess.run(
                    cmd, shell=True, cwd=cwd, 
                    input=input_text.encode() if input_text else None,
                    capture_output=True, text=True, timeout=timeout
                )
            else:
                # 재시작 시도 - checkpoint 파일 확인
                checkpoint_files = []
                for output_file in expected_output:
                    base_name = output_file.replace('.gro', '').replace('.xtc', '')
                    checkpoint_file = os.path.join(cwd, f"{base_name}.cpt")
                    if os.path.exists(checkpoint_file):
                        checkpoint_files.append(checkpoint_file)
                
                if not checkpoint_files:
                    log(f"재시작 시도 {attempt}: checkpoint 파일을 찾을 수 없음")
                    break
                
                # checkpoint에서 재시작하는 명령어 구성
                restart_cmd = cmd + " -cpi"
                log(f"MD 재시작 시도 {attempt + 1}/{max_retries + 1} (checkpoint에서)")
                log(f"재시작 명령어: {restart_cmd}")
                
                result = subprocess.run(
                    restart_cmd, shell=True, cwd=cwd,
                    capture_output=True, text=True, timeout=timeout
                )
            
            # 실행 결과 확인
            if result.returncode == 0:
                # 출력 파일 존재 확인
                all_files_exist = True
                for output_file in expected_output:
                    full_path = os.path.join(cwd, output_file)
                    if not os.path.exists(full_path) or os.path.getsize(full_path) == 0:
                        all_files_exist = False
                        break
                
                if all_files_exist:
                    log(f"MD 실행 성공 (시도 {attempt + 1})")
                    return True
                else:
                    log(f"MD 실행 후 출력 파일 확인 실패 (시도 {attempt + 1})")
            else:
                log(f"MD 실행 실패 (시도 {attempt + 1}): Return code {result.returncode}")
                log(f"Error: {result.stderr}")
            
        except subprocess.TimeoutExpired:
            log(f"MD 실행 시간 초과 (시도 {attempt + 1}/{max_retries + 1})")
            if is_long_md:
                log(f"Long MD 시간 초과 - 다음 시도에서 checkpoint 복구 시도")
        except Exception as e:
            log(f"MD 실행 중 예외 발생 (시도 {attempt + 1}): {e}")
        
        attempt += 1
        
        if attempt <= max_retries:
            log(f"다음 시도까지 5초 대기...")
            time.sleep(5)
    
    log(f"MD 실행 최종 실패: {max_retries + 1}회 시도 모두 실패")
    return False

def create_mdp_files(work_dir, long_md=False):
    """MDP 파일들 생성 - SD integrator 및 랜덤 시드 적용"""
    
    try:
        em_settings = MDP_SETTINGS["em"]
        nvt_settings = MDP_SETTINGS["nvt"] 
        npt_settings = MDP_SETTINGS["npt"]
        npt2_settings = MDP_SETTINGS["npt2"]
        md_settings = MDP_SETTINGS["md"]
        output_freq = OUTPUT_FREQUENCY
        
        if long_md and "long_md" in MDP_SETTINGS:
            md_settings = MDP_SETTINGS["long_md"]
    except (NameError, KeyError):
        em_settings = {"integrator": "steep", "nsteps": 50000, "emtol": 1000.0, "emstep": 0.01}
        nvt_settings = {"integrator": "sd", "dt": 0.002, "nsteps": 25000, "temperature": 300}
        npt_settings = {"integrator": "sd", "dt": 0.002, "nsteps": 25000, "temperature": 300, "pressure": 1.0}
        npt2_settings = {"integrator": "sd", "dt": 0.002, "nsteps": 25000, "temperature": 300, "pressure": 1.0}
        md_settings = {"integrator": "sd", "dt": 0.002, "temperature": 300, "pressure": 1.0}
        output_freq = {"energy": 5000, "log": 5000, "trajectory": 5000}
    
    # EM MDP
    em_mdp = f"""integrator = {em_settings["integrator"]}
nsteps = {em_settings["nsteps"]}
emtol = {em_settings["emtol"]}
emstep = {em_settings["emstep"]}
nstlist = 1
cutoff-scheme = Verlet
ns_type = grid
coulombtype = PME
rcoulomb = 1.0
rvdw = 1.0
pbc = xyz
"""
    
    # NVT MDP - SD integrator with random seed
    nvt_mdp = f"""integrator = sd
dt = {nvt_settings["dt"]}
nsteps = {nvt_settings["nsteps"]}
nstenergy = {output_freq["energy"]//10}
nstlog = {output_freq["log"]//10}
nstxout-compressed = {output_freq["trajectory"]//10}
constraints = h-bonds
constraint_algorithm = lincs
cutoff-scheme = Verlet
ns_type = grid
nstlist = 10
rcoulomb = 1.0
rvdw = 1.0
DispCorr = EnerPres
coulombtype = PME
tc-grps = System
tau_t = 0.1
ref_t = {nvt_settings["temperature"]}
bd-fric = 0
ld-seed = -1
pcoupl = no
pbc = xyz
gen_vel = yes
gen_temp = {nvt_settings["temperature"]}
gen_seed = -1
"""
    
    # NPT MDP - SD integrator with random seed
    npt_mdp = f"""define = -DPOSRES
integrator = sd
dt = {npt_settings["dt"]}
nsteps = {npt_settings["nsteps"]}
nstenergy = {output_freq["energy"]//10}
nstlog = {output_freq["log"]//10}
nstxout-compressed = {output_freq["trajectory"]//10}
continuation = yes
constraints = h-bonds
constraint_algorithm = lincs
cutoff-scheme = Verlet
ns_type = grid
nstlist = 10
rcoulomb = 1.0
rvdw = 1.0
DispCorr = EnerPres
coulombtype = PME
tc-grps = System
tau_t = 0.1
ref_t = {npt_settings["temperature"]}
bd-fric = 0
ld-seed = -1
pcoupl = C-rescale
pcoupltype = isotropic
tau_p = 2.0
ref_p = {npt_settings["pressure"]}
compressibility = 4.5e-5
pbc = xyz
gen_vel = no
"""

# NPT2 MDP - SD integrator with random seed
    npt2_mdp = f"""integrator = sd
dt = {npt2_settings["dt"]}
nsteps = {npt2_settings["nsteps"]}
nstenergy = {output_freq["energy"]//10}
nstlog = {output_freq["log"]//10}
nstxout-compressed = {output_freq["trajectory"]//10}
continuation = yes
constraints = h-bonds
constraint_algorithm = lincs
cutoff-scheme = Verlet
ns_type = grid
nstlist = 10
rcoulomb = 1.0
rvdw = 1.0
DispCorr = EnerPres
coulombtype = PME
tc-grps = System
tau_t = 0.1
ref_t = {npt_settings["temperature"]}
bd-fric = 0
ld-seed = -1
pcoupl = Parrinello-Rahman
pcoupltype = isotropic
tau_p = 2.0
ref_p = {npt_settings["pressure"]}
compressibility = 4.5e-5
pbc = xyz
gen_vel = no
"""
    
    # MD MDP - SD integrator with random seed
    simulation_time = LONG_MD_TIME_NS if long_md else SIMULATION_TIME_NS
    nsteps = int(simulation_time * 1000 / md_settings["dt"])
    
    md_mdp = f"""integrator = sd
dt = {md_settings["dt"]}
nsteps = {nsteps}
nstenergy = {output_freq["energy"]}
nstlog = {output_freq["log"]}
nstxout-compressed = {output_freq["trajectory"]}
tc-grps = System
tau_t = 0.1
ref_t = {md_settings["temperature"]}
bd-fric = 0
ld-seed = -1
pcoupl = C-rescale
pcoupltype = isotropic
tau_p = 2.0
ref_p = {md_settings["pressure"]}
compressibility = 4.5e-5
constraints = h-bonds
constraint_algorithm = LINCS
cutoff-scheme = Verlet
nstlist = 40
ns_type = grid
coulombtype = PME
rcoulomb = 1.0
rvdw = 1.0
DispCorr = EnerPres
pbc = xyz
"""
    
    # 파일들 저장
    for name, content in [("em.mdp", em_mdp), ("nvt.mdp", nvt_mdp), 
                         ("npt.mdp", npt_mdp), ("npt2.mdp", npt2_mdp), ("md.mdp", md_mdp)]:
        with open(os.path.join(work_dir, name), "w") as f:
            f.write(content)

def get_chains_by_size(structure):
    """구조의 chain들을 크기 순으로 정렬하여 반환"""
    chain_info = []
    for model in structure:
        for chain in model:
            atom_count = len(list(chain.get_atoms()))
            chain_info.append((chain.id, atom_count, chain))
    
    # 크기 순 정렬 (작은 것부터)
    chain_info.sort(key=lambda x: x[1])
    return chain_info[:1]

def calculate_rmsd_between_structures(pdb1, pdb2):
    """두 PDB 구조 간의 RMSD 계산 (단백질만)"""
    try:
        from Bio.PDB import PDBParser, Superimposer
        import numpy as np
        
        parser = PDBParser(QUIET=True)
        structure1 = parser.get_structure("struct1", pdb1)
        structure2 = parser.get_structure("struct2", pdb2)

        # Chain ID를 기준으로 매칭 (순서 무관)
        chains1_info = get_chains_by_size(structure1)  # [(H, 100, chainH), (L, 300, chainL)]
        chains2_info = get_chains_by_size(structure2)  # [(A, 100, chainA), (B, 300, chainB)]

        if len(chains1_info) != len(chains2_info):
            log(f"Chain 개수 불일치: {len(chains1_info)} vs {len(chains2_info)}")
            return float('inf')
        
        if len(chains1_info) < 2:
            log("비교할 chain이 충분하지 않음")
            return float('inf')
        
        log(f"Reference chains: {[(info[0], info[1]) for info in chains1_info]}")
        log(f"Target chains: {[(info[0], info[1]) for info in chains2_info]}")
        
        # 크기 순서대로 매칭 (ligand끼리, receptor끼리)
        atoms1 = []
        atoms2 = []
        
        for (id1, size1, chain1), (id2, size2, chain2) in zip(chains1_info, chains2_info):
            log(f"매칭: {id1}({size1} atoms) <-> {id2}({size2} atoms)")
            
            # 각 chain의 residue들을 정렬
            residues1 = sorted(chain1.get_residues(), key=lambda r: r.id)
            residues2 = sorted(chain2.get_residues(), key=lambda r: r.id)
            
            # Residue ID 기준 매칭 (번호 기준)
            res_dict1 = {r.id[1]: r for r in residues1}  # residue number만 사용
            res_dict2 = {r.id[1]: r for r in residues2}
            
            common_res_nums = set(res_dict1.keys()) & set(res_dict2.keys())
            
            for res_num in sorted(common_res_nums):
                residue1 = res_dict1[res_num]
                residue2 = res_dict2[res_num]
                
                # CA 원자가 둘 다 있는 경우만
                if 'CA' in residue1 and 'CA' in residue2:
                    atoms1.append(residue1['CA'])
                    atoms2.append(residue2['CA'])
        
        if len(atoms1) != len(atoms2) or len(atoms1) < 3:
            log(f"RMSD 계산 실패: 대응되는 CA 원자 부족 ({len(atoms1)} vs {len(atoms2)})")
            return float('inf')
        
        # Superimposer를 사용하여 RMSD 계산
        super_imposer = Superimposer()
        super_imposer.set_atoms(atoms1, atoms2)
        rmsd = super_imposer.rms
        
        log(f"크기 기준 RMSD 계산 성공: {rmsd:.3f}Å (CA 원자 {len(atoms1)}개)")
        return float(rmsd)
        
    except Exception as e:
        log(f"크기 기준 RMSD 계산 중 오류: {e}")
        import traceback
        log(f"상세 오류: {traceback.format_exc()}")
        return float('inf')

def calculate_combined_score(distance, rmsd, distance_weight=0.6, rmsd_weight=0.4):
    """거리와 RMSD를 조합한 점수 계산 (낮을수록 좋음)"""
    try:
        # 정규화를 위한 기준값들
        max_reasonable_distance = 50.0  # 50Å
        max_reasonable_rmsd = 10.0      # 10Å
        
        # 0-1로 정규화
        norm_distance = min(distance / max_reasonable_distance, 1.0)
        norm_rmsd = min(rmsd / max_reasonable_rmsd, 1.0)
        
        # 가중 평균 계산
        combined_score = 1 - (distance_weight * norm_distance + rmsd_weight * norm_rmsd)
        
        return combined_score
        
    except Exception as e:
        log(f"점수 계산 중 오류: {e}")
        return float('inf')


def extract_distances_and_top_structures_with_rmsd(tpr_file, xtc_file, binding_site_residues, 
                                                  work_dir, reference_pdb, save_top_n=5):
    """궤적에서 거리와 RMSD를 고려한 Top N 구조 저장"""
    try:
        # 궤적을 개별 PDB 파일로 변환
        cmd = f"echo 'Protein' | gmx trjconv -s {tpr_file} -f {xtc_file} -o trajectory.pdb -sep"
        if not run_command_with_output_check(cmd, work_dir, expected_output="trajectory0.pdb"):
            log("궤적 변환 실패")
            return [], []
        
        structures_with_scores = []  # (frame_num, distance, rmsd, combined_score, pdb_file) 튜플 리스트
        frame_num = 0
        
        log("프레임별 거리 및 RMSD 계산 시작...")
        
        while True:
            frame_pdb = os.path.join(work_dir, f"trajectory{frame_num}.pdb")
            if not os.path.exists(frame_pdb):
                break
            
            if os.path.getsize(frame_pdb) == 0:
                log(f"빈 프레임 파일: trajectory{frame_num}.pdb")
                break
            
            # 거리 계산
            distance = calculate_distance_binding_site(frame_pdb, binding_site_residues)
            if distance == float('inf'):
                frame_num += 1
                continue
            
            # RMSD 계산 (reference 구조 대비)
            rmsd = calculate_rmsd_between_structures(reference_pdb, frame_pdb)
            if rmsd == float('inf'):
                log(f"Frame {frame_num}: RMSD 계산 실패, 거리만 사용")
                rmsd = 0.0  # RMSD 계산 실패시 거리만 고려
            
            # 조합 점수 계산
            combined_score = calculate_combined_score(distance, rmsd)
            
            structures_with_scores.append((frame_num, distance, rmsd, combined_score, frame_pdb))
            
            if frame_num % 50 == 0:  # 진행상황 로그
                log(f"Frame {frame_num}: 거리={distance:.2f}Å, RMSD={rmsd:.2f}Å, 점수={combined_score:.4f}")
            
            frame_num += 1
        
        log(f"총 {len(structures_with_scores)}개 프레임에서 거리 및 RMSD 추출")
        
        # Top N 구조 찾기 (조합 점수 기준으로 정렬)
        top_structures = []
        if structures_with_scores:
            # 조합 점수 순으로 정렬 (낮을수록 좋음)
            sorted_structures = sorted(structures_with_scores, key=lambda x: x[3])
            top_n = min(save_top_n, len(sorted_structures))
            
            log(f"Top {top_n} 구조 저장 중...")
            
            for i in range(top_n):
                frame_num, distance, rmsd, combined_score, frame_pdb = sorted_structures[i]
                
                # Top 구조 저장
                top_structure_name = f"top_{i+1}_frame_{frame_num}_dist_{distance:.2f}A_rmsd_{rmsd:.2f}A_score_{combined_score:.4f}.pdb"
                top_structure_path = os.path.join(work_dir, top_structure_name)
                
                try:
                    shutil.copy(frame_pdb, top_structure_path)
                    
                    structure_info = {
                        'rank': i + 1,
                        'frame': frame_num,
                        'distance': distance,
                        'rmsd': rmsd,
                        'combined_score': combined_score,
                        'filename': top_structure_name,
                        'path': top_structure_path,
                        'distance_rank': None,  # 나중에 계산
                        'rmsd_rank': None       # 나중에 계산
                    }
                    top_structures.append(structure_info)
                    
                    log(f"Top {i+1}: Frame {frame_num}, 거리 {distance:.2f}Å, RMSD {rmsd:.2f}Å, 점수 {combined_score:.4f}")
                    
                except Exception as e:
                    log(f"Top 구조 저장 실패 (rank {i+1}): {e}")
            
            # 개별 랭킹 정보 추가 (참고용)
            distance_sorted = sorted(structures_with_scores, key=lambda x: x[1])
            rmsd_sorted = sorted(structures_with_scores, key=lambda x: x[2])
            
            for struct_info in top_structures:
                frame_num = struct_info['frame']
                
                # 거리 랭킹 찾기
                for rank, (f_num, _, _, _, _) in enumerate(distance_sorted, 1):
                    if f_num == frame_num:
                        struct_info['distance_rank'] = rank
                        break
                
                # RMSD 랭킹 찾기
                for rank, (f_num, _, _, _, _) in enumerate(rmsd_sorted, 1):
                    if f_num == frame_num:
                        struct_info['rmsd_rank'] = rank
                        break
        
        # 임시 frame 파일들 정리
        frame_num = 0
        while True:
            frame_pdb = os.path.join(work_dir, f"trajectory{frame_num}.pdb")
            if not os.path.exists(frame_pdb):
                break
            try:
                os.remove(frame_pdb)
            except:
                pass
            frame_num += 1
        
        # 거리만 따로 추출 (기존 호환성을 위해)
        distances = [item[1] for item in structures_with_scores]
        
        # 결과 요약 로그
        if top_structures:
            log(f"Top 구조 선정 완료:")
            log(f"  최고 조합 점수: {top_structures[0]['combined_score']:.4f}")
            log(f"  최고 거리: {min([s['distance'] for s in top_structures]):.2f}Å")
            log(f"  최고 RMSD: {min([s['rmsd'] for s in top_structures]):.2f}Å")
        
        return distances, top_structures
        
    except Exception as e:
        log(f"궤적 분석 오류: {e}")
        return [], []
    
def extract_distances_from_trajectory(tpr_file, xtc_file, binding_site_residues, work_dir):
    """궤적에서 거리 추출"""
    try:
        cmd = f"echo 'Protein' | gmx trjconv -s {tpr_file} -f {xtc_file} -o trajectory.pdb -sep"
        if not run_command_with_output_check(cmd, work_dir, expected_output="trajectory0.pdb"):
            log("궤적 변환 실패")
            return []
        
        distances = []
        frame_num = 0
        
        while True:
            frame_pdb = os.path.join(work_dir, f"trajectory{frame_num}.pdb")
            if not os.path.exists(frame_pdb):
                break
            
            if os.path.getsize(frame_pdb) == 0:
                log(f"빈 프레임 파일: trajectory{frame_num}.pdb")
                break
            
            distance = calculate_distance_binding_site(frame_pdb, binding_site_residues)
            if distance != float('inf'):
                distances.append(distance)
            
            os.remove(frame_pdb)
            frame_num += 1
        
        log(f"이 {len(distances)}개 프레임에서 거리 추출")
        return distances
        
    except Exception as e:
        log(f"궤적 분석 오류: {e}")
        return []

def calculate_slope(distances):
    """기울기 계산"""
    if len(distances) < 2:
        return 0.0
    
    n = len(distances)
    x = list(range(n))
    x_mean = sum(x) / n
    y_mean = sum(distances) / n
    
    numerator = sum((x[i] - x_mean) * (distances[i] - y_mean) for i in range(n))
    denominator = sum((x[i] - x_mean) ** 2 for i in range(n))
    
    return numerator / denominator if denominator != 0 else 0.0

def run_gromacs_pipeline(work_dir, input_pdb, long_md=False):
    """GROMACS 파이프라인 실행"""
    stages = []
    
    # 1. pdb2gmx
    log("pdb2gmx 실행")
    cmd = f"echo '1\\n1' | gmx pdb2gmx -f {input_pdb} -o complex.gro -p topol.top \
          -water {WATER_MODEL} -ff {FORCE_FIELD} -ignh"
    success = run_command_with_output_check(cmd, work_dir, expected_output=["complex.gro", "topol.top"])
    stages.append({"stage": "pdb2gmx", "success": success})
    if not success:
        return stages
    
    # 2. editconf
    log("editconf 실행")
    cmd = f"gmx editconf -f complex.gro -o box.gro -c -d {BOX_DISTANCE} -bt cubic"
    success = run_command_with_output_check(cmd, work_dir, expected_output="box.gro")
    stages.append({"stage": "editconf", "success": success})
    if not success:
        return stages
    
    # 3. solvate
    log("solvate 실행")
    cmd = "gmx solvate -cp box.gro -cs spc216.gro -o solv.gro -p topol.top"
    success = run_command_with_output_check(cmd, work_dir, expected_output="solv.gro")
    stages.append({"stage": "solvate", "success": success})
    if not success:
        return stages
    
    # 4. ions grompp
    log("ions grompp 실행")
    with open(os.path.join(work_dir, "ions.mdp"), "w") as f:
        f.write("""integrator = steep
emtol = 1000.0
emstep = 0.01
nsteps = 50000
nstlist = 1
cutoff-scheme = Verlet
ns_type = grid
coulombtype = cutoff
rcoulomb = 1.0
rvdw = 1.0
pbc = xyz
""")
    cmd = f"gmx grompp -f ions.mdp -c solv.gro -p topol.top -o ions.tpr -maxwarn {MAX_WARNINGS}"
    success = run_command_with_output_check(cmd, work_dir, expected_output="ions.tpr")
    stages.append({"stage": "ions_grompp", "success": success})
    if not success:
        return stages
    
    # 5. genion
    log("genion 실행")
    cmd = "echo 'SOL' | gmx genion -s ions.tpr -o solv_ions.gro -p topol.top -pname NA -nname CL -neutral"
    success = run_command_with_output_check(cmd, work_dir, expected_output="solv_ions.gro")
    stages.append({"stage": "genion", "success": success})
    if not success:
        return stages
    
    # MDP 파일들 생성
    create_mdp_files(work_dir, long_md)
    
    # 6. EM
    log("EM 실행")
    cmd = f"gmx grompp -f em.mdp -c solv_ions.gro -p topol.top -o em.tpr -maxwarn {MAX_WARNINGS}"
    success = run_command_with_output_check(cmd, work_dir, expected_output="em.tpr")
    if success:
        cmd = f"mpirun --allow-run-as-root -np {MPI_RANKS} gmx_mpi mdrun -v -deffnm em -ntomp {NTOMP} \
              -nb gpu -gpu_id {GPU_ID}"
        success = run_command_with_output_check(cmd, work_dir, expected_output=["em.gro", "em.edr"])
    stages.append({"stage": "em", "success": success})
    if not success:
        return stages
    
    # 7. NVT
    log("NVT 실행")
    cmd = f"gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr -maxwarn {MAX_WARNINGS}"
    success = run_command_with_output_check(cmd, work_dir, expected_output="nvt.tpr")
    if success:
        cmd = f"mpirun --allow-run-as-root -np {MPI_RANKS} gmx_mpi mdrun -v -deffnm nvt -ntomp {NTOMP} \
              -nb gpu -gpu_id {GPU_ID} -npme 1 -pme gpu -bonded gpu"
        success = run_command_with_output_check(cmd, work_dir, expected_output=["nvt.gro", "nvt.cpt"])
    stages.append({"stage": "nvt", "success": success})
    if not success:
        return stages
    
    # 8. NPT
    log("NPT 실행")
    cmd = f"gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr -maxwarn {MAX_WARNINGS}"
    success = run_command_with_output_check(cmd, work_dir, expected_output="npt.tpr")
    if success:
        cmd = f"mpirun --allow-run-as-root -np {MPI_RANKS} gmx_mpi mdrun -v -deffnm npt -ntomp {NTOMP} \
              -nb gpu -gpu_id {GPU_ID} -npme 1 -pme gpu -bonded gpu"
        success = run_command_with_output_check(cmd, work_dir, expected_output=["npt.gro", "npt.cpt"])
    stages.append({"stage": "npt", "success": success})
    if not success:
        return stages
    
    # 9. MD
    md_label = "긴 MD" if long_md else "MD"
    log(f"{md_label} 실행")
    cmd = f"gmx grompp -f md.mdp -c npt.gro -p topol.top -o md.tpr -maxwarn {MAX_WARNINGS}"
    success = run_command_with_output_check(cmd, work_dir, expected_output="md.tpr")
    if success:
        timeout = TIMEOUT_LONG_MD if long_md else 3600
        max_retries = 2 if long_md else 0  # Long MD만 재시작 시도
        cmd = f"mpirun --allow-run-as-root -np {MPI_RANKS} gmx_mpi mdrun -v -deffnm md -ntomp {NTOMP} \
        -nb gpu -gpu_id {GPU_ID} -npme 1 -pme gpu -bonded gpu"
        if long_md:
            log(f"{md_label} - checkpoint 복구 기능 활성화 (최대 {max_retries}회 재시작)")
            success = run_mdrun_with_checkpoint_recovery(
                cmd, work_dir, expected_output=["md.gro", "md.xtc"], timeout=timeout, max_retries=max_retries, is_long_md=True
            )
        else:
            success = run_command_with_output_check(cmd, work_dir, expected_output=["md.gro", "md.xtc"], timeout=timeout)
    
    stages.append({"stage": md_label, "success": success})
    
    return stages

# ===== 체인 복원 함수들 =====

def get_original_chain_order(input_pdb):
    """원본 PDB 파일에서 체인 순서 추출"""
    try:
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("structure", input_pdb)
        
        chain_order = []
        for model in structure:
            for chain in model:
                if chain.id not in chain_order:
                    chain_order.append(chain.id)
        
        log(f"원본 체인 순서: {chain_order}")
        return chain_order
        
    except Exception as e:
        log(f"원본 체인 순서 추출 실패: {e}")
        return []

def extract_chain_id_from_filename(filename):
    """파일명에서 체인 ID 추출"""
    import re
    match = re.search(r'topol_Protein_chain_([A-Z])\.itp', os.path.basename(filename))
    return match.group(1) if match else None

def count_atoms_in_topology(topology_file):
    """토폴로지 파일에서 원자 수 카운트"""
    try:
        atom_count = 0
        in_atoms_section = False
        
        with open(topology_file, 'r') as f:
            for line in f:
                line = line.strip()
                
                # [ atoms ] 섹션 시작
                if line == "[ atoms ]":
                    in_atoms_section = True
                    continue
                
                # 다른 섹션 시작하면 atoms 섹션 종료
                if line.startswith("[") and in_atoms_section:
                    break
                
                # atoms 섹션 내에서 원자 라인 카운트
                if in_atoms_section and line and not line.startswith(";"):
                    parts = line.split()
                    if len(parts) >= 5:  # 최소한의 원자 정보가 있는 라인
                        atom_count += 1
        
        return atom_count
        
    except Exception as e:
        log(f"토폴로지 파일 읽기 실패 {topology_file}: {e}")
        return 0

def assign_chains_by_atom_ranges(input_pdb, output_pdb, chain_atom_counts, original_chain_order):
    """원자 범위 기반 체인 할당 - 원본 체인 순서 유지"""
    try:
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("structure", input_pdb)
        if structure is None:
            log(f"PDB 구조가 None: {input_pdb}")
            raise
        
        
        # 원본 순서에 있는 체인들만 필터링
        sorted_chains = [chain_id for chain_id in original_chain_order 
                        if chain_id in chain_atom_counts]
        # 혹시 빠진 체인이 있으면 추가
        missing_chains = [chain_id for chain_id in chain_atom_counts.keys() 
                        if chain_id not in sorted_chains]
        sorted_chains.extend(sorted(missing_chains))
    
        
        log(f"사용할 체인 순서: {sorted_chains}")
        
        # 각 체인의 시작/끝 원자 인덱스 계산
        chain_ranges = {}
        current_start = 1
        
        for chain_id in sorted_chains:
            atom_count = chain_atom_counts[chain_id]
            chain_ranges[chain_id] = (current_start, current_start + atom_count - 1)
            current_start += atom_count
        
        log(f"체인별 원자 범위: {chain_ranges}")
        
        # PDB 구조에서 원자 할당
        atom_index = 0
        
        for model in structure:
            for chain in model:
                chain_atoms = list(chain.get_atoms())
                
                for atom in chain_atoms:
                    atom_index += 1
                    
                    # 현재 원자가 속할 체인 찾기
                    for target_chain_id, (start_idx, end_idx) in chain_ranges.items():
                        if start_idx <= atom_index <= end_idx:
                            chain.id = target_chain_id
                            break
                    else:
                        # 범위를 벗어나면 마지막 체인에 할당
                        chain.id = sorted_chains[-1] if sorted_chains else 'A'
                        log(f"경고: 원자 {atom_index}가 범위를 벗어남, 체인 {chain.id}에 할당")
        
        # 수정된 구조 저장
        io = PDBIO()
        io.set_structure(structure)
        io.save(output_pdb)
        
        log(f"토폴로지 기반 체인 복원 완료: {sorted_chains}")
        return True
        
    except Exception as e:
        log(f"체인 할당 실패: {e}")
        return False

def restore_chain_info_from_topology(work_dir, converted_pdb, output_pdb, original_chain_order):
    """토폴로지 파일 기반 체인 정보 복원 (원본 순서 유지)"""
    try:
        log("=== 토폴로지 기반 체인 정보 복원 시작 ===")
        
        # 1. 토폴로지 파일들 스캔 (정렬하지 않음)
        topology_files = glob.glob(os.path.join(work_dir, "topol_Protein_chain_*.itp"))
        
        if not topology_files:
            log("체인별 토폴로지 파일을 찾을 수 없음 - 원본 파일 복사")
            shutil.copy(converted_pdb, output_pdb)
            return False
        
        log(f"발견된 토폴로지 파일들: {[os.path.basename(f) for f in topology_files]}")
        
        # 2. 각 체인의 원자 수 파악
        chain_atom_counts = {}
        total_atoms = 0
        
        for topo_file in topology_files:
            chain_id = extract_chain_id_from_filename(topo_file)
            atom_count = count_atoms_in_topology(topo_file)
            
            if chain_id and atom_count > 0:
                chain_atom_counts[chain_id] = atom_count
                total_atoms += atom_count
                log(f"체인 {chain_id}: {atom_count} 원자")
        
        if not chain_atom_counts:
            log("토폴로지에서 체인 정보 추출 실패")
            shutil.copy(converted_pdb, output_pdb)
            return False
        
        log(f"이 {len(chain_atom_counts)}개 체인, {total_atoms} 원자")
        
        # 3. PDB 파일에서 체인 할당 (원본 순서 유지)
        success = assign_chains_by_atom_ranges(converted_pdb, output_pdb, chain_atom_counts, original_chain_order)
        
        if success:
            log("토폴로지 기반 체인 정보 복원 완료")
        else:
            log("토폴로지 기반 체인 정보 복원 실패 - 원본 파일 복사")
            shutil.copy(converted_pdb, output_pdb)
        
        return success
        
    except Exception as e:
        log(f"토폴로지 기반 체인 복원 실패: {e}")
        try:
            shutil.copy(converted_pdb, output_pdb)
        except:
            pass
        return False

def restore_original_chain_ids(gromacs_pdb, output_pdb, work_dir):
    """토폴로지 기반 체인 복원 (원본 순서 유지)"""
    
    try:
        # 원본 체인 순서 추출
        original_chain_order = None
        original_pdb=os.path.join(work_dir, "input.pdb")
        original_chain_order = get_original_chain_order(original_pdb)
        
        # 토폴로지 파일 기반 복원 시도
        success = restore_chain_info_from_topology(work_dir, gromacs_pdb, output_pdb, original_chain_order)
        
        if success:
            log("토폴로지 기반 체인 복원 성공")
            return True
        else:
            raise
        
    except Exception as e:
        log(f"체인 복원 실패: {e}")
        import traceback
        log(f"상세 오류: {traceback.format_exc()}")
        try:
            shutil.copy(gromacs_pdb, output_pdb)
        except:
            pass
        return False

def create_next_iteration_structure(work_dir, output_pdb, long_md=False):
    """다음 iteration을 위한 구조 생성 (protein group 선택, 체인 정보 보존)"""
    try:
        tpr_file = os.path.join(work_dir, "md.tpr")
        gro_file = os.path.join(work_dir, "md.gro")
        
        # GROMACS로 protein만 추출
        temp_pdb = os.path.join(work_dir, "protein_only.pdb")
        cmd = f"echo 'Protein' | gmx trjconv -s {tpr_file} -f {gro_file} -o {temp_pdb}"
        
        if not run_command_with_output_check(cmd, work_dir, expected_output="protein_only.pdb"):
            log("protein 추출 실패")
            return False
        
        # 원본 chain ID 복원 (긴 MD이거나 설정이 활성화된 경우)
        if ENABLE_CHAIN_RESTORATION or long_md:
            success = restore_original_chain_ids(temp_pdb, output_pdb, work_dir)
            if success:
                restoration_reason = "긴 MD" if long_md else "설정 활성화"
                log(f"다음 iteration용 구조 생성 완료 (체인 복원됨 - {restoration_reason})")
            else:
                log("체인 복원 실패, GROMACS 결과 사용")
                shutil.copy(temp_pdb, output_pdb)
        else:
            shutil.copy(temp_pdb, output_pdb)
            log("다음 iteration용 구조 생성 완료 (체인 복원 비활성화)")
        
        # 임시 파일 정리
        if os.path.exists(temp_pdb):
            os.remove(temp_pdb)
        
        return True
        
    except Exception as e:
        log(f"다음 iteration용 구조 생성 실패: {e}")
        return False


def print_top_structures_summary(top_structures):
    """Top 구조들의 상세 정보 출력"""
    if not top_structures:
        log("저장된 Top 구조가 없습니다.")
        return
    
    log(f"\n{'='*60}")
    log("Top 구조 상세 정보")
    log(f"{'='*60}")
    log(f"{'순위':<4} {'Frame':<6} {'거리(Å)':<8} {'RMSD(Å)':<9} {'조합점수':<10} {'거리순위':<8} {'RMSD순위':<8}")
    log("-" * 60)
    
    for struct in top_structures:
        log(f"{struct['rank']:<4} {struct['frame']:<6} {struct['distance']:<8.2f} "
            f"{struct['rmsd']:<9.2f} {struct['combined_score']:<10.4f} "
            f"{struct.get('distance_rank', 'N/A'):<8} {struct.get('rmsd_rank', 'N/A'):<8}")
    
    log(f"\n통계 요약:")
    log(f"  평균 거리: {np.mean([s['distance'] for s in top_structures]):.2f}Å")
    log(f"  평균 RMSD: {np.mean([s['rmsd'] for s in top_structures]):.2f}Å")
    log(f"  평균 조합점수: {np.mean([s['combined_score'] for s in top_structures]):.4f}")

# ===== 시뮬레이션 실행 함수들 =====
# run_attempt 함수 수정 부분
def run_attempt(work_dir, input_pdb, attempt_num, binding_site_residues, reference_pdb, long_md=False):
    """단일 attempt 실행"""
    log(f"Attempt {attempt_num} 시작 {'(긴 MD)' if long_md else ''}")
    
    # 작업 디렉토리 준비
    attempt_dir = os.path.join(work_dir, f"attempt_{attempt_num}")
    if os.path.exists(attempt_dir):
        shutil.rmtree(attempt_dir)
    os.makedirs(attempt_dir)
    
    # 입력 PDB 복사
    shutil.copy(input_pdb, os.path.join(attempt_dir, "input.pdb"))
    
    # GROMACS 파이프라인 실행
    stages = run_gromacs_pipeline(attempt_dir, "input.pdb", long_md)
    
    # 마지막 단계가 성공했는지 확인
    if not stages or not stages[-1]["success"]:
        log(f"Attempt {attempt_num} 실패: GROMACS 파이프라인 오류")
        return {
            "attempt": attempt_num,
            "success": False,
            "stages": stages,
            "reason": "gromacs_pipeline_failed",
            "long_md_executed": long_md
        }
    
    # 궤적 분석
    tpr_file = os.path.join(attempt_dir, "md.tpr")
    xtc_file = os.path.join(attempt_dir, "md.xtc")
    
    # Long MD인 경우 RMSD도 고려한 Top 구조 저장, 일반 MD인 경우 기존 방식
    if long_md:
        # Reference 구조로 input.pdb 사용
        distances, top_structures = extract_distances_and_top_structures_with_rmsd(
            tpr_file, xtc_file, binding_site_residues, attempt_dir, reference_pdb, save_top_n=SAVE_TOP_N
        )
        log(f"Long MD Top 구조 {len(top_structures)}개 저장됨 (거리+RMSD 기준)")
        
        # Top 구조 상세 정보 출력
        if top_structures:
            print_top_structures_summary(top_structures)
        
    else:
        distances = extract_distances_from_trajectory(tpr_file, xtc_file, binding_site_residues, attempt_dir)
        top_structures = []
    
    if len(distances) < 2:
        log(f"Attempt {attempt_num} 실패: 거리 데이터 부족")
        return {
            "attempt": attempt_num,
            "success": False,
            "stages": stages,
            "distances": distances,
            "reason": "insufficient_distance_data",
            "long_md_executed": long_md,
            "top_structures": top_structures
        }
    
    # 기울기 계산
    slope = calculate_slope(distances)
    min_distance = min(distances)
    accepted = True if long_md else (slope < SLOPE_THRESHOLD)
    
    log(f"Attempt {attempt_num} - 기울기: {slope:.6f}, 최소거리: {min_distance:.2f}Å, 채택: {accepted}")
    
    if accepted:
        # 다음 iteration용 구조 생성 (Long MD가 아닌 경우만)
        if not long_md:
            next_structure = os.path.join(work_dir, "next_structure.pdb")
            structure_success = create_next_iteration_structure(attempt_dir, next_structure, long_md)
            
            if structure_success:
                log(f"Attempt {attempt_num} 성공! 다음 iteration용 구조 생성 완료")
            else:
                log(f"Attempt {attempt_num} 성공! 하지만 구조 생성 실패")
        else:
            log(f"Long MD Attempt {attempt_num} 성공! Top {len(top_structures)}개 구조 저장됨")
            
            # Long MD 결과에 대한 추가 정보
            if top_structures:
                best_struct = top_structures[0]  # 조합 점수 기준 최고 구조
                log(f"최고 구조 정보:")
                log(f"  Frame: {best_struct['frame']}")
                log(f"  거리: {best_struct['distance']:.2f}Å")
                log(f"  RMSD: {best_struct['rmsd']:.2f}Å") 
                log(f"  조합점수: {best_struct['combined_score']:.4f}")
    
    return {
        "attempt": attempt_num,
        "success": accepted,
        "stages": stages,
        "distances": distances,
        "slope": slope,
        "initial_distance": distances[0],
        "final_distance": distances[-1],
        "min_distance": min_distance,
        "long_md_executed": long_md,
        "close_contact_detected": min_distance <= CLOSE_DISTANCE_THRESHOLD,
        "top_structures": top_structures
    }

def run_iteration(work_dir, input_pdb, iteration_num, binding_site_residues, reference_pdb, long_md=False):
    """단일 iteration 실행 - binding_site_residues 매개변수 추가됨"""
    log(f"=== Iteration {iteration_num} 시작 {'(긴 MD)' if long_md else ''} ===")

    for attempt in range(1, MAX_ATTEMPTS + 1):
        result = run_attempt(work_dir, input_pdb, attempt, binding_site_residues, reference_pdb, long_md)
        
        # JSON에 attempt 결과 저장
        attempt_file = os.path.join(work_dir, f"iteration_{iteration_num}_attempt_{attempt}.json")
        with open(attempt_file, "w") as f:
            json.dump(result, f, indent=2, default=str)
        
        if result["success"]:
            log(f"Iteration {iteration_num} 성공! (Attempt {attempt})")
            iteration_result = {
                "iteration": iteration_num,
                "success": True,
                "attempts_used": attempt,
                "final_result": result,
                "long_md": long_md,
                "close_contact_in_iteration": result.get("close_contact_detected", False),
                "final_with_top_structures": long_md  # Long MD인 경우 최종 종료
            }

            if long_md:
                log(f"Iteration {iteration_num} Long MD 완료! 시뮬레이션 최종 종료 (Attempt {attempt})")
                log(f"Top {len(result['top_structures'])}개 구조가 최종 결과로 저장됨")
            else:
                log(f"Iteration {iteration_num} 성공! (Attempt {attempt})")
            
            return iteration_result
    
    log(f"Iteration {iteration_num} 실패: {MAX_ATTEMPTS}번 시도 모두 실패")
    return {
        "iteration": iteration_num,
        "success": False,
        "attempts_used": MAX_ATTEMPTS,
        "final_result": None,
        "long_md": long_md,
        "close_contact_in_iteration": False,
        "final_with_top_structures": False
    }

def run_structure_simulation_with_gpu(structure_info, binding_site_residues, gpu_queue, results_queue, process_id, reference_pdb=None):
    """GPU 할당된 단일 구조 시뮬레이션 실행"""
    structure_pdb, structure_name, output_dir = structure_info
    
    # GPU 할당받기
    assigned_gpu = gpu_queue.get()
    
    try:
        log(f"[Process {process_id}] 구조 {structure_name} GPU {assigned_gpu}에 할당됨")
        
        # 현재 프로세스의 GPU 환경 설정
        original_gpu_id = globals().get('GPU_ID', '0')
        globals()['GPU_ID'] = assigned_gpu
        
        # 구조별 출력 디렉토리 생성
        struct_dir = os.path.join(output_dir, f"structure_{structure_name}")
        if os.path.exists(struct_dir):
            shutil.rmtree(struct_dir)
        os.makedirs(struct_dir)
        
        current_pdb = structure_pdb
        structure_results = {
            "structure_name": structure_name,
            "structure_file": os.path.basename(structure_pdb),
            "assigned_gpu": assigned_gpu,
            "process_id": process_id,
            "start_time": datetime.now().isoformat(),
            "iterations": [],
            "final_top_structures": []  # 최종 Top 구조들
        }
        
        iteration = 0
        need_long_md = False
        first_dir = None
        simulation_completed = False  # 시뮬레이션 완료 플래그

        while iteration < MAX_ITERATIONS and not simulation_completed:
            iteration += 1
            
            # iteration 디렉토리 생성
            iter_dir = os.path.join(struct_dir, f'iteration_{iteration}')
            if os.path.exists(iter_dir):
                shutil.rmtree(iter_dir)
            os.makedirs(iter_dir)
            
            if first_dir is None:
                first_dir = os.path.join(iter_dir, 'attempt_1')
            
            # iteration 실행
            iteration_result = run_iteration(iter_dir, current_pdb, iteration, binding_site_residues, reference_pdb, need_long_md)
            structure_results["iterations"].append(iteration_result)
            
            # iteration 결과 JSON 저장
            iteration_file = os.path.join(struct_dir, f"iteration_{iteration}_summary.json")
            with open(iteration_file, "w") as f:
                json.dump(iteration_result, f, indent=2, default=str)
            
            if iteration_result["success"]:
                # Long MD 완료시 최종 종료
                if iteration_result.get("final_with_top_structures", False):
                    log(f"[Process {process_id}] 구조 {structure_name}: Long MD 완료로 시뮬레이션 최종 종료")
                    
                    # Top 구조들을 전체 결과에 저장
                    if iteration_result["final_result"].get("top_structures"):
                        structure_results["final_top_structures"] = iteration_result["final_result"]["top_structures"]
                        
                        # Top 구조들을 구조별 디렉토리로 복사
                        for i, top_struct in enumerate(structure_results["final_top_structures"]):
                            source_path = top_struct["path"]
                            dest_name = f"final_top_{i+1}_{top_struct['filename']}"
                            dest_path = os.path.join(struct_dir, dest_name)
                            
                            success = restore_original_chain_ids(source_path, source_path, first_dir)
                            if success:
                                log(f"[Process {process_id}] 구조 {structure_name} 최종 구조 {dest_name}: 저장 완료 (체인 복원됨)")
                            else:
                                log(f"[Process {process_id}] 구조 {structure_name} 최종 구조 {dest_name}: 체인 복원 실패, 변환된 구조 사용")

                            try:
                                shutil.copy(source_path, dest_path)
                                # 경로 업데이트
                                top_struct["final_path"] = dest_path
                                log(f"[Process {process_id}] Top {i+1} 구조 복사: {dest_name}")
                            except Exception as e:
                                log(f"[Process {process_id}] Top 구조 복사 실패: {e}")
                    
                    simulation_completed = True
                else:
                    # 일반 iteration 성공 - 다음 iteration용 PDB 업데이트
                    next_structure = os.path.join(iter_dir, "next_structure.pdb")
                    if ENABLE_CHAIN_RESTORATION:
                        success = restore_original_chain_ids(next_structure, next_structure, first_dir)
                        if success:
                            log(f"[Process {process_id}] 구조 {structure_name}: 다음 구조 저장 완료 (체인 복원됨)")
                        else:
                            log(f"[Process {process_id}] 구조 {structure_name}: 체인 복원 실패, 변환된 구조 사용")
                    
                    if os.path.exists(next_structure):
                        current_pdb = next_structure
                    
                    # 근접 접촉 검사
                    if ENABLE_LONG_MD and iteration_result.get("close_contact_in_iteration", False):
                        if not need_long_md:
                            need_long_md = True
                            log(f"[Process {process_id}] 구조 {structure_name}: 근접 접촉 감지! 긴 MD 예정")
                    else:
                        need_long_md = False
                if iteration+1==MAX_ITERATIONS-1:
                    need_long_md=True 
            else:
                # iteration 실패 시 처음부터 다시 시작
                log(f"[Process {process_id}] 구조 {structure_name}: Iteration {iteration} 실패 - 재시작")
                current_pdb = structure_pdb
                iteration = 0
                need_long_md = False
        
        # 구조 시뮬레이션 완료
        structure_results["end_time"] = datetime.now().isoformat()
        structure_results["total_iterations"] = len(structure_results["iterations"])
        structure_results["successful_iterations"] = len([r for r in structure_results["iterations"] if r["success"]])
        structure_results["long_md_executed"] = any(r.get("long_md", False) for r in structure_results["iterations"])
        structure_results["simulation_completed"] = simulation_completed
        
        # 구조별 결과 저장
        structure_result_file = os.path.join(struct_dir, "structure_results.json")
        with open(structure_result_file, "w") as f:
            json.dump(structure_results, f, indent=2, default=str)
        
        completion_msg = "Long MD로 완료" if simulation_completed else f"{structure_results['successful_iterations']}/{structure_results['total_iterations']} 성공"
        log(f"[Process {process_id}] 구조 {structure_name} (GPU {assigned_gpu}) 완료: {completion_msg}")
        
        if structure_results["final_top_structures"]:
            log(f"[Process {process_id}] 최종 Top 구조 {len(structure_results['final_top_structures'])}개 저장됨")
        
        # 결과를 큐에 넣기
        results_queue.put(structure_results)
        
    except Exception as e:
        log(f"[Process {process_id}] 구조 {structure_name} (GPU {assigned_gpu}) 오류: {e}")
        import traceback
        log(f"상세 오류: {traceback.format_exc()}")
        
        # 오류 발생 시에도 기본 결과 반환
        error_result = {
            "structure_name": structure_name,
            "structure_file": os.path.basename(structure_pdb),
            "assigned_gpu": assigned_gpu,
            "process_id": process_id,
            "start_time": datetime.now().isoformat(),
            "end_time": datetime.now().isoformat(),
            "iterations": [],
            "total_iterations": 0,
            "successful_iterations": 0,
            "final_top_structures": [],
            "error": str(e)
        }
        results_queue.put(error_result)
        
    finally:
        # GPU 반납
        gpu_queue.put(assigned_gpu)
        log(f"[Process {process_id}] GPU {assigned_gpu} 반납됨")

def run_parallel_gpu_simulation(structure_pool, output_dir, binding_site_residues, reference_pdb=None):
    """GPU 병렬 시뮬레이션 실행"""
    log("=== GPU 병렬 시뮬레이션 시작 ===")
    
    # GPU ID 파싱
    available_gpus = parse_gpu_ids(GPU_ID)
    max_processes = len(available_gpus)
    
    # 멀티프로세싱 매니저 생성
    manager = Manager()
    gpu_queue = manager.Queue()
    results_queue = manager.Queue()
    
    # GPU 큐에 사용 가능한 GPU들 추가
    for gpu in available_gpus:
        gpu_queue.put(gpu)
    
    log(f"병렬 처리: {max_processes}개 GPU 동시 사용")
    log(f"처리할 구조 수: {len(structure_pool)}")
    
    # 구조 정보 리스트 준비
    structure_infos = []
    for i, structure_pdb in enumerate(structure_pool):
        structure_name = f"struct_{i:02d}"
        structure_infos.append((structure_pdb, structure_name, output_dir))
    
    # 프로세스 리스트
    processes = []
    all_structure_results = []
    completed_count = 0
    
    # 구조별로 프로세스 생성 및 실행
    for i, info in enumerate(structure_infos):
        # 최대 프로세스 수만큼만 동시 실행
        if len(processes) >= max_processes:
            # 완료된 프로세스 대기
            for p in processes[:]:
                if not p.is_alive():
                    p.join()
                    processes.remove(p)
                    
                    # 결과 수집
                    if not results_queue.empty():
                        result = results_queue.get()
                        all_structure_results.append(result)
                        completed_count += 1
                        log(f"[{completed_count}/{len(structure_infos)}] 구조 {result['structure_name']} 완료: {result.get('successful_iterations', 0)}회 성공")
            
            # 아직도 최대치면 잠시 대기
            while len(processes) >= max_processes:
                time.sleep(1)
                for p in processes[:]:
                    if not p.is_alive():
                        p.join()
                        processes.remove(p)
                        
                        # 결과 수집
                        if not results_queue.empty():
                            result = results_queue.get()
                            all_structure_results.append(result)
                            completed_count += 1
                            log(f"[{completed_count}/{len(structure_infos)}] 구조 {result['structure_name']} 완료: {result.get('successful_iterations', 0)}회 성공")
        
        # 새 프로세스 시작
        process = Process(
            target=run_structure_simulation_with_gpu,
            args=(info, binding_site_residues, gpu_queue, results_queue, i, reference_pdb)
        )
        process.start()
        processes.append(process)
        log(f"구조 {info[1]} 프로세스 시작 (PID: {process.pid})")
    
    # 모든 프로세스 완료 대기
    for process in processes:
        process.join()
    
    # 남은 결과 수집
    while not results_queue.empty():
        result = results_queue.get()
        all_structure_results.append(result)
        completed_count += 1
        log(f"[{completed_count}/{len(structure_infos)}] 구조 {result['structure_name']} 완료: {result.get('successful_iterations', 0)}회 성공")
    
    log("=== 모든 GPU 병렬 시뮬레이션 완료 ===")
    
    # 결과를 원래 순서대로 정렬
    all_structure_results.sort(key=lambda x: x['structure_name'])
    
    return all_structure_results

# ===== 메인 함수 =====
def main():
    if len(sys.argv) < 4:
        print("사용법: python simple_sumd.py <input_pdb> <chain1> <chain2> [output_dir]")
        sys.exit(1)
    
    input_pdb = sys.argv[1]
    chain1 = sys.argv[2]
    chain2 = sys.argv[3]
    output_dir = sys.argv[4] if len(sys.argv) > 4 else "sumd_output"
    
    # 출력 디렉토리 준비
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)
    os.makedirs(output_dir)
    
    log("=== Simple SuMD 시작 (모든 구조에 대한 완전한 시뮬레이션) ===")
    log(f"입력 PDB: {input_pdb}")
    log(f"체인: {chain1} - {chain2}")
    log(f"출력 디렉토리: {output_dir}")
    log(f"최대 iterations: {MAX_ITERATIONS}")
    log(f"다방향 이격: {'활성화' if ENABLE_MULTI_DIRECTION_SEPARATION else '비활성화'}")
    log(f"회전 변형: {'활성화' if ENABLE_ROTATIONAL_VARIANTS else '비활성화'}")
    
    # 원본 chain 정보 저장
    # original_chains = save_original_chain_info(input_pdb, chain1, chain2)
    # if original_chains is None:
    #     log("원본 chain 정보 저장 실패")
    #     sys.exit(1)
    
    # 타겟 체인만 추출
    target_pdb = os.path.join(output_dir, "target_chains.pdb")
    if not extract_target_chains_pdb(input_pdb, target_pdb, chain1, chain2):
        log("타겟 체인 추출 실패")
        sys.exit(1)
    
    # Ligand/Receptor 식별 및 Binding site 정의
    ligand_chain, receptor_chain = identify_ligand_receptor(target_pdb)
    if ligand_chain is None or receptor_chain is None:
        log("Ligand/Receptor 식별 실패")
        sys.exit(1)
    
    binding_site_residues = define_binding_site(target_pdb, ligand_chain, receptor_chain, distance_cutoff=BINDING_SITE_CUTOFF)
    if binding_site_residues is None:
        log("Binding site 정의 실패")
        sys.exit(1)
    
    # 구조 풀 생성
    log("=== 구조 풀 생성 ===")
    structure_pool = create_initial_structure_pool(target_pdb, chain1, chain2, output_dir)
    
    log(f"생성된 구조 수: {len(structure_pool)}")
    
    # 각 구조에 대해 완전한 시뮬레이션 수행
    all_results = {
        "start_time": datetime.now().isoformat(),
        "input_pdb": input_pdb,
        "chains": [chain1, chain2],
        "ligand_chain": ligand_chain,
        "receptor_chain": receptor_chain,
        "binding_site_residues": binding_site_residues,
        "total_structures": len(structure_pool),
        "structure_results": [],
        "final_top_structures_summary": {
            "total_structures_with_tops": 0,
            "best_overall_distance": float('inf'),
            "best_overall_structure": None,
            "all_top_structures": []
        }
    }
    
    # 시뮬레이션 설정 정보 로깅
    log("=== 시뮬레이션 설정 정보 ===")
    log(f"입력 PDB: {all_results['input_pdb']}")
    log(f"체인: {all_results['chains'][0]} - {all_results['chains'][1]}")
    log(f"Ligand: {all_results['ligand_chain']} | Receptor: {all_results['receptor_chain']}")
    log(f"Binding site 잔기 수: {len(all_results['binding_site_residues'])}")
    log(f"총 구조 수: {all_results['total_structures']}")
    log(f"시작 시간: {all_results['start_time']}")
    log("=" * 50)
    
    # GPU 병렬 시뮬레이션 실행
    log("=== GPU 병렬 시뮬레이션 시작 ===")

    # 구조 정보 리스트 준비
    structure_infos = []
    for i, structure_pdb in enumerate(structure_pool):
        structure_name = f"struct_{i:02d}"
        structure_infos.append((structure_pdb, structure_name, output_dir))
    
    # 병렬 실행
    all_structure_results = []
    # completed_count = 0
    
    all_structure_results = run_parallel_gpu_simulation(structure_pool, output_dir, binding_site_residues, target_pdb)
    all_structure_results.sort(key=lambda x: x['structure_name'])
    all_results["structure_results"] = all_structure_results
    
    log("=== 모든 GPU 병렬 시뮬레이션 완료 ===")
    
    # 전체 결과 종합
    all_results["end_time"] = datetime.now().isoformat()
    all_results["total_successful_structures"] = len([r for r in all_results["structure_results"] 
                                                     if r["successful_iterations"] > 0])
    
    # Top 구조 정보 수집 및 분석 (개선됨)
    log("=== Top 구조 정보 수집 및 분석 ===")
    all_top_structures = []
    structures_with_tops = 0
    best_combined_score = float('inf')
    best_distance = float('inf')
    best_rmsd = float('inf')
    best_structure_info = None
    
    for struct_result in all_results["structure_results"]:
        if struct_result.get("final_top_structures"):
            structures_with_tops += 1
            structure_name = struct_result["structure_name"]
            
            for top_struct in struct_result["final_top_structures"]:
                # 구조 정보에 추가 메타데이터 포함
                enhanced_top = {
                    **top_struct,
                    "source_structure": structure_name,
                    "source_iterations": struct_result.get("successful_iterations", 0)
                }
                all_top_structures.append(enhanced_top)
                
                # 전체 최고 조합 점수 구조 찾기
                combined_score = top_struct.get('combined_score', float('inf'))
                if combined_score < best_combined_score:
                    best_combined_score = combined_score
                    best_structure_info = {
                        "structure_name": structure_name,
                        "top_info": enhanced_top,
                        "structure_dir": os.path.join(output_dir, f"structure_{structure_name}"),
                        "selection_criteria": "combined_score"
                    }
                
                # 개별 기준 최고값들도 추적
                if top_struct["distance"] < best_distance:
                    best_distance = top_struct["distance"]
                if top_struct.get("rmsd", float('inf')) < best_rmsd:
                    best_rmsd = top_struct.get("rmsd", float('inf'))
    
    # Top 구조 요약 정보 업데이트 (개선됨)
    all_results["final_top_structures_summary"] = {
        "total_structures_with_tops": structures_with_tops,
        "best_overall_combined_score": best_combined_score if best_combined_score != float('inf') else None,
        "best_overall_distance": best_distance if best_distance != float('inf') else None,
        "best_overall_rmsd": best_rmsd if best_rmsd != float('inf') else None,
        "best_overall_structure": best_structure_info,
        "total_top_structures": len(all_top_structures),
        "scoring_method": "distance_rmsd_combined",
        "top_structures_by_combined_score": sorted(all_top_structures, key=lambda x: x.get('combined_score', float('inf')))[:10],
        "top_structures_by_distance": sorted(all_top_structures, key=lambda x: x['distance'])[:5],
        "top_structures_by_rmsd": sorted([s for s in all_top_structures if 'rmsd' in s], key=lambda x: x['rmsd'])[:5]
    }
    
    # 최고 성능 구조 및 Top 구조 처리 (개선됨)
    if best_structure_info:
        log(f"전체 최고 구조 발견: {best_structure_info['structure_name']}")
        log(f"선정 기준: 조합점수 (거리+RMSD)")
        log(f"조합점수: {best_combined_score:.4f}")
        log(f"거리: {best_structure_info['top_info']['distance']:.2f}Å")
        if 'rmsd' in best_structure_info['top_info']:
            log(f"RMSD: {best_structure_info['top_info']['rmsd']:.2f}Å")
        log(f"Frame: {best_structure_info['top_info']['frame']}")
        
        # 전체 최고 구조를 메인 디렉토리에 복사
        best_top_source = best_structure_info['top_info'].get('final_path')
        if not best_top_source or not os.path.exists(best_top_source):
            best_top_source = best_structure_info['top_info'].get('path')
        
        if best_top_source and os.path.exists(best_top_source):
            best_final_path = os.path.join(output_dir, "best_overall_structure.pdb")
            shutil.copy(best_top_source, best_final_path)
            log(f"전체 최고 구조 복사: {best_final_path}")
        
        # 다양한 기준별 상위 구조들을 메인 디렉토리에 복사
        log("다양한 기준별 Top 구조 복사 중...")
        
        # 조합점수 기준 상위 5개
        top_combined = all_results["final_top_structures_summary"]["top_structures_by_combined_score"][:5]
        for i, top_struct in enumerate(top_combined):
            source_path = top_struct.get('final_path') or top_struct.get('path')
            if source_path and os.path.exists(source_path):
                combined_score = top_struct.get('combined_score', 0)
                dest_name = f"top_combined_{i+1}_{top_struct['source_structure']}_frame_{top_struct['frame']}_score_{combined_score:.4f}.pdb"
                dest_path = os.path.join(output_dir, dest_name)
                
                try:
                    shutil.copy(source_path, dest_path)
                    log(f"조합점수 Top {i+1}: {dest_name} (점수: {combined_score:.4f})")
                except Exception as e:
                    log(f"조합점수 Top {i+1} 복사 실패: {e}")
        
        # 거리 기준 상위 3개
        top_distance = all_results["final_top_structures_summary"]["top_structures_by_distance"][:3]
        for i, top_struct in enumerate(top_distance):
            source_path = top_struct.get('final_path') or top_struct.get('path')
            if source_path and os.path.exists(source_path):
                dest_name = f"top_distance_{i+1}_{top_struct['source_structure']}_frame_{top_struct['frame']}_dist_{top_struct['distance']:.2f}A.pdb"
                dest_path = os.path.join(output_dir, dest_name)
                
                try:
                    shutil.copy(source_path, dest_path)
                    log(f"거리 Top {i+1}: {dest_name} ({top_struct['distance']:.2f}Å)")
                except Exception as e:
                    log(f"거리 Top {i+1} 복사 실패: {e}")
        
        # RMSD 기준 상위 3개 (RMSD가 있는 구조만)
        top_rmsd = all_results["final_top_structures_summary"]["top_structures_by_rmsd"][:3]
        for i, top_struct in enumerate(top_rmsd):
            source_path = top_struct.get('final_path') or top_struct.get('path')
            if source_path and os.path.exists(source_path):
                dest_name = f"top_rmsd_{i+1}_{top_struct['source_structure']}_frame_{top_struct['frame']}_rmsd_{top_struct['rmsd']:.2f}A.pdb"
                dest_path = os.path.join(output_dir, dest_name)
                
                try:
                    shutil.copy(source_path, dest_path)
                    log(f"RMSD Top {i+1}: {dest_name} ({top_struct['rmsd']:.2f}Å)")
                except Exception as e:
                    log(f"RMSD Top {i+1} 복사 실패: {e}")
    
    # 최종 결과 로그 (개선됨)
    log("=== Simple SuMD 완료 ===")
    log(f"처리된 구조 수: {len(structure_pool)}")
    log(f"성공한 구조 수: {all_results['total_successful_structures']}")
    log(f"전체 성공률: {all_results['total_successful_structures']}/{len(structure_pool)}")
    
    # Top 구조 결과 요약 (개선됨)
    if all_results["final_top_structures_summary"]["total_structures_with_tops"] > 0:
        log(f"Top 구조 보유 구조 수: {all_results['final_top_structures_summary']['total_structures_with_tops']}")
        log(f"총 Top 구조 수: {all_results['final_top_structures_summary']['total_top_structures']}")
        
        if all_results["final_top_structures_summary"]["best_overall_combined_score"]:
            log(f"전체 최고 조합점수: {all_results['final_top_structures_summary']['best_overall_combined_score']:.4f}")
        if all_results["final_top_structures_summary"]["best_overall_distance"]:
            log(f"전체 최고 거리: {all_results['final_top_structures_summary']['best_overall_distance']:.2f}Å")
        if all_results["final_top_structures_summary"]["best_overall_rmsd"]:
            log(f"전체 최고 RMSD: {all_results['final_top_structures_summary']['best_overall_rmsd']:.2f}Å")
    
    log(f"최종 결과: {final_result_file}")
    log("주요 결과 파일:")
    log(f"  - best_overall_structure.pdb (전체 최고 조합점수 구조)")
    log(f"  - top_combined_*.pdb (조합점수 기준 상위 구조들)")
    log(f"  - top_distance_*.pdb (거리 기준 상위 구조들)")
    log(f"  - top_rmsd_*.pdb (RMSD 기준 상위 구조들)")

if __name__ == "__main__":
    main()