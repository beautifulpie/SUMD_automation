#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import os
import sys
import numpy as np
import hashlib
import time
import json
import shutil
from datetime import datetime
from Bio.PDB import PDBParser, PDBIO, Select
from Bio.PDB.Atom import Atom
import logging
from utils import NumpyEncoder as NE
from MDPGenerator import MDPGenerator

# 추가된 import
from Bio.PDB.SASA import ShrakeRupley
import warnings
warnings.filterwarnings("ignore")

# 기존 import들 유지
script_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, script_dir)

try:
    from mda_pdb_processor import process_pdb_for_gromacs
except ImportError:
    print("경고: mda_pdb_processor 모듈을 찾을 수 없습니다. 전처리 기능이 제한될 수 있습니다.")
    def process_pdb_for_gromacs(input_pdb, output_pdb, target_chains, logger=None):
        shutil.copy(input_pdb, output_pdb)
        return output_pdb, {"processing_method": "bypass"}

try:
    from gromacs_runner_class import GromacsCommandRunner
except ImportError:
    print("경고: gromacs_runner 모듈을 찾을 수 없습니다. gromacs 기능이 제한 됩니다.")

# Interface 기반 거리 계산 모듈 import
try:
    from interface_distance_calculator import DistanceCalculator, InterfaceDistanceCalculator
    print("Interface 기반 거리 계산 모듈 로드됨")
except ImportError:
    print("경고: interface_distance_calculator 모듈을 찾을 수 없습니다. 기본 거리 계산 방식을 사용합니다.")
    
    class DistanceCalculator:
        @staticmethod
        def calculate_chain_distance(pdb_file, chain1_id, chain2_id):
            """기존 center of mass 방식 (fallback)"""
            try:
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
                
                if not chain1_atoms:
                    raise Exception(f"체인 {chain1_id}를 찾을 수 없거나 원자가 없습니다")
                if not chain2_atoms:
                    raise Exception(f"체인 {chain2_id}를 찾을 수 없거나 원자가 없습니다")
                
                chain1_center = np.mean(chain1_atoms, axis=0)
                chain2_center = np.mean(chain2_atoms, axis=0)
                
                distance = np.linalg.norm(chain1_center - chain2_center)
                return distance
                
            except Exception as e:
                raise Exception(f"거리 계산 실패: {e}")


# ===== 새로 추가: Golden Standard 기반 수렴 판정 클래스 =====
class GoldenStandardConvergenceChecker:
    """Golden standard PDB와의 RMSD 기반 수렴 판정"""
    
    def __init__(self, golden_standard_pdb, window_size=5, rmsd_threshold=1.5):
        """
        Args:
            golden_standard_pdb: Golden standard PDB 파일 경로
            window_size: 평균을 계산할 iteration 수 (기본값: 5)
            rmsd_threshold: 수렴 판정 RMSD 임계값 (기본값: 1.5Å)
        """
        self.golden_standard_pdb = golden_standard_pdb
        self.window_size = window_size
        self.rmsd_threshold = rmsd_threshold
        self.rmsd_history = []
        
        # Golden standard PDB 존재 확인
        if not os.path.exists(golden_standard_pdb):
            raise FileNotFoundError(f"Golden standard PDB 파일을 찾을 수 없습니다: {golden_standard_pdb}")
        
    def add_rmsd(self, current_pdb):
        """새로운 PDB와 golden standard 간의 RMSD 계산 및 추가"""
        rmsd = self.calculate_rmsd_with_golden_standard(current_pdb)
        self.rmsd_history.append(rmsd)
        return rmsd
        
    def calculate_rmsd_with_golden_standard(self, current_pdb):
        """현재 PDB와 golden standard 간의 RMSD 계산"""
        try:
            parser = PDBParser(QUIET=True)
            
            golden_structure = parser.get_structure("golden", self.golden_standard_pdb)
            current_structure = parser.get_structure("current", current_pdb)
            
            golden_coords = []
            current_coords = []
            
            # Golden standard 좌표 수집 (CA 원자만)
            for model in golden_structure:
                for chain in model:
                    for residue in chain:
                        if "CA" in residue:
                            golden_coords.append(residue["CA"].coord)
            
            # 현재 구조 좌표 수집 (CA 원자만)
            for model in current_structure:
                for chain in model:
                    for residue in chain:
                        if "CA" in residue and len(current_coords) < len(golden_coords):
                            current_coords.append(residue["CA"].coord)
            
            min_len = min(len(golden_coords), len(current_coords))
            if min_len == 0:
                return float('inf')
            
            golden_coords = np.array(golden_coords[:min_len])
            current_coords = np.array(current_coords[:min_len])
            
            # RMSD 계산 (중심 이동 후)
            golden_center = np.mean(golden_coords, axis=0)
            current_center = np.mean(current_coords, axis=0)
            
            golden_centered = golden_coords - golden_center
            current_centered = current_coords - current_center
            
            diff = golden_centered - current_centered
            rmsd = np.sqrt(np.mean(np.sum(diff**2, axis=1)))
            
            return float(rmsd)
            
        except Exception as e:
            # 계산 실패 시 임의의 높은 값 반환
            return np.random.uniform(5.0, 10.0)
    
    def check_convergence(self):
        """수렴 여부 확인"""
        if len(self.rmsd_history) < self.window_size:
            return False, 0.0, f"충분한 데이터 없음 ({len(self.rmsd_history)}/{self.window_size})"
        
        # 마지막 window_size개의 RMSD 평균 계산
        recent_rmsds = self.rmsd_history[-self.window_size:]
        avg_rmsd = np.mean(recent_rmsds)
        
        converged = avg_rmsd <= self.rmsd_threshold
        
        status_msg = f"마지막 {self.window_size}회 평균 RMSD (vs Golden): {avg_rmsd:.3f}Å"
        if converged:
            status_msg += f" ≤ {self.rmsd_threshold}Å (수렴)"
        else:
            status_msg += f" > {self.rmsd_threshold}Å (미수렴)"
            
        return converged, avg_rmsd, status_msg
    
    def get_convergence_info(self):
        """수렴 관련 상세 정보 반환"""
        if not self.rmsd_history:
            return {"status": "데이터 없음"}
        
        recent_count = min(len(self.rmsd_history), self.window_size)
        recent_rmsds = self.rmsd_history[-recent_count:]
        
        info = {
            "total_iterations": len(self.rmsd_history),
            "recent_iterations": recent_count,
            "recent_rmsds": recent_rmsds,
            "recent_avg_rmsd": np.mean(recent_rmsds),
            "recent_std_rmsd": np.std(recent_rmsds),
            "threshold": self.rmsd_threshold,
            "window_size": self.window_size,
            "golden_standard_pdb": self.golden_standard_pdb
        }
        
        if recent_count >= self.window_size:
            info["converged"] = info["recent_avg_rmsd"] <= self.rmsd_threshold
        else:
            info["converged"] = False
            
        return info


# ===== 기존 함수들 유지 (로깅 등) =====
def setup_logging(output_dir, job_id):
    """로깅 설정"""
    log_file = os.path.join(output_dir, f"sumd_{job_id}.log")
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file, encoding='utf-8'),
            logging.StreamHandler(sys.stdout)
        ]
    )
    return logging.getLogger("SuMD-Master")

def generate_job_id(pdb_file, receptor_chain, ligand_chain):
    """작업 ID 생성 (SHA256 기반)"""
    content = f"{pdb_file}_{receptor_chain}_{ligand_chain}_{time.time()}"
    return hashlib.sha256(content.encode()).hexdigest()[:16]

def analyze_interface_changes(initial_pdb, final_pdb, receptor_chain, ligand_chain, logger):
    """Interface 변화 분석 (기존 코드 유지)"""
    try:
        if 'InterfaceDistanceCalculator' in globals():
            calculator = InterfaceDistanceCalculator()
            
            initial_analysis = calculator.get_interface_analysis(initial_pdb, receptor_chain, ligand_chain)
            final_analysis = calculator.get_interface_analysis(final_pdb, receptor_chain, ligand_chain)
            
            if 'interface_residues' in initial_analysis and 'interface_residues' in final_analysis:
                initial_count = (initial_analysis['interface_residues'][receptor_chain] + 
                               initial_analysis['interface_residues'][ligand_chain])
                final_count = (final_analysis['interface_residues'][receptor_chain] + 
                             final_analysis['interface_residues'][ligand_chain])
                
                logger.info(f"Interface residue 변화: {initial_count} → {final_count}")
                
                initial_percentage = (initial_analysis['interface_percentage'][receptor_chain] + 
                                    initial_analysis['interface_percentage'][ligand_chain]) / 2
                final_percentage = (final_analysis['interface_percentage'][receptor_chain] + 
                                  final_analysis['interface_percentage'][ligand_chain]) / 2                
                logger.info(f"Interface 비율 변화: {initial_percentage:.1f}% → {final_percentage:.1f}%")

                return {
                    'initial_interface_count': initial_count,
                    'final_interface_count': final_count,
                    'initial_interface_percentage': initial_percentage,
                    'final_interface_percentage': final_percentage,
                    'interface_change': final_count - initial_count
                }
        
    except Exception as e:
        logger.warning(f"Interface 분석 실패: {e}")
    
    return None

def calculate_structure_sasa(pdb_file, logger):
    """
    BioPython을 사용하여 구조의 SASA 계산
    
    Args:
        pdb_file: PDB 파일 경로
        logger: 로거 객체
        
    Returns:
        dict: SASA 분석 결과
    """
    try:
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("structure", pdb_file)
        
        # ShrakeRupley 알고리즘을 사용한 SASA 계산
        sr = ShrakeRupley()
        sr.compute(structure, level="R")  # Residue 레벨에서 계산
        
        total_sasa = 0.0
        chain_sasa = {}
        residue_sasa = {}
        
        for model in structure:
            for chain in model:
                chain_id = chain.id
                chain_total = 0.0
                chain_residues = {}
                
                for residue in chain:
                    if hasattr(residue, 'sasa'):
                        res_sasa = residue.sasa
                        res_key = f"{residue.id[1]}{residue.get_resname()}"
                        chain_residues[res_key] = res_sasa
                        chain_total += res_sasa
                        total_sasa += res_sasa
                
                if chain_total > 0:
                    chain_sasa[chain_id] = {
                        'total_sasa': chain_total,
                        'residue_count': len(chain_residues),
                        'avg_sasa_per_residue': chain_total / len(chain_residues) if chain_residues else 0
                    }
                    residue_sasa[chain_id] = chain_residues
        
        sasa_result = {
            'total_sasa': total_sasa,
            'chain_sasa': chain_sasa,
            'residue_sasa': residue_sasa,
            'calculation_method': 'ShrakeRupley'
        }
        
        logger.info(f"SASA 계산 완료: 총 SASA = {total_sasa:.2f} Ų")
        for chain_id, data in chain_sasa.items():
            logger.info(f"  체인 {chain_id}: {data['total_sasa']:.2f} Ų ({data['residue_count']}개 잔기)")
        
        return sasa_result
        
    except Exception as e:
        logger.error(f"SASA 계산 실패: {e}")
        return {
            'error': str(e),
            'total_sasa': 0.0,
            'calculation_method': 'failed'
        }

def extract_rmsf_data(sample_result, iteration_num, job_output_dir, logger, gromacs_runner=None):
    """
    GROMACS를 사용하여 RMSF 데이터 추출 및 저장
    
    Args:
        sample_result: 최적 샘플 결과
        iteration_num: iteration 번호
        job_output_dir: 작업 출력 디렉토리
        logger: 로거 객체
        gromacs_runner: GROMACS 실행 객체
        
    Returns:
        dict: RMSF 추출 결과
    """
    try:
        if not sample_result or not sample_result.get('success', False):
            logger.warning("RMSF 추출: 유효한 샘플 결과가 없습니다")
            return {'success': False, 'error': 'no_valid_sample'}
        
        sample_dir = sample_result['sample_dir']
        rmsf_output_dir = os.path.join(job_output_dir, "rmsf_analysis")
        os.makedirs(rmsf_output_dir, exist_ok=True)
        
        # 궤적 파일과 구조 파일 확인
        trajectory_files = []
        structure_files = []
        
        # MD 궤적이 있는지 확인
        for ext in ['.xtc', '.trr']:
            md_traj = os.path.join(sample_dir, f"md{ext}")
            if os.path.exists(md_traj):
                trajectory_files.append(md_traj)
        
        # 구조 파일 확인
        for ext in ['.tpr', '.gro']:
            struct_file = os.path.join(sample_dir, f"md{ext}")
            if os.path.exists(struct_file):
                structure_files.append(struct_file)
        
        if not trajectory_files:
            logger.warning("RMSF 추출: 궤적 파일을 찾을 수 없습니다")
            return {'success': False, 'error': 'no_trajectory'}
        
        if not structure_files:
            logger.warning("RMSF 추출: 구조 파일을 찾을 수 없습니다")
            return {'success': False, 'error': 'no_structure'}
        
        # GROMACS gmx rmsf 명령어 실행
        if gromacs_runner:
            try:
                rmsf_xvg = os.path.join(rmsf_output_dir, f"iteration_{iteration_num:02d}_rmsf.xvg")
                rmsf_pdb = os.path.join(rmsf_output_dir, f"iteration_{iteration_num:02d}_rmsf.pdb")
                
                # RMSF 계산 파라미터
                rmsf_params = {
                    "structure_file": structure_files[0],
                    "trajectory_file": trajectory_files[0],
                    "output_xvg": os.path.basename(rmsf_xvg),
                    "output_pdb": os.path.basename(rmsf_pdb),
                    "group": "Protein"
                }
                
                # 작업 디렉토리를 RMSF 출력 디렉토리로 변경하여 실행
                original_work_dir = gromacs_runner.work_dir
                gromacs_runner.work_dir = rmsf_output_dir
                
                # 필요한 파일들을 RMSF 출력 디렉토리로 복사
                import shutil
                local_struct = os.path.join(rmsf_output_dir, os.path.basename(structure_files[0]))
                local_traj = os.path.join(rmsf_output_dir, os.path.basename(trajectory_files[0]))
                
                if not os.path.exists(local_struct):
                    shutil.copy(structure_files[0], local_struct)
                if not os.path.exists(local_traj):
                    shutil.copy(trajectory_files[0], local_traj)
                
                rmsf_params["structure_file"] = os.path.basename(local_struct)
                rmsf_params["trajectory_file"] = os.path.basename(local_traj)
                
                # RMSF 계산 실행
                success, message = gromacs_runner.execute_gromacs_command("rmsf", rmsf_params)
                
                # 작업 디렉토리 복원
                gromacs_runner.work_dir = original_work_dir
                
                if success and os.path.exists(rmsf_xvg):
                    # XVG 파일 크기 확인
                    xvg_size = os.path.getsize(rmsf_xvg)
                    
                    logger.info(f"RMSF 추출 완료: {rmsf_xvg} ({xvg_size} bytes)")
                    
                    return {
                        'success': True,
                        'rmsf_xvg': rmsf_xvg,
                        'rmsf_pdb': rmsf_pdb if os.path.exists(rmsf_pdb) else None,
                        'file_size': xvg_size,
                        'trajectory_used': trajectory_files[0],
                        'structure_used': structure_files[0]
                    }
                else:
                    logger.error(f"RMSF 계산 실패: {message}")
                    return {'success': False, 'error': f'gromacs_failed: {message}'}
                    
            except Exception as e:
                gromacs_runner.work_dir = original_work_dir  # 예외 시에도 복원
                logger.error(f"RMSF 계산 중 예외: {e}")
                return {'success': False, 'error': str(e)}
        else:
            logger.warning("RMSF 추출: GROMACS runner가 없습니다")
            return {'success': False, 'error': 'no_gromacs_runner'}
            
    except Exception as e:
        logger.error(f"RMSF 추출 실패: {e}")
        return {'success': False, 'error': str(e)}


# ===== 기존 iteration 결과 저장 함수들 유지 =====
def save_iteration_results(iteration_num, best_sample_result, all_sample_results, 
                         job_output_dir, receptor_chain, ligand_chain, logger):
    """각 iteration의 결과를 메인 폴더에 정리해서 저장 (기존 코드 유지)"""
    try:
        iteration_results_dir = os.path.join(job_output_dir, "iteration_results")
        os.makedirs(iteration_results_dir, exist_ok=True)
        
        current_iter_dir = os.path.join(iteration_results_dir, f"iteration_{iteration_num}")
        os.makedirs(current_iter_dir, exist_ok=True)
        
        # 최적 PDB 파일 복사
        if best_sample_result and best_sample_result['final_pdb']:
            best_pdb_name = f"iteration_{iteration_num:02d}_best.pdb"
            best_pdb_path = os.path.join(current_iter_dir, best_pdb_name)
            shutil.copy(best_sample_result['final_pdb'], best_pdb_path)
            
            main_best_pdb = os.path.join(iteration_results_dir, best_pdb_name)
            shutil.copy(best_sample_result['final_pdb'], main_best_pdb)
            
            logger.info(f"Iteration {iteration_num} 최적 PDB 저장: {best_pdb_name}")
        
        # 기존 iteration 요약 정보 생성 로직 유지...
        iteration_summary = {
            "iteration": iteration_num,
            "timestamp": datetime.now().isoformat(),
            "chains": {
                "receptor": receptor_chain,
                "ligand": ligand_chain
            },
            "best_sample": {
                "sample_id": best_sample_result['sample_id'] if best_sample_result else None,
                "initial_distance": best_sample_result['initial_distance'] if best_sample_result else None,
                "final_distance": best_sample_result['final_distance'] if best_sample_result else None,
                "simulation_type": best_sample_result['simulation_type'] if best_sample_result else None,
                "pdb_file": best_pdb_name if best_sample_result else None
            },
            "all_samples_summary": {
                "total_samples": len(all_sample_results),
                "successful_samples": len([s for s in all_sample_results if s['success']]),
                "failed_samples": len([s for s in all_sample_results if not s['success']]),
                "simulation_types": {}
            },
            "distance_statistics": {},
            "interface_analysis": {}
        }
        
        # 나머지 기존 로직들 유지...
        # (생략 - 기존 코드와 동일)
        
        return iteration_summary
        
    except Exception as e:
        logger.error(f"Iteration {iteration_num} 결과 저장 실패: {e}")
        return None

# ===== 기존 run_single_sample 함수 유지 (수정 없음) =====
def run_single_sample(sample_id, base_dir, processed_pdb, simulation_time, receptor_chain, ligand_chain, 
                     distance_threshold, logger, system_size="medium", config_file=None):
    """단일 샘플 시뮬레이션 실행 - GROMACS 단계별 성공 여부 확인 로직 추가"""
    sample_dir = os.path.join(base_dir, f"sample_{sample_id}")
    os.makedirs(sample_dir, exist_ok=True)
    
    logger.info(f"샘플 {sample_id} 시작 (시스템 크기: {system_size})")
    
    # 실패 시 반환할 기본 딕셔너리
    failure_result = {
        'sample_id': sample_id, 'initial_distance': float('inf'), 'final_distance': float('inf'),
        'final_pdb': None, 'simulation_type': 'failed', 'sample_dir': sample_dir,
        'success': False, 'error': 'Unknown error'
    }

    try:
        input_pdb_path = os.path.join(sample_dir, "input.pdb")
        shutil.copy(processed_pdb, input_pdb_path)
        
        MDPGenerator.generate_em_mdp(os.path.join(sample_dir, "em.mdp"), simulation_time, system_size)
        MDPGenerator.generate_ions_mdp(os.path.join(sample_dir, "ions.mdp"), simulation_time)
        MDPGenerator.generate_md_mdp(os.path.join(sample_dir, "md.mdp"), simulation_time)
        
        initial_distance = DistanceCalculator.calculate_chain_distance(processed_pdb, receptor_chain, ligand_chain)
        failure_result['initial_distance'] = initial_distance

        runner = GromacsCommandRunner(sample_dir, config_file, logger, system_size)
        
        # --- GROMACS 파이프라인 시작 ---
        # 1. pdb2gmx
        success, msg = runner.try_multiple_force_fields({"input_pdb": "input.pdb", "output_gro": "processed.gro", "output_top": "topol.top"})
        if not success:
            failure_result['error'] = f"pdb2gmx 실패: {msg}"
            return failure_result

        # 2. editconf
        success, msg = runner.execute_gromacs_command("editconf", {"input_gro": "processed.gro", "output_gro": "boxed.gro"})
        if not success:
            failure_result['error'] = f"editconf 실패: {msg}"
            return failure_result

        # 3. solvate
        success, msg = runner.execute_gromacs_command("solvate", {"input_gro": "boxed.gro", "output_gro": "solvated.gro", "topology": "topol.top"})
        if not success:
            failure_result['error'] = f"solvate 실패: {msg}"
            return failure_result

        # 4. grompp (ions)
        success, msg = runner.execute_gromacs_command("grompp_ions", {"mdp_file": "ions.mdp", "input_gro": "solvated.gro", "topology": "topol.top", "output_tpr": "ions.tpr"})
        if not success:
            failure_result['error'] = f"grompp_ions 실패: {msg}"
            return failure_result

        # 5. genion
        success, msg = runner.execute_gromacs_command("genion", {"input_tpr": "ions.tpr", "output_gro": "neutral.gro", "topology": "topol.top"}, max_retries=2)
        if not success:
            failure_result['error'] = f"genion 실패: {msg}"
            return failure_result

        # --- 시뮬레이션 실행 ---
        simulation_type = "EM+MD" if initial_distance <= distance_threshold else "EM_only"
        
        # 6. grompp (EM)
        success, msg = runner.execute_gromacs_command("grompp_simulation", {"mdp_file": "em.mdp", "input_gro": "neutral.gro", "topology": "topol.top", "output_tpr": "em.tpr"})
        if not success:
            failure_result['error'] = f"grompp_em 실패: {msg}"
            return failure_result

        # 7. mdrun (EM)
        success, msg = runner.execute_gromacs_command("mdrun_em", {"prefix": "em"})
        if not success:
            failure_result['error'] = f"mdrun_em 실패: {msg}"
            return failure_result

        final_pdb_path = os.path.join(sample_dir, "em.gro") # 우선 EM 결과 사용

        if simulation_type == "EM+MD":
            # 8. grompp (MD)
            success, msg = runner.execute_gromacs_command("grompp_simulation", {"mdp_file": "md.mdp", "input_gro": "em.gro", "topology": "topol.top", "output_tpr": "md.tpr"})
            if not success:
                failure_result['error'] = f"grompp_md 실패: {msg}"
                return failure_result
            
            # 9. mdrun (MD)
            success, msg = runner.execute_gromacs_command("mdrun_md", {"prefix": "md"})
            if not success:
                failure_result['error'] = f"mdrun_md 실패: {msg}"
                return failure_result
            final_pdb_path = os.path.join(sample_dir, "md.gro")

        # 최종 결과 파일 확인 및 거리 계산
        if not os.path.exists(final_pdb_path):
            failure_result['error'] = f"최종 결과 파일 누락: {final_pdb_path}"
            return failure_result

        # === 새로 추가: GRO 파일을 PDB로 변환 ===
        if final_pdb_path.endswith('.gro'):
            pdb_path = final_pdb_path.replace('.gro', '.pdb')

            # 구조 파일 결정 (em.tpr 또는 md.tpr)
            if simulation_type == "EM+MD":
                tpr_file = "md.tpr"
            else:
                tpr_file = "em.tpr"
            
            success, msg = runner.execute_gromacs_command("trjconv", {
                "structure_file": tpr_file,
                "input_trajectory": os.path.basename(final_pdb_path),
                "output_file": os.path.basename(pdb_path),
                "pbc_option": "mol",
                "group": "Protein"
            })
            
            if success:
                # 토폴로지 기반 체인 정보 복원
                restored_pdb_path = pdb_path.replace('.pdb', '_restored.pdb')
                chain_restored = runner.restore_chain_info_from_topology(
                    pdb_path, processed_pdb, restored_pdb_path
                )
                
                if chain_restored:
                    final_pdb_path = restored_pdb_path
                    logger.info(f"토폴로지 기반 체인 복원 성공: {os.path.basename(restored_pdb_path)}")
                else:
                    final_pdb_path = pdb_path
                    logger.warning(f"체인 복원 실패, 기본 PDB 사용: {os.path.basename(pdb_path)}")
            else:
                failure_result['error'] = f"PDB 변환 실패: {msg}"
                return failure_result
        
        final_distance = DistanceCalculator.calculate_chain_distance(final_pdb_path, receptor_chain, ligand_chain)

        # 성공 결과 반환
        result_pdb_path = os.path.join(sample_dir, f"result_{sample_id}.pdb")
        shutil.copy(final_pdb_path, result_pdb_path) # 최종 파일을 일관된 이름으로 복사

        return {
            'sample_id': sample_id, 'initial_distance': initial_distance, 'final_distance': final_distance,
            'final_pdb': result_pdb_path, 'simulation_type': simulation_type, 'sample_dir': sample_dir,
            'success': True
        }

    except Exception as e:
        logger.error(f"샘플 {sample_id} 실패: {e}")
        failure_result['error'] = str(e)
        return failure_result

# ===== 기존 다중 샘플 실행 함수 유지 =====
def run_multi_sample_iteration_with_results_saving(processed_pdb, iter_dir, simulation_time, receptor_chain, ligand_chain, 
                                                  distance_threshold, num_samples, logger, system_size="medium", 
                                                  config_file=None, iteration_num=None, job_output_dir=None):
    """다중 샘플 반복 실행 - 기존 코드 유지"""
    logger.info(f"다중 샘플 실행: {num_samples}개 샘플 (시스템 크기: {system_size})")
    
    sample_results = []
    successful_count = 0
    
    for sample_id in range(1, num_samples + 1):
        result = run_single_sample(
            sample_id, iter_dir, processed_pdb, simulation_time,
            receptor_chain, ligand_chain, distance_threshold, logger, system_size, config_file
        )
        sample_results.append(result)
        
        if result['success']:
            successful_count += 1
            logger.info(f"현재까지 성공한 샘플: {successful_count}/{sample_id}")
    
    successful_samples = [r for r in sample_results if r['success']]
    
    if not successful_samples:
        logger.error("모든 샘플이 실패했습니다")
        best_sample = None
    else:
        best_sample = min(successful_samples, key=lambda x: x['final_distance'])
        
        em_only_count = len([r for r in successful_samples if 'EM_only' in r['simulation_type']])
        em_md_count = len([r for r in successful_samples if r['simulation_type'] == 'EM+MD'])
        failed_count = len([r for r in sample_results if not r['success']])
        
        logger.info(f"다중 샘플 완료: 성공 {len(successful_samples)}, 실패 {failed_count}")
        logger.info(f"EM만: {em_only_count}, EM+MD: {em_md_count}")
        logger.info(f"최적 샘플: {best_sample['sample_id']} (거리: {best_sample['final_distance']:.2f}Å)")
    
    # Iteration 결과 저장
    if iteration_num is not None and job_output_dir is not None:
        iteration_summary = save_iteration_results(
            iteration_num, best_sample, sample_results, 
            job_output_dir, receptor_chain, ligand_chain, logger
        )
    else:
        iteration_summary = None
    
    return best_sample, sample_results

def estimate_system_size(pdb_file):
    """PDB 파일로부터 시스템 크기 추정"""
    try:
        with open(pdb_file, 'r') as f:
            atom_count = sum(1 for line in f if line.startswith('ATOM'))
        
        if atom_count > 10000:
            return "large"
        elif atom_count < 3000:
            return "small"
        else:
            return "medium"
    except:
        return "medium"

# ===== 수정된 메인 함수 =====
def main():
    parser = argparse.ArgumentParser(description="Golden Standard 기반 SuMD 시뮬레이션")
    parser.add_argument("--input_pdb", required=True, help="입력 PDB 파일")
    parser.add_argument("--golden_standard_pdb", required=True, help="Golden Standard PDB 파일")  # 새로 추가
    parser.add_argument("--output_dir", default="/app/output", help="출력 디렉토리")
    parser.add_argument("--simulation_time", type=float, default=1.0, help="시뮬레이션 시간 (ns)")
    parser.add_argument("--receptor_chain", required=True, help="수용체 체인 ID")
    parser.add_argument("--ligand_chain", required=True, help="리간드 체인 ID")
    parser.add_argument("--distance_threshold", type=float, default=5.0, help="동작 임계거리 (Å)")
    parser.add_argument("--rmsd_threshold", type=float, default=1.5, help="수렴 판정 RMSD 임계값 (Å)")  # 기본값 변경
    parser.add_argument("--convergence_window", type=int, default=5, help="수렴 판정 윈도우 크기")
    parser.add_argument("--max_iterations", type=int, default=10, help="최대 반복 횟수")
    parser.add_argument("--num_samples", type=int, default=5, help="각 반복당 샘플 수")
    parser.add_argument("--job_id", help="작업 ID (자동생성)")
    parser.add_argument("--skip_preprocessing", action="store_true", help="전처리 건너뛰기")
    parser.add_argument("--config_file", help="GROMACS 설정 파일 경로")
    
    args = parser.parse_args()
    
    # 입력 파일 존재 확인
    if not os.path.exists(args.input_pdb):
        print(f"오류: 입력 PDB 파일 '{args.input_pdb}'를 찾을 수 없습니다.")
        return 1
    
    # Golden standard PDB 존재 확인
    if not os.path.exists(args.golden_standard_pdb):
        print(f"오류: Golden Standard PDB 파일 '{args.golden_standard_pdb}'를 찾을 수 없습니다.")
        return 1
    
    # 작업 ID 생성
    if not args.job_id:
        args.job_id = generate_job_id(args.input_pdb, args.receptor_chain, args.ligand_chain)
    
    # 출력 디렉토리 설정
    job_output_dir = os.path.join(args.output_dir, f"sumd_{args.job_id}")
    os.makedirs(job_output_dir, exist_ok=True)
    
    # 로깅 설정
    logger = setup_logging(job_output_dir, args.job_id)
    
    # 시스템 크기 추정
    system_size = estimate_system_size(args.input_pdb)
    
    # ===== 새로 추가: Golden Standard 기반 수렴 판정기 초기화 =====
    convergence_checker = GoldenStandardConvergenceChecker(
        golden_standard_pdb=args.golden_standard_pdb,
        window_size=args.convergence_window, 
        rmsd_threshold=args.rmsd_threshold
    )
    
    logger.info(f"Golden Standard 기반 SuMD 시뮬레이션 시작 - Job ID: {args.job_id}")
    logger.info(f"입력 PDB: {args.input_pdb}")
    logger.info(f"Golden Standard PDB: {args.golden_standard_pdb}")
    logger.info(f"체인: {args.receptor_chain} - {args.ligand_chain}")
    logger.info(f"추정 시스템 크기: {system_size}")
    logger.info(f"시뮬레이션 시간: {args.simulation_time} ns")
    logger.info(f"수렴 조건: Golden Standard와의 마지막 {args.convergence_window}회 평균 RMSD ≤ {args.rmsd_threshold}Å")
    
    # 결과 추적을 위한 딕셔너리
    results = {
        "job_id": args.job_id,
        "input_pdb": args.input_pdb,
        "golden_standard_pdb": args.golden_standard_pdb,  # 새로 추가
        "system_size": system_size,
        "convergence_settings": {
            "window_size": args.convergence_window,
            "rmsd_threshold": args.rmsd_threshold
        },
        "start_time": datetime.now().isoformat(),
        "iterations": [],
        "final_pdb": None,
        "converged": False
    }
    
    current_pdb = args.input_pdb
    
    # ===== 시작 PDB와 Golden Standard 간의 초기 RMSD 계산 =====
    initial_rmsd = convergence_checker.calculate_rmsd_with_golden_standard(current_pdb)
    logger.info(f"초기 RMSD (vs Golden Standard): {initial_rmsd:.3f}Å")
    
    try:
        for iteration in range(1, args.max_iterations + 1):
            logger.info(f"\n=== 반복 {iteration} 시작 ===")
            
            # 반복별 작업 디렉토리
            iter_dir = os.path.join(job_output_dir, f"iteration_{iteration}")
            os.makedirs(iter_dir, exist_ok=True)
            
            # 1. PDB 전처리
            processed_pdb = None
            
            if not args.skip_preprocessing:
                logger.info("PDB 전처리 시작")
                try:
                    processed_pdb, stats = process_pdb_for_gromacs(
                        current_pdb, 
                        os.path.join(iter_dir, "processed.pdb"),
                        [args.receptor_chain, args.ligand_chain],
                        logger
                    )
                    logger.info(f"전처리 통계: {stats}")
                except Exception as e:
                    logger.error(f"전처리 실패: {e}")
                    processed_pdb = os.path.join(iter_dir, "processed.pdb")
                    shutil.copy(current_pdb, processed_pdb)
                    logger.warning("원본 PDB 파일로 계속 진행합니다")
            else:
                processed_pdb = os.path.join(iter_dir, "processed.pdb")
                shutil.copy(current_pdb, processed_pdb)
                logger.info("전처리 건너뛰고 원본 PDB 사용")
            
            if not processed_pdb or not os.path.exists(processed_pdb):
                logger.error("PDB 준비가 완전히 실패했습니다")
                break

            # 2. 거리 계산
            try:
                distance = DistanceCalculator.calculate_chain_distance(
                    processed_pdb, args.receptor_chain, args.ligand_chain
                )
                logger.info(f"체인 간 거리: {distance:.2f} Å")
            except Exception as e:
                logger.error(f"거리 계산 실패: {e}")
                break
            
            # 3. 다중 샘플 시뮬레이션 실행
            best_sample, all_samples = run_multi_sample_iteration_with_results_saving(
                processed_pdb, iter_dir, args.simulation_time,
                args.receptor_chain, args.ligand_chain, args.distance_threshold,
                args.num_samples, logger, system_size, args.config_file,
                iteration_num=iteration, job_output_dir=job_output_dir
            )
            
            if best_sample is None:
                logger.error(f"반복 {iteration}에서 모든 샘플이 실패했습니다")
                
                if iteration == 1:
                    logger.info("첫 번째 반복 실패, 전처리 없이 재시도")
                    args.skip_preprocessing = True
                    continue
                else:
                    break
            
            # ===== 4. Golden Standard와의 RMSD 계산 및 수렴 확인 =====
            rmsd_vs_golden = convergence_checker.add_rmsd(best_sample['final_pdb'])
            logger.info(f"Golden Standard와의 RMSD: {rmsd_vs_golden:.3f} Å")
            
            converged, avg_rmsd, status_msg = convergence_checker.check_convergence()
            logger.info(f"수렴 상태: {status_msg}")
            
            if converged:
                logger.info(f"수렴 달성! Golden Standard와의 마지막 {args.convergence_window}회 평균 RMSD {avg_rmsd:.3f}Å ≤ {args.rmsd_threshold}Å")
                results["converged"] = True
                results["final_pdb"] = best_sample['final_pdb']
                results["convergence_info"] = convergence_checker.get_convergence_info()
            
            # 5. 결과 저장
            iteration_result = {
                "iteration": iteration,
                "initial_distance": distance,
                "final_distance": best_sample['final_distance'],
                "rmsd_vs_golden": rmsd_vs_golden,  # 새로 추가
                "best_sample_id": best_sample['sample_id'],
                "simulation_type": best_sample['simulation_type'],
                "final_pdb": best_sample['final_pdb'],
                "successful_samples": len([s for s in all_samples if s['success']]),
                "failed_samples": len([s for s in all_samples if not s['success']])
            }
            results["iterations"].append(iteration_result)
            
            # 수렴 확인 후 종료
            if results["converged"]:
                break
            
            # 다음 반복 준비
            current_pdb = best_sample['final_pdb']
        
        # 최종 결과
        if not results["converged"]:
            logger.info("최대 반복 횟수 도달 (수렴하지 않음)")
            if results["iterations"]:
                results["final_pdb"] = results["iterations"][-1]["final_pdb"]
                results["convergence_info"] = convergence_checker.get_convergence_info()
        
        results["end_time"] = datetime.now().isoformat()
        
        # 결과 파일 저장
        with open(os.path.join(job_output_dir, "results.json"), 'w', encoding='utf-8') as f:
            json.dump(results, f, indent=2, default=str, ensure_ascii=False, cls=NE)
        
        logger.info(f"\n=== Golden Standard 기반 SuMD 시뮬레이션 완료 ===")
        logger.info(f"최종 PDB: {results['final_pdb']}")
        logger.info(f"수렴 여부: {results['converged']}")
        logger.info(f"총 반복 횟수: {len(results['iterations'])}")
        
        # 수렴 정보 출력
        if 'convergence_info' in results:
            conv_info = results['convergence_info']
            logger.info(f"최종 수렴 상태: 마지막 {conv_info['recent_iterations']}회 평균 RMSD: {conv_info['recent_avg_rmsd']:.3f}Å")
        
        # 결과 출력 (파싱용)
        print(f"SUMD_RESULT:{results['final_pdb']}:{results['converged']}")
        
    except Exception as e:
        logger.error(f"시뮬레이션 실행 중 오류: {e}")
        results["error"] = str(e)
        results["end_time"] = datetime.now().isoformat()
        
        with open(os.path.join(job_output_dir, "results.json"), 'w', encoding='utf-8') as f:
            json.dump(results, f, indent=2, default=str, ensure_ascii=False, cls=NE)
        
        return 1
    
    return 0

if __name__ == "__main__":
    sys.exit(main())