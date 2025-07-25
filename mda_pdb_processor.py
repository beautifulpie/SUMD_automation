#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import numpy as np
import logging
from Bio.PDB import PDBParser, PDBIO, Select
from Bio.PDB.Atom import Atom
from typing import List, Optional, Tuple

class AdvancedPDBProcessor:
    """MDAnalysis를 중심으로 한 고급 PDB 전처리 클래스"""
    
    def __init__(self, logger=None):
        self.logger = logger or logging.getLogger(__name__)
        self.standard_amino_acids = {
            'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY',
            'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER',
            'THR', 'TRP', 'TYR', 'VAL'
        }
        
        # 각 아미노산의 필수 백본 원자
        self.required_backbone = {'N', 'CA', 'C', 'O'}
        
        # 중요한 사이드체인 원자들 (누락되면 문제가 되는 것들)
        self.critical_sidechain = {
            'GLN': {'CB', 'CG'},
            'GLU': {'CB', 'CG'},
            'ARG': {'CB', 'CG', 'CD'},
            'LYS': {'CB', 'CG', 'CD'},
            'ASN': {'CB', 'CG'},
            'ASP': {'CB', 'CG'},
            'CYS': {'CB', 'SG'},
            'MET': {'CB', 'CG', 'SD'},
            'HIS': {'CB', 'CG'},
            'TRP': {'CB', 'CG'},
            'PHE': {'CB', 'CG'},
            'TYR': {'CB', 'CG'},
            'LEU': {'CB', 'CG'},
            'ILE': {'CB', 'CG1'},
            'VAL': {'CB', 'CG1'},
            'THR': {'CB', 'OG1'},
            'SER': {'CB', 'OG'},
            'PRO': {'CB', 'CG', 'CD'}
        }
    
    def process_pdb_pipeline(self, input_pdb: str, output_pdb: str, 
                           target_chains: List[str]) -> Tuple[str, dict]:
        """
        완전한 PDB 전처리 파이프라인
        
        Args:
            input_pdb: 입력 PDB 파일 경로
            output_pdb: 출력 PDB 파일 경로  
            target_chains: 대상 체인 리스트 (예: ['A', 'C'])
            
        Returns:
            tuple: (처리된 파일 경로, 처리 통계)
        """
        self.logger.info(f"=== PDB 전처리 파이프라인 시작 ===")
        self.logger.info(f"입력 파일: {input_pdb}")
        self.logger.info(f"대상 체인: {target_chains}")
        
        stats = {
            'original_atoms': 0,
            'original_residues': 0,
            'final_atoms': 0,
            'final_residues': 0,
            'removed_residues': [],
            'processed_cysteines': 0,
            'processing_method': None
        }
        
        # 1단계: 기본 전처리 (MDAnalysis 우선)
        temp_pdb1 = f"{output_pdb}.temp1"
        processed_pdb = self._process_with_mdanalysis(input_pdb, temp_pdb1, target_chains, stats)
        
        if not processed_pdb:
            self.logger.warning("MDAnalysis 실패, Bio.PDB로 fallback")
            processed_pdb = self._process_with_biopdb(input_pdb, temp_pdb1, target_chains, stats)
        
        # 2단계: 리간드 제거
        temp_pdb2 = f"{output_pdb}.temp2"
        self._remove_ligands(processed_pdb, temp_pdb2, stats)
        
        # 3단계: 시스테인 정비
        temp_pdb3 = f"{output_pdb}.temp3"
        self._process_cysteines(temp_pdb2, temp_pdb3, stats)
        
        # 4단계: 최종 검증 및 정리
        self._final_validation(temp_pdb3, output_pdb, stats)
        
        # 임시 파일 정리
        for temp_file in [temp_pdb1, temp_pdb2, temp_pdb3]:
            if os.path.exists(temp_file):
                os.remove(temp_file)
        
        self.logger.info(f"=== PDB 전처리 완료 ===")
        self.logger.info(f"최종 파일: {output_pdb}")
        self.logger.info(f"원자 수: {stats['original_atoms']} → {stats['final_atoms']}")
        self.logger.info(f"잔기 수: {stats['original_residues']} → {stats['final_residues']}")
        self.logger.info(f"제거된 잔기: {len(stats['removed_residues'])}개")
        self.logger.info(f"처리된 시스테인: {stats['processed_cysteines']}개")
        
        return output_pdb, stats
    
    def _process_with_mdanalysis(self, input_pdb: str, output_pdb: str, 
                               target_chains: List[str], stats: dict) -> Optional[str]:
        """MDAnalysis를 사용한 PDB 전처리"""
        try:
            import MDAnalysis as mda
            
            self.logger.info("MDAnalysis로 PDB 전처리 시작")
            
            # Universe 생성
            u = mda.Universe(input_pdb)
            stats['original_atoms'] = len(u.atoms)
            stats['original_residues'] = len(u.residues)
            
            # 사용 가능한 체인 확인
            available_chains = np.unique([atom.chainID for atom in u.atoms if atom.chainID.strip()])
            self.logger.info(f"사용 가능한 체인: {list(available_chains)}")
            
            # 체인 유효성 검증
            missing_chains = [chain for chain in target_chains if chain not in available_chains]
            if missing_chains:
                raise ValueError(f"지정된 체인 {missing_chains}이 PDB에 없습니다. 사용 가능한 체인: {list(available_chains)}")
            
            # 체인 선택
            chain_selection = " or ".join([f"chainID {chain}" for chain in target_chains])
            protein_atoms = u.select_atoms(f"protein and ({chain_selection})")
            
            if len(protein_atoms) == 0:
                raise ValueError("선택된 단백질 원자가 없습니다")
            
            self.logger.info(f"선택된 원자 수: {len(protein_atoms)}")
            
            # 잔기별 완전성 검사
            complete_residues = []
            incomplete_residues = []
            
            residue_count = 0
            for residue in protein_atoms.residues:
                residue_count += 1
                resname = residue.resname
                
                # 디버깅 정보는 debug 레벨로 변경하고 주기적으로만 출력
                if residue_count % 50 == 0 or residue_count <= 10:  # 처음 10개와 50개마다 출력
                    self.logger.debug(f"처리 중: {residue.chainID}:{residue.resid}{resname} ({residue_count}번째)")

                # 표준 아미노산 확인
                if resname not in self.standard_amino_acids:
                    incomplete_residues.append(f"{residue.chainID}:{residue.resid}{resname} (non-standard)")
                    continue
                
                # 백본 원자 확인
                backbone_atoms = residue.atoms.select_atoms("name N C CA O")
                if len(backbone_atoms) < 4:
                    incomplete_residues.append(f"{residue.chainID}:{residue.resid}{resname} (missing backbone)")
                    continue
                
                # 중요한 사이드체인 원자 확인
                if resname in self.critical_sidechain:
                    critical_atoms = self.critical_sidechain[resname]
                    present_atoms = set(residue.atoms.names)
                    missing_critical = critical_atoms - present_atoms
                    
                    if missing_critical:
                        incomplete_residues.append(f"{residue.chainID}:{residue.resid}{resname} (missing {missing_critical})")
                        continue
                
                complete_residues.append(residue.resid)
            
            stats['removed_residues'] = incomplete_residues
            self.logger.info(f"완전한 잔기: {len(complete_residues)}개")
            self.logger.info(f"제거된 잔기: {len(incomplete_residues)}개")
            
            # 제거된 잔기가 많은 경우 일부만 로그로 출력
            if len(incomplete_residues) > 0:
                if len(incomplete_residues) <= 10:
                    self.logger.warning(f"제거된 잔기: {incomplete_residues}")
                else:
                    self.logger.warning(f"제거된 잔기 (처음 5개): {incomplete_residues[:5]}")
                    self.logger.warning(f"제거된 잔기 (마지막 5개): {incomplete_residues[-5:]}")
            
            if not complete_residues:
                raise ValueError("완전한 잔기가 없습니다")
            
            # 완전한 잔기만 선택
            resid_selection = " or ".join([f"resid {resid}" for resid in complete_residues])
            final_selection = u.select_atoms(f"protein and ({chain_selection}) and ({resid_selection})")
            
            # PDB 저장
            final_selection.write(output_pdb)
            
            stats['processing_method'] = 'MDAnalysis'
            self.logger.info(f"MDAnalysis 전처리 완료: {len(final_selection)} 원자")
            
            return output_pdb
            
        except ImportError:
            self.logger.warning("MDAnalysis가 설치되지 않음")
            return None
        except Exception as e:
            self.logger.error(f"MDAnalysis 처리 실패: {e}")
            return None
    
    def _process_with_biopdb(self, input_pdb: str, output_pdb: str, 
                           target_chains: List[str], stats: dict) -> str:
        """Bio.PDB를 사용한 fallback 전처리"""
        self.logger.info("Bio.PDB로 PDB 전처리 시작")
        
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("structure", input_pdb)
        
        # 통계 수집
        total_atoms = sum(len(list(chain.get_atoms())) for model in structure for chain in model)
        total_residues = sum(len(list(chain.get_residues())) for model in structure for chain in model)
        stats['original_atoms'] = total_atoms
        stats['original_residues'] = total_residues
        
        # 사용 가능한 체인 확인
        available_chains = []
        for model in structure:
            for chain in model:
                if chain.id.strip():
                    available_chains.append(chain.id)
        
        unique_chains = sorted(set(available_chains))
        self.logger.info(f"사용 가능한 체인: {unique_chains}")
        
        # 체인 유효성 검증
        missing_chains = [chain for chain in target_chains if chain not in unique_chains]
        if missing_chains:
            raise ValueError(f"지정된 체인 {missing_chains}이 PDB에 없습니다. 사용 가능한 체인: {unique_chains}")
        
        # 선택기 클래스
        class StrictProteinSelect(Select):
            def __init__(self, target_chains, processor):
                self.target_chains = target_chains
                self.processor = processor
                self.rejected_residues = []
                self.processed_count = 0
            
            def accept_chain(self, chain):
                return chain.id in self.target_chains
            
            def accept_residue(self, residue):
                self.processed_count += 1
                resname = residue.get_resname()
                
                # 디버깅 정보는 주기적으로만 출력
                if self.processed_count % 50 == 0 or self.processed_count <= 10:
                    self.processor.logger.debug(f"Bio.PDB 처리 중: {residue.parent.id}:{residue.id[1]}{resname} ({self.processed_count}번째)")
                
                # 표준 아미노산 확인
                if resname not in self.processor.standard_amino_acids:
                    self.rejected_residues.append(f"{residue.parent.id}:{residue.id[1]}{resname} (non-standard)")
                    return False
                
                # 백본 원자 확인
                present_atoms = {atom.get_name() for atom in residue}
                if not self.processor.required_backbone.issubset(present_atoms):
                    self.rejected_residues.append(f"{residue.parent.id}:{residue.id[1]}{resname} (missing backbone)")
                    return False
                
                # 중요한 사이드체인 확인
                if resname in self.processor.critical_sidechain:
                    critical_atoms = self.processor.critical_sidechain[resname]
                    missing_critical = critical_atoms - present_atoms
                    if missing_critical:
                        self.rejected_residues.append(f"{residue.parent.id}:{residue.id[1]}{resname} (missing {missing_critical})")
                        return False
                
                return True
        
        # 선택 및 저장
        selector = StrictProteinSelect(target_chains, self)
        io = PDBIO()
        io.set_structure(structure)
        io.save(output_pdb, selector)
        
        stats['removed_residues'] = selector.rejected_residues
        stats['processing_method'] = 'Bio.PDB'
        
        self.logger.info(f"Bio.PDB 전처리 완료, 제거된 잔기: {len(selector.rejected_residues)}개")
        
        # 제거된 잔기가 많은 경우 일부만 로그로 출력
        if len(selector.rejected_residues) > 0:
            if len(selector.rejected_residues) <= 10:
                self.logger.warning(f"제거된 잔기: {selector.rejected_residues}")
            else:
                self.logger.warning(f"제거된 잔기 (처음 5개): {selector.rejected_residues[:5]}")
                self.logger.warning(f"제거된 잔기 (마지막 5개): {selector.rejected_residues[-5:]}")
        
        return output_pdb
    
    def _remove_ligands(self, input_pdb: str, output_pdb: str, stats: dict):
        """리간드 및 용매 분자 제거"""
        self.logger.info("리간드 제거 시작")
        
        try:
            import MDAnalysis as mda
            
            u = mda.Universe(input_pdb)
            
            # 단백질만 선택 (리간드, 용매, 이온 제외)
            protein_only = u.select_atoms("protein")
            
            if len(protein_only) == 0:
                raise ValueError("단백질 원자가 없습니다")
            
            # 제거된 원자 수 계산
            removed_atoms = len(u.atoms) - len(protein_only)
            
            protein_only.write(output_pdb)
            
            self.logger.info(f"리간드 제거 완료: {removed_atoms}개 원자 제거")
            
        except ImportError:
            # MDAnalysis가 없으면 Bio.PDB 사용
            self._remove_ligands_with_biopdb(input_pdb, output_pdb)
    
    def _remove_ligands_with_biopdb(self, input_pdb: str, output_pdb: str):
        """Bio.PDB를 사용한 리간드 제거"""
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("structure", input_pdb)
        
        class ProteinOnlySelect(Select):
            def accept_residue(self, residue):
                # 단백질 잔기만 허용
                return residue.get_resname() in {
                    'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY',
                    'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER',
                    'THR', 'TRP', 'TYR', 'VAL'
                }
        
        io = PDBIO()
        io.set_structure(structure)
        io.save(output_pdb, ProteinOnlySelect())
        
        self.logger.info("Bio.PDB로 리간드 제거 완료")
    
    def _process_cysteines(self, input_pdb: str, output_pdb: str, stats: dict):
        """시스테인 처리 - HG 원자 추가"""
        self.logger.info("시스테인 처리 시작")
        
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("structure", input_pdb)
        
        processed_count = 0
        
        for model in structure:
            for chain in model:
                for residue in chain:
                    if residue.get_resname() == "CYS":
                        # HG 원자가 없고 SG, CA 원자가 있는 경우
                        if "HG" not in residue and "SG" in residue and "CA" in residue:
                            sg_atom = residue["SG"]
                            ca_atom = residue["CA"]
                            
                            # SG에서 CA로의 벡터 계산
                            vector = sg_atom.coord - ca_atom.coord
                            norm = np.linalg.norm(vector)
                            
                            if norm > 0:
                                unit_vector = vector / norm
                                bond_length = 1.33  # S-H 결합 길이 (Å)
                                hg_coord = sg_atom.coord + unit_vector * bond_length
                                
                                # HG 원자 생성
                                hg_atom = Atom(
                                    name="HG", 
                                    coord=hg_coord,
                                    bfactor=sg_atom.get_bfactor(),
                                    occupancy=sg_atom.get_occupancy(),
                                    altloc=sg_atom.get_altloc(),
                                    fullname=" HG ",
                                    serial_number=0, 
                                    element="H"
                                )
                                
                                residue.add(hg_atom)
                                processed_count += 1
        
        # 구조 저장
        io = PDBIO()
        io.set_structure(structure)
        io.save(output_pdb)
        
        stats['processed_cysteines'] = processed_count
        self.logger.info(f"시스테인 처리 완료: {processed_count}개 시스테인에 HG 원자 추가")
    
    def _final_validation(self, input_pdb: str, output_pdb: str, stats: dict):
        """최종 검증 및 통계 업데이트"""
        self.logger.info("최종 검증 시작")
        
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("structure", input_pdb)
        
        # 최종 통계 계산
        final_atoms = sum(len(list(chain.get_atoms())) for model in structure for chain in model)
        final_residues = sum(len(list(chain.get_residues())) for model in structure for chain in model)
        
        stats['final_atoms'] = final_atoms
        stats['final_residues'] = final_residues
        
        # 체인별 통계
        chain_stats = {}
        for model in structure:
            for chain in model:
                chain_id = chain.id
                chain_atoms = len(list(chain.get_atoms()))
                chain_residues = len(list(chain.get_residues()))
                chain_stats[chain_id] = {
                    'atoms': chain_atoms,
                    'residues': chain_residues
                }
        
        self.logger.info("체인별 통계:")
        for chain_id, stat in chain_stats.items():
            self.logger.info(f"  체인 {chain_id}: {stat['residues']} 잔기, {stat['atoms']} 원자")
        
        # 최종 파일 복사
        import shutil
        shutil.copy(input_pdb, output_pdb)
        
        self.logger.info("최종 검증 완료")

# 사용 예시 함수
def process_pdb_for_gromacs(input_pdb: str, output_pdb: str, 
                          target_chains: List[str], logger=None) -> Tuple[str, dict]:
    """
    GROMACS 시뮬레이션을 위한 PDB 전처리
    
    Args:
        input_pdb: 입력 PDB 파일 경로
        output_pdb: 출력 PDB 파일 경로
        target_chains: 대상 체인 리스트 (예: ['A', 'C'])
        logger: 로거 객체 (선택사항)
    
    Returns:
        tuple: (처리된 파일 경로, 처리 통계)
    """
    processor = AdvancedPDBProcessor(logger=logger)
    return processor.process_pdb_pipeline(input_pdb, output_pdb, target_chains)

# 테스트 함수
def test_pdb_processing():
    """PDB 처리 테스트"""
    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger("PDB_Processor")
    
    input_pdb = "test_input.pdb"
    output_pdb = "test_processed.pdb"
    target_chains = ['A', 'C']
    
    try:
        result_pdb, stats = process_pdb_for_gromacs(
            input_pdb, output_pdb, target_chains, logger
        )
        print(f"\n=== 처리 완료 ===")
        print(f"결과 파일: {result_pdb}")
        print(f"처리 방법: {stats['processing_method']}")
        print(f"원자 수 변화: {stats['original_atoms']} → {stats['final_atoms']}")
        print(f"잔기 수 변화: {stats['original_residues']} → {stats['final_residues']}")
        print(f"제거된 잔기: {len(stats['removed_residues'])}개")
        print(f"처리된 시스테인: {stats['processed_cysteines']}개")
        
    except Exception as e:
        print(f"처리 실패: {e}")

if __name__ == "__main__":
    test_pdb_processing()