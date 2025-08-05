#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
간단한 Chain 추출기

외부 의존성 없이 BioPython만으로 체인 추출
/app/SUMD_automation/pdb_preprocessing/ 디렉토리에 위치
"""

import os
import sys
from Bio.PDB import PDBParser, PDBIO, Select
from typing import List, Tuple, Dict
import logging

# 상위 디렉토리의 모듈들에 접근 (필요시)
script_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(script_dir)  # /app/SUMD_automation
sys.path.insert(0, parent_dir)

class SimpleChainExtractor:
    """간단한 체인 추출 클래스"""
    
    def __init__(self, logger=None):
        self.logger = logger or logging.getLogger(__name__)
        self.parser = PDBParser(QUIET=True)
    
    def extract_chains(self, input_pdb: str, output_pdb: str, 
                  chains_to_extract: List[str]) -> Tuple[bool, Dict]:
        """
        지정된 체인들만 추출
        
        Args:
            input_pdb: 입력 PDB 파일
            output_pdb: 출력 PDB 파일 (이미 고유 이름으로 전달됨)
            chains_to_extract: 추출할 체인 ID 리스트
            
        Returns:
            tuple: (성공 여부, 통계 정보)
        """
        
        try:
            # PDB 구조 로드
            structure = self.parser.get_structure("structure", input_pdb)
            
            # 사용 가능한 체인 확인
            available_chains = []
            chain_info = {}
            
            for model in structure:
                for chain in model:
                    chain_id = chain.id
                    if chain_id not in available_chains:
                        available_chains.append(chain_id)
                        
                        # 체인 정보 수집
                        residues = list(chain.get_residues())
                        atoms = list(chain.get_atoms())
                        
                        chain_info[chain_id] = {
                            'residue_count': len(residues),
                            'atom_count': len(atoms)
                        }
            
            self.logger.info(f"사용 가능한 체인: {available_chains}")
            self.logger.info(f"추출 대상 체인: {chains_to_extract}")
            for chain_id, info in chain_info.items():
                self.logger.info(f"  체인 {chain_id}: {info['residue_count']} 잔기, {info['atom_count']} 원자")
            
            # 요청된 체인이 있는지 확인
            missing_chains = [c for c in chains_to_extract if c not in available_chains]
            if missing_chains:
                error_msg = f"누락된 체인: {missing_chains}, 사용 가능한 체인: {available_chains}"
                self.logger.error(error_msg)
                return False, {'error': error_msg}
            
            # 체인 선택기 클래스
            class ChainSelector(Select):
                def __init__(self, target_chains, logger):
                    self.target_chains = set(target_chains)
                    self.logger = logger
                    self.extracted_chains = set()
                    self.extracted_residues = 0
                    self.extracted_atoms = 0
                    self.skipped_residues = 0
                
                def accept_chain(self, chain):
                    accepted = chain.id in self.target_chains
                    if accepted:
                        self.extracted_chains.add(chain.id)
                        self.logger.debug(f"체인 {chain.id} 선택됨")
                    return accepted
                
                def accept_residue(self, residue):
                    # 기본적인 잔기 필터링
                    resname = residue.get_resname().strip()
                    
                    # 표준 아미노산만 허용 (기본 필터링)
                    standard_aa = {
                        'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY',
                        'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER',
                        'THR', 'TRP', 'TYR', 'VAL'
                    }
                    
                    if resname not in standard_aa:
                        self.skipped_residues += 1
                        self.logger.debug(f"비표준 잔기 제외: {residue.parent.id}:{residue.id[1]}{resname}")
                        return False
                    
                    # 백본 원자 확인
                    present_atoms = {atom.get_name().strip() for atom in residue}
                    backbone = {'N', 'CA', 'C', 'O'}
                    
                    if not backbone.issubset(present_atoms):
                        missing_backbone = backbone - present_atoms
                        self.skipped_residues += 1
                        self.logger.debug(f"백본 원자 누락으로 제외: {residue.parent.id}:{residue.id[1]}{resname} (누락: {missing_backbone})")
                        return False
                    
                    self.extracted_residues += 1
                    return True
                
                def accept_atom(self, atom):
                    self.extracted_atoms += 1
                    return True
            
            # 추출 실행
            selector = ChainSelector(chains_to_extract, self.logger)
            io = PDBIO()
            io.set_structure(structure)
            io.save(output_pdb, selector)
            
            # 결과 검증
            if not os.path.exists(output_pdb):
                return False, {'error': 'output_file_not_created'}
            
            # 추출된 파일 검증
            try:
                extracted_structure = self.parser.get_structure("extracted", output_pdb)
                
                final_chain_info = {}
                total_residues = 0
                total_atoms = 0
                
                for model in extracted_structure:
                    for chain in model:
                        chain_id = chain.id
                        residues = list(chain.get_residues())
                        atoms = list(chain.get_atoms())
                        
                        final_chain_info[chain_id] = {
                            'residue_count': len(residues),
                            'atom_count': len(atoms)
                        }
                        
                        total_residues += len(residues)
                        total_atoms += len(atoms)
                
                # 성공 통계
                stats = {
                    'success': True,
                    'input_file': input_pdb,
                    'output_file': output_pdb,
                    'requested_chains': chains_to_extract,
                    'extracted_chains': list(selector.extracted_chains),
                    'total_residues': total_residues,
                    'total_atoms': total_atoms,
                    'skipped_residues': selector.skipped_residues,
                    'chain_details': final_chain_info
                }
                
                self.logger.info(f"체인 추출 완료:")
                self.logger.info(f"  요청: {chains_to_extract}")
                self.logger.info(f"  추출: {list(selector.extracted_chains)}")
                self.logger.info(f"  총 잔기: {total_residues} (제외: {selector.skipped_residues})")
                self.logger.info(f"  총 원자: {total_atoms}")
                
                for chain_id, info in final_chain_info.items():
                    self.logger.info(f"  체인 {chain_id}: {info['residue_count']} 잔기, {info['atom_count']} 원자")
                
                return True, stats
                
            except Exception as e:
                return False, {'error': f'validation_failed: {e}'}
            
        except Exception as e:
            self.logger.error(f"체인 추출 실패: {e}")
            return False, {'error': str(e)}
    
    def get_chain_info(self, pdb_file: str) -> Dict:
        """PDB 파일의 체인 정보 조회"""
        
        try:
            structure = self.parser.get_structure("structure", pdb_file)
            
            chain_info = {}
            
            for model in structure:
                for chain in model:
                    chain_id = chain.id
                    residues = list(chain.get_residues())
                    atoms = list(chain.get_atoms())
                    
                    # 잔기 종류별 통계
                    residue_types = {}
                    for residue in residues:
                        resname = residue.get_resname().strip()
                        if resname in residue_types:
                            residue_types[resname] += 1
                        else:
                            residue_types[resname] = 1
                    
                    chain_info[chain_id] = {
                        'residue_count': len(residues),
                        'atom_count': len(atoms),
                        'first_residue': residues[0].id[1] if residues else None,
                        'last_residue': residues[-1].id[1] if residues else None,
                        'residue_types': residue_types
                    }
            
            return chain_info
            
        except Exception as e:
            self.logger.error(f"체인 정보 조회 실패: {e}")
            return {}
    
    def validate_pdb_file(self, pdb_file: str) -> Dict:
        """PDB 파일 기본 검증"""
        
        validation = {
            'valid': False,
            'file_exists': False,
            'file_size': 0,
            'chains': [],
            'total_residues': 0,
            'total_atoms': 0,
            'errors': []
        }
        
        try:
            # 파일 존재 확인
            if not os.path.exists(pdb_file):
                validation['errors'].append("File does not exist")
                return validation
            
            validation['file_exists'] = True
            validation['file_size'] = os.path.getsize(pdb_file)
            
            if validation['file_size'] == 0:
                validation['errors'].append("File is empty")
                return validation
            
            # PDB 구조 로드
            structure = self.parser.get_structure("validation", pdb_file)
            
            chains = []
            total_residues = 0
            total_atoms = 0
            
            for model in structure:
                for chain in model:
                    chains.append(chain.id)
                    residues = list(chain.get_residues())
                    atoms = list(chain.get_atoms())
                    total_residues += len(residues)
                    total_atoms += len(atoms)
            
            validation['chains'] = list(set(chains))  # 중복 제거
            validation['total_residues'] = total_residues
            validation['total_atoms'] = total_atoms
            
            if not chains:
                validation['errors'].append("No chains found")
            elif total_residues == 0:
                validation['errors'].append("No residues found")
            elif total_atoms == 0:
                validation['errors'].append("No atoms found")
            else:
                validation['valid'] = True
            
            return validation
            
        except Exception as e:
            validation['errors'].append(f"Parsing error: {e}")
            return validation


def extract_receptor_ligand_chains(input_pdb: str, output_pdb: str, 
                                 receptor_chain: str, ligand_chain: str,
                                 logger=None) -> Tuple[bool, Dict]:
    """
    Receptor와 Ligand 체인 추출 (간단한 래퍼 함수)
    
    Args:
        input_pdb: 입력 PDB 파일
        output_pdb: 출력 PDB 파일  
        receptor_chain: 수용체 체인 ID
        ligand_chain: 리간드 체인 ID
        logger: 로거
        
    Returns:
        tuple: (성공 여부, 결과 정보)
    """
    
    extractor = SimpleChainExtractor(logger)
    chains_to_extract = [receptor_chain, ligand_chain]
    
    success, result = extractor.extract_chains(input_pdb, output_pdb, chains_to_extract)
    
    if success:
        # 추가 정보 포함
        result['main_receptor'] = receptor_chain
        result['main_ligand'] = ligand_chain
        result['extraction_method'] = 'SimpleChainExtractor'
    
    return success, result


def main():
    """테스트용 메인 함수"""
    import argparse
    
    parser = argparse.ArgumentParser(description="간단한 체인 추출기")
    parser.add_argument("--input", required=True, help="입력 PDB 파일")
    parser.add_argument("--output", required=True, help="출력 PDB 파일")
    parser.add_argument("--chains", required=True, help="추출할 체인 (쉼표로 구분, 예: A,B)")
    parser.add_argument("--info", action="store_true", help="체인 정보만 출력")
    parser.add_argument("--validate", action="store_true", help="PDB 파일 검증")
    parser.add_argument("--verbose", "-v", action="store_true", help="상세 로그")
    
    args = parser.parse_args()
    
    # 로깅 설정
    level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(level=level, format='%(levelname)s: %(message)s')
    logger = logging.getLogger(__name__)
    
    if not os.path.exists(args.input):
        print(f"❌ 입력 파일을 찾을 수 없습니다: {args.input}")
        return 1
    
    extractor = SimpleChainExtractor(logger)
    
    # PDB 파일 검증
    if args.validate:
        print("PDB 파일 검증 중...")
        validation = extractor.validate_pdb_file(args.input)
        
        print(f"파일 유효성: {'✅' if validation['valid'] else '❌'}")
        print(f"파일 크기: {validation['file_size']} bytes")
        print(f"체인: {validation['chains']}")
        print(f"총 잔기: {validation['total_residues']}")
        print(f"총 원자: {validation['total_atoms']}")
        
        if validation['errors']:
            print(f"오류: {validation['errors']}")
        
        if not validation['valid']:
            return 1
    
    # 체인 정보 출력
    if args.info:
        print("체인 정보 조회 중...")
        chain_info = extractor.get_chain_info(args.input)
        
        for chain_id, info in chain_info.items():
            print(f"\n체인 {chain_id}:")
            print(f"  잔기 수: {info['residue_count']}")
            print(f"  원자 수: {info['atom_count']}")
            print(f"  잔기 범위: {info['first_residue']}-{info['last_residue']}")
            
            if info['residue_types']:
                print(f"  잔기 종류: {dict(list(info['residue_types'].items())[:5])}")
        
        return 0
    
    # 체인 추출
    chains = [c.strip() for c in args.chains.split(',')]
    print(f"체인 추출: {chains}")
    
    success, result = extractor.extract_chains(args.input, args.output, chains)
    
    if success:
        print(f"✅ 체인 추출 성공: {args.output}")
        print(f"추출된 체인: {result['extracted_chains']}")
        print(f"총 잔기: {result['total_residues']}")
        print(f"총 원자: {result['total_atoms']}")
        return 0
    else:
        print(f"❌ 체인 추출 실패: {result['error']}")
        return 1


if __name__ == "__main__":
    exit(main())
