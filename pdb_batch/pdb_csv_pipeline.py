#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CSV 기반 PDB 처리 파이프라인

CSV 파일에서 PDB 정보를 읽어서:
1. PDB 파일 다운로드
2. Receptor/Ligand chain 추출 (simple_chain_extractor 사용)
3. PDB 수정 (GROMACS 호환성)
4. pdb2gmx 실행

파일 구조:
- 기존 파일들: /app/SUMD_automation/
- 새 파일들: /app/SUMD_automation/pdb_preprocessing/
"""

import os
import sys
import pandas as pd
import requests
import logging
import json
import time
from typing import Dict, List, Tuple, Optional
from pathlib import Path
import shutil

# 상위 디렉토리 (/app/SUMD_automation)의 모듈들에 접근
script_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(script_dir)  # /app/SUMD_automation
sys.path.insert(0, parent_dir)
sys.path.insert(0, script_dir)

# 기존 모듈들 import (상위 디렉토리에서)
try:
    from gromacs_runner_class import GromacsCommandRunner
    GROMACS_RUNNER_AVAILABLE = True
    print("✅ gromacs_runner_class 모듈 로드 성공")
except ImportError as e:
    GROMACS_RUNNER_AVAILABLE = False
    print(f"❌ gromacs_runner_class 모듈 로드 실패: {e}")

try:
    from interface_distance_calculator import DistanceCalculator
    DISTANCE_CALCULATOR_AVAILABLE = True
    print("✅ interface_distance_calculator 모듈 로드 성공")
except ImportError as e:
    DISTANCE_CALCULATOR_AVAILABLE = False
    print(f"❌ interface_distance_calculator 모듈 로드 실패: {e}")

try:
    from mda_pdb_processor import process_pdb_for_gromacs
    PDB_PROCESSOR_AVAILABLE = True
    print("✅ mda_pdb_processor 모듈 로드 성공")
except ImportError as e:
    PDB_PROCESSOR_AVAILABLE = False
    print(f"❌ mda_pdb_processor 모듈 로드 실패: {e}")

# 현재 디렉토리의 모듈들 import
try:
    from simple_chain_extractor import SimpleChainExtractor
    CHAIN_EXTRACTOR_AVAILABLE = True
    print("✅ simple_chain_extractor 모듈 로드 성공")
except ImportError as e:
    CHAIN_EXTRACTOR_AVAILABLE = False
    print(f"❌ simple_chain_extractor 모듈 로드 실패: {e}")


class PDBProcessingPipeline:
    """CSV 기반 PDB 처리 파이프라인"""
    
    def __init__(self, output_dir: str = "pdb_processing_results", 
                 config_file: str = None, logger=None):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)
        
        self.logger = logger or self._setup_logger()
        self.config_file = config_file
        
        # 통계 추적
        self.stats = {
            'total_entries': 0,
            'successful_downloads': 0,
            'successful_extractions': 0,
            'successful_pdb2gmx': 0,
            'failed_entries': [],
            'processing_times': {}
        }
        
        self.logger.info(f"PDB 처리 파이프라인 초기화: {self.output_dir}")
        self._log_module_availability()
    
    def _setup_logger(self) -> logging.Logger:
        """로거 설정"""
        logger = logging.getLogger("PDB_Pipeline")
        if not logger.handlers:
            handler = logging.StreamHandler()
            formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
            handler.setFormatter(formatter)
            logger.addHandler(handler)
            logger.setLevel(logging.INFO)
        return logger
    
    def _log_module_availability(self):
        """모듈 사용 가능 여부 로깅"""
        self.logger.info("=== 모듈 사용 가능 여부 ===")
        self.logger.info(f"GROMACS Runner: {'✅' if GROMACS_RUNNER_AVAILABLE else '❌'}")
        self.logger.info(f"Chain Extractor: {'✅' if CHAIN_EXTRACTOR_AVAILABLE else '❌'}")
        self.logger.info(f"PDB Processor: {'✅' if PDB_PROCESSOR_AVAILABLE else '❌'}")
        self.logger.info(f"Distance Calculator: {'✅' if DISTANCE_CALCULATOR_AVAILABLE else '❌'}")
    
    def process_csv(self, csv_file: str, max_entries: int = None) -> Dict:
        """CSV 파일 전체 처리"""
        self.logger.info(f"=== CSV 기반 PDB 처리 시작 ===")
        self.logger.info(f"CSV 파일: {csv_file}")
        
        start_time = time.time()
        
        try:
            # CSV 읽기
            df = pd.read_csv(csv_file)
            self.logger.info(f"CSV 로드 완료: {len(df)} 항목")
            
            # 필수 컬럼 확인
            required_columns = ['PDB', 'Peptide Chain', 'Receptor Chain']
            missing_columns = [col for col in required_columns if col not in df.columns]
            if missing_columns:
                raise ValueError(f"필수 컬럼 누락: {missing_columns}")
            
            # 처리할 항목 수 제한
            if max_entries:
                df = df.head(max_entries)
                self.logger.info(f"처리 항목 제한: {len(df)} 항목")
            
            self.stats['total_entries'] = len(df)
            
            # 각 항목 처리
            results = []
            
            for idx, row in df.iterrows():
                try:
                    self.logger.info(f"\n--- 처리 중 ({idx+1}/{len(df)}): {row['PDB']} ---")
                    result = self.process_single_entry(row)
                    results.append(result)
                    
                    if result['success']:
                        self.logger.info(f"✅ {row['PDB']} 처리 완료")
                    else:
                        self.logger.error(f"❌ {row['PDB']} 처리 실패: {result.get('error', 'Unknown')}")
                        
                except Exception as e:
                    self.logger.error(f"❌ {row['PDB']} 처리 중 예외: {e}")
                    results.append({
                        'pdb_id': row['PDB'],
                        'success': False,
                        'error': str(e)
                    })
            
            # 최종 통계
            end_time = time.time()
            total_time = end_time - start_time
            
            successful = len([r for r in results if r['success']])
            failed = len([r for r in results if not r['success']])
            
            self.logger.info(f"\n=== 처리 완료 ===")
            self.logger.info(f"총 항목: {len(df)}")
            self.logger.info(f"성공: {successful}")
            self.logger.info(f"실패: {failed}")
            self.logger.info(f"총 처리 시간: {total_time:.1f}초")
            
            # 결과 저장
            summary = {
                'total_entries': len(df),
                'successful': successful,
                'failed': failed,
                'total_time': total_time,
                'results': results,
                'stats': self.stats
            }
            
            summary_file = self.output_dir / "processing_summary.json"
            with open(summary_file, 'w') as f:
                json.dump(summary, f, indent=2, default=str)
            
            self.logger.info(f"처리 요약 저장: {summary_file}")
            
            return summary
            
        except Exception as e:
            self.logger.error(f"CSV 처리 실패: {e}")
            raise
    
    def process_single_entry(self, row: pd.Series) -> Dict:
        """단일 항목 처리"""
        pdb_id = row['PDB'].lower()
        peptide_chain = str(row['Peptide Chain']).strip()
        receptor_chain = str(row['Receptor Chain']).strip()
        
        entry_start_time = time.time()
        
        result = {
            'pdb_id': pdb_id,
            'peptide_chain': peptide_chain,
            'receptor_chain': receptor_chain,
            'success': False,
            'downloaded_pdb': None,
            'extracted_pdb': None,
            'processed_pdb': None,
            'pdb2gmx_success': False
        }
        
        try:
            # 1. PDB 다운로드
            self.logger.info(f"1. PDB 다운로드: {pdb_id}")
            downloaded_pdb = self.download_pdb(pdb_id)
            if not downloaded_pdb:
                result['error'] = "PDB 다운로드 실패"
                return result
            
            result['downloaded_pdb'] = str(downloaded_pdb)
            self.stats['successful_downloads'] += 1
            
            # 2. Chain 추출
            self.logger.info(f"2. Chain 추출: {receptor_chain} (receptor), {peptide_chain} (ligand)")
            extracted_pdb = self.extract_chains(downloaded_pdb, receptor_chain, peptide_chain)
            if not extracted_pdb:
                result['error'] = "Chain 추출 실패"
                return result
            
            result['extracted_pdb'] = str(extracted_pdb)
            self.stats['successful_extractions'] += 1
            
            # 3. PDB 전처리 (GROMACS 호환성)
            self.logger.info(f"3. PDB 전처리 (GROMACS 호환성)")
            processed_pdb = self.process_pdb_for_gromacs_compatibility(
                extracted_pdb, [receptor_chain, peptide_chain]
            )
            result['processed_pdb'] = str(processed_pdb)
            
            # 4. pdb2gmx 실행
            self.logger.info(f"4. pdb2gmx 실행")
            pdb2gmx_success = self.run_pdb2gmx(processed_pdb, pdb_id)
            result['pdb2gmx_success'] = pdb2gmx_success
            
            if pdb2gmx_success:
                self.stats['successful_pdb2gmx'] += 1
                result['success'] = True
            else:
                result['error'] = "pdb2gmx 실행 실패"
            
            # 처리 시간 기록
            processing_time = time.time() - entry_start_time
            self.stats['processing_times'][pdb_id] = processing_time
            result['processing_time'] = processing_time
            
            return result
            
        except Exception as e:
            result['error'] = str(e)
            self.stats['failed_entries'].append({
                'pdb_id': pdb_id,
                'error': str(e)
            })
            return result
    
    def download_pdb(self, pdb_id: str) -> Optional[Path]:
        """PDB 파일 다운로드"""
        pdb_dir = self.output_dir / "downloaded_pdbs"
        pdb_dir.mkdir(exist_ok=True)
        
        pdb_file = pdb_dir / f"{pdb_id}.pdb"
        
        # 이미 다운로드된 경우 건너뛰기
        if pdb_file.exists():
            self.logger.info(f"이미 다운로드됨: {pdb_file}")
            return pdb_file
        
        try:
            # RCSB PDB에서 다운로드
            url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
            
            self.logger.info(f"다운로드 시도: {url}")
            response = requests.get(url, timeout=30)
            response.raise_for_status()
            
            # PDB 내용 확인
            content = response.text
            if len(content) < 100 or "HEADER" not in content:
                raise ValueError("유효하지 않은 PDB 내용")
            
            # 파일 저장
            with open(pdb_file, 'w') as f:
                f.write(content)
            
            self.logger.info(f"다운로드 완료: {pdb_file}")
            return pdb_file
            
        except Exception as e:
            self.logger.error(f"PDB 다운로드 실패 {pdb_id}: {e}")
            return None
    
    def extract_chains(self, pdb_file: Path, receptor_chain: str, 
                      peptide_chain: str) -> Optional[Path]:
        """Receptor/Ligand chain 추출"""
        
        # 출력 디렉토리
        extracted_dir = self.output_dir / "extracted_chains"
        extracted_dir.mkdir(exist_ok=True)
        
        pdb_id = pdb_file.stem
        output_pdb = extracted_dir / f"{pdb_id}_extracted.pdb"
        
        # 이미 추출된 경우 건너뛰기
        if output_pdb.exists():
            self.logger.info(f"이미 추출됨: {output_pdb}")
            return output_pdb
        
        try:
            if not CHAIN_EXTRACTOR_AVAILABLE:
                self.logger.error("Chain 추출기를 사용할 수 없습니다")
                return None
            
            # SimpleChainExtractor 사용
            extractor = SimpleChainExtractor(logger=self.logger)
            chains_to_extract = [receptor_chain, peptide_chain]
            
            success, result = extractor.extract_chains(
                str(pdb_file), str(output_pdb), chains_to_extract
            )
            
            if success:
                self.logger.info(f"Chain 추출 완료: {output_pdb}")
                return output_pdb
            else:
                self.logger.error(f"Chain 추출 실패: {result.get('error', 'Unknown')}")
                return None
                
        except Exception as e:
            self.logger.error(f"Chain 추출 중 예외: {e}")
            return None
    
    def process_pdb_for_gromacs_compatibility(self, pdb_file: Path, 
                                            target_chains: List[str]) -> Path:
        """GROMACS 호환성을 위한 PDB 전처리"""
        
        processed_dir = self.output_dir / "processed_pdbs"
        processed_dir.mkdir(exist_ok=True)
        
        pdb_id = pdb_file.stem
        processed_pdb = processed_dir / f"{pdb_id}_processed.pdb"
        
        # 이미 처리된 경우 건너뛰기
        if processed_pdb.exists():
            self.logger.info(f"이미 처리됨: {processed_pdb}")
            return processed_pdb
        
        try:
            # 방법 1: mda_pdb_processor 사용 (권장)
            if PDB_PROCESSOR_AVAILABLE:
                self.logger.info("mda_pdb_processor 사용")
                result_pdb, stats = process_pdb_for_gromacs(
                    str(pdb_file), str(processed_pdb), target_chains, self.logger
                )
                self.logger.info(f"PDB 전처리 완료: {stats}")
                return processed_pdb
            
            # 방법 2: 기본 전처리 (fallback)
            else:
                self.logger.info("기본 PDB 전처리 수행")
                self._basic_pdb_processing(pdb_file, processed_pdb, target_chains)
                return processed_pdb
            
        except Exception as e:
            self.logger.error(f"PDB 전처리 실패: {e}")
            # 실패시 원본 파일 복사
            shutil.copy(pdb_file, processed_pdb)
            self.logger.warning("원본 파일로 대체")
            return processed_pdb
    
    def _basic_pdb_processing(self, input_pdb: Path, output_pdb: Path, 
                            target_chains: List[str]):
        """기본 PDB 전처리 (fallback)"""
        
        from Bio.PDB import PDBParser, PDBIO, Select
        
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("structure", input_pdb)
        
        class BasicProcessingSelect(Select):
            def __init__(self, target_chains):
                self.target_chains = set(target_chains)
                self.removed_residues = 0
            
            def accept_chain(self, chain):
                return chain.id in self.target_chains
            
            def accept_residue(self, residue):
                resname = residue.get_resname().strip()
                
                # 표준 아미노산만 허용
                standard_aa = {
                    'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY',
                    'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER',
                    'THR', 'TRP', 'TYR', 'VAL'
                }
                
                if resname not in standard_aa:
                    self.removed_residues += 1
                    return False
                
                # 백본 원자 확인
                present_atoms = {atom.get_name().strip() for atom in residue}
                backbone = {'N', 'CA', 'C', 'O'}
                
                if not backbone.issubset(present_atoms):
                    self.removed_residues += 1
                    return False
                
                return True
        
        selector = BasicProcessingSelect(target_chains)
        io = PDBIO()
        io.set_structure(structure)
        io.save(output_pdb, selector)
        
        self.logger.info(f"기본 전처리 완료: {selector.removed_residues}개 잔기 제거")
    
    def run_pdb2gmx(self, pdb_file: Path, pdb_id: str) -> bool:
        """pdb2gmx 실행"""
        
        if not GROMACS_RUNNER_AVAILABLE:
            self.logger.error("GROMACS Runner를 사용할 수 없습니다")
            return False
        
        try:
            # 작업 디렉토리
            work_dir = self.output_dir / "gromacs_results" / pdb_id
            work_dir.mkdir(parents=True, exist_ok=True)
            
            # 입력 PDB 복사
            input_pdb = work_dir / "input.pdb"
            shutil.copy(pdb_file, input_pdb)
            
            # GROMACS Runner 초기화
            runner = GromacsCommandRunner(
                str(work_dir), 
                self.config_file, 
                self.logger, 
                "medium"
            )
            
            # pdb2gmx 실행
            parameters = {
                "input_pdb": "input.pdb",
                "output_gro": "processed.gro",
                "output_top": "topol.top"
            }
            
            success, message = runner.try_multiple_force_fields(parameters)
            
            if success:
                self.logger.info(f"pdb2gmx 성공: {message}")
                
                # 결과 파일 확인
                required_files = ["processed.gro", "topol.top"]
                missing_files = []
                
                for file_name in required_files:
                    file_path = work_dir / file_name
                    if not file_path.exists():
                        missing_files.append(file_name)
                
                if missing_files:
                    self.logger.error(f"결과 파일 누락: {missing_files}")
                    return False
                
                return True
            else:
                self.logger.error(f"pdb2gmx 실패: {message}")
                return False
                
        except Exception as e:
            self.logger.error(f"pdb2gmx 실행 중 예외: {e}")
            return False


def main():
    """메인 함수 - 명령행 인터페이스"""
    import argparse
    
    parser = argparse.ArgumentParser(
        description="CSV 기반 PDB 처리 파이프라인",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
사용 예시:
  # 기본 사용 (전체 CSV 처리)
  python pdb_csv_pipeline.py --csv data.csv --output results
  
  # 처리 항목 수 제한 (테스트용)
  python pdb_csv_pipeline.py --csv data.csv --output results --max-entries 5
  
  # GROMACS 설정 파일 지정
  python pdb_csv_pipeline.py --csv data.csv --output results --config ../gromacs_commands_config.json
        """
    )
    
    parser.add_argument("--csv", required=True, help="입력 CSV 파일")
    parser.add_argument("--output", default="pdb_processing_results", help="출력 디렉토리")
    parser.add_argument("--max-entries", type=int, help="처리할 최대 항목 수")
    parser.add_argument("--config", help="GROMACS 설정 파일 경로")
    parser.add_argument("--verbose", "-v", action="store_true", help="상세 로그 출력")
    
    args = parser.parse_args()
    
    # 로깅 레벨 설정
    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)
    
    # 입력 파일 확인
    if not os.path.exists(args.csv):
        print(f"❌ 입력 CSV 파일을 찾을 수 없습니다: {args.csv}")
        return 1
    
    try:
        # 파이프라인 실행
        pipeline = PDBProcessingPipeline(
            output_dir=args.output,
            config_file=args.config
        )
        
        summary = pipeline.process_csv(args.csv, args.max_entries)
        
        # 결과 출력
        print(f"\n=== 처리 완료 ===")
        print(f"총 항목: {summary['total_entries']}")
        print(f"성공: {summary['successful']}")
        print(f"실패: {summary['failed']}")
        print(f"성공률: {summary['successful']/summary['total_entries']*100:.1f}%")
        print(f"총 처리 시간: {summary['total_time']:.1f}초")
        print(f"결과 디렉토리: {args.output}")
        
        # 실패한 항목들 출력
        failed_results = [r for r in summary['results'] if not r['success']]
        if failed_results:
            print(f"\n실패한 항목들:")
            for result in failed_results[:10]:  # 처음 10개만 출력
                print(f"  - {result['pdb_id']}: {result.get('error', 'Unknown')}")
            
            if len(failed_results) > 10:
                print(f"  ... 및 {len(failed_results) - 10}개 더")
        
        return 0
        
    except Exception as e:
        print(f"❌ 파이프라인 실행 실패: {e}")
        return 1


if __name__ == "__main__":
    exit(main())
