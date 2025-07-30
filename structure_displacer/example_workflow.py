#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
50Å 거리 기반 구조 변형 및 SuMD 시뮬레이션 전체 워크플로우

새로운 변형 방식에 맞춘 완전한 워크플로우:
1. Golden Standard PDB에서 50Å 거리로 변형된 구조들 생성
2. 변형된 구조들의 품질 분석 (거리 달성도, Clash, 경로 차단)
3. 최적의 유효한 구조 선택
4. 선택된 구조로 SuMD 시뮬레이션 실행

사용법:
    python example_workflow.py --golden_standard native_complex.pdb --receptor_chains A,C --ligand_chains B,D --num_variants 5
"""

import os
import sys
import argparse
import logging
import json
import subprocess
from datetime import datetime
from typing import List, Dict, Tuple, Optional

# 새로운 모듈들 import
try:
    from structure_displacer import MultipleDisplacer, DisplacementConfig
    from displacement_utils import Advanced50AAnalyzer, collect_displaced_files
except ImportError:
    print("오류: structure_displacer.py와 displacement_utils.py가 필요합니다.")
    sys.exit(1)


class Advanced50AWorkflowManager:
    """50Å 거리 기반 SuMD 워크플로우 관리자"""
    
    def __init__(self, config_dict: Dict = None, logger=None):
        self.logger = logger or logging.getLogger(__name__)
        self.config = config_dict or {}
        
        # 기본 설정
        self.default_config = {
            'num_variants': 5,
            'displacement_config': {
                'target_distance': 50.0,
                'clash_threshold': 2.0,
                'max_attempts': 200,
                'path_check_threshold': 5.0
            },
            'sumd_config': {
                'simulation_time': 2.0,
                'distance_threshold': 5.0,
                'rmsd_threshold': 1.5,
                'max_iterations': 10,
                'num_samples': 5
            },
            'analysis_config': {
                'create_plots': True,
                'select_best_valid_only': True,  # 유효한 구조들 중에서만 선택
                'min_quality_score': 0.7        # 최소 품질 점수 요구
            }
        }
        
        # 설정 병합
        self.merged_config = self._merge_configs(self.default_config, self.config)
    
    def _merge_configs(self, default: Dict, user: Dict) -> Dict:
        """설정 딕셔너리 병합"""
        merged = default.copy()
        for key, value in user.items():
            if isinstance(value, dict) and key in merged:
                merged[key].update(value)
            else:
                merged[key] = value
        return merged
    
    def step1_generate_50A_displaced_structures(self, golden_standard: str, receptor_chains: List[str], ligand_chains: List[str],
                                              output_dir: str) -> List[Tuple[str, Dict]]:
        """
        단계 1: 50Å 거리 기반 변형된 구조들 생성
        
        Args:
            golden_standard: Golden Standard PDB 파일
            receptor_chains: 수용체 체인 ID 리스트 (고정됨)
            ligand_chains: 리간드 체인 ID 리스트 (이동함)
            output_dir: 출력 디렉토리
            
        Returns:
            List[Tuple[str, Dict]]: 생성된 구조 파일들과 정보
        """
        self.logger.info("=== 단계 1: 50Å 거리 기반 변형된 구조 생성 ===")
        
        if not receptor_chains or not ligand_chains:
            raise ValueError("수용체 체인과 리간드 체인이 모두 필요합니다")
        
        # 변형 설정 생성
        displacement_config = DisplacementConfig()
        config_dict = self.merged_config['displacement_config']
        
        displacement_config.actual_displacement_distance = 60.0  # 고정 거리 모드용 (사용안함)
        displacement_config.target_distance_for_validation = 50.0  # 목표 표면간 거리
        displacement_config.use_adaptive_displacement = True  # 적응형 거리 계산 사용
        displacement_config.clash_threshold = config_dict['clash_threshold']
        displacement_config.max_attempts = config_dict['max_attempts']
        displacement_config.path_check_threshold = config_dict['path_check_threshold']
        
        # 변형 구조 생성기
        displacer = MultipleDisplacer(displacement_config, self.logger)
        
        # 변형 구조 생성
        displaced_dir = os.path.join(output_dir, "displaced_structures_50A")
        results = displacer.generate_multiple_variants(
            golden_standard, 
            receptor_chains,
            ligand_chains,
            self.merged_config['num_variants'],
            displaced_dir
        )
        
        self.logger.info(f"단계 1 완료: {len(results)}개 변형 구조 생성")
        
        # 생성 통계
        total_attempts = 0
        successful_structures = 0
        
        for _, info in results:
            if info['final_result']['success']:
                successful_structures += 1
                total_attempts += info['final_result']['total_attempts']
        
        avg_attempts = total_attempts / successful_structures if successful_structures > 0 else 0
        
        self.logger.info(f"생성 통계: 성공 {successful_structures}/{self.merged_config['num_variants']}개, "
                        f"평균 시도 횟수 {avg_attempts:.1f}회")
        
        return results
    
    def step2_analyze_50A_displacement_quality(self, golden_standard: str, displaced_results: List[Tuple[str, Dict]],
                                             receptor_chains: List[str], ligand_chains: List[str], output_dir: str) -> Tuple[str, Dict]:
        """
        단계 2: 50Å 변형 품질 분석
        
        Args:
            golden_standard: Golden Standard PDB 파일
            displaced_results: 변형된 구조 결과들
            receptor_chains: 수용체 체인 ID 리스트
            ligand_chains: 리간드 체인 ID 리스트
            output_dir: 출력 디렉토리
            
        Returns:
            Tuple[str, Dict]: (최적 구조 파일, 분석 결과)
        """
        self.logger.info("=== 단계 2: 50Å 변형 품질 분석 ===")
        
        # 변형된 파일들 추출 (성공한 것들만)
        displaced_files = []
        for result in displaced_results:
            output_file, info = result
            if info['final_result']['success']:
                displaced_files.append(output_file)
        
        if not displaced_files:
            raise Exception("성공적으로 생성된 변형 구조가 없습니다.")
        
        self.logger.info(f"분석할 변형 구조: {len(displaced_files)}개")
        
        # 분석기 생성
        target_distance = self.merged_config['displacement_config']['target_distance']
        analyzer = Advanced50AAnalyzer(target_distance, self.logger)
        
        # 품질 분석 실행 (호환성을 위해 첫 번째 수용체와 첫 번째 리간드 사용)
        analysis_dir = os.path.join(output_dir, "quality_analysis_50A")
        target_chains = [receptor_chains[0], ligand_chains[0]]  # 분석용 대표 체인
        df = analyzer.analyze_batch_structures(
            golden_standard, displaced_files, target_chains, analysis_dir
        )
        
        if df.empty:
            raise Exception("품질 분석 결과가 없습니다.")
        
        # 시각화 생성
        if self.merged_config['analysis_config']['create_plots']:
            analyzer.create_analysis_plots(df, analysis_dir)
        
        # 최적 구조 선택 로직
        select_best_valid_only = self.merged_config['analysis_config']['select_best_valid_only']
        min_quality_score = self.merged_config['analysis_config']['min_quality_score']
        
        # 선택 기준에 따라 후보 필터링
        if select_best_valid_only:
            candidates = df[df['is_valid'] == True].copy()
            if candidates.empty:
                self.logger.warning("유효한 구조가 없습니다. 전체 구조에서 선택합니다.")
                candidates = df.copy()
        else:
            candidates = df.copy()
        
        # 최소 품질 점수 필터링
        quality_filtered = candidates[candidates['quality_score'] >= min_quality_score].copy()
        if quality_filtered.empty:
            self.logger.warning(f"품질 점수 {min_quality_score} 이상인 구조가 없습니다. 기준을 완화합니다.")
            quality_filtered = candidates.copy()
        
        # 최적 구조 선택 (품질 점수 기준)
        if quality_filtered.empty:
            raise Exception("선택할 수 있는 구조가 없습니다.")
        
        best_structure = quality_filtered.iloc[0]
        best_file = best_structure['filepath']
        
        # 분석 요약 생성
        analysis_summary = {
            'selection_criteria': {
                'select_valid_only': select_best_valid_only,
                'min_quality_score': min_quality_score,
                'target_distance': target_distance
            },
            'best_structure': {
                'file': best_file,
                'filename': best_structure['filename'],
                'rank': int(best_structure.get('rank', 1)),
                'quality_score': float(best_structure['quality_score']),
                'is_valid': bool(best_structure['is_valid']),
                'min_distance': float(best_structure['min_distance']),
                'distance_error': float(best_structure['distance_error']),
                'has_clash': bool(best_structure['has_clash']),
                'is_path_obstructed': bool(best_structure['is_path_obstructed']),
                'rmsd_vs_golden': float(best_structure.get('rmsd_vs_golden', float('inf')))
            },
            'analysis_stats': {
                'total_analyzed': len(df),
                'valid_structures': len(df[df['is_valid'] == True]),
                'invalid_structures': len(df[df['is_valid'] == False]),
                'validation_rate': float(len(df[df['is_valid'] == True]) / len(df) * 100),
                'distance_achievement_rate': float(len(df[df['distance_acceptable']]) / len(df) * 100),
                'clash_free_rate': float(len(df[~df['has_clash']]) / len(df) * 100),
                'path_clear_rate': float(len(df[~df['is_path_obstructed']]) / len(df) * 100),
                'mean_quality_score': float(df['quality_score'].mean()),
                'mean_distance': float(df['min_distance'].mean()),
                'distance_std': float(df['min_distance'].std())
            }
        }
        
        # 상위 구조들 정보 (최대 5개)
        top_structures = quality_filtered.head(5)
        analysis_summary['top_structures'] = []
        
        for idx, row in top_structures.iterrows():
            analysis_summary['top_structures'].append({
                'rank': int(row.get('rank', idx + 1)),
                'filename': row['filename'],
                'quality_score': float(row['quality_score']),
                'is_valid': bool(row['is_valid']),
                'min_distance': float(row['min_distance']),
                'distance_error': float(row['distance_error'])
            })
        
        self.logger.info(f"단계 2 완료: 최적 구조 선택")
        self.logger.info(f"선택된 구조: {best_structure['filename']} "
                        f"(품질점수: {best_structure['quality_score']:.3f}, "
                        f"유효성: {best_structure['is_valid']}, "
                        f"거리: {best_structure['min_distance']:.2f}Å)")
        
        return best_file, analysis_summary
    
    def step3_run_sumd_simulation(self, displaced_structure: str, golden_standard: str,
                                receptor_chains: List[str], ligand_chains: List[str], output_dir: str) -> Dict:
        """
        단계 3: SuMD 시뮬레이션 실행
        
        Args:
            displaced_structure: 변형된 구조 파일 (시작점)
            golden_standard: Golden Standard 파일 (목표점)
            receptor_chains: 수용체 체인 ID 리스트
            ligand_chains: 리간드 체인 ID 리스트
            output_dir: 출력 디렉토리
            
        Returns:
            Dict: 시뮬레이션 결과
        """
        self.logger.info("=== 단계 3: SuMD 시뮬레이션 실행 ===")
        
        if not receptor_chains or not ligand_chains:
            raise ValueError("수용체 체인과 리간드 체인이 모두 필요합니다")
        
        sumd_config = self.merged_config['sumd_config']
        sumd_output_dir = os.path.join(output_dir, "sumd_simulation_50A")
        
        # SuMD 스크립트 경로 찾기
        sumd_script_candidates = [
            "golden_standard_sumd_master.py",
            "robust_sumd_master.py", 
            "/app/scripts/golden_standard_sumd_master.py",
            "/app/scripts/robust_sumd_master.py"
        ]
        
        sumd_script = None
        for candidate in sumd_script_candidates:
            if os.path.exists(candidate):
                sumd_script = candidate
                break
        
        if sumd_script is None:
            raise Exception("SuMD 스크립트를 찾을 수 없습니다.")
        
        self.logger.info(f"SuMD 스크립트 사용: {sumd_script}")
        
        # SuMD 명령어 구성
        cmd = [
            "python3", sumd_script,
            "--input_pdb", displaced_structure,
            "--output_dir", sumd_output_dir,
            "--simulation_time", str(sumd_config['simulation_time']),
            "--receptor_chain", receptor_chains[0],  # 첫 번째 수용체 체인 사용
            "--ligand_chain", ligand_chains[0],      # 첫 번째 리간드 체인 사용
            "--distance_threshold", str(sumd_config['distance_threshold']),
            "--rmsd_threshold", str(sumd_config['rmsd_threshold']),
            "--max_iterations", str(sumd_config['max_iterations']),
            "--num_samples", str(sumd_config['num_samples'])
        ]
        
        # Golden Standard 스크립트인 경우 추가 파라미터
        if "golden_standard" in sumd_script:
            cmd.extend(["--golden_standard_pdb", golden_standard])
        
        self.logger.info(f"SuMD 실행 명령어: {' '.join(cmd)}")
        
        simulation_result = {'success': False}
        
        try:
            # SuMD 실행
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=7200)  # 2시간 타임아웃
            
            if result.returncode == 0:
                # 성공 결과 파싱
                output_lines = result.stdout.split('\n')
                sumd_result_line = None
                
                for line in output_lines:
                    if line.startswith("SUMD_RESULT:"):
                        sumd_result_line = line
                        break
                
                if sumd_result_line:
                    parts = sumd_result_line.split(':')
                    if len(parts) >= 3:
                        final_pdb = parts[1]
                        converged = parts[2] == 'True'
                        
                        simulation_result = {
                            'success': True,
                            'final_pdb': final_pdb,
                            'converged': converged,
                            'output_dir': sumd_output_dir,
                            'command_used': cmd,
                            'stdout_excerpt': result.stdout[-1000:],  # 마지막 1000자만 저장
                            'execution_info': {
                                'return_code': result.returncode,
                                'script_used': sumd_script,
                                'golden_standard_mode': "golden_standard" in sumd_script
                            }
                        }
                        
                        self.logger.info(f"SuMD 시뮬레이션 완료: 수렴={converged}")
                        if final_pdb:
                            self.logger.info(f"최종 구조: {os.path.basename(final_pdb)}")
                        
                    else:
                        raise Exception("SuMD 결과 파싱 실패: 잘못된 형식")
                else:
                    raise Exception("SuMD 결과를 찾을 수 없습니다.")
            else:
                raise Exception(f"SuMD 실행 실패 (코드: {result.returncode})\n{result.stderr}")
        
        except subprocess.TimeoutExpired:
            simulation_result = {
                'success': False,
                'error': "SuMD 실행 시간 초과 (2시간)",
                'timeout': True
            }
            self.logger.error("SuMD 실행 시간 초과")
        except Exception as e:
            simulation_result = {
                'success': False,
                'error': str(e),
                'command_used': cmd,
                'stdout': result.stdout if 'result' in locals() else '',
                'stderr': result.stderr if 'result' in locals() else ''
            }
            self.logger.error(f"SuMD 시뮬레이션 실패: {e}")
        
        return simulation_result
    
    def run_complete_50A_workflow(self, golden_standard: str, receptor_chains: List[str], ligand_chains: List[str],
                                output_dir: str) -> Dict:
        """
        전체 50Å 기반 워크플로우 실행
        
        Args:
            golden_standard: Golden Standard PDB 파일
            receptor_chains: 수용체 체인 ID 리스트
            ligand_chains: 리간드 체인 ID 리스트
            output_dir: 출력 디렉토리
            
        Returns:
            Dict: 전체 워크플로우 결과
        """
        os.makedirs(output_dir, exist_ok=True)
        
        workflow_result = {
            'workflow_type': '50A_distance_based',
            'start_time': datetime.now().isoformat(),
            'golden_standard': golden_standard,
            'receptor_chains': receptor_chains,
            'ligand_chains': ligand_chains,
            'config': self.merged_config,
            'steps': {}
        }
        
        try:
            # 단계 1: 50Å 변형 구조 생성
            displaced_results = self.step1_generate_50A_displaced_structures(
                golden_standard, receptor_chains, ligand_chains, output_dir
            )
            
            successful_structures = [r for r in displaced_results if r[1]['final_result']['success']]
            
            workflow_result['steps']['displacement'] = {
                'success': len(successful_structures) > 0,
                'total_attempted': len(displaced_results),
                'successful_structures': len(successful_structures),
                'success_rate': len(successful_structures) / len(displaced_results) * 100 if displaced_results else 0,
                'structures': [result[0] for result in successful_structures],
                'target_distance': self.merged_config['displacement_config']['target_distance']
            }
            
            if len(successful_structures) == 0:
                raise Exception("변형된 구조 생성에 완전히 실패했습니다.")
            
            # 단계 2: 품질 분석
            best_structure, analysis_summary = self.step2_analyze_50A_displacement_quality(
                golden_standard, displaced_results, receptor_chains, ligand_chains, output_dir
            )
            workflow_result['steps']['analysis'] = {
                'success': True,
                'best_structure': best_structure,
                'analysis_summary': analysis_summary
            }
            
            # 단계 3: SuMD 시뮬레이션
            simulation_result = self.step3_run_sumd_simulation(
                best_structure, golden_standard, receptor_chains, ligand_chains, output_dir
            )
            workflow_result['steps']['simulation'] = simulation_result
            
            workflow_result['overall_success'] = simulation_result['success']
            
        except Exception as e:
            workflow_result['overall_success'] = False
            workflow_result['error'] = str(e)
            self.logger.error(f"워크플로우 실행 실패: {e}")
        
        workflow_result['end_time'] = datetime.now().isoformat()
        
        # 결과 저장
        result_file = os.path.join(output_dir, "workflow_result_50A.json")
        with open(result_file, 'w') as f:
            json.dump(workflow_result, f, indent=2, default=str)
        
        return workflow_result


def main():
    """메인 함수"""
    parser = argparse.ArgumentParser(
        description="50Å 거리 기반 구조 변형 및 SuMD 시뮬레이션 전체 워크플로우",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
50Å 거리 기반 워크플로우:
1. Golden Standard PDB에서 50Å 거리로 변형된 구조들 생성
2. 변형된 구조들의 품질 분석 (거리 달성도, Clash, 경로 차단)
3. 최적의 유효한 구조 선택  
4. 선택된 구조로 SuMD 시뮬레이션 실행

특징:
- 두 체인 간 최소 거리가 정확히 50Å이 되도록 이동
- Clash 검사 및 경로 차단 검사
- 유효한 구조들 중에서 최고 품질 선택

사용 예시:
  python example_workflow.py --golden_standard native_complex.pdb --receptor_chains A,C --ligand_chains B,D --num_variants 5
        """
    )
    
    # 필수 인수
    parser.add_argument("--golden_standard", required=True,
                       help="Golden Standard PDB 파일")
    parser.add_argument("--receptor_chains", required=True,
                       help="수용체 체인 ID (쉼표로 구분, 고정됨, 예: A,C)")
    parser.add_argument("--ligand_chains", required=True,
                       help="리간드 체인 ID (쉼표로 구분, 이동함, 예: B,D)")
    
    # 워크플로우 설정
    parser.add_argument("--num_variants", type=int, default=5,
                       help="생성할 변형 구조 수 (기본값: 5)")
    parser.add_argument("--output_dir", default="workflow_output_50A",
                       help="출력 디렉토리 (기본값: workflow_output_50A)")
    
    # 변형 설정
    parser.add_argument("--target_distance", type=float, default=50.0,
                       help="목표 거리 (Å, 기본값: 50.0)")
    parser.add_argument("--clash_threshold", type=float, default=2.0,
                       help="Clash 판정 거리 (Å, 기본값: 2.0)")
    parser.add_argument("--max_attempts", type=int, default=200,
                       help="변형 생성 최대 시도 횟수 (기본값: 200)")
    
    # SuMD 설정
    parser.add_argument("--simulation_time", type=float, default=2.0,
                       help="SuMD 시뮬레이션 시간 (ns, 기본값: 2.0)")
    parser.add_argument("--rmsd_threshold", type=float, default=1.5,
                       help="SuMD 수렴 RMSD 임계값 (Å, 기본값: 1.5)")
    parser.add_argument("--max_iterations", type=int, default=10,
                       help="SuMD 최대 반복 횟수 (기본값: 10)")
    
    # 분석 설정
    parser.add_argument("--min_quality_score", type=float, default=0.7,
                       help="최소 품질 점수 요구 (기본값: 0.7)")
    parser.add_argument("--allow_invalid", action="store_true",
                       help="무효한 구조도 선택 허용")
    
    # 기타
    parser.add_argument("--verbose", "-v", action="store_true", help="상세 로그 출력")
    parser.add_argument("--skip_sumd", action="store_true", help="SuMD 시뮬레이션 건너뛰기")
    
    args = parser.parse_args()
    
    # 로깅 설정
    level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )
    logger = logging.getLogger("50A-WorkflowManager")
    
    # 입력 파일 확인
    if not os.path.exists(args.golden_standard):
        logger.error(f"Golden Standard 파일을 찾을 수 없습니다: {args.golden_standard}")
        return 1
    
    # 체인 파싱
    receptor_chains = [chain.strip().upper() for chain in args.receptor_chains.split(',')]
    ligand_chains = [chain.strip().upper() for chain in args.ligand_chains.split(',')]
    
    if not receptor_chains or not ligand_chains:
        logger.error("수용체 체인과 리간드 체인이 모두 필요합니다")
        return 1
    
    logger.info(f"수용체 체인 (고정): {receptor_chains}")
    logger.info(f"리간드 체인 (이동): {ligand_chains}")
    
    # 출력 디렉토리 자동 생성
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir, exist_ok=True)
        logger.info(f"출력 디렉토리 생성: {args.output_dir}")
    else:
        logger.info(f"기존 출력 디렉토리 사용: {args.output_dir}")
        
    # 설정 준비
    config = {
        'num_variants': args.num_variants,
        'displacement_config': {
            'target_distance': args.target_distance,
            'clash_threshold': args.clash_threshold,
            'max_attempts': args.max_attempts
        },
        'sumd_config': {
            'simulation_time': args.simulation_time,
            'rmsd_threshold': args.rmsd_threshold,
            'max_iterations': args.max_iterations
        },
        'analysis_config': {
            'select_best_valid_only': not args.allow_invalid,
            'min_quality_score': args.min_quality_score
        }
    }
    
    logger.info(f"워크플로우 설정:")
    logger.info(f"  목표 거리: {args.target_distance}Å")
    logger.info(f"  변형 구조 수: {args.num_variants}개")
    logger.info(f"  유효한 구조만 선택: {not args.allow_invalid}")
    logger.info(f"  최소 품질 점수: {args.min_quality_score}")
    
    try:
        # 워크플로우 관리자 생성
        workflow_manager = Advanced50AWorkflowManager(config, logger)
        
        # 전체 워크플로우 실행
        if args.skip_sumd:
            logger.info("SuMD 시뮬레이션을 건너뜁니다.")
            
            # 단계 1, 2만 실행
            displaced_results = workflow_manager.step1_generate_50A_displaced_structures(
                args.golden_standard, receptor_chains, ligand_chains, args.output_dir
            )
            
            best_structure, analysis_summary = workflow_manager.step2_analyze_50A_displacement_quality(
                args.golden_standard, displaced_results, receptor_chains, ligand_chains, args.output_dir
            )
            
            result = {
                'displacement_completed': True,
                'analysis_completed': True,
                'best_structure': best_structure,
                'analysis_summary': analysis_summary,
                'sumd_skipped': True
            }
        else:
            result = workflow_manager.run_complete_50A_workflow(
                args.golden_standard, receptor_chains, ligand_chains, args.output_dir
            )
        
        # 결과 요약 출력
        print(f"\n=== 50Å 기반 워크플로우 실행 결과 ===")
        
        if 'steps' in result:
            # 단계별 결과
            if 'displacement' in result['steps']:
                disp_result = result['steps']['displacement']
                print(f"1. 50Å 구조 변형: {'성공' if disp_result['success'] else '실패'}")
                if disp_result['success']:
                    print(f"   성공률: {disp_result['success_rate']:.1f}% ({disp_result['successful_structures']}/{disp_result['total_attempted']}개)")
                    print(f"   목표 거리: {disp_result['target_distance']}Å")
            
            if 'analysis' in result['steps']:
                analysis_result = result['steps']['analysis']
                print(f"2. 품질 분석: {'성공' if analysis_result['success'] else '실패'}")
                if analysis_result['success']:
                    best = analysis_result['analysis_summary']['best_structure']
                    stats = analysis_result['analysis_summary']['analysis_stats']
                    print(f"   최적 구조: {best['filename']}")
                    print(f"   품질 점수: {best['quality_score']:.3f}")
                    print(f"   유효성: {'유효' if best['is_valid'] else '무효'}")
                    print(f"   거리: {best['min_distance']:.2f}Å (오차: {best['distance_error']:.2f}Å)")
                    print(f"   전체 유효율: {stats['validation_rate']:.1f}%")
            
            if 'simulation' in result['steps']:
                sim_result = result['steps']['simulation']
                print(f"3. SuMD 시뮬레이션: {'성공' if sim_result['success'] else '실패'}")
                if sim_result['success']:
                    print(f"   수렴: {'예' if sim_result['converged'] else '아니오'}")
                    if sim_result.get('final_pdb'):
                        print(f"   최종 구조: {os.path.basename(sim_result['final_pdb'])}")
                elif 'error' in sim_result:
                    print(f"   오류: {sim_result['error']}")
        
        if args.skip_sumd and 'best_structure' in result:
            print(f"선택된 최적 구조: {os.path.basename(result['best_structure'])}")
        
        print(f"\n결과 디렉토리: {args.output_dir}")
        print(f"  - 변형된 구조들: {args.output_dir}/displaced_structures_50A/")
        print(f"  - 품질 분석 결과: {args.output_dir}/quality_analysis_50A/")
        
        if not args.skip_sumd:
            print(f"  - SuMD 시뮬레이션: {args.output_dir}/sumd_simulation_50A/")
        
        print(f"  - 전체 결과: {args.output_dir}/workflow_result_50A.json")
        
        return 0
        
    except Exception as e:
        logger.error(f"워크플로우 실행 실패: {e}")
        return 1


if __name__ == "__main__":
    exit(main())
