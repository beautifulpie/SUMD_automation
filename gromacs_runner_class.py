#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import json
import subprocess
import logging
import multiprocessing
from typing import Dict, List, Optional, Tuple, Any

class GromacsCommandRunner:
    """GROMACS 명령어 실행을 위한 통합 클래스 - 정리된 버전"""
    
    def __init__(self, work_dir: str, config_file: str = None, logger=None, system_size: str = "medium"):
        self.work_dir = work_dir
        self.logger = logger or logging.getLogger(__name__)
        self.system_size = system_size
        
        # 설정 파일 로드
        if config_file and os.path.exists(config_file):
            self.config = self._load_config(config_file)
        else:
            # 기본 설정 (JSON 파일이 없는 경우)
            self.config = self._get_default_config()
        
        # MPI 환경 설정 (컨테이너 대응)
        self._setup_mpi_environment()
        
        self.logger.info(f"GROMACS Runner 초기화: 작업 디렉토리={work_dir}, 시스템 크기={system_size}")
    
    def _setup_mpi_environment(self):
        """MPI 환경 설정 (컨테이너에서 root 실행 허용)"""
        mpi_settings = self.config.get("mpi_settings", {})
        container_env = mpi_settings.get("container_environment", {})
        
        if container_env.get("allow_root", False):
            env_vars = container_env.get("environment_variables", {})
            for key, value in env_vars.items():
                os.environ[key] = str(value)
                self.logger.debug(f"MPI 환경변수 설정: {key}={value}")
        
        # GPU 성능 최적화 환경변수
        perf_settings = mpi_settings.get("performance_tuning", {})
        if perf_settings.get("gpu_direct_communication", False):
            env_vars = perf_settings.get("environment_variables", {})
            for key, value in env_vars.items():
                os.environ[key] = str(value)
                self.logger.debug(f"GPU 최적화 환경변수 설정: {key}={value}")
    
    def _load_config(self, config_file: str) -> Dict:
        """설정 파일 로드"""
        try:
            with open(config_file, 'r', encoding='utf-8') as f:
                return json.load(f)
        except Exception as e:
            self.logger.error(f"설정 파일 로드 실패: {e}")
            return self._get_default_config()
    
    def _get_default_config(self) -> Dict:
        """기본 설정 반환 (JSON 파일이 없는 경우)"""
        return {
            "gromacs_commands": {
                "pdb2gmx": {
                    "template": "gmx pdb2gmx -f {input_pdb} -o {output_gro} -p {output_top} -water {water_model} -ff {force_field} -ignh {extra_flags}",
                    "defaults": {
                        "water_model": "tip3p",
                        "force_field": "charmm36-jul2022",
                        "extra_flags": "-maxwarn 10"
                    },
                    "force_field_options": ["charmm36-jul2022", "amber99sb-ildn", "gromos54a7", "oplsaa"]
                },
                "mdrun_em": {
                    "template": "mpirun --allow-run-as-root -np {mpi_ranks} gmx_mpi mdrun -v {prefix} -ntomp {ntomp} -nb {nb_mode} -gpu_id {gpu_id}",
                    "defaults": {
                        "mpi_ranks": "4",
                        "ntomp": "2",
                        "nb_mode": "gpu", 
                        "gpu_id": "3"
                    }
                },
                "mdrun_md": {
                    "template": "mpirun --allow-run-as-root -np {mpi_ranks} gmx_mpi mdrun -v {prefix} -ntomp {ntomp} -nb {nb_mode} -gpu_id {gpu_id} -pme {pme_mode}",
                    "defaults": {
                        "mpi_ranks": "4",
                        "ntomp": "2",
                        "nb_mode": "gpu",
                        "gpu_id": "3",
                        "pme_mode": "auto"
                    }
                }
            }
        }
    
    def _substitute_parameters(self, template: str, parameters: Dict[str, Any]) -> str:
        """템플릿에 파라미터 치환"""
        try:
            return template.format(**parameters)
        except KeyError as e:
            raise ValueError(f"필수 파라미터 누락: {e}")
    
    def _get_timeout(self, command_type: str) -> int:
        """시스템 크기에 따른 타임아웃 반환"""
        timeout_settings = self.config.get("timeout_settings", {})
        system_timeouts = timeout_settings.get(self.system_size, {})
        
        # 명령어 타입에 따른 기본 타임아웃
        default_timeouts = {
            "pdb2gmx": 300,
            "editconf": 60,
            "solvate": 300,
            "grompp": 120,
            "genion": 120,
            "mdrun_em": 1800,
            "mdrun_md": 3600,
            "mdrun_mpi": 2400,
            "trjconv": 120
        }
        
        return system_timeouts.get(command_type, default_timeouts.get(command_type, 300))
    
    def _check_files_exist(self, file_paths: List[str]) -> Tuple[bool, List[str]]:
        """파일 존재 여부 확인"""
        missing_files = []
        for file_path in file_paths:
            full_path = os.path.join(self.work_dir, file_path)
            if not os.path.exists(full_path):
                missing_files.append(file_path)
        
        return len(missing_files) == 0, missing_files
    
    def run_command(self, cmd_list: List[str], timeout: int = None, input_str: str = None, 
                   check_files: List[str] = None) -> Tuple[bool, str]:
        """단일 명령어 실행"""
        cmd_str = ' '.join(cmd_list)
        self.logger.info(f"명령어 실행: {cmd_str}")
        
        try:
            # input 문자열 처리 개선
            input_bytes = None
            if input_str:
                # 개행 문자 자동 추가 (GROMACS가 기대하는 형식)
                if not input_str.endswith('\n'):
                    input_str += '\n'
                input_bytes = input_str.encode('utf-8')
                self.logger.debug(f"입력 데이터: {repr(input_str)}")
            
            # subprocess 실행 (중요: check=False로 설정)
            result = subprocess.run(
                cmd_list,
                cwd=self.work_dir,
                input=input_bytes,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                timeout=timeout,
                check=False  # 중요: GROMACS 경고를 오류로 보지 않음
            )
            
            # stderr은 항상 디코딩하여 로그 및 오류 메시지에 활용
            stderr_text = result.stderr.decode('utf-8', errors='ignore')

            # 로그 출력
            if result.stdout:
                stdout_text = result.stdout.decode('utf-8', errors='ignore')
                self.logger.debug(f"STDOUT: {stdout_text[:500]}")
            
            if result.stderr:
                self.logger.debug(f"STDERR: {stderr_text[:500]}")
            
            # 성공 여부 판정: 파일 존재 여부로 판단 (returncode 무시)
            if check_files:
                files_exist, missing = self._check_files_exist(check_files)
                if files_exist:
                    self.logger.info(f"명령어 성공: 출력 파일 확인됨")
                    return True, "성공"
                else:
                    # 실패 시 GROMACS의 stderr 내용을 함께 반환하여 원인 파악을 돕는다.
                    error_msg = f"출력 파일 누락: {missing}. GROMACS 오류: {stderr_text[:200].strip()}"
                    self.logger.error(error_msg)
                    return False, error_msg
            else:
                # 파일 체크가 없는 경우 returncode로만 판정
                if result.returncode == 0:
                    return True, "성공"
                else:
                    error_msg = f"명령어 실행 실패 (코드: {result.returncode}). GROMACS 오류: {stderr_text[:200].strip()}"
                    self.logger.error(error_msg)
                    return False, error_msg
                    
        except subprocess.TimeoutExpired:
            error_msg = f"명령어 실행 시간 초과: {cmd_str}"
            self.logger.error(error_msg)
            return False, error_msg
        except Exception as e:
            error_msg = f"명령어 실행 중 예외: {e}"
            self.logger.error(error_msg)
            return False, error_msg
    
    def execute_gromacs_command(self, command_name: str, parameters: Dict[str, Any], 
                               max_retries: int = 1) -> Tuple[bool, str]:
        """GROMACS 명령어 실행 (재시도 포함)"""
        
        if command_name not in self.config["gromacs_commands"]:
            return False, f"알 수 없는 명령어: {command_name}"
        
        cmd_config = self.config["gromacs_commands"][command_name]
        
        # 기본값 적용
        full_params = cmd_config.get("defaults", {}).copy()
        
        # 시스템 크기별 오버라이드 적용
        if self.system_size == "large" and "large_system_overrides" in cmd_config:
            full_params.update(cmd_config["large_system_overrides"])
        
        # 사용자 파라미터 적용
        full_params.update(parameters)
        
        # 명령어 생성
        try:
            cmd_template = cmd_config["template"]
            cmd_str = self._substitute_parameters(cmd_template, full_params)
            cmd_list = cmd_str.split()
        except Exception as e:
            return False, f"명령어 생성 실패: {e}"
        
        # input 처리
        input_str = None
        if cmd_config.get("input_required", False):
            try:
                input_template = cmd_config.get("input_template", "")
                if input_template:
                    input_str = self._substitute_parameters(input_template, full_params)
                    self.logger.debug(f"입력 문자열 생성: {input_str}")
            except Exception as e:
                return False, f"입력 문자열 생성 실패: {e}"
        
        # 출력 파일 목록 생성
        output_files = []
        if "output_files" in cmd_config:
            for file_template in cmd_config["output_files"]:
                try:
                    output_file = self._substitute_parameters(file_template, full_params)
                    output_files.append(output_file)
                except:
                    pass  # 파라미터 치환 실패 시 무시
        
        # 타임아웃 설정
        timeout = self._get_timeout(command_name)
        
        # 재시도 로직
        for attempt in range(max_retries + 1):
            if attempt > 0:
                self.logger.info(f"{command_name} 재시도 {attempt}")
                # 재시도 시 더 관대한 설정
                if "maxwarn" in full_params:
                    full_params["maxwarn"] = str(int(full_params["maxwarn"]) + 5)
                    cmd_str = self._substitute_parameters(cmd_template, full_params)
                    cmd_list = cmd_str.split()
                
                timeout = timeout * 2  # 타임아웃 2배 증가
            
            # 명령어 실행
            success, message = self.run_command(
                cmd_list, 
                timeout=timeout, 
                input_str=input_str,
                check_files=output_files if output_files else None
            )
            
            if success:
                return True, message
            elif attempt < max_retries:
                self.logger.warning(f"{command_name} 시도 {attempt + 1} 실패: {message}")
            else:
                return False, f"{command_name} 최종 실패: {message}"
        
        return False, "예상치 못한 오류"
    
    def try_multiple_force_fields(self, base_parameters: Dict[str, Any]) -> Tuple[bool, str]:
        """여러 포스필드로 pdb2gmx 시도"""
        
        pdb2gmx_config = self.config["gromacs_commands"].get("pdb2gmx", {})
        force_fields = pdb2gmx_config.get("force_field_options", ["charmm36-jul2022"])
        
        for ff in force_fields:
            self.logger.info(f"포스필드 시도: {ff}")
            
            params = base_parameters.copy()
            params["force_field"] = ff
            
            success, message = self.execute_gromacs_command("pdb2gmx", params)
            
            if success:
                self.logger.info(f"포스필드 성공: {ff}")
                return True, f"성공한 포스필드: {ff}"
            else:
                self.logger.warning(f"포스필드 실패: {ff} - {message}")
        
        return False, "모든 포스필드 시도 실패"

# 사용 예시
def test_gromacs_runner():
    """GROMACS Runner 테스트"""
    import tempfile
    
    # 임시 디렉토리 생성
    with tempfile.TemporaryDirectory() as temp_dir:
        # 로거 설정
        logging.basicConfig(level=logging.INFO)
        logger = logging.getLogger("Test")
        
        # Runner 생성
        runner = GromacsCommandRunner(temp_dir, logger=logger, system_size="medium")
        
        # 단일 명령어 테스트 (MPI + GPU)
        params = {
            "prefix": "em",
            "mpi_ranks": "4",
            "ntomp": "2",
            "gpu_id": "3"
        }
        
        success, message = runner.execute_gromacs_command("mdrun_em", params)
        print(f"mdrun_em 결과: {success}, {message}")

if __name__ == "__main__":
    test_gromacs_runner()