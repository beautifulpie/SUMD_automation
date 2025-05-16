#!/bin/bash
# test_environment.sh - SuMD Gromacs 테스트 환경 설정 및 실행 스크립트

set -e  # 오류 발생 시 스크립트 중단

# 실행 권한 부여
chmod +x scripts/*.py scripts/*.sh

# 테스트 데이터 다운로드
echo "테스트 데이터 다운로드 중..."
mkdir -p test_data
wget -O test_data/example.pdb "https://files.rcsb.org/download/1BRS.pdb" || echo "테스트 PDB 파일 다운로드 실패. 직접 PDB 파일을 test_data 디렉토리에 추가하세요."

# Docker 이미지 빌드 및 실행
echo "Docker 이미지 빌드 및 실행 중..."
docker-compose up -d --build

echo "테스트 환경 구축 완료!"
echo "컨테이너에 접속하려면 다음 명령어를 사용하세요:"
echo "docker exec -it sumd-gromacs bash"
echo ""
echo "시스테인 처리가 포함된 테스트 실행 예시:"
echo "docker exec -it sumd-gromacs /app/scripts/run_sumd_gromacs.sh /app/test_data/example.pdb \"45 46 47\" \"123 124 125\" /app/output 0.5 true"
