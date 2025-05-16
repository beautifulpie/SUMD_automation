FROM kcaladram/gromacs-new:latest

# 기본 패키지 설치
RUN apt-get update && apt-get upgrade -y && \
    apt-get install -y nano python3-pip git wget && \
    pip3 install numpy biopython pandas matplotlib tqdm

# 작업 디렉토리 생성
WORKDIR /app

# 필요한 디렉토리 생성
RUN mkdir -p /app/mdp_templates /app/scripts /app/test_data /app/output

# MDP 템플릿 복사
COPY mdp_templates/ /app/mdp_templates/

# 스크립트 복사
COPY scripts/ /app/scripts/

# 실행 권한 부여
RUN chmod +x /app/scripts/*.py /app/scripts/*.sh

# 환경 변수 설정
ENV PATH="/app/scripts:${PATH}"

# 기본 명령어 설정
CMD ["/bin/bash"]
