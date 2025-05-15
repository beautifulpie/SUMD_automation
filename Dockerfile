FROM gromacs-gpu

# 비대화형 모드 설정
ENV DEBIAN_FRONTEND=noninteractive

# 필요한 추가 패키지 설치
RUN apt-get update && apt-get install -y \
    python3 \
    python3-pip \
    python3-numpy \
    python3-matplotlib \
    python3-biopython \
    git \
    wget \
    curl \
    vim \
    && rm -rf /var/lib/apt/lists/*

# Python 패키지 설치
RUN pip3 install --no-cache-dir tqdm

# 작업 디렉토리 설정
WORKDIR /sumd

# GROMACS 환경 활성화를 위한 스크립트 생성
RUN echo '#!/bin/bash\n\
if [ -f "/usr/local/gromacs/bin/GMXRC" ]; then\n\
    source /usr/local/gromacs/bin/GMXRC\n\
elif [ -f "/usr/bin/GMXRC" ]; then\n\
    source /usr/bin/GMXRC\n\
else\n\
    echo "경고: GROMACS 환경 설정 파일을 찾을 수 없습니다."\n\
fi\n\
exec "$@"' > /sumd/entrypoint.sh && \
    chmod +x /sumd/entrypoint.sh

# SuMD 스크립트 복사
COPY SuMD_script.py /sumd/
COPY Running_SuMD_Automation.py /sumd/
COPY README.md /sumd/

# 필요한 디렉토리 생성
RUN mkdir -p /sumd/data /sumd/results /sumd/forcefield

# 실행 권한 부여
RUN chmod +x /sumd/SuMD_script.py /sumd/Running_SuMD_Automation.py

# 환경 변수 설정
ENV PYTHONPATH="/sumd:${PYTHONPATH:-}"
ENV PATH="/sumd:${PATH}"

# PDB 파일 처리를 위한 유틸리티 스크립트 생성
RUN echo '#!/bin/bash\n\
# 원본 PDB 파일\n\
input_pdb="${1:-input.pdb}"\n\
output_pdb="${2:-protein_only.pdb}"\n\
\n\
# 표준 아미노산만 유지\n\
grep "^ATOM\\|^TER\\|^END\\|^CRYST1\\|^HEADER\\|^TITLE" "$input_pdb" > "$output_pdb"\n\
\n\
echo "단백질만 포함된 PDB 파일이 생성되었습니다: $output_pdb"\n\
' > /sumd/clean_pdb.sh && \
    chmod +x /sumd/clean_pdb.sh

# 엔트리포인트 설정
ENTRYPOINT ["/sumd/entrypoint.sh"]
CMD ["bash"]