FROM python:3.13.7-slim
LABEL authors="Aeneas Tews"

WORKDIR /mesa_designer

RUN apt update && apt upgrade -y && apt install -y \
    curl \
    git \
    dos2unix \
    && rm -rf /var/lib/apt/lists/*

RUN git clone https://github.com/igem-munich/MESA-Designer.git .

RUN dos2unix ./docker_entrypoint.sh

RUN chmod +x ./docker_entrypoint.sh

RUN pip install -r requirements.txt

RUN python3 ./setup.py

EXPOSE 8501
EXPOSE 8000

HEALTHCHECK CMD curl --fail http://localhost:8501/_stcore/health

ENTRYPOINT ["./docker_entrypoint.sh"]