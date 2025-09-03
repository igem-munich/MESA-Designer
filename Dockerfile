FROM ubuntu:24.04
LABEL authors="Aeneas Tews"

WORKDIR /mesa_designer

RUN apt-get update && apt-get install -y \
    build-essential \
    software-properties-common \
    python3.12 \
    python3-pip \
    curl \
    git \
    && rm -rf /var/lib/apt/lists/*

RUN git clone https://github.com/igem-munich/MESA-Designer.git .

RUN pip install -r requirements.txt --break-system-packages

RUN python3 ./setup.py

EXPOSE 8501

HEALTHCHECK CMD curl --fail http://localhost:8501/_stcore/health

ENTRYPOINT ["streamlit", "run", "app/main.py", "--server.port=8501", "--server.address=0.0.0.0"]