FROM python:3.12-slim
LABEL authors="Aeneas Tews"

WORKDIR /mesa_designer

RUN apt-get update && apt-get install -y \
    curl \
    git \
    && rm -rf /var/lib/apt/lists/*

RUN git clone https://github.com/igem-munich/MESA-Designer.git .

RUN pip install -r requirements.txt

RUN python3 ./setup.py

EXPOSE 8501

HEALTHCHECK CMD curl --fail http://localhost:8501/_stcore/health

ENTRYPOINT ["streamlit", "run", "app/main.py", "--server.port=8501", "--server.address=0.0.0.0"]