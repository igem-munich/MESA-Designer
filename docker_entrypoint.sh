#!/bin/sh
streamlit run app/main.py --server.port=8501 --server.address=0.0.0.0 &
uvicorn api.main:app --host 0.0.0.0