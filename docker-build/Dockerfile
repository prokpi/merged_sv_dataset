FROM python:3.9-slim
ENV DEBIAN_FRONTEND=noninteractive
WORKDIR /app

COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt
