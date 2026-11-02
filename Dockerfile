FROM python:3.12-slim

WORKDIR /app

COPY requirements.txt .

RUN pip install --no-cache-dir -r requirements.txt

COPY src/ ./src/

COPY references/ ./references/

RUN cat references/gencode.part-* > references/gencode.v50.annotation.gtf.gz \
    && rm references/gencode.part-*

ENTRYPOINT ["python", "-m", "src.cli"]
