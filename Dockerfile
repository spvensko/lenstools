FROM ubuntu:xenial

RUN apt-get update \
  && apt-get install -y bash \
    python3 \
    python3-biopython \
    python3-pyvcf \
    python3-scipy \
  && rm -rf /var/lib/apt/lists/*

COPY lenstools.py lentools.py
