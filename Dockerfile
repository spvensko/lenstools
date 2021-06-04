FROM ubuntu:xenial

RUN apt-get update \
  && apt-get install -y bash \
    python3 \
    python3-biopython \
    python3-pyvcf \
    python3-scipy \
  && rm -rf /var/lib/apt/lists/*

COPY lenstools.py /opt/lenstools/lenstools.py
COPY lenstools_w_final_report.py /opt/lenstools/lenstools_w_final_report.py
