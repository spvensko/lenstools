FROM ubuntu:16.04

ENV PATH /opt/conda/bin:$PATH

RUN apt-get update --fix-missing &&  \
    apt-get install -qy wget curl git bzip2 ca-certificates procps zlib1g-dev \
    make build-essential cmake libncurses-dev ncurses-dev g++ gcc \
    nfs-common pigz bedtools gawk fuse mdadm time \
    libbz2-dev lzma-dev liblzma-dev libglib2.0-0 libxext6 libsm6 libxrender1 \
    syslog-ng libssl-dev libtool autoconf automake \
    libcurl4-openssl-dev libffi-dev libblas-dev liblapack-dev \
    libatlas-base-dev libroot-math-mathmore-dev

COPY lenstools.py /opt/lenstools/lenstools.py
COPY lenstools_w_final_report.py /opt/lenstools/lenstools_w_final_report.py
RUN chmod +x /opt/lenstools/*

RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda2-4.5.11-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh && \
    ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
    echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc && \
    echo "conda activate base" >> ~/.bashrc

ENV HTSLIB_LIBRARY_DIR /usr/local/lib
ENV HTSLIB_INCLUDE_DIR /usr/local/include
ENV LD_LIBRARY_PATH /usr/local/lib

RUN conda update conda
RUN conda create -n run-env bioconda::pyvcf bioconda::biopython anaconda::numpy anaconda::scipy bioconda::pysam anaconda::pandas
RUN echo "source activate run-env" > ~/.bashrc
ENV PATH /opt/conda/envs/run-env/bin:$PATH
