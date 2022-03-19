# hash:sha256:2fe5c07ee137bfbf1b8db2f03e93e888e86bee55c5fa6dfd616bac3bc4582682
FROM registry.acmecorp-edge.codeocean.com/codeocean/miniconda3:4.7.10-cuda10.1-cudnn7-ubuntu18.04

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        curl \
        libjpeg-dev \
    && rm -rf /var/lib/apt/lists/*

RUN conda install -y --channel=conda-forge \
        jupyterlab==3.2.0 \
    && conda clean -ya

RUN pip install -U --no-cache-dir \
    matplotlib==3.5.1 \
    mygene==3.2.2 \
    networkx==2.6.3 \
    pandas==1.3.5 \
    scipy==1.7.3 \
    seaborn==0.11.2 \
    sklearn==0.0

COPY postInstall /
RUN /postInstall
