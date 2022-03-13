# hash:sha256:567765db95850bdd92aff5a69251663fe019d5dd6ec5bc56f1be9a26948e22c7
FROM registry.acmecorp-edge.codeocean.com/codeocean/miniconda3:4.7.10-cuda10.1-cudnn7-ubuntu18.04

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        curl=7.58.0-2ubuntu3.16 \
        libjpeg-dev=8c-2ubuntu8 \
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
