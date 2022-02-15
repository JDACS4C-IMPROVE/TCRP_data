# hash:sha256:ff61e49be6c4d0a28e23761eb875b88ac420f0bfa3a1995c0115aced56a0069d
FROM registry.acmecorp-edge.codeocean.com/codeocean/miniconda3:4.7.10-cuda10.1-cudnn7-ubuntu18.04

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        libjpeg-dev=8c-2ubuntu8 \
    && rm -rf /var/lib/apt/lists/*

RUN conda install -y --channel=conda-forge \
        jupyterlab==3.2.0 \
    && conda clean -ya

RUN pip install -U --no-cache-dir \
    matplotlib \
    mygene \
    networkx \
    pandas \
    scipy==1.7.3 \
    seaborn \
    sklearn==0.0

COPY postInstall /
RUN /postInstall
