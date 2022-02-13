# hash:sha256:38496b5da59012becc98baf194985b2ae5e851dd4d83cf36dfab0d57f86454c7
FROM registry.acmecorp-edge.codeocean.com/codeocean/miniconda3:4.7.10-cuda10.1-cudnn7-ubuntu18.04

ARG DEBIAN_FRONTEND=noninteractive

RUN conda install -y \
        jupyterlab==3.2.1 \
    && conda clean -ya

RUN pip install -U --no-cache-dir \
    scipy

COPY postInstall /
RUN /postInstall
