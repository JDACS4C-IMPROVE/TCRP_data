# hash:sha256:d056c783cb53f922cff1854c35b0fb21c41e5e9b658313b86d0738da253039e3
FROM registry.acmecorp-edge.codeocean.com/codeocean/miniconda3:4.7.10-cuda10.1-cudnn7-ubuntu18.04

ARG DEBIAN_FRONTEND=noninteractive

RUN conda install -y --channel=conda-forge \
        jupyterlab==3.2.0 \
    && conda clean -ya

RUN pip install -U --no-cache-dir \
    scipy==1.7.3 \
    sklearn==0.0

COPY postInstall /
RUN /postInstall
