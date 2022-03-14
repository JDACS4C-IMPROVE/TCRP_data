# hash:sha256:349b4a8fff50ff90539870f1249f2baeceb2fce0781c70209a78c91856b31e82
FROM registry.acmecorp-edge.codeocean.com/codeocean/miniconda3:4.7.10-cuda10.1-cudnn7-ubuntu18.04

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        curl=7.58.0-2ubuntu3.16 \
        libjpeg-dev=1.5.2-0ubuntu5.18.04.4 \
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

RUN mkdir /.code-server \
	&& cd /.code-server \
	&& curl -L "https://github.com/coder/code-server/releases/download/v4.0.2/code-server-4.0.2-linux-amd64.tar.gz" -o code-server.tar.gz \
	&& tar -xvf code-server.tar.gz \
	&& rm code-server.tar.gz \
	&& ln -s /.code-server/code-server-4.0.2-linux-amd64/bin/code-server  /usr/bin/code-server

COPY postInstall /
RUN /postInstall
