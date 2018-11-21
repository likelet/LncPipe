FROM nfcore/base
MAINTAINER Qi Zhao <zhaoqi@sysucc.org.cn>
LABEL authors="zhaoqi@sysucc.org.cn" \
    description="Docker image containing all requirements for the nfcore/lncpipe pipeline"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/nf-core-lncpipe-1.0dev/bin:$PATH