FROM mambaorg/micromamba:2.9.0-ubuntu22.04

COPY --chown=$MAMBA_USER:$MAMBA_USER envs/environment.yml /tmp/environment.yml
RUN micromamba install --yes --name base --file /tmp/environment.yml && \
    micromamba clean --all --yes

COPY --chown=$MAMBA_USER:$MAMBA_USER R /opt/celltype/R
COPY --chown=$MAMBA_USER:$MAMBA_USER scripts /opt/celltype/scripts
ENV PATH=/opt/conda/bin:$PATH
ENV R_LIBS_USER=/opt/conda/lib/R/library

ARG MAMBA_DOCKERFILE_ACTIVATE=1
RUN Rscript -e 'stopifnot(as.character(packageVersion("dtangle")) == "2.0.10")' && \
    Rscript -e 'stopifnot(as.character(packageVersion("TCA")) == "1.2.1")'

WORKDIR /workspace
