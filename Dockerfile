FROM mambaorg/micromamba:latest
#Metadata
LABEL \
	image.name="panta" \
	image.version="1.0.1" \
	image.description="" \
	maintainer="amromics.org" 

# Set Environment Variables
#ENV PATH=$PATH:/opt/conda/bin
USER root
RUN apt update -y && apt install -y --no-install-recommends \
    time \
    && \
    apt-get clean && \
    rm -rf /var/cache/apt/* /var/lib/apt/lists/*
USER mambauser
# Install amromics via micromamba
WORKDIR /tmp

RUN micromamba create -y -c conda-forge -c defaults --name panta python=3.10 git && \
    eval "$(micromamba shell hook --shell bash)" && \
	micromamba activate panta && \
    git clone https://github.com/amromics/panta.git && \
	cd panta && \
    micromamba install -y -c conda-forge -c bioconda -c anaconda -c defaults  --file requirements.txt && \
    pip install . && \
    micromamba list && \
	micromamba clean -afy

WORKDIR /tmp/panta
