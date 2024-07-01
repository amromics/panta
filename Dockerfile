FROM mambaorg/micromamba:latest
#Metadata
LABEL \
	image.name="panta" \
	image.version="1.0.1" \
	image.description="" \
	maintainer="amromics.org" 

# Set Environment Variables
ENV ENV_NAME=panta
USER root
RUN apt update -y && apt install -y --no-install-recommends \
    time \
    && \
    apt-get clean && \
    rm -rf /var/cache/apt/* /var/lib/apt/lists/*
#USER mambauser
# Install amromics via micromamba
WORKDIR /tmp

RUN micromamba create -y -c conda-forge -c defaults --name ${ENV_NAME} python=3.10 git && \
    eval "$(micromamba shell hook --shell bash)" && \
	  micromamba activate ${ENV_NAME} && \
    git clone https://github.com/amromics/panta.git && \
	  cd panta && \
    micromamba install -y -c conda-forge -c bioconda -c anaconda -c defaults  --file requirements.txt && \
    pip install . && \
    cd .. && rm -rf panta && \
    micromamba list && \
	  micromamba clean -afy && \
    rm -rf .cache && \
    find /opt/conda/envs/panta/ -follow -type f -name '*.a' -delete && \
    find /opt/conda/envs/panta/ -follow -type f -name '*.pyc' -delete && \
    find /opt/conda/envs/panta/ -follow -type f -name '*.js.map' -delete

# Container Startup
ENTRYPOINT ["/usr/local/bin/_entrypoint.sh"]
# Container Default Command
CMD [ "/bin/bash", "-c", "panta --help" ]
