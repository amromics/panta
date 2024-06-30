FROM continuumio/miniconda3:4.12.0
#Metadata
LABEL \
	image.name="panta" \
	image.version="1.0.1" \
	image.description="" \
	maintainer="amromics.org" 

# Set Environment Variables
#ENV PATH=$PATH:/opt/conda/bin
RUN apt-get update
RUN apt-get install -y parallel make cmake wget git locales

# Install amromics via micromamba
WORKDIR /tmp
RUN conda install -y -c conda-forge python=3.10 mamba
RUN git clone https://github.com/amromics/panta.git && \
	cd panta && \
    mamba install -y -c conda-forge -c bioconda -c anaconda -c defaults   --file requirements.txt && \
    pip install . 

WORKDIR /tmp/panta
