FROM continuumio/miniconda3

WORKDIR /pan-genome

COPY . /pan-genome

RUN apt-get update &&  apt-get install -y time
RUN conda update --all
RUN conda install -y -c conda-forge python=3.7 mamba
RUN mamba install -y -c conda-forge -c bioconda -c anaconda -c defaults  --file /pan-genome/requirements.txt

ENV PATH="/pan-genome:${PATH}" PYTHONPATH="/pan-genome:$PYTHONPATH"
