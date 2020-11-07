# amromics-visualization
## Introduction
## Installation
### Using Docker image
The easiest way to run Amromics-viz is to use the predefined docker image.
### Install from source
#### Set up pipeline
```bash
conda install -c conda-forge -c bioconda -c defaults prokka pysam samtools mlst abricate roary
```
Download and extract Parsnp
```bash
wget https://github.com/marbl/parsnp/releases/download/v1.2/parsnp-Linux64-v1.2.tar.gz && \
    tar -xvf parsnp-Linux64-v1.2.tar.gz && rm parsnp-Linux64-v1.2.tar.gz
```
Add parsnp to PATH 
Install requirement packs using pip
```bash
pip install -r requirements.txt
```
#### Set up Amromics-viz web application
Install NodeJS (v10.20 or later)
```bash
sudo apt install nodejs
```
Setup web application using npm
```bash
cd amromics-vis
```
```bash
npm install
```
## Usage
### Run from docker image

### Run from source
#### To run pipeline
##### Prepare input file
- Check sample input file in /samples/set1.tsv
- Note:
  + Gram column should be empty
  + Metadata is empty or in format: key1:value1;key2:value2;...  Ex: Geographic Location:Houston,USA;Insert Date:8/8/2017;Host Name:Human, Homo sapiens;ampicillin:Resistant;aztreonam:Resistant;ciprofloxacin:Resistant;gentamicin:Susceptible;tetracycline:Susceptible
##### Run pipeline.py
python pipeline.py pa -i samples/set1.tsv -t 2 -m 8 -e amromics-vis/static/data 
#### To visualize
```bash
cd amromics-vis
```
```bash
 npm run dev
```
