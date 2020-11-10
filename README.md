# amromics-visualization
## Introduction
## Installation
### Using Docker image (comming soon)
The easiest way to run Amromics-viz is to use the predefined docker image.
### Install from source
#### Set up pipeline
```bash
conda install -c conda-forge -c bioconda -c defaults prokka pysam samtools mlst abricate roary
```
- Download and extract Parsnp
```bash
wget https://github.com/marbl/parsnp/releases/download/v1.2/parsnp-Linux64-v1.2.tar.gz && \
    tar -xvf parsnp-Linux64-v1.2.tar.gz && rm parsnp-Linux64-v1.2.tar.gz
```
- Add parsnp to PATH 
- Install requirement packs using pip
```bash
pip install -r requirements.txt
```
#### Set up Amromics-viz web application
- Install NodeJS (v10.20 or later)
```bash
sudo apt install nodejs
```
- Install live-server
```bash
npm install -g live-server
```
- Setup and build web application using npm (option for developing)
```bash
cd amromics-vis
```
```bash
npm install
```
```bash
npm run build --modern
```
## Usage
### Run from docker image (Comming soon)

### Run from source
#### To run pipeline
##### Prepare input file
- Check sample input file in data/samples/set1.tsv
- Note:
  + Gram column should be empty
  + Metadata is empty or in format: key1:value1;key2:value2;...  Ex: Geographic Location:Houston,USA;Insert Date:8/8/2017;Host Name:Human, Homo sapiens;ampicillin:Resistant;aztreonam:Resistant;ciprofloxacin:Resistant;gentamicin:Susceptible;tetracycline:Susceptible
##### Run pipeline and visualize tools
```bash
python amromics-viz.py --id sample1 --input data/samples/set1.tsv -t 2 -m 16
```
```bash
--id : Collection ID, a string without space or special character
--input: the file listing samples of collection
--t : threads number
--m: memory (in GB)
```
The web application is auto opened on URL localhost:3000 (or other port if this port is occupated)
