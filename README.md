# amromics-visualization

## Installation

### Install with conda

The easiest method is to create a conda environment with all the neccesary dependencies.

1. Download and install the appropriate anaconda from https://repo.anaconda.com/archive/
2. Create an environment to run amromics-viz:
```bash
conda create -c bioconda -c conda-forge -c anaconda --name amromics-viz python=3.7 ipykernel numpy pandas biopython prokka pysam samtools mlst abricate snippy tqdm shovill roary parsnp nodejs
```
3. Setup nodejs
```bash
soure activate amromics-viz
npm install -g live-server
```
4. (Optional) Setup and build web application using npm
```bash
cd amromics-vis
npm install
npm run build --modern
```

## Usage

#### To run pipeline
##### Prepare input file
- Check sample input file in data/samples/set1.tsv
- Note:
  + Gram column should be empty
  + Metadata is empty or in format: key1:value1;key2:value2;...  Ex: Geographic Location:Houston,USA;Insert Date:8/8/2017;Host Name:Human, Homo sapiens;ampicillin:Resistant;aztreonam:Resistant;ciprofloxacin:Resistant;gentamicin:Susceptible;tetracycline:Susceptible
##### Run pipeline and export visualization data to web app
```bash
python amromics-viz.py --id sample1 --input data/samples/set1.tsv -t 2 -m 16
```
```bash
--id : Collection ID, a string without space or special character
--input: the file listing samples of collection
--t : threads number
--m: memory (in GB)
```
##### Run web app
```bash
cd web-app && live-server --port=3000  --entry-file=index.html
```
```bash
The web application is auto opened on URL localhost:3000 (or other port if this port is occupated)
