# Amromics-viz

## Welcome to Amromics-viz

**Amromics-viz** is an analysis and visualization package for antimicrobial 
resistance (AMR) genome study. The core of Amromics-viz is a Python pipeline 
that is run via a simple terminal command line but performs a comprehensive 
multifaceted analysis of AMR genome data. The pipeline analysis results are 
represented and visualized via a web application. The web application also 
provides efficient data management as well as supports exporting of results 
and visualization displays ready for insertion in publications.

((To-do: Is there any requirements/assumptions regarding users? For example: 
do they need any special knowledge/background to be able to use the package?))

### Computational requirements: 

Amromics-viz requires Python ((To-do: Include version number)) to compile 
and run. ((To-do: Include other requirements))

## Installation

1. Download and install the appropriate anaconda from 
https://repo.anaconda.com/archive/
2. Create a conda environment with all the neccesary dependencies to run 
amromics-viz:
```bash
conda create -c bioconda -c conda-forge -c anaconda --name amromics-viz python=3.7 ipykernel numpy pandas biopython prokka pysam samtools mlst abricate snippy tqdm shovill roary parsnp nodejs
```
3. Setup nodejs
```bash
source activate amromics-viz
npm install -g live-server
```
4. Clone **amromics-vis** GitHub repository to your local computer:
Change the current working directory to the location where you want the 
cloned directory.
```bash
git clone https://github.com/amromics-org/amromics-vis.git
```
5. (Optional) Setup and build web application using npm 
((To-do: Document what this is for.))
```bash
cd amromics-vis
npm install
npm run build --modern
```

## Usage

Everytime you start using Amromics-viz, make sure you have activated
the *amromics-viz* conda environment. If not yet, activate the environment
with the following command:
```bash
source activate amromics-viz
```

Amromics-viz comprises two components: a web application and an analysis 
pipeline. 

### To run the web application
Change the current working directory to the **amromics-vis** cloned directory 
in the step 4 of Intallation above.
```bash
cd web-app && live-server --port=3000  --entry-file=index.html
```

The web application is auto opened on the URL **localhost:3000** (or another 
port if this port is occupied). ((To-do: Modify this to report to user if
this port is occupied and instruct user to choose a different port.))


### To run the pipeline
#### Prepare input file
- Data file inputted for analysis needs to be in *.tsv* format 
((To-do: Check if .tsv format is required)) and follows specific requirements. 
Please check the sample input file *data/samples/set1.tsv* for an example.
- Note:
  + Column names need to be as follow:
    - sample_id	
    - sample_name	
    - input_type	
    - files	
    - genus	
    - species	
    - strain	
    - gram	
    - metadata
  + *gram* column should be empty. ((To-do: Delete gram column?))
  + *metadata* is empty or in the format: key1:value1;key2:value2;...  
  For example: Geographic Location:Houston,USA;Insert Date:8/8/2017;Host Name:Human, Homo sapiens;ampicillin:Resistant;aztreonam:Resistant;ciprofloxacin:Resistant;gentamicin:Susceptible;tetracycline:Susceptible
#### Run pipeline and export visualization data to web application
```bash
python amromics-viz.py --id sample1 --input data/samples/set1.tsv -t 2 -m 16
```
```bash
--id : Collection ID, a string without space or special character
--input: the file listing samples of collection
-t : threads number
-m: memory (in GB)
```

## Authors and Contributors

Amromics-viz is currently maintained by Minh Duc Cao, Hoang Anh Nguyen and Le Duc Quang. The following people (in alphatical order) have contributed to the development of Amromics-viz, including ideas, algorithms, implementation, documentation and feedback: ((To-do: 1-Add contributors if any. 2-Add team contact info))

## License

Amromics-viz is released under the accompanying BSD-like license.

((To-do: 1-Review license type. 2. Do we need a separate LICENSE.md file? 3-Add additional requirements for users regarding reference/acknowledgement/authorship etc. in publications))


