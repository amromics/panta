
# AMR-Viz


[![CRAN](https://www.r-pkg.org/badges/version-ago/AMR)](https://cran.r-project.org/package=AMR)
[![CRANlogs](https://cranlogs.r-pkg.org/badges/grand-total/AMR)](https://cran.r-project.org/package=AMR)
![R-code-check](https://github.com/msberends/AMR/workflows/R-code-check/badge.svg?branch=master)
[![CodeFactor](https://www.codefactor.io/repository/github/msberends/amr/badge)](https://www.codefactor.io/repository/github/msberends/amr)
[![Codecov](https://codecov.io/gh/msberends/AMR/branch/master/graph/badge.svg)](https://codecov.io/gh/msberends/AMR?branch=master)


## Welcome to AMR-Viz

**AMR-Viz** is a package for genomics analysis antimicrobial resistant bacteria. 
The core of AMR-Viz is a pipeline that bundles the current best practice for 
multiple aspects of AMR analyses. The pipeline analysis results are 
represented and visualized via a web application. The web application also 
provides efficient data management.
 
AMR-Viz is written in python, and its web back-end is implemented with nodejs. 
It includes the followings dependencies:
 * blast (known to work with 2.10.1+)
 * samtools (1.11)
 * trimmomatic (0.39)
 * spades (3.14.1)
 * shovill (1.1.0)
 * prokka (1.14.6)
 * mlst (2.19.6)
 * abricate (1.0.1 | Database: vfdb ecoli_vf ecoh card megares resfinder argannot ncbi plasmidfinder)
 * roary (3.13.0)
 * parsnp (1.5.3)
 * iqtree (2.0.3)

## Installation

The simplest method is installed via conda:

0. Download and install the appropriate conda, such as anaconda from 
   https://repo.anaconda.com/archive/
   
1. Create a conda environment with all the necessary dependencies: 
```bash
conda create -c bioconda -c conda-forge -c anaconda --name amromics-viz \
      python=3.7 \
      ipykernel \
      numpy \
      pandas \
      tqdm \
      biopython \
      pysam \
      prokka \
      samtools \
      mlst \
      abricate \
      shovill \
      roary \
      iqtree \
      parsnp \
      nodejs
```
2. Setup nodejs
```bash
source activate amromics-viz
npm install -g live-server
```
3. Clone **amromics-vis** GitHub repository to your local computer:
Change the current working directory to the location where you want the 
cloned directory.
```bash
git clone https://github.com/amromics-org/amromics-vis.git
```
4. (Optional) Setup and build web application using npm 

```bash
cd amromics-vis
npm install
npm run build --modern
```

## Usage

Everytime you start using Amromics-viz, make sure you have activated
the *amromics-viz* conda environment. This can beIf not yet, activate the environment
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
./amrviz.py start [-p 3000] [--webapp-dir web-app]
```

The web application is auto opened on the URL **localhost:3000** (or another 
port if this port is occupied). 

### To run the pipeline


### Examples

We prepare a dataset for amrviz which from public data base. To download the raw data,
```bash
cd examples/Kp24/raw
./download.sh
cd ../../
```
The following command will run that 24 samples through the pipeline, and import the results
to the web-app for visualization:

```bash
./amrviz.py pa --time-log k24_time.log  -t 7 -m 25 -c KpClinicalGRBZ -i examples/Kp24/Kp24.tsv --work-dir data/work --webapp-dir web-app  -n "Collection of 24 clinical isolates from Greek and Brazil"
```



<!--

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
-t : threads number ((To-do: Change this to number of CPUs?))
-m: memory (in GB)
```

## Authors and Contributors

Amromics-viz is currently maintained by Minh Duc Cao, Hoang Anh Nguyen and Le Duc Quang. The following people (in alphatical order) have contributed to the development of Amromics-viz, including ideas, algorithms, implementation, documentation and feedback: ((To-do: 1-Add contributors if any. 2-Add team contact info))

## License

Amromics-viz is released under the accompanying BSD-like license.

((To-do: 1-Review license type. 2. Do we need a separate LICENSE.md file? 3-Add additional requirements for users regarding reference/acknowledgement/authorship etc. in publications))


-->
