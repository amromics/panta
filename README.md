# AMomics   


## Welcome to AMRomics

**AMRomics** is a software package to analyze microbial genomics data collections.
It provides a pipeline that bundles the current best practice for 
multiple aspects of AMR analyses. The pipeline analysis results can be 
represented and visualized via a web application. The web application also 
provides efficient data management.
 
AMRomics is written in python, it includes the followings dependencies:
 * blast (known to work with 2.10.1+)
 * samtools (1.11)
 * trimmomatic (0.39)
 * spades (3.15.2)
 * shovill (1.1.0)
 * prokka (1.14.6)
 * mlst (2.19.6)
 * abricate (1.0.1 | Database: vfdb ecoli_vf ecoh card megares resfinder argannot ncbi plasmidfinder)
 * roary (3.13.0) 
 * iqtree (2.1.2)

## Installation

The simplest method is installed via conda:

0. Download and install the appropriate conda, such as anaconda from 
   https://repo.anaconda.com/archive/
   
1. Create a conda environment with all the necessary dependencies: From the repository directory run

```bash
conda create -y -c conda-forge -c bioconda -c anaconda -c etetoolkit -c defaults --name amromics --file requirements.txt
```
2. Optional: install amromics library into conda environment
```bash
source activate amromics
pip install . --use-feature=in-tree-build 
``
## Usage


```bash
source activate amromics
```
### To run the pipeline


### Examples

We prepare a collection of dataset for public data base. To download the raw data,
```bash
cd examples/Kp24/raw
./download.sh
cd ../../
```
The following command will run that 24 samples through the pipeline, and import the results
to the web-app for visualization:

```bash
amr-analysis.py pg --time-log k24_time.log  -t 7 -m 25 -c KpClinicalGRBZ -i examples/Kp24/Kp24.tsv --work-dir data/work  -n "Collection of 24 clinical isolates from Greek and Brazil"
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
