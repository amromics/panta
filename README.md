# amromics-visualization
## To run visualization
### Prerequisites
- NodeJS
- npm
### Run visualization
cd amromics-vis
npm run dev
## To run pipeline
### Prerequisites
- Python
- Spades
- Skesa
- Prokka
- abricate
- mlst
- Roary
- Parsnp
- Python pack: Biopython, json, csv
### Prepare input file
- Check sample input file in /samples/set1.tsv
- Note:
  + Gram column should be empty
  + Metadata is empty or in format: key1:value1;key2:value2;...  Ex: Geographic Location:Houston,USA;Insert Date:8/8/2017;Host Name:Human, Homo sapiens;ampicillin:Resistant;aztreonam:Resistant;ciprofloxacin:Resistant;gentamicin:Susceptible;tetracycline:Susceptible
### Run pipeline.py
python pipeline.py pa -i samples/set1.tsv -t 2 -m 8 -e amromics-vis/static/data 
