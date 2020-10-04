# amromics-visualization
## To run visualization
### Prerequisites
- NodeJS
- npm
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
- Check sample input file in /samples/set1.csv
### Run pipeline.py
python pipeline.py pa -i samples/set1.csv -t 2 -m 8 -e amromics-vis/static/data 
