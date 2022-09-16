The Makefiles included here download the most recent versions of the raw datasets.
For the analyses in the manuscript and in the notebooks, we used a dataset that was downloaded on November 15, 2021. More recent datasets can contain a slightly different set of proteins and annotations, which can influence the evaluation results and the dataset analysis.
If you want to reproduce or validate the results from the manuscript, follow the section *Reproduce results from manuscript* below.

## Setup:

1. Install miniconda
2. Recreate conda environment:
```
conda env create --file environment.yml
```
3. Activate conda environment: 
```
conda activate subpred
```
4. Install code as python package into environment: 
```
pip install -e .
```
5. Download raw transporter data: 
```
make raw_data
```
6. Create BLAST databases (Needs >100GB of space and takes several hours): 
```
make blast_databases
```

## Reproduce results from manuscript:

1. Install miniconda
2. Recreate conda environment:
```
conda env create --file environment.yml
```
3. Activate conda environment: 
```
conda activate subpred
```
4. Install code as python package into environment: 
```
pip install -e .
```
5. Download data_full.tar from https://cloud.hiz-saarland.de/s/sGTyGApAqdgAQiB
6. Rename existing data folder:
```
mv data data_bak
```
7. Extract tar archive:
```
make raw_data_manuscript
```
8. Create BLAST databases (Needs >100GB of space and takes several hours):
    - This step is optional, as the previous step extracts pre-computed PSSMs for all proteins to *data/intermediate/blast*
  
```
make blast_databases
```
