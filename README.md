Note: Currently work in progress. The code is being migrated from an older version of the subpred package. The old versions of the code and the notebooks which were used for the manuscript can be found on the branch *subpred2*. 

Setup:

1. Install miniconda
2. Recreate conda environment:
```
conda env create --file environment.yml
```
3. Activate conda environment: 
```
conda activate subpred
```
4. Install code as python package: 
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


The Makefiles included here download the most recent versions of the raw datasets. 
For the analysis in the manuscript, we used a dataset downloaded on the 15th of November of 2021.
If you want to reproduce or validate the results, we will send you a link to the data we used. 
