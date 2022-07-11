Note: The code is being migrated from an older version of the subpred package. The old versions of the code and the notebooks, which were used to calculate the results included in the manuscript, can be found on the branch *subpred2*. These old notebooks are also located in the folder *notebooks_old* on the *main* branch.

The Makefiles included here download the most recent versions of the raw datasets.
For the analyses in the manuscript and in the notebooks, we used a dataset that was downloaded on November 15, 2021. More recent datasets can contain a slightly different set of proteins and annotations, which can influence the evaluation results and the dataset analysis.
If you want to reproduce or validate the results, we will send you a link to the exact datasets we used. 

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
