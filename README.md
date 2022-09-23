# subpred2

## Note

*This is an older version of the project, and contains the notebooks and code used to calculate the results included in the first revision of the Manuscript.*

*Tested on Ubuntu 18.04*

## How to reproduce results:

1. Install miniconda
2. Recreate conda environment:
```
conda env create --file environment.yml
```
3. Activate conda environment: 
```
conda activate subpred
```
4. Download **data_full.tar** from https://cloud.hiz-saarland.de/s/sGTyGApAqdgAQiB and place in repo folder
5. Extract tar archives:
```
make raw_data
```
6. Create BLAST databases (Needs >100GB of space and 1-2 hours): 
    - This step is optional, as the tar file also contains pre-computed PSSMs for each protein in the dataset in **data/intermediate/blast**
```
make blast_databases
```




