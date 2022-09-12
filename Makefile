# .PHONY: clean lint format requirements blastdb_uniref50 blastdb_uniref90 dataset expression_data

#################################################################################
# CONVENTIONS                                                                   #
#################################################################################

# SUBSTRATES = "Amino-acid transport;Sugar transport"
# DATASETNAME = "swissprot_transporters"
# TODO organism, other keywords, parameters etc.
# Or: JSON file with parameters?

#################################################################################
# COMMANDS                                                                      #
#################################################################################

# "human": 		9606
# "athaliana":	3702
# "ecoli": 		83333
# "yeast": 		559292

#################################################################################
# Setup                                                                         #
#################################################################################

## Install packages required on Ubuntu 20.04 LTS on WSL
setup_ubuntu:
	sudo apt install build-essential
	wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh
	bash ~/miniconda.sh -p ~/miniconda3
	rm ~/miniconda.sh
	@echo Reload shell to put conda in path: source ~/.bashrc

## Install packages.
requirements:
	conda update -n base -c defaults conda
	conda install -n base -c conda-forge mamba
	mamba env create --file environment.yml
	
# conda install numpy pandas scikit-learn imbalanced-learn xgboost matplotlib seaborn yellowbrick black flake8 ipykernel joblib scikit-learn-intelex qnorm fastcluster


#################################################################################
# Raw data                                                                      #
#################################################################################

## Extract raw data
raw_data: data_full.tar
	tar xvf data_full.tar
	tar xf data/intermediate/blast.tar.xz --directory=data/intermediate/blast
	rm data/intermediate/blast.tar.xz

#################################################################################
# Raw data: BLAST databases                                                     #
#################################################################################

## Extract and init blast dbs for recalculating PSSM feature. >100GB needed
blast_databases: blastdb_extract blastdb_uniref50 blastdb_uniref90

# blastdb_extract: subpred_data_uniref.tar
#	tar xvf subpred_data_uniref.tar
#	rm subpred_data_uniref.tar

blastdb_uniref50: 
	@echo Creating BLASTDB...
	$(MAKE) -C data/raw/uniref/uniref50 uniref50.fasta.pdb

blastdb_uniref90: 
	@echo Creating BLASTDB...
	$(MAKE) -C data/raw/uniref/uniref90 uniref90.fasta.pdb

#################################################################################
# Preprocessing                                                                 #
#################################################################################

## Preprocess raw data (data/intermediate)
preprocessing: preprocessing_expression_data preprocessing_go_data preprocessing_gene_positions

preprocessing_expression_data:
	# python3 src/preprocessing/athaliana_exp.py --aggregate -e data/raw/gene_expression/athaliana/E-TABM-17-processed-data-1343490029.txt.gz -g data/raw/gene_expression/athaliana/GPL198.annot.gz -m data/raw/uniprot_swissprot/ARATH_3702_idmapping.dat.gz -a data/raw/gene_expression/athaliana/E-TABM-17.sdrf.txt -o data/intermediate/gene_expression/athaliana/e_tabm_17.tsv
	python3 src/preprocessing/athaliana_exp.py --aggregate -e data/raw/gene_expression/athaliana/E-TABM-17-processed-data-1343490029.txt.gz -g data/raw/gene_expression/athaliana/GPL198.annot.gz -m data/raw/uniprot_swissprot/ARATH_3702_idmapping.dat.gz -a data/raw/gene_expression/athaliana/E-TABM-17.sdrf.txt --ecotype Columbia-0 -o data/intermediate/gene_expression/athaliana/athaliana_columbia.tsv
	python3 src/preprocessing/athaliana_exp.py --aggregate -e data/raw/gene_expression/athaliana/E-TABM-17-processed-data-1343490029.txt.gz -g data/raw/gene_expression/athaliana/GPL198.annot.gz -m data/raw/uniprot_swissprot/ARATH_3702_idmapping.dat.gz -a data/raw/gene_expression/athaliana/E-TABM-17.sdrf.txt --ecotype Columbia-0 --organism-parts flower -o data/intermediate/gene_expression/athaliana/athaliana_columbia_flower.tsv
	python3 src/preprocessing/athaliana_exp.py --aggregate -e data/raw/gene_expression/athaliana/E-TABM-17-processed-data-1343490029.txt.gz -g data/raw/gene_expression/athaliana/GPL198.annot.gz -m data/raw/uniprot_swissprot/ARATH_3702_idmapping.dat.gz -a data/raw/gene_expression/athaliana/E-TABM-17.sdrf.txt --ecotype Columbia-0 --organism-parts root -o data/intermediate/gene_expression/athaliana/athaliana_columbia_root.tsv
	python3 src/preprocessing/athaliana_exp.py --aggregate -e data/raw/gene_expression/athaliana/E-TABM-17-processed-data-1343490029.txt.gz -g data/raw/gene_expression/athaliana/GPL198.annot.gz -m data/raw/uniprot_swissprot/ARATH_3702_idmapping.dat.gz -a data/raw/gene_expression/athaliana/E-TABM-17.sdrf.txt --ecotype Columbia-0 --organism-parts seed -o data/intermediate/gene_expression/athaliana/athaliana_columbia_seed.tsv
	python3 src/preprocessing/athaliana_exp.py --aggregate -e data/raw/gene_expression/athaliana/E-TABM-17-processed-data-1343490029.txt.gz -g data/raw/gene_expression/athaliana/GPL198.annot.gz -m data/raw/uniprot_swissprot/ARATH_3702_idmapping.dat.gz -a data/raw/gene_expression/athaliana/E-TABM-17.sdrf.txt --ecotype Columbia-0 --organism-parts "rosette leaf" -o data/intermediate/gene_expression/athaliana/athaliana_columbia_rosette_leaf.tsv
	python3 src/preprocessing/athaliana_exp.py --aggregate -e data/raw/gene_expression/athaliana/E-TABM-17-processed-data-1343490029.txt.gz -g data/raw/gene_expression/athaliana/GPL198.annot.gz -m data/raw/uniprot_swissprot/ARATH_3702_idmapping.dat.gz -a data/raw/gene_expression/athaliana/E-TABM-17.sdrf.txt --ecotype Columbia-0 --organism-parts "rosette leaf" root flower seed -o data/intermediate/gene_expression/athaliana/athaliana_columbia_top4.tsv
	# TODO no more median aggregation, check if normalized already
	# python3 src/preprocessing/preprocessing_geo.py --remove-zero-var -i data/raw/gene_expression/ecoli/GDS680.soft.gz -m data/raw/uniprot_swissprot/ECOLI_83333_idmapping.dat.gz -o data/intermediate/gene_expression/ecoli/GDS680.tsv
	# python3 src/preprocessing/preprocessing_geo.py --remove-zero-var --quantile-normalize -i data/raw/gene_expression/ecoli/GDS680.soft.gz -m data/raw/uniprot_swissprot/ECOLI_83333_idmapping.dat.gz -o data/intermediate/gene_expression/ecoli/GDS680_norm.tsv
	# python3 src/preprocessing/preprocessing_geo.py --remove-zero-var --exp 2 -i data/raw/gene_expression/yeast/GDS91.soft.gz -m data/raw/uniprot_swissprot/YEAST_559292_idmapping.dat.gz -o data/intermediate/gene_expression/yeast/GDS91.tsv
	# python3 src/preprocessing/preprocessing_geo.py --remove-zero-var --exp 2 --quantile-normalize -i data/raw/gene_expression/yeast/GDS91.soft.gz -m data/raw/uniprot_swissprot/YEAST_559292_idmapping.dat.gz -o data/intermediate/gene_expression/yeast/GDS91_norm.tsv
	# python3 src/preprocessing/preprocessing_gtex_parallel.py -i data/raw/gene_expression/human -o data/intermediate/gene_expression/human

preprocessing_go_data:
	python3 src/preprocessing/gene_ontology.py -i data/raw/gene_ontology/goa_human.gaf.gz -o data/intermediate/gene_ontology/goa_human.tsv 
	python3 src/preprocessing/gene_ontology.py -i data/raw/gene_ontology/18.E_coli_MG1655.goa -o data/intermediate/gene_ontology/goa_ecoli.tsv 
	python3 src/preprocessing/gene_ontology.py -i data/raw/gene_ontology/goa_arabidopsis.gaf.gz -o data/intermediate/gene_ontology/goa_athaliana.tsv
	python3 src/preprocessing/gene_ontology.py -i data/raw/gene_ontology/goa_yeast.gaf.gz -o data/intermediate/gene_ontology/goa_yeast.tsv

preprocessing_gene_positions:
	python3 src/preprocessing/gene_positions.py -i data/raw/gene_positions/Homo_sapiens.GRCh38.104.gff3.gz -o data/intermediate/gene_positions/gene_positions_human.tsv -m data/raw/uniprot_swissprot/HUMAN_9606_idmapping.dat.gz -c Ensembl
	python3 src/preprocessing/gene_positions.py -i data/raw/gene_positions/Arabidopsis_thaliana.TAIR10.51.gff3.gz -o data/intermediate/gene_positions/gene_positions_athaliana.tsv -m data/raw/uniprot_swissprot/ARATH_3702_idmapping.dat.gz -c EnsemblGenome
	python3 src/preprocessing/gene_positions.py -i data/raw/gene_positions/Escherichia_coli_str_k_12_substr_mg1655_gca_000005845.ASM584v2.49.gff3.gz -o data/intermediate/gene_positions/gene_positions_ecoli.tsv -m data/raw/uniprot_swissprot/ECOLI_83333_idmapping.dat.gz -c EnsemblGenome
	python3 src/preprocessing/gene_positions.py -i data/raw/gene_positions/Saccharomyces_cerevisiae.R64-1-1.104.gff3.gz -o data/intermediate/gene_positions/gene_positions_yeast.tsv -m data/raw/uniprot_swissprot/YEAST_559292_idmapping.dat.gz -c EnsemblGenome


## Remove all dataset and feature data
clean:
	rm data/datasets/*.tsv
	rm data/datasets/*.fasta
	rm data/datasets/*.clstr
	find data/features/ -name *.tsv -delete

#################################################################################
# Self Documenting Commands                                                     #
#################################################################################

.DEFAULT_GOAL := help

# Inspired by <http://marmelab.com/blog/2016/02/29/auto-documented-makefile.html>
# sed script explained:
# /^##/:
# 	* save line in hold space
# 	* purge line
# 	* Loop:
# 		* append newline + line to hold space
# 		* go to next line
# 		* if line starts with doc comment, strip comment character off and loop
# 	* remove target prerequisites
# 	* append hold space (+ newline) to line
# 	* replace newline plus comments by `---`
# 	* print line
# Separate expressions are necessary because labels cannot be delimited by
# semicolon; see <http://stackoverflow.com/a/11799865/1968>
.PHONY: help
help:
	@echo "$$(tput bold)Available rules:$$(tput sgr0)"
	@echo
	@sed -n -e "/^## / { \
		h; \
		s/.*//; \
		:doc" \
		-e "H; \
		n; \
		s/^## //; \
		t doc" \
		-e "s/:.*//; \
		G; \
		s/\\n## /---/; \
		s/\\n/ /g; \
		p; \
	}" ${MAKEFILE_LIST} \
	| LC_ALL='C' sort --ignore-case \
	| awk -F '---' \
		-v ncol=$$(tput cols) \
		-v indent=19 \
		-v col_on="$$(tput setaf 6)" \
		-v col_off="$$(tput sgr0)" \
	'{ \
		printf "%s%*s%s ", col_on, -indent, $$1, col_off; \
		n = split($$2, words, " "); \
		line_length = ncol - indent; \
		for (i = 1; i <= n; i++) { \
			line_length -= length(words[i]) + 1; \
			if (line_length <= 0) { \
				line_length = ncol - indent - length(words[i]) - 1; \
				printf "\n%*s ", -indent, " "; \
			} \
			printf "%s ", words[i]; \
		} \
		printf "\n"; \
	}' \
	| more $(shell test $(shell uname) = Darwin && echo '--no-init --raw-control-chars')
