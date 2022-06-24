.PHONY: setup_ubuntu requirements package raw_data extract_pssms blast_databases blastdb_uniref50 blastdb_uniref90

#################################################################################
# Conventions                                                                   #
#################################################################################

# "human": 		9606
# "athaliana":		3702
# "ecoli": 		83333
# "yeast": 		559292

#################################################################################
# Setup                                                                         #
#################################################################################

## Install packages required on Ubuntu 20.04 LTS WSL
setup_ubuntu:
	sudo apt update && sudo apt upgrade -y
	sudo apt install build-essential
	wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh
	bash ~/miniconda.sh -p ~/miniconda3
	rm ~/miniconda.sh
	@echo Reload shell to put conda in path: source ~/.bashrc

## Install packages.
requirements:
	conda update -n base -c defaults conda
	conda env create --file environment.yml

## Install code as python package to use it in notebooks
package:
	pip install -e .

#################################################################################
# Raw data                                                                      #
#################################################################################

## Download raw data
raw_data:
	curl "https://www.uniprot.org/uniprot/?query=reviewed:yes&format=tab&columns=id,genes,protein%20names,organism,organism-id,keyword-id,keywords,go-id,go,database(TCDB),existence,sequence,fragment&sort=score" > data/raw/swissprot/sp_data.tsv

#################################################################################
# Raw data: BLAST                           		                            #
#################################################################################

## Extract pssms from archive
extract_pssms: 
	mkdir -p data/intermediate/blast
	tar xf data/intermediate/blast.tar.xz --directory=data/intermediate/blast
	rm data/intermediate/blast.tar.xz

## Init blast dbs for creating additional PSSMs. >100GB needed
blast_databases: blastdb_uniref50 blastdb_uniref90

blastdb_uniref50: 
	@echo Creating BLASTDB...
	$(MAKE) -C data/raw/uniref/uniref50 uniref50.fasta.pdb

blastdb_uniref90: 
	@echo Creating BLASTDB...
	$(MAKE) -C data/raw/uniref/uniref90 uniref90.fasta.pdb


#################################################################################
# Self Documenting Commands                                                     #
#################################################################################

.DEFAULT_GOAL := help

# Copied from https://github.com/drivendata/cookiecutter-data-science
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
