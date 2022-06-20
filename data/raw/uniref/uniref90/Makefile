uniref90.fasta.gz:
	@echo Downloading Uniref file...
	wget -c ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/uniref90.fasta.gz -O uniref90.fasta.gz.tmp
	mv uniref90.fasta.gz.tmp uniref90.fasta.gz

uniref90.fasta: uniref90.fasta.gz
	@echo Unzipping fasta file...
	gunzip -c uniref90.fasta.gz > uniref90.fasta.tmp
	mv uniref90.fasta.tmp uniref90.fasta
	rm uniref90.fasta.gz

# uniref90.fasta:
# 	@echo Unzipping fasta file...
# 	xz -d -c -T0 uniref90.fasta.xz > uniref90.fasta.tmp
# 	mv uniref90.fasta.tmp uniref90.fasta
# 	rm uniref90.fasta.xz

uniref90.fasta.pdb: uniref90.fasta
	@echo Creating BLASTDB...
	makeblastdb -in uniref90.fasta -parse_seqids -dbtype prot