uniref50.fasta.gz:
	@echo Downloading Uniref file...
	wget -c ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref50/uniref50.fasta.gz -O uniref50.fasta.gz.tmp
	mv uniref50.fasta.gz.tmp uniref50.fasta.gz

uniref50.fasta: uniref50.fasta.gz
	@echo Unzipping fasta file...
	gunzip -c uniref50.fasta.gz > uniref50.fasta.tmp
	mv uniref50.fasta.tmp uniref50.fasta
	rm uniref50.fasta.gz

# uniref50.fasta:
# 	@echo Unzipping fasta file...
# 	xz -d -c -T0 uniref50.fasta.xz > uniref50.fasta.tmp
# 	mv uniref50.fasta.tmp uniref50.fasta
# 	rm uniref50.fasta.xz

uniref50.fasta.pdb: uniref50.fasta
	@echo Creating BLASTDB...
	makeblastdb -in uniref50.fasta -parse_seqids -dbtype prot