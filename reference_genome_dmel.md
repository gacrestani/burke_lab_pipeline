## Drosophila melanogaster genome
# Downloaded from https://flybase.org/genomes/Drosophila_melanogaster/dmel_r6.51_FB2023_02/fasta/index.html
# As of February 2026, the Burke Lab uses the release 6.51 (from 2023-04-19). This is not the most current one, but this is the one that our 'known sites' .vcf file used as a reference.

# How to get files to work with the pipeline
# FASTA
wget https://flybase.org/genomes/Drosophila_melanogaster/dmel_r6.51_FB2023_02/fasta/dmel-all-chromosome-r6.51.fasta.gz
gzip -d dmel-all-chromosome-r6.51.fasta.gz

# Create the .dict dictionary:
gatk CreateSequenceDictionary -R dmel-all-chromosome-r6.51.fasta

# Index the .fasta file
samtools samtools faidx dmel-all-chromosome-r6.51.fasta
bwa index dmel-all-chromosome-r6.51.fasta
# ----

# VCF ----
# Now download the known sites from DGRP2 (.vcf file)
wget https://resources.aertslab.org/DGRP2/NCSU/final/dm6/DGRP2.source_NCSU.dm6.final.vcf.gz

# Create index files
tabix -p vcf DGRP2.source_NCSU.dm6.final.vcf.gz

# Done!
