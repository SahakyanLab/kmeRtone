#!/bin/bash

mkdir data

# Download and prepare genome GRCh38

mkdir data/GRCh38

wget ftp://ftp.ensembl.org/pub/release-99/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.*.fa.gz

gunzip Homo_sapiens.GRCh38.dna.chromosome.*.fa.gz

gzip -9 Homo_sapiens.GRCh38.dna.chromosome.*.fa

for f in Homo_sapiens.GRCh38.dna.chromosome.*.fa.gz; do echo $f | awk -F. '{print "$5.$6.$7"}'; done

mv chr*.fa.gz data/GRCh38


# Download and prepare genome GRCh37

mkdir data/GRCh37

wget ftp://ftp.ensembl.org/pub/grch37/current/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.chromosome.*.fa.gz

gunzip Homo_sapiens.GRCh37.dna.chromosome.*.fa.gz

gzip -9 Homo_sapiens.GRCh37.dna.chromosome.*.fa

mv chr*.fa.gz data/GRCh37


# Download test bedfile UV damage

wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2585nnn/GSM2585711/suppl/GSM2585711%5FG%2D6%2D4A22%2E1%2Ecu%2Ebo%2Ehg19%2EcoToBa%2EcoToBe%2EunSo%2EcoBeToSiFr%2EslBeb6%2EcoToFiRa10%2EsoBe%2EcoBeToFa%2EgePyDi%2EsoBe%2Ebed%2Egz

mv *.bed.gz data/