#!/bin/bash

mkdir data

chr=$(seq 22; echo MT X Y)

# Download and prepare genome GRCh38
mkdir data/GRCh38
for i in $chr; do
curl -o chr"$i".fa.gz ftp://ftp.ensembl.org/pub/release-99/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome."$i".fa.gz ;
done
gunzip *.fa.gz
gzip -9 *.fa
mv chrMT.fa.gz chrM.fa.gz
mv chr*.fa.gz data/GRCh38

# Download and prepare genome GRCh37
mkdir data/GRCh37
for i in $chr; do
curl -o chr"$i".fa.gz ftp://ftp.ensembl.org/pub/grch37/current/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.chromosome."$i".fa.gz ;
done
gunzip *.fa.gz
gzip -9 *.fa
mv chrMT.fa.gz chrM.fa.gz
mv chr*.fa.gz data/GRCh37

# Download bedfile UV damage for testing
curl -o table.csv.gz ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2585nnn/GSM2585711/suppl/GSM2585711%5FG%2D6%2D4A22%2E1%2Ecu%2Ebo%2Ehg19%2EcoToBa%2EcoToBe%2EunSo%2EcoBeToSiFr%2EslBeb6%2EcoToFiRa10%2EsoBe%2EcoBeToFa%2EgePyDi%2EsoBe%2Ebed%2Egz
gunzip table.csv
echo "chromosome,start,end,strand" > temp
awk '{print $1 "," $2+5 "," $3-4 "," $6}' table.csv >> temp
mv temp table.csv
gzip -9 table.csv
mv table.csv.gz data/
