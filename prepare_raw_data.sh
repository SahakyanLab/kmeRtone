#!/bin/bash

mkdir data

# Download and prepare genome GRCh38

mkdir data/GRCh38

wget ftp://ftp.ensembl.org/pub/release-99/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.*.fa.gz

gunzip Homo_sapiens.GRCh38.dna.chromosome.*.fa.gz

gzip -9 Homo_sapiens.GRCh38.dna.chromosome.*.fa

for f in Homo_sapiens.GRCh38.dna.chromosome.*.fa.gz; do echo $f | mv $f `awk -F. '{print "chr" $5 "." $6 "." $7}'`; done

mv chrMT.fa.gz chrM.fa.gz

mv chr*.fa.gz data/GRCh38


# Download and prepare genome GRCh37

mkdir data/GRCh37

wget ftp://ftp.ensembl.org/pub/grch37/current/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.chromosome.*.fa.gz

gunzip Homo_sapiens.GRCh37.dna.chromosome.*.fa.gz

gzip -9 Homo_sapiens.GRCh37.dna.chromosome.*.fa

for f in Homo_sapiens.GRCh37.dna.chromosome.*.fa.gz; do echo $f | mv $f `awk -F. '{print "chr" $5 "." $6 "." $7}'`; done

mv chrMT.fa.gz chrM.fa.gz

mv chr*.fa.gz data/GRCh37


# Download bedfile UV damage for testing

wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2585nnn/GSM2585711/suppl/GSM2585711%5FG%2D6%2D4A22%2E1%2Ecu%2Ebo%2Ehg19%2EcoToBa%2EcoToBe%2EunSo%2EcoBeToSiFr%2EslBeb6%2EcoToFiRa10%2EsoBe%2EcoBeToFa%2EgePyDi%2EsoBe%2Ebed%2Egz

gunzip GSM2585711_G-6-4A22.1.cu.bo.hg19.coToBa.coToBe.unSo.coBeToSiFr.slBeb6.coToFiRa10.soBe.coBeToFa.gePyDi.soBe.bed.gz

echo "chromosome,start,end,strand" > table.csv
awk '{print $1 "," $2+5 "," $3-4 "," $6}' GSM2585711_G-6-4A22.1.cu.bo.hg19.coToBa.coToBe.unSo.coBeToSiFr.slBeb6.coToFiRa10.soBe.coBeToFa.gePyDi.soBe.bed >> table.csv

rm -f *.bed.gz
gzip -9 table.csv.gz
mv table.csv.gz data/