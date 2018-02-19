#!/usr/bin/env bash

# ref genomes: CFT073: GCA_000007445.1, MG1655: GCA_000005845.2
# curl -s ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/assembly_summary.txt| cat|grep "MG1655"|cut -f1,4,11
genbank=ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/
targetTaxa=Escherichia_coli
outDir=coliGenomes
mkdir -p $outDir

#genomes=`curl -ls  ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/$targetTaxa/latest_assembly_ versions/ | head

#for g in $genomes
 #do mkdir $outDir/$g
 #wget ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/$targetTaxa/latest_assembly_ve rsions/$g/* -P $outDir/$g
 #done


mkdir data/CFT073
cd data/CFT073

files=$(curl -ls ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/Escherichia_coli/latest_assembly_versions/GCA_000007445.1_ASM744v1/)

for f in $files:
do
curl ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/Escherichia_coli/latest_assembly_versions/GCA_000007445.1_ASM744v1/$f -O
done

cd ..
mkdir MG1655
cd MG1655
curl -ls ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/Escherichia_coli/latest_assembly_versions/|grep "GCA_000005845.2"
files=$(curl -ls ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/Escherichia_coli/latest_assembly_versions/GCA_000005845.2_ASM584v2/)
for f in $files:
do
curl ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/Escherichia_coli/latest_assembly_versions/GCA_000005845.2_ASM584v2/$f -O
done
