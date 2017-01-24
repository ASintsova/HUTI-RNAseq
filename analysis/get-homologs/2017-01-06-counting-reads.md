# 2017-01-06-RNAseq-DE-analysis

## INPUT: ANNOTATIONS
### Working with get-homologs output
* get-homologs run by Ali
	* Questions: all the genes found are always on HM#_1, what are the other sequences?, why we do not find anything on them?
	
* folder of fasta files, name of the file: name of the homolog, contains fasta for each actual protein from all the genomes where it was present. The info line (">") contains name of genome, location, strand, name. 
	* confirm that this is correct?

	```
	
	```
	
* For each file name insert name of the file (i.e. unique homolog identifier) as first line after ">>"

* In shell:

```
for i in *fna
do FILE=$(ls $i)
echo ">>"${FILE%.*} > ${FILE%.*}"_name.fna"
cat $FILE >> ${FILE%.*}"_name.fna"
rm $FILE
done

```
* Concatenate all the files: `cat *fna`
* run **get_gffs_from_homolgs.py**:

```
# -*- coding: utf-8 -*-
"""
Created on Fri Jan  6 10:50:55 2017

@author: annasintsova
"""
import sys
filename_in = str(sys.argv[1])
filename_out = str(sys.argv[2])

def getGFFs (filename_in, filename_out):
 
    gff = open(filename_out, "w+")
    with open (filename_in, 'r') as fh:
        for line in fh:
            if line.startswith(">>"):
                gene_id = line.rstrip()[2:]
                
                
            elif line.startswith(">"):
                next_line = line
                description = next_line.split("|")
                genome = description[2]+ "_1"
                gene_name = description[3]
                start = description[5].split(":")[1].split("-")[0]
                end = description[5].split(":")[1].split("-")[1]
                strand = description[5].split(":")[2].split("^^")[0].strip()
                prokka = description[0].split(":")[1]
                if strand == '1':
                    strand = '+'
                elif strand == '-1':
                    strand = '-'
                else:
                    strand = '.'

                gff.write('{}\tcustom\tCDS\t{}\t{}\t.\t{}\t.\tgene_id "{}"; gene_name "{}"; PROKKA_ID "{}"\n'.format
                      (genome, start, end, strand, gene_id, gene_name, prokka))     
                
getGFFs(filename_in, filename_out)


```

* Sort these by first column (i.e. genomes) and separate into separate gffs for each genome
	- `while read p; do grep "^HM"$p"_1" test.gff > $p"_out"; done < strain.numbers`

* Final Pbs script: get_gffs.pbs

```



```


## INPUT: ALIGNMENTS:
* Bam files sorted by positions

## COUNTING READS
* Use htseq-count with Ali's parameters to count reads

* Trial PBS Script (- preamble): **MODIFY**

```
DATA=/scratch/hmobley_fluxod/annasint/HUTI_RNAseq/Data
ALIGN=/scratch/hmobley_fluxod/annasint/HUTI_RNAseq/Data/Alignments/BAM_by_Position_New_Samples
ANNOT=//scratch/hmobley_fluxod/annasint/HUTI_RNAseq/Data/Annotations/get-homologs-gff

htseq-count -t CDS $ALIGN/UR17.sam $ANNOT/HM17_gh.gff > HM17.out

```
* Creating a loop to go through all alignment files:

```
for i in $ALIGN/*.bam; do NUM=$(ls $i|sed -e s/[^0-9]//g); echo "HM"$NUM"_gh.gff"; done

```
* Right now have issues with htseq-count: not reads are being counted
	* troubleshooting: is bam file not sorted properly?:
		* resort by name see if that helps
	* is gff file not being read in correctly?
	* problem identified: genome names not identical in bam and gff file
		* solution: on further examination learned that all the genes found by get homologues are always on HM*NUM*_1 sequence (Question then: the rest are things not asembled well?? )
		* have to re-do gff files
* Merge count data with pandas

## CORE GENOME ANALYSIS 
* Import final counts into R
* Look for complete cases
* See 2017-01-20-counts-analysis.R

	
