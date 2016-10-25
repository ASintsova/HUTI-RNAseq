# Comparative genomics and transcriptomics of *eut* operon

##*eutR* as an example

1. _Extracting eutR genomic location from HUTI strains and CFT073_
	* Download CFT073 sequence from NCBI: c2971(yfeG).  Location:2832730 to 2831678 (Rev Comp)
> 
ATGAAAAAGACCCGTACAGCCAATTTGCACCATCTTTATCATGAACCCTTACCCGAAAACCTGAAGCTCA
CGCCGAAGGTCGAAGTGGATAATGTTCATCAACGACAGACAACGGACGTCTATGAACATGCTTTGACAAT
TACCGCCTGGCAGCAGATTTACGATCAGCTGCATCCGGGCAAGTTTCATGGTGAATTTACGGAAATTCTA
CTCGATGATATTCAGGTTTTTCGTGAATACACTGGTCTGGCGCTGCGTCAGTCGTGCCTGGTCTGGCCGA
ACTCGTTCTGGTTTGGCATTCCGGCGACGCGCGGTGAGCAGGGATTTATCGGTTCGCAATGTCTGGGAAG
TGCAGAAATTGCGACGCGCCCTGGTGGAACCGAGTTTGAATTAAGTACGCCGGACGATTACACGATCCTT
GGCGTAGTGCTTTCTGAAGATGTCATCACTCGTCAGGCTAACTTTTTGCATAATCCGGATCGGGTATTAC
ATATGCTGCGTAGCCAGTCGGCGCTGGAAGTGAAAGAGCAGCATAAAGCCGCGCTGTGGGGCTTTGTCCA
ACAGGCGCTGGCGACGTTTTGCGAGAACCCGGAAAATCTCCATCAGCCTGCAGTGCGAAAAGTGCTGGGG
GATAATTTGCTAATGGCGATGGGGGCGATGCTGGAAGAAGCGCAGCCGATGATGACGGCGGAAAGCATCA
GTCATCAGAGTTACCGTCGATTGCTTTCCCGCGCCCGTGAATATGTGCTGGAAAACATGTCCGAACCGGT
GACGGTGCTGGACTTGTGTAATCAACTGCATGTCAGCCGCCGCACGCTACAAAACGCGTTTCACGCCATT
TTGGGCATTGGTCCGAACGCGTGGCTGAAACGTATTCGCCTGAACGCCGTGCGCCGCGAACTCATAAGTC
CGTGGTCACAAAGCACAACGGTAAAAGACGCCGCCATGCAGTGGGGATTCTGGCATCTTGGGCAATTTGC
CACGGATTATCAGCAACTGTTTGCCGAGAAGCCGTCACTGACGCTGCATCAGCGGATGCGGGAGTGGGGG
TGA
	
	* Create balstdb from 20 *E.coli* genomes:
	
	```
	cat CFT073.fasta *HM* > ../EUT/HUTI_20.fasta  
	
   makeblastdb -in HUTI_20.fasta -title HUTI_ref -dbtype nucl -out databases/HUTI_ref
   ```  
   * Blast EutR against the databse (need the following fields: sseqid: what genome this is coming from, s. start, s. end, evalue, qseqid: query seqid)

   ```
   blastn -db databases/HUTI_ref -query eutR.fasta -evalue 1e-3 -word_size 11 -outfmt "6 sseqid sstart send evalue qseqid" > eutR2.test 
   ```
   * This will give us genomic location for eutR in each of the sequences. Now will need to create a separate file for each genomic sequence. 
   ```
   grep "HM1_" eutR2.test > HM1_eutR
   ```
   * First created txt file that contains #s corresponding to strains
   * Then:
   ```
   while read p;do grep "HM"$p"_" eutR2.test > "HM"$p"_eutR"; done < strain.numbers
   ```
   * move these files into separate directory, for each run blast_to_gtf.py
```
   for i in HM*; do  python blast_to_gtf.py $i; done
```

	* Move these annotation files onto local system to run R


```
csvfile <- file.path("sample_info.csv")
(sample_info <- read.csv(csvfile,row.names=1))

strain_nums <- c(1, 3, 6,7,14, 17, 43, 54, 56, 57,66, 68, 86)
#strain_nums <- c(1,3)
res <- c()
for (i in strain_nums){
    filenames <- paste0(levels(sample_info$Condition), i, "_SE_aln_sort.bam")
    file.exists(filenames)
    
    library("Rsamtools")
    bamfiles <- BamFileList(filenames, yieldSize=2000000)
    library("GenomicFeatures")
    gtffile <- paste0("HM", i, "_eutR.gtf")
    file.exists(gtffile)
    
    (txdb <- makeTxDbFromGFF(gtffile, format="gtf", circ_seqs=character()))
    (ebg <- exonsBy(txdb, by="tx"))
    library("GenomicAlignments")
    library("BiocParallel")
    se <- summarizeOverlaps(features=ebg, reads=bamfiles,
                            mode="Union",
                            singleEnd=TRUE,
                            ignore.strand=FALSE )
    res <- append(res, assay(se))
    }
final_eutR <- matrix(res, ncol=2, byrow= TRUE) 
```
* These will be raw counts: since only one gene silly to calculate RPKM? 


_Expression of Aerobactin operon_

* Now want to expand to the whole operon: if can work for 4 genes, should work for 2000, right?

* python to get genomic location for said genes out of gbk file: crude but all I got right now

```
import re
def getAllGenesLocs (filename):
	
	x = []
	z = ''
	y = ''
	with open(filename) as fh:
		while True:
			z = fh.readline()
			if len (z) == 0:#if reached the end of the file want to break
				break
			if ' gene' in z:
				y = fh.readline().strip()
				x.append((z.strip(), y))
			else:
				pass
	return x


def getYourGenesLocs(geneList, allgenelocs, fn):
	gen_locs = open(fn, 'w+')
	for gene in geneList:
		for ag in allgenelocs:
			if gene in ag[1]:
				if 'complement' in ag[0]:
					loc = ag[0].split('(')					
					print "gi|26111730|gb|AE014075.1|" +  "\t" + loc[1].split('..')[0] + '\t' + loc[1].split('..')[1][:-1]+ '\t' + str(gene)
					gen_locs.write("gi|26111730|gb|AE014075.1|"+ "\t" + loc[1].split('..')[0] + '\t' + loc[1].split('..')[1][:-1]+ '\t' + str(gene)+"\n")
				else:	
					loc = ag[0].split('            ')
					print "gi|26111730|gb|AE014075.1|" + "\t" + loc[1].split('..')[0] + '\t' + loc[1].split('..')[1]+ '\t' + str(gene)
					gen_locs.write("gi|26111730|gb|AE014075.1|" + "\t" + loc[1].split('..')[0] + '\t' + loc[1].split('..')[1]+ '\t' + str(gene)+"\n")
				
				
	
	
		
#print (getAllGenesLocs("test.gbk"))

print (getYourGenesLocs(['iucA', 'iutC', 'iucC', 'iucD', 'iutA', 'fitA'],getAllGenesLocs("CFT073.gbk"), "test2.bed"))
```

* use bedtools getfasta to get the sequences out of genome:

```
bedtools getfasta -name  -fi CFT073.fasta -bed test2.bed -fo iut.test  
```

* Now back to the beginning:

	* blast:
	* delete all empty files
	* ```find . -size 0 -delete```
	* 




 

