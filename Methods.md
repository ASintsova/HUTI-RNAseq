## Genetically diverse uropathogenic Escherichia coli adopt a common transcriptional program in patients with urinary tract infections.



### Study Design


### RNA Sequencing


### RNA Sequencing Data Processing

**Trimmomatic v. 0.36**

`java jar <Trimmomatic path>/trimmomatic-0.36.jar SE <fastq_file> <trimmed_fastq_file> ILLUMINACLIP:<Trimmomatic-0.36 path>/adapters/TruSeq3-SE.fa:2:30:10:8:true SLIDINGWINDOW:4:15 MINLEN:40 HEADCROP:0`


**bowtie2 v. 2.3.4**

* Build index for each genome with `bowtie2-build <genome_file> <genome_index>`
* Align with `bowtie2 -x <genome_index> <fastq_file> -S <sam_file>`


**samtools v 1.5**

* convert to bam with `samtools view -b -o <bam_file> -@ 4 <sam_file>`
* sort with `samtools sort -o <sorted_bam> -@4 <bam_file>`
* index with `samtools index -@ 4 <sorted_bam>`


**HTseq v. 0.9.1**

`htseq-count -f bam -r pos -m union -s yes -t CDS -i ID <bam_file> <gff_file> > <count_file>`


### Quality Control

* Link to Jupyter notebook


### Identification of core genome

[find_core_genome.py](find_core_genome.py): Run get_homologues to identify core genome using directory of gbk files

**get_homologues v. 20170807**

<genome_folder> contains 14 clinical UPEC genomes isolated in this study, as well as *E. coli* K-12 MG1655 genome.

BDBHalg: `get_homologues.pl <genome_folder> -t 0 -C 50 -S 90 -c -e`
COGalg: `get_homologues.pl <genome_folder> -G -t 0 -C 50 -S 90 -A -c -e`
OMCLalg: `get_homologues.pl <genome_folder> -M -t 0 -C 50 -S 90 -A -c -e`


`compare_clusters.pl -o <core_genome_out_dir>  -d <BDBHalg_output>,<COGalg_output>,<OMCLalg_output>,  -t 15 -m`

[Ortholog table](link) was parsed from <core_genome_out_dir> using `run_gethomologues.py::cross_ref_from_core_genome`


### Identification of virulence factors

[find_virulence_factors.py](find_virulence_factors.py): Using multifasta with extracted virulence factor sequences from reference genomes find corresponding PROKKA IDs for each of the clinical genomes

**blast 2.7.1+**

Parameters: evalue < 0.000001 and identity > 80 and coverage > 90

    * To pull out fasta sequences for virulence factors:




