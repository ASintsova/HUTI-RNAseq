#LS-BSR

## Running LS-BSR  

- Requires quite a few modules/dependencies:
	- in flux: 
	
	```
	module load usearch 
	module unload python
	module load python/2.7.3 
	module load biopython 
	module load ls-bsr 
	module load prodigal  
	``` 

	- Click [here](https://github.com/jasonsahl/LS-BSR/blob/master/LS_BSR_manual.pdf) for more details
- can either supply a non-redundant gene list or use the de novo gene prediction method
- Inputs for LS-BSR: directory of your genomes, a fasta file of query genes  

	- 	something like: `python /../../ls_bsr.py -d genomes -c usearch `  

	
	
## Analyzing bsr_matrix in R

- Limit analysis to CFT073 and 13 (not HM60) strains that we did most recent RNAseq on
- Load in the matrix:
`bsr_mat_PG = read.table('~/Desktop/bsr_matrix_values_annot.txt', sep = "\t", row.names = 1, header = TRUE, quote = "")`
- Core genome:  
`sum(rowSums(bsr_mat_PG > 0.4) == 14)` -> 4089 
`sum(rowSums(bsr_mat_PG > 0.9) == 14)` -> 2407
`sum(rowSums(bsr_mat_PG > 0.95) == 14)` -> 2167 
- Set arbitrary cut off at 0.9
`core <- bsr_mat[rowSums(bsr_mat >0.9)==14,]`
- Ok, this is not pretty but I'm too tired to figure this out:
`write.table(rownames(core), "core_ids2.txt", quote = FALSE, row.names = FALSE)`
For whaterver reason first line is 'x' -> delete this, put file on flux
- using **select_seqs_by_IDs.py**, provided by LS-BSR to extract all sequences for all core genes from consensus.fasta
- using blast as before to define genomic regions in each of the genomes as described in **Comparative genomics and transcriptomics of eut operon.md** 
- This time using .pbs script

```

####  PBS preamble

#PBS -N blast

#PBS -M annasint@umich.edu
#PBS -m bea
#PBS -j oe

#PBS -l nodes=1:ppn=4,pmem=4gb,walltime=12:00:00
#PBS -V

#PBS -A hmobley_fluxod
#PBS -l qos=flux
#PBS -q fluxod

####  End PBS preamble
####  Commands follow this line

# This script requires ncbi-blast modules loaded

if [ -s "$PBS_NODEFILE" ] ; then
    echo "Running on"
    cat $PBS_NODEFILE
fi

if [ -d "$PBS_O_WORKDIR" ] ; then
    cd $PBS_O_WORKDIR
    echo "Running from $PBS_O_WORKDIR"
fi


blastn -db databases/HUTI_ref -query core_genes.fasta  -evalue 1e-3 -word_size 11 -outfmt "6 sseqid sstart send evalue qse$


```
- separate by genomes
``` 
while read p;do grep "HM"$p"_" core.genes.blast.result > core_genes_by_genome/"HM"$p"_genes"; done < strain.numbers
```

- for i in HM*; do  python ../blast_to_dataframe.py $i; done
- move file to desktop
- run ReadCount.R  
- 