#### 2018-01-16
##### Summary: ASM absrtact submitted. Trying to reconcile genes between all the different genomes. Using CFT073 and MG1655 for reference

* Writing bash script to retrieve genomes from ncbi


#### 2018-01-17
##### Summary: Had an awful experience downloading genomes through ncbi ftp. Setting up get-homologs run

* See `getting_ncbi_genomes.sh`

* putting all genomes into one folder `genomes`, delete this when done with analysis

* Writing a script for get_homologues based on [this tutorial](http://digital.csic.es/bitstream/10261/146411/1/pangenome_workshop09032017.html)
* Softwares is [here](https://github.com/eead-csic-compbio/get_homologues)

* The pfam thing is not working for me. Resubmit without it. 
* Don't know how to stop this script from running. Nvm. The usual way. 
* It suggests you do the -D thing on the cluster, so might want to try again later. 
* Need to install get_homologues on cluster

#### 2018-01-18
##### Summary: get_homologues analysis cntd.

* The run took a long time - seems like it finished yesterday around 8. Started around 4.
* Uses protein_id instead of locus_tags for references. Try to parse genbanks with biopython

* [Parsing gbk](https://warwick.ac.uk/fac/sci/moac/people/students/peter_cock/python/genbank/)


#### 2018-01-23
##### Summary: Working on a script to get counts and orthologs for a specific clinical genome

#### 2018-01-24
##### Summary: Working on script to put all data together

#### 2018-01-25
##### Summary: Finished the script that would put together all the counts (or counts from whatever samples you want)

* Want you to clean the code on the weekend, think of best way to structure, add argparse, and flexibility with regards to the reference. 
* Need to extend up - include running get_homologues as an option

#### 2018-01-26
#### Summary: Have all the count data, putting it all together. As always, fingers crossed. 

* Downloading all the counts

```bash
scp annasint@flux-xfer.arc-ts.umich.edu:/scratch/hmobley_fluxod/annasint/HUTI-RNAseq/data/counts/*counts /Users/annasintsova/git_repos/HUTI-RNAseq/data/counts
```

* Running clean_counts.py
* Problem with HM43 samples, after some investigation figured out that the gff has FASTA genome attached at the end.. WHY???? Easy enough fix. 

```bash
sed -e '/##FASTA/,$d' data/annot/HM43_with_FASTA.gff > data/annot/HM43.gff
```
* now have to re-run htseq on HM43 samples :'(


```bash
bam=data/bt2_align/bam/sorted/HM43*bam
ref=data/annot/HM43.gff
python analysis/pipeline/lib/pipeline.py -a count -i $bam -o data/counts -ref $ref
```

#### 2018-01-29
##### Summary: Finally have all the counts, putting them together.

* Downloaded all the counts into `git_repos/HUTI-RNAseq/analysis/DE/data/`:

```bash

scp annasint@flux-xfer.arc-ts.umich.edu:/scratch/hmobley_fluxod/annasint/HUTI-RNAseq/data/counts/*sorted_counts .

```


* Seemingly successful

#### 2018-01-30
##### Summary: At this point I am satisfied with the OCGs I found, in the future want to come back and turn this into a small pipeline

* Given a set of genomes and their gene expression counts (raw RNAseq reads), return a table of all orthologs, with counts and RPKMs

* As is right now will need gbks for get_homologues, htseq-count output, gffs for RPKM calculation, flagstats for total read counts, alternatively just take in a file with the total read numbers

* This could be an interesting tool for gene expression meta analysis

* What are the caveats you have to think about? 
* Can in theory infer total # of reads from summing all the counts for all the features 


* To be of use needs to start from the beginning: take in short reads, gbks, gffs - produce counts and meta-analysis

* To be of further use can't rely on get_homologues alone, if want to do more genomes

* Also needs to perform downstream analysis

* What other visualizations besides PCA can you provide

* Need to develop downstream analysis

* Collaborate with Ali?

* Problem with HTSeq: based on python2, need to migrate everything to python3

#### 2018-01-31
##### Summary: Played around to get CEACAM3 RPKMs, so right now it's set up to use human_config and using combineCustomCounts instead of combineCounts

* Overall, I think combineCustomCounts is an improvement, incorporate that into main pipeline



#### 2018-02-08
##### Summary: Need to come up with metrics to evaluate get_homologues performance

* want to call get_homologues from python

```bash

pipe = subprocess.Popen(["perl", "uireplace.pl", var])#just an example, obviously adjust accordingly

```
* Looking back at the tutorial:


* There are [other pangenome](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4411478/table/t0005/) tools to consider

1.  Want to produce clusters of all sizes **(-t 0)**
2.  Started working on python script that would cycle thorough a bunch of conditions. Haven't figured out the naming conventions get_homologues using. Doing a test run, see what output directories will be called. 
3. Doesn't look like -I option worked


* After running a range of conditions compare distributions of cluster sizes. 


#### 2018-02-10
##### Summary: Code refactoring

* Main script should be responsible for:
    - Calling get_homologues with right parameters if needed
    - Creating a crossreference matrix
    - If counts are provided, putting counts and calculated RPKMs in as well
    
#### 2018-02-20
##### Summary: Re-worked the script, still can't figure out how to run these in parallel

* Re-worked the script so that there are csv files produced after each run, otherwise it's too frustrating
* Fixed some minor bugs (looking at `size` not `taxa` in the cluster.lists)
* Also now preserving some of the essential outputs
* Testing on 3 genomes if works will set up a bunch for overnight run
*  Getting an error in my clusterSizeDistribution, but can still do the runGH for conditions of interest

* glob.glob is not working. 


#### 2018-02-21
##### Summary: Fixed another bug, made so can submit multiple jobs at the same time to flux

* Re-worked version of this is now in `multistrain_rnaseq`, see that directory for more detailed documentation

#### 2018-02-23
##### Summary: Refactoring code that adds counts and RPKMS to OGC matrix

* Needs to be tested


#### 2018-02-26
##### Summary: Testing code from Feb 23

* Pulling data from flux into `git_repos/HUTI-RNAseq/data/get_homologs_output`
*  `scp` has an r argument that allows transferring whole directories

* Running crossRefOGC.py: success
* Running parsing_gbk.py for CFT073 and MG1655:success

* **NOTE**: all of this will need to be tested. Ex. of test to run: look up CFT073 and MG1655 in KEGG database and see if they are the same thing
* Running getGeneID to add counts and RPKM: success
* Exploratory data analysis






