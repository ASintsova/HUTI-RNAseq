#### 2017-08-11
##### Summary: Learning get_homologues

### First attempt at making strain specific model

1. **get_homologues**
- Following [this tutorial](http://digital.csic.es/bitstream/10261/146411/1/pangenome_workshop09032017.html), and [this guide](http://161.111.227.80/compbio/soft/manual_get_homologues.pdf)

- Installing get_homologues

```
tar xfz get_homologues-macosx-20170807.tar
cd get_homologues-macosx-20170807
perl install.pl

```
- Seemed to have worked ok. Also installed optional databases
- Can use 3 different clustering algs: BDBH, OMCL, and COGtriangles

- prefered input gbk files, but can also use fasta

- Following one of the example protocols from the tutorial:

```

export GETHOMS=~/tools/get_homologues-macosx-20170807
# set data paths
export gbk_dir=~/git_repos/HUTI-RNAseq/analysis/fba/2017-08-11-making-of-the-model

# algorithm: BDBH
# alignment coverage: 75% [-C 75]
# E-value: 1e-05 [-E 1e-05]
# cluster size: clusters with 1+ representative protein from each proteome analyzed [-t number_of_proteomes]
# run mode: local [-m local])


cd $gbk_dir
perl $GETHOMS/get_homologues.pl -d data -n 2

export blast_dir=$gbk_dir/data_homologues

# 1.3 explore contents by file extension names
ls | cut -d\. -f2 | sort | uniq -c

#creates orthologoues clusters, other options available, i.e. do not allow paralogs, Pfam annotation, etc.

## 3 compute robust core and pan-genomes
# 3.1 produce pan-genome clusters of all sizes (-t 0) with both the COG and OrthoMCL algorithms
# [Note that Average Peptide Identities (-A) pan-genome growth simulation are also produced (-c),
#  please check the produced logfiles]
# If Average Nucleotide Identities are desired then option -a should be used, enforcing the
# alignment if nucleotide sequences
# If option -X is set DIAMOND would be called instead of BLASTP, which reduces computing time
# when many input genomes are to be analyzed


cd $gbk_dir
$GETHOMS/get_homologues.pl -d data -M -D -t 0 -A -c &> log.get_homologues_pIncAC_MDt0c

$GETHOMS/get_homologues.pl -d data -G -D -t 0 -A -c &> log.get_homologues_pIncAC_GDt0c

$GETHOMS/compare_clusters.pl -d EscherichiacolistrainSCEC2plasmidpSCEC2NC022377_f0_0taxa_algCOG_Pfam_e0_,\
EscherichiacolistrainSCEC2plasmidpSCEC2NC022377_f0_0taxa_algOMCL_Pfam_e0_,\
EscherichiacolistrainSCEC2plasmidpSCEC2NC022377_f0_alltaxa_algBDBH_Pfam_e0_\
 -o intersect_core_BCM_Dt12 -t 12 -m

```



#### 2017-08-14
##### Summary: Playing with my first model, running get_homologues on HM01

### First attempt at making strain specific model cntd


- Going to start by running get_homologues BDBH alg on MG1655 and HM01, then use the cluster groups and map it on existing MG1655 model

- Download MG 1655 gbk file. Downloaded from Patric - very easy to use!

- Using Ali-generated gbk file for HM01, with PROKKA gene ids

- Running with all default params for now

```
export GETHOMS=~/tools/get_homologues-macosx-20170807
export gbk_dir=~/git_repos/HUTI-RNAseq/analysis/fba/2017-08-11-making-of-the-model/data
cd $gbk_dir
perl $GETHOMS/get_homologues.pl -d hm01 -n 2

```
- start: 12:23, end 12:26

- looking at the output:


```

>ID:NP_414543.1 |[Escherichia coli str. K-12 substr. MG1655]|K-12|thrA|2463|NC_000913(4639675):337-2799:1 ^ASAP:ABE-0000008,ECOCYC:EG10998,EcoGene:EG10998,GI:16127996,GeneID:945803,UniProtKB/Swiss-Prot:P00561^ Escherichia coli str. K-12 substr. MG1655, complete genome.|neighbours:ID:NP_414542.1(1),ID:NP_414544.1(1)|neighbour_genes:thrL,thrB| | aligned:1-820 (820)
```

- doesn't parse out gene_locus, but can use ASAP ID for now, might need to parse gene locus after from gbk file

- wrote python script mapping prokkas to ASAP ids, for reaction mapping need to covert to locus tags, a.k.a. b numbers


- plan of action would be: read in E.coli pangenomic model, for each reaction update gene reaction rule -> this should incorporate deleting reactions not found in HM01 and keeping track of the deleted ones.

- if its blank to begin with, leave it blank, however if it isn't blank, but not homologs found delete reaction


#### 2017-08-15
##### Summary: Building first model with CobraPy

- MG1655 model was downloaded from [this paper](https://bmcsystbiol.biomedcentral.com/articles/10.1186/1752-0509-5-182)

- This is going to be base, but first hoping to work out a script with sample model that comes with cobrapy package

- have a semi working script, want to go back, make sure I'm dealing with paralogs ok

```
grep -c '>' *faa|grep -v ':2'|wc -l
# overall had 78 clusters genes that had paralogs
```

- These can contain both multiple MG1655 genes, and multiple hm01 genes,
thus for each prokka number in the file want to have all possible ASAP number_of_proteomes

- modify script accordingly, center everything around prokka numbers
- mapOCG.py seems to work well, including paralogs, side benfit - finally can reliably map bnums onto HM genomes! Check all my previous work. This wasn't so hard after all


- worked out a script to update gene reaction rule with new gene_ids, might have to work around some bracket issues, but in general this works!

- next step: read in mg1655 model, and see if can update with HM01 info. See how many reactions are removed. Not worrying about adding extra at this time.


#### 2017-08-16
##### Summary: Contd working on Cobrapy script for strain specific model

### Building a model contd

- mg1655 model in `data` folder
- had to install another package `libsbml`, because the sbml file seems to be in outdated format
- successfully loaded the model ;D
- had to change some syntax, as this models gene reaction rules are formated slightly differently
- script not working
- nevermind, something about not saving before running
- everything is poised for my first model creation
- fingers crossed

- ran

```
python lib/mapOGCs.py data/hm01_homologues/mg1655_f0_alltaxa_algBDBH_e0_.cluster_list data/hm01_homologues/mg1655_f0_alltaxa_algBDBH_e0_ data/MG1655_asap_to_bnum.csv data/hm01/hm01_ids_to_mg1655_bnums.txt

```

- then

```
python lib/buildingModel.py data/hm01/hm01_ids_to_mg1655_bnums.txt data/mg1655_model.XML data/hm01/hm01_model1.xml

```


- analyze created model after lunch
- this updated all the reactions and removed quite a few (118) - might have to put some back if are essential. Did not do anything to metabolites. Also did not remove the bnum genes - need to do that.

- how do I know if the model is 'gapless'

#### 2017-08-30
##### Summary: Contd working on Cobrapy script for strain specific model. Trying to avoid gaps. 

### Building a model contd

- adding removal of all bnum genes from the model, should not make a big difference in long term, but cleaning up.

- Cannot find a biomass reaction in the base model I'm using. Add one if I can find it? Use a different base that already has this incorporated? Why is it missing?

- I think idea behind creating a gapless model is trying to run your biomass optimization after deletion of each reaction, and if flux is inhibited put it back in.


#### 2017-09-01
##### Summary: Generated model not workable due to too many reactions being removed. Figure out a way around this.

### Building a model contd

- Looking for biomass reaction function. Found one in a different model of MG1655. See what the differences are, probably just use the other model as a base? Will this be a problem when incorporating multiple models???


- Looking for the latest/most updated MG1655 model. Going to stick with just incorporating one model for now. This latest (is it latest?) seems to have two different bioMass reactions. What is the difference?

- trying to run old script on new base model

```
python lib/buildingModel.py data/hm01/hm01_ids_to_mg1655_bnums.txt data/mg1655_2011_model.xml data/hm01/2017-09-01-hm01-model-1.xml

```

- as expected, deleted way too many reactions, and consequently flux via model is infeasible

- see if deleting disables flux

- success, but now significantly diminished growth rate, possible make a different cut off, not just > 0 flux through biomass reaction.


- all right onto mapping confidence scores. First get those pesky gh ids onto prokkas ... fml...



#### 2017-09-13
##### Summary: Evan has no confidence in what I am doing, wants me to use Model Seed

### Building a model cntd

- Last week met with Evan. Not surprisingly he was not very confident in my ability to create a decent model, so instead of what I am doing he wants me to use modelSeed.
- To do so need to re-annotate genomes with RAST then build a model. Going to do so via PATRIC. How am I going to integrate gene expression data into this??
- Keeps failing I'm not sure what this is about, too dumb right now to figure it out



#### 2017-12-05
##### Summary: Trying to revisit this model after a long break. Switched to playing with PATRIC models as a warm up. 

### Building model: starting again.

- Although I'm pretty confident I will have to come up with my own model at some point - for fun if for nothing else, right now playing around with PATRIC and modelSEED.
- [modelSEED original paper](https://www.nature.com/articles/nbt.1672) 

- uploaded all fasta files into PATRIC
- ran the annotation tool, all queued right now, probably will have to build the models tomorrow. 
- going to go read that paper now


#### 2017-12-13
##### Summary: Found an awesome resource with detailed workflow that will allow me to build strain specific models. 

### Building Strains Specific Model 2.0

- Overall I was on the right track. A few python shortcuts to make life easier that I did not know about. 
- As I suspected the hard but important part will be to incorporate multiple models into my reference instead of using E.coli alone. 
- Really excite to start playing with code etc. 