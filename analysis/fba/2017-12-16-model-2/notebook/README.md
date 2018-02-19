#### 2017-12-16
##### Summary: Setting up for Get_homologs run. Building infrastructure.

Downloaded gbf files from Ali's directory on flux. Going to create a directory for each strain and work on each model separately.

```
gbf_files=/Users/annasintsova/git_repos/HUTI-RNAseq/analysis/fba/2017-12-16-model-2/data
for f in $gbf_files/*; do  p_f=${f%%.*}; p_f=${p_f##*/};mkdir $p_f; mv $f $p_f; done

```

#### 2017-12-18
##### Summary: Actually running get_homologs. Default params as before. Modified mapOGC.py. Running on all genomes.

* Need to copy ref genome into each of the folders
```
for d in */; do cp ../../*making*/data/hm01/mg1655.gbf "$d";done
```

* Going to run get_homologues with defaults, except also set `-S`, which is minimal identity to 70% for now.

```
export GETHOMS=~/tools/get_homologues-macosx-20170807
for d in */; do gbk_dir=~/git_repos/HUTI-RNAseq/analysis/fba/2017-12-16-model-2/data/$d; perl $GETHOMS/get_homologues.pl -d $gbk_dir -S 70; done

```

* Running mapping script

```
cluster_list=mg1655_f0_alltaxa_algBDBH_e0_S70_.cluster_list
hom_directory=mg1655_f0_alltaxa_algBDBH_e0_S70_
b_num=/Users/annasintsova/git_repos/HUTI-RNAseq/analysis/fba/2017-12-16-model-2/data/ref/MG1655_asap_to_bnum.csv
map=/Users/annasintsova/git_repos/HUTI-RNAseq/analysis/fba/2017-12-16-model-2/lib/mapOGCs.py
for d in *homologues/
  do
  cL=$d$cluster_list
  hom_dir=$d$hom_directory
  name=${d%%_*}"_OGC.tab"
  python $map $cL $hom_dir $b_num $d$name
  done

```
* Success


#### 2017-12-19
##### Summary: Working on create_Strain_Model script, using chapter from Metabolic textbook as reference.

* Actually I'm an idiot. This is much easier if run get-homs on all strains at the same time and then work with the matrix. Going to try this now.
* Actually maybe not. Working with what I got. Get all the bnumbers from OGC files - these are all the genes that were found in my strain. Check logic on this one. Cross reference with all the genes in the model - delete whateveer is not found in the strain's genome with `cobra.manipulation.delete_model_genes`

#### 2018-01-03
##### Summary: When creating a model check if gene being deleted is 'essential'

* Added function to `create_Strain_Model.py` that check whether each gene to be deleted is 'essential', also can now see which 'essential' genes were not found in the clinical strain
* Function arbitrary decides that if deleting gene results in log reduction in growth it is essentials
* Played around more with gapfilling - removed all the exchange reactions from the model that is being gapfilled, and tried to gapfill it with all those exchange reactions with bounds set to -1000 and 0. Took a long long time, but still no go. Abandoning this line of investigation for now.
* Need to set up meeting with Evan for sometime next week to show him what I have double_gene_deletion

#### 2018-01-04
##### Summary: Testing model created yesterday. Created models for all 13 strains.

* Seemingly working: both base and HM01 model show reduction in growth rate under anaerobic conditions. But base model much more so. Overall though HM01 is somewhat crippled in it's growth rate compared to the base model. Might want to play with the ration (systematically) to see what makes most sense

* Few issues to address: change model name, do I need to trim orphan metabolites? Are there orphan metabolites? The way the model files are being read in is not great: each metabolite name has LPAREN and RPAREN in the name... Should remove parenthesis from the model? :S

* First things first: let's see why HM01 is missing 4 'essential' genes

* In the future need to come up with automated way to test growth on different conditions

* Modularized `create_Strain_Model.py` and created strain specific models for 13 clinical strains


* Working on script to automate lab notebook updates

#### 2018-01-05
##### Summary: Created a function that tests all 13 models for growth on 13 different carbon sources in aerobic and anaerobic conditions

* See `lib/Testing_strain_models.ipynb`, and  results are saved as csv files in the analysis folder

* Most strains/models look ok, something is off about HM7 and HM66. Need to go back to the list of deleted genes and see whats up.

* Other things to consider: develop scripts that quickly identify what genes were deleted from each model, and then possibly re-run blast to make sure those are absent. See if there are any patterns between different clinical strains.


#### 2018-01-11
##### Summary: Ran iron simulations. NO difference between models.

* Probably not surprising and good is that not differences between models.
* There are differences between aerobic and anaerobic conditions in terms both of how much iron is needed for at least minimal growth and need for ferrous vs ferric iron. Don't really know what is present during infection. 



#### 2018-01-15
##### Summary: Brainstorming paper outline

* Made a list of figures I want by meeting in Feb: 

    - Tree
    - PCA
    - heatmap of presence/absence of genes
    - heatmap of metabolic potential
    - heatmpat showing differences in expresion of metabolic pathway
    - zoom in on a couple of specific pathways
    
* Need to submit ASM abstract tomorrow.

* Want to take a closer at the genes being deleted, fill like something is not quite right there - same genes being deleted, but they are all related to LPS biosynthesis - maybe just too divergent?
* Easy script to go through all those files and find genes in common: count how many times each gene appears. Dictionary like.
* Interesting? all of them missing an lps biosynthetic cluster??
