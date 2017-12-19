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
