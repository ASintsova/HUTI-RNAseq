#!/usr/bin/env bash

# Set get_homologs path
export GETHOMS=~/tools/get_homologues-macosx-20170807

# Set data path

export top_dir=~/git_repos/HUTI-RNAseq/analysis/DE/2018-01-16-sync-genomes/data

export gbk_dir=$top_dir/genomes

# Put together all the gbk files I want to analyze. Moving to step 3 in the pipeline: generate robust core and pan-genomes

cd $top_dir
# D: require equal Pfam domain composition, might not want to use this. Might want to drop C to 70/
# removed it for now

$GETHOMS/get_homologues.pl -d $gbk_dir -n 6 &> log.get_homologues_16_BDBH_C75
$GETHOMS/get_homologues.pl -d $gbk_dir -G -t 0 -A -c &> log.get_homologues_16_Gt0c
$GETHOMS/get_homologues.pl -d $gbk_dir -M -t 0 -A -c &> log.get_homologues_16_Mt0c


# Generating single-copy orthologous gene clusters using the intersection of BDBH, COG, and OrthoMCL

export blast_dir=$gbk_dir"_homologues"
cd $blast_dir

$GETHOMS/compare_clusters.pl -d MG1655_f0_0taxa_algCOG_e0_,\
MG1655_f0_0taxa_algOMCL_e0_,\
MG1655_f0_alltaxa_algBDBH_e0_/\
 -o intersect_core_BCM_t16 -t 16

open intersect_core_BCM_t16/venn_t16.pdf

cd intersect_core_BCM_t16  && ls && ls *faa | wc

# Confirm that all clusters have only one sequence from each genome
grep '>' *faa | cut -d\| -f6 | sort | uniq -c

# Create robust pangenome
cd $blast_dir
$GETHOMS/compare_clusters.pl -d MG1655_f0_0taxa_algCOG_e0_,\
MG1655_f0_0taxa_algOMCL_e0_,\
 -o intersect_pan_CM_t0 -t 0 \
-m -T

open intersect_pan_CM_t0/venn_t0.pdf

cut -f 1-25 intersect_pan_CM_t0/pangenome_matrix_t0.tab

# There is a parsimony absence/presence tree too look at
