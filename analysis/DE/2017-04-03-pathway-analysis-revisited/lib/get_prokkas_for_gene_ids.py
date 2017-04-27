# -*- coding: utf-8 -*-
"""
author: annasintsova
takes in gff file that has info for all the genomes and list of gene_ids that
you are interested in (i.e. all the core genes) and returns a file with
genome specific prokka number (so total of 13 numbers) for each gene_ids
columns will be sorted in alphabetical order

"""


gff_file = "/Users/annasintsova/git_repos/HUTI-RNAseq/data/annotations/kegg/new_upec_homologs.gff"

gene_ids_txt = "data/gene_ids.txt"

output_txt = "data/gene_ids_with_prokkas.txt"


def getProkkaForGeneID(gff_file, gene_ids_txt, output_txt):

    gi = open(gene_ids_txt)
    gene_ids = []
    for line in gi:
        gene_ids.append(line.rstrip())
    out = open (output_txt, "w+")

    title = 'gene_id\t'
    p_d = {'HM1':'a','HM3':'b',  'HM7':'d','HM14':'d', 'HM17':'e', 'HM43':'f',
                   'HM54':'g', 'HM56':'h', 'HM66':'i','HM57':'j', 'HM68':'k', 'HM86':'l','HM6':'c'}

    for key in sorted(p_d.keys()):
            title += str(key) + '\t'
    out.write(title.rstrip("\t")+ "\n")

    for gene in gene_ids:

        s = gene + '\t'
        fh = open(gff_file)
        prokka_dict = {'HM1':'a','HM3':'b',  'HM7':'d','HM14':'d', 'HM17':'e', 'HM43':'f',
                   'HM54':'g', 'HM56':'h', 'HM66':'i','HM57':'j', 'HM68':'k', 'HM86':'l','HM6':'c'}


        for line in fh:
            if gene in line:
                words = line.rstrip().split("\t")
                genome = words[0].split("_")[0]
                prk = words[8].split("PROKKA_ID=")[1]
                prokka_dict[genome] = prk

        for key in sorted(prokka_dict.keys()):

            s += prokka_dict[key] + '\t'
        out.write(s.rstrip("\t") + "\n")

        print (s.rstrip("\t"))


getProkkaForGeneID(gff_file, gene_ids_txt, output_txt)
