import argparse
from collections import OrderedDict
import datetime as dt
import pandas as pd
import os
import re
import urllib.request
from termcolor import cprint

class Gene(object):

    def __init__(self, gene_id):
        self.gene_id = gene_id  # This could be either locus tag or gene name
        self.entry = ""
        self.name = ""
        self.definition = ""
        self.pathways = []
        self.aa_seq = ""
        self.nt_seq = ""
        self.pathway_string = ""

    def get_full_info(self):
        return "\t".join([self.entry, self.name, self.definition,
                          self.pathway_string]).rstrip('\t')

    def get_short_info(self):
        return "\t".join([self.entry, self.name, self.definition])

    def get_aa_seq(self):
        return ">{}|{}\n{}\n".format(self.entry, self.name, self.aa_seq)

    def get_nt_seq(self):
        return ">{}|{}\n{}\n".format(self.entry, self.name, self.nt_seq)

    def __str__(self):
        return "GENE:\n{}\t{}\t{}\nPATHWAYS:\n{}".format(self.entry, self.name,
                                                         self.definition,
                                                         "\n".join(self.pathways))

    def print_fancy(self):
        gene_string = "{}\t{}\t{}".format(self.entry, self.name, self.definition).expandtabs(tabsize=4)
        x = 'test'
        max_p = max([len(p) for p in self.pathways]) if self.pathways else 0
        max_len = max(max_p, len(gene_string))
        border = "_"*(max_len+1) if max_len else ""
        print(type(border))
        print("")
        cprint(border, "white")
        cprint("GENE:", "red", attrs=["bold"])
        print(gene_string)
        cprint(border, "white")
        cprint("PATHWAYS:", "blue", attrs=["bold"])
        print("\n".join(self.pathways))
        cprint(border, "white")
        print("")


    def kegg_get(self):
        try:
            url = "http://rest.kegg.jp/get/{}".format(self.gene_id)

            with urllib.request.urlopen(url) as f:
                lines = f.read().decode('utf-8').splitlines()

            for i in range(len(lines)):
                field = lines[i].split()[0]
                if field == 'ENTRY':
                    self.entry = lines[i].split()[1]

                elif field == "NAME":
                    self.name = lines[i].split("NAME        ")[1]

                elif field == "DEFINITION":
                    p = lines[i].split("DEFINITION  ")[1]
                    self.definition = re.sub(r"^\(.*?\)", "", p).strip()

                elif field == "PATHWAY":
                    p = lines[i].split("PATHWAY     ")[1]
                    count = 1
                    while lines[i + count].startswith("           ") and len(lines[i + count]) > 0:
                        p += lines[i + count]
                        count += 1
                    self.pathways = p.split("            ")
                    self.pathway_string = ", ".join(self.pathways).rstrip(", ")

                elif field == "AASEQ":
                    p = ''
                    count = 1
                    while lines[i + count].startswith("           ") and len(lines[i + count]) > 0:
                        p += lines[i + count]
                        count += 1
                    self.aa_seq = p.replace("            ", '')

                elif field == "NTSEQ":
                    p = ''
                    count = 1
                    while lines[i + count].startswith("           ") and len(lines[i + count]) > 0:
                        p += lines[i + count]
                        count += 1
                    self.nt_seq = p.replace("            ", '')

        except urllib.error.HTTPError:
            print("Bad gene identifier")


class GeneSet(object):

    def __init__(self, gene_id_list, genome="eco", out_dir="."):
        self.genome = genome
        self.kegg_ids = ["{}:{}".format(genome, gene_id) for gene_id in gene_id_list]
        self.name = ""
        self.gene_set = []
        self.out_dir = out_dir

    def get_info_for_each_gene(self):
        for gene_id in self.kegg_ids:
            fg = Gene(gene_id)
            fg.kegg_get()
            self.gene_set.append(fg)

    def write_nt_seq(self):
        today = dt.datetime.today().strftime("%Y_%m_%d")
        filename = os.path.join(self.out_dir,  "{}_{}_gene_set_nt_seq.fasta".format(today, self.genome))
        with open(filename, "w") as fo:
            for gene in self.gene_set:
                fo.write(gene.get_nt_seq())

    def write_aa_seq(self):
        today = dt.datetime.today().strftime("%Y_%m_%d")
        filename = os.path.join(self.out_dir,  "{}_{}_gene_set_nt_seq.fasta".format(today, self.genome))
        with open(filename, "w") as fo:
            for gene in self.gene_set:
                fo.write(gene.get_aa_seq())

    def get_info_df(self):
        genes = []
        labels = ["Entry", "Name", "Function", "Pathways"]
        for gene in self.gene_set:
            gene_info = (gene.entry, gene.name, gene.definition, gene.pathway_string)
            genes.append(gene_info)
        return pd.DataFrame.from_records(genes, index="Entry", columns=labels)

    def __str__(self):
        s = ""
        for g in self.gene_set:
            s += g.get_short_info() + "\n"
        return s.strip()


if __name__ == "__main__":

    gene_list = ["b0356", "b1001"]
    gs = GeneSet(gene_list)
    print(gs.kegg_ids)
    gs.get_info_for_each_gene()
    y = gs.gene_set
    for gene in y:
        print(gene.name)
    print(gs.get_info_df())
