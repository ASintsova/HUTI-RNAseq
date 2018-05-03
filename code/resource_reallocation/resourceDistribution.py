

import urllib.request
from collections import OrderedDict
import datetime as dt
import os


# 1. Find crp regulated genes
def processGeneregulationFile(gene_regulation_file):

    """

    Building up data structure from the RegulonDB file

    :param gene_regulation_file:
    :return: list of tfs, list of genes regulated by tfs, tf:{actv:[gene list], repr:[gene list]}, gene:{activ:[tf list]}

    """
    regulators = []
    all_regulated = []
    regulons = {}
    genes_to_tf = {}
    with open(gene_regulation_file, "r") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            else:
                regulator = line.split()[1]
                regulated = line.split()[7]
                modulation = line.split()[8]
                if regulator not in regulators:
                    regulators.append(regulator)
                if regulated not in all_regulated:
                    all_regulated.append(regulated)

                if not regulator in regulons.keys():
                    regulons[regulator] = {modulation: [regulated]}
                elif modulation not in regulons[regulator].keys():
                    regulons[regulator][modulation] = [regulated]
                elif regulated not in regulons[regulator][modulation]:
                    regulons[regulator][modulation] = regulons[regulator][modulation] + [regulated]

                if not regulated in genes_to_tf.keys():
                    genes_to_tf[regulated] = {modulation: [regulator]}
                elif modulation not in genes_to_tf[regulated].keys():
                    genes_to_tf[regulated][modulation] = [regulator]
                elif regulator not in genes_to_tf[regulated][modulation]:
                    genes_to_tf[regulated][modulation] = genes_to_tf[regulated][modulation] + [regulator]

    return regulators, all_regulated, regulons, genes_to_tf


def findRegulon(gene_regulation_file, tf, out_dir, modulation = ''):

    """

    :param gene_regulation_file:
    :param tf:
    :param out_dir:
    :param modulation: options activator repressor dual unknown
    :return:
    """
    regulon = []
    regulators, all_regulated, regulons, genes_to_tf = processGeneregulationFile(gene_regulation_file)
    if modulation:
        regulated = regulons[tf][modulation]
        regulon = [(r, modulation) for r in regulated]
    else:
        mods = regulons[tf].keys()
        for mod in mods:
            regulon += [(r, mod) for r in regulons[tf][mod]]
    out_file = "{}_{}_{}.csv".format(dt.datetime.today().strftime("%Y-%m-%d"),tf,
                                 modulation if modulation else 'all')
    with open(os.path.join(out_dir, out_file), "w") as fo:
        for tu in regulon:
            fo.write("{},{}\n".format(tu[0], tu[1]))

    genes = ["eco:"+tu[0] for tu in regulon]
    prefix = os.path.join(out_dir, out_file.split(".")[0])

    getgeneInfo(genes, prefix)
    return regulon

################################################################################################

# This should be imported from code/methods/, but I can't figure out how to do that

def getgeneInfo(gene_call_list, output_prefix):

    for gene in gene_call_list:
        try:
            url = "http://rest.kegg.jp/get/{}".format(gene)
            with urllib.request.urlopen(url) as f:
                lines = f.read().decode('utf-8').splitlines()

            gene_info = OrderedDict([('locus_tag',''),
                                     ('name',''),
                                     ('definition',''),
                                     ('pathway',''),
                                     ('aa_seq','')]) # can add more if needed

            for i in range(len(lines)):
                field = lines[i].split()[0]
                if field == 'ENTRY':
                    gene_info['locus_tag'] = lines[i].split()[1]

                elif field == "NAME":
                    p = lines[i].split("NAME        ")[1]
                    gene_info['name'] = p

                elif field == "DEFINITION":

                    p = lines[i].split("DEFINITION  ")[1]
                    print(p)
                    gene_info['definition'] = p

                elif field == "PATHWAY":
                    p = lines[i].split("PATHWAY     ")[1]
                    count = 1
                    while lines[i+count].startswith("           ") and len(lines[i+count]) > 0:
                        p += lines[i+count]
                        count += 1
                    gene_info['pathway'] = ", ".join(p.split("            ")).rstrip(', ')

                elif field == "AASEQ":
                    p = ''
                    count = 1
                    while lines[i+count].startswith("           ") and len(lines[i+count]) > 0:
                        p += lines[i+count]
                        count +=1
                    gene_info['aa_seq'] = p.replace("            ", '')

            info_str = "\t".join([v for k,v in gene_info.items() if k in ["locus_tag", "name", "definition","pathway"]]).rstrip('\t')

            aa_str = ">{}\n{}".format(gene, gene_info['aa_seq'])

            # Want to write to file right away, in case connection breaks
            with open(output_prefix+"_info.tab", "a") as of:
                of.write(info_str+"\n")
            with open(output_prefix+".fasta", "a") as of:
                of.write(aa_str+"\n")
        except:
            with open(output_prefix+"_info.tab", "a") as of:
                of.write("{}\tERROR\n".format(gene.split(":")[1]))
            continue


#############################################################################################


def getgenes(genome, gene_list_txt, frmt = 'csv'):
    """
    :param genome: string, KEGG genome designation, ex. ecc for CFT073, eco for MG1655
    :param gene_list_txt: list of genes for which info/seq is needed, have to Gene entry name or accession (from KEGG API)
    :return: list KEGG API compatible gene requests, ex. genome:geneid
    """
    genes = []
    if frmt == 'csv':
        with open (gene_list_txt) as gl:
            genes = [line.rstrip().split(',')[0] for line in gl]
    return ["{}:{}".format(genome, gene) for gene in genes]



if __name__ == "__main__":

    gene_regulation_file = "/Users/annasintsova/git_repos/HUTI-RNAseq/data/RegulonDB/generegulation_tmp.txt"
    out_dir = "/Users/annasintsova/git_repos/HUTI-RNAseq/results/resource_allocation"
    regulon = findRegulon(gene_regulation_file, "crp", out_dir, modulation='activator')
    print(len(regulon))

