import pandas as pd
import sys
import sqlite3

sys.path.append('/Users/annasintsova/git_repos/HUTI-RNAseq/code/methods')
import keggAPI
import helpers

def get_bnums(external_db_link_file):
    eck_to_bnums = {}
    with open(external_db_link_file, "r") as fh:
        for line in fh:
            if line.startswith(("#")):
                continue
            else:
                ids = line.strip().split()
                if ids[1] == "GID000000044":
                    eck_to_bnums[ids[0]] = ids[2]
    return eck_to_bnums


def make_gene_sets(gene_reg_file, t):
    gene_sets = {}
    with open(gene_reg_file , "r") as fh:
        for line in fh:
            if line.startswith(("#")):
                continue
            else:
                words = line.strip().split()
                regulator = words[0]

                regulated = words[6]
                kind = words[8]
                if kind == t:
                    if regulator not in gene_sets.keys():
                        gene_sets[regulator] = [regulated]
                    else:
                        gene_sets[regulator].append(regulated)
                    # todo refactor
    return gene_sets


def write_gmt_gene_sets(gene_reg_file, external_db_link_file,
                        out_file, t="activator"):
    gene_sets = make_gene_sets(gene_reg_file, t)
    eck_to_bnums = get_bnums(external_db_link_file)
    with open(out_file, "w") as fo:
        for key, val in gene_sets.items():
            #print(key)
            #print(val)
            try:
                key_bnum = eck_to_bnums[key]
            except:
                key_bnum = key
            val_bnums = ""
            for v in val:
                try:
                    v_bnum = eck_to_bnums[v]
                except:
                    v_bnum = v

                val_bnums += "{}\t".format(v_bnum)
            fo.write("{}\t{}\t{}\n".format(key, key_bnum, val_bnums.strip("\t")))


def get_tf_name(gsea_out_csv, external_db_link_file, output_prefix):
    df = pd.read_csv(gsea_out_csv, index_col=0)
    ids = df.index
    print(ids)
    eck_to_bnums = get_bnums(external_db_link_file)

    bnums = [eck_to_bnums[eck] for eck in ids]
    gnames = keggAPI.get_gene_names("eco", bnums)
    keggAPI.get_gene_info(gnames, output_prefix)


def final_tables(gsea_out_csv, kegg_info_tab, external_db_link_file, out_file):
    df = pd.read_csv(gsea_out_csv, index_col=0)
    ids = df.index
    eck_to_bnums = get_bnums(external_db_link_file)
    bnums = [eck_to_bnums[eck] for eck in ids]
    print(bnums)
    info_tab = pd.read_table(kegg_info_tab, sep="\t", names=["bnum", "Name", "function", "pathway"]).set_index(["bnum"])
    names = [info_tab.loc[b].Name for b in bnums]
    df.index = names
    df.to_csv(out_file)

######################################################################################################



if __name__ == "__main__":

    external_link_file = "/Users/annasintsova/git_repos/HUTI-RNAseq/data/" \
                         "regulondb-9.4/object_external_db_link.txt"
    gene_regulation_file = "/Users/annasintsova/git_repos/HUTI-RNAseq/data/" \
                           "regulondb-9.4/generegulation_tmp.txt"
    out_f = "/Users/annasintsova/git_repos/HUTI-RNAseq/results/" \
             "regulon_analysis/2018-07-17-t-test-a-ur-edited.csv"
    #links = get_bnums(external_link_file)
    # write_gmt_gene_sets(gene_regulation_file, external_link_file,
    #                     o_file, t="repressor")
    csv = "/Users/annasintsova/git_repos/HUTI-RNAseq/results/regulon_analysis/2018-07-17-t-test-a-ur.csv"
    kegg_info = "/Users/annasintsova/git_repos/HUTI-RNAseq/results/regulon_analysis/act_ur_info.tab"


    final_tables(csv, kegg_info, external_link_file, out_f)