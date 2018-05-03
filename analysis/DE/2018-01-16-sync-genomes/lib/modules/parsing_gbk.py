from Bio import SeqIO
import pandas as pd

# def cft073(parsing_func):
#
#     gbk ="/Users/annasintsova/git_repos/HUTI-RNAseq/analysis/DE/" \
#          "2018-01-16-sync-genomes/data/" \
#          "CFT073/GCA_000007445.1_ASM744v1_genomic.gbff"
#     genome = "CFT073"
#     def wrapper(csv_file):
#         return parsing_func(gbk, genome, csv_file)
#     return wrapper
#
# def MG1655(parsing_func):
#     gbk = "/Users/annasintsova/git_repos/HUTI-RNAseq/analysis/DE/" \
#                  "2018-01-16-sync-genomes/data/" \
#                   "MG1655/GCA_000005845.2_ASM584v2_genomic.gbff"
#     genome = "MG1655"
#     def wrapper(csv_file):
#         return parsing_func(gbk, genome, csv_file)
#     return wrapper

#@cft073
#@MG1655



def replaceProteinIdWithLocusTags(gbk_file, genome,table, feature_type = 'CDS',
                                  from_qual = "protein_id",
                                  to_qual = "locus_tag"):
    # Need to add logging
    # Parses genbank file
    answer = {}
    for gb_record in SeqIO.parse(gbk_file, "genbank"):
        for feature in gb_record.features:
            if feature.type == feature_type:
                if from_qual in feature.qualifiers:
                    # There should only be one locus_tag per feature, but there
                    # are usually several db_xref entries
                    if len(feature.qualifiers[from_qual]) > 1:
                        print("WARNING DUPLICATE ENTRIES FOR {}, skipping".format(feature))
                    for value in feature.qualifiers[from_qual]:
                        answer[value] = feature.qualifiers[to_qual][0]
                else:
                    print("No such qualifier in the {}".format(feature))
        #print(answer["AAN78540.1"])
        #return answer
    # Replaces values in the table file

    pgm = pd.read_csv(table, index_col=0)
    pgm.replace(to_replace={genome:answer}, inplace=True)
    #out_file = table.split(".")[0] + "_{}_edited.csv".format(genome) #when sure everything is running properly just replace the original file
    pgm.to_csv(table)



if __name__ == "__main__":
    table = "/Users/annasintsova/git_repos/HUTI-RNAseq/data/get_homologs_output/C50_S90_e0_/run_C50_S90_e0__pan_C50_S90/pangenome_matrix_t0_crossRef.csv"
    print("CFT073")
    gbk = "/Users/annasintsova/git_repos/HUTI-RNAseq/analysis/DE/" \
                  "2018-01-16-sync-genomes/data/" \
              "CFT073/GCA_000007445.1_ASM744v1_genomic.gbff"
    genome = "CFT073"

    replaceProteinIdWithLocusTags(gbk, genome, table)
    print("MG1655")

    gbk = "/Users/annasintsova/git_repos/HUTI-RNAseq/analysis/DE/" \
                          "2018-01-16-sync-genomes/data/" \
                       "MG1655/GCA_000005845.2_ASM584v2_genomic.gbff"
    genome = "MG1655"

    replaceProteinIdWithLocusTags(gbk, genome, table)
