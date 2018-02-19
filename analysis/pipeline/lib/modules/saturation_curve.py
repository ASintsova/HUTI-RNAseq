"""

testing script found on the internet, first need file with gene_id, length of gene, full length reads


"""

counts = "/Users/annasintsova/git_repos/HUTI-RNAseq/analysis/DE/data/HM01_UR_trimmed_sorted_counts"
gff = "/Users/annasintsova/git_repos/HUTI-RNAseq/analysis/DE/2018-01-16-sync-genomes/data/gff_files/HM01.gff"
subsample = "test_out/subsample.test"

# Read through counts, find prokka in gff, find length, add to new file


def prepareSubSampleFile(counts, gff, subsample):

    with open(subsample, "w") as ss:
        ss.write("pbid\tlength\tfl_count\n")
        with open(counts, "r") as cc:
            for line in cc:
                prokka = line.split()[0].strip()
                count = line.split()[1].strip()
                with open(gff, "r") as ref:
                    for ll in ref:
                        if prokka in ll:
                            start = int(ll.split("\t")[3])
                            end = int(ll.split("\t")[4])
                            gene_length = str(end - start)
                            ss.write("{}\t{}\t{}\n".format(prokka, gene_length, count))


prepareSubSampleFile(counts, gff, subsample)