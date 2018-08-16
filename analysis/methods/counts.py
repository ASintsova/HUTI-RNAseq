import os
import datetime as dt




def calculate_rpkms(count_files_dir, flagstat_directory, gff_folder, output_dir):

    count_files = [os.path.join(count_files_dir, f) for f in os.listdir(count_files_dir)]
    flagstat_summary = process_flagstat(flagstat_directory)
    rpkm_dict = {}


    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    for cf in count_files:
        if not cf.endswith("counts"):
            continue
        print(cf)
        genome = os.path.basename(cf).split("_")[0]

        gff_file = os.path.join(gff_folder, "{}.gff".format(genome))
        gff = process_gff(gff_file)

        N = int(flagstat_summary[genome][1]) # mapped reads
        with open(cf, "r") as fh:
            for line in fh:
                if line.startswith("__"):
                    continue
                gene = line.split()[0].rstrip()
                print(genome)
                print(gff[gene])
                count = int(line.split()[1].rstrip())
                locus_tag = gff[gene][0]
                gene_length = int(gff[gene][1])
                rpkm = calculate_rpkm(count, N, gene_length)
                rpkm_dict[locus_tag] = rpkm

            suffix = os.path.basename(cf).split("_trimmed")[0]
            rpkm_fn = os.path.join(output_dir,
                    ("{}_{}_rpkm.csv".format(dt.datetime.today().strftime("%Y-%m-%d"),
                                            suffix)))
            with open(rpkm_fn, "w") as fo:
                for key, val in rpkm_dict.items():
                    fo.write("{},{}\n".format(key, val))

def process_gff(gff_file):
    gff = {}
    with open(gff_file, "r") as fh:
        for line in fh:
            line = line.strip()
            if line.startswith("#") or line.startswith(" ") or len(line) == 0:
                continue
            elif line.split('\t')[2] != "CDS":
                continue
            else:
                start = int(line.split("\t")[3].strip())
                end = int(line.split("\t")[4].strip())
                gene_length = end - start

                ID = line.split("\t")[-1].split(";")[0].strip("gene_id=")

                gff[ID] = (ID, gene_length)
    return gff

def process_flagstat(flagstat_directory):

    fstats = [os.path.join(flagstat_directory, f) for f in os.listdir(flagstat_directory)]
    flag_summary = {}
    for fstat in fstats:
        if "flagstat" not in fstat:
            continue
        with open(fstat, "r") as fh:
            sample = os.path.basename(fstat).split('_')[0]
            summary = fh.readlines()
            flag_summary[sample] = (summary[0].split()[0], summary[4].split()[0])

    return flag_summary




def calculate_rpkm(count, N, gene_length):
    return round((10 ** 9 * count) / (N * gene_length), 2)







if __name__ == "__main__":

    flagstat_dir = "/Users/annasintsova/git_repos/HUTI-RNAseq/data/counts/flagstat"
    counts_path = "/Users/annasintsova/git_repos/HUTI-RNAseq/data/counts"
    gff_folder = "/Users/annasintsova/git_repos/HUTI-RNAseq/data/annotations/gff_files"
    output_dir = "/Users/annasintsova/git_repos/HUTI-RNAseq/data/counts/rpkms"
    calculate_rpkms(counts_path, flagstat_dir, gff_folder, output_dir)

