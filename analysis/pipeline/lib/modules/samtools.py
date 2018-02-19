import configparser
import os
import logging

def sam2bam(files, out_dir, config, sort=True):
    script = ''
    samtools_bin = config.get("samtools", "bin")

    for fh in files:
        bam = os.path.join(out_dir, (os.path.basename(fh).split(".")[0] + ".bam"))
        fl_st = os.path.join(out_dir, (os.path.basename(fh).split(".")[0] + "_flagstat.txt"))
        script += "{0} view -b -o {1} -@ 4 {2}\n" \
                  "{0} flagstat -@ 4 {1} > {3}\n".format(samtools_bin, bam, fh, fl_st)
        if sort:
            sorted_name = bam.split(".")[0] + "_sorted.bam"
            script += "{0} sort -o {1} -@ 4 {2}\n" \
                      "{0} index -@ 4 {1}\n".format(samtools_bin, sorted_name, bam)
    return script


def countReads(files, out_dir, ref, config,
                form = "bam", order = "pos",
                mode = "union", stranded = "yes",
                feature = "CDS"):
    # htseq_bin = ""
    script = ""
    if os.path.isdir(ref):
        annot_files = [os.path.join(ref, r) for r in os.listdir(ref)]
    else:
        annot_files = [ref]
    for fh in files:
        prefix = os.path.basename(fh).split("_")[0]
        if len(annot_files) > 1:
            annot = [an for an in annot_files if prefix in an][0]
        else:
            annot = annot_files[0]
        count_file = os.path.join(out_dir, (os.path.basename(fh).split(".")[0] + "_counts"))
        script += "htseq-count -f {0} -r {1} -m {2}" \
                  " -s {3} -t {4} {5} {6} > {7}\n".format(form, order, mode, stranded, feature, fh, annot, count_file)

    return script
