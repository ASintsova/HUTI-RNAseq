import configparser
import os

def sam2bam(files, out_dir, sort=True):
    script = ''
    samtools_bin = "/home/annasint/bin/samtools-1.2/samtools"
    for fh in files:
        bam = os.path.join(out_dir, (os.path.basename(fh).split(".")[0] + ".bam"))
        fl_st = os.path.join(out_dir, (os.path.basename(fh).split(".")[0] + "_flagstat.txt"))
        script += "{0} view -Sb {1} > {2}\n{0} flagstat {2} > {3}\n".format(samtools_bin, fh, bam, fl_st)
        if sort:
            sorted_name = bam.split(".")[0] + "_sorted.bam"
            script += "{0} sort -o {1} -O bam {2}\n".format(samtools_bin, sorted_name, bam)
    #keep_logging('SAM to BAM Conversion', 'SAM to BAM Conversion', logger, 'info')
    #keep_logging(cmd, cmd, logger, 'debug')
    return script


def countReads(files, out_dir, ref,
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
        prefix = os.path.basename(fh.split("_"))[0]
        if len(annot_files) > 1:
            annot = [an for an in annot_files if prefix in an][0]
        else:
            annot = annot_files[0]
        count_file = os.path.join(out_dir, (os.path.basename(fh).split(".")[0] + "_counts"))
        script += "htseq-count -f {0} -r {1} -m {2}" \
                  " -s {3} -t {4} {5} {6} > {7}".format(form, order, mode, stranded, feature, fh, annot, count_file)

    return script
