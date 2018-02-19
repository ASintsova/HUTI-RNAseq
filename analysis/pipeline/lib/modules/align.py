import os
import logging
import configparser

def fastqAligner(files, reference, out_dir,  analysis_name, config,multi_ref=False, ix=True):
    if config.get("aligner", "name") == "bowtie2":
        return BowTieAlign(files, reference, out_dir, analysis_name, config, multi_ref, ix)
    else:
        return "Work in progress..."


def BowTieAlign(files, reference, out_dir, analysis_name, config, multi_ref=False, ix=True):
    log = logging.getLogger("{}.align.BowTieAlign".format(analysis_name))
    align_bin = config.get("aligner", "bin")

    # Generate an index.
    log.info("Generating an index")
    script = "cd {}\n".format(out_dir)
    ge_indices = []
    if multi_ref:  # build all the indices
        genomes = [os.path.join(reference, r) for r in os.listdir(reference)]  # reference needs to be path to folder

        for ge in genomes:
            bt2_base = os.path.basename(ge).split(".")[0] + "_index"
            if ix:
                script += "{}/bowtie2-build {} {}\n".format(align_bin, ge, bt2_base)
            ge_indices.append(bt2_base)
    else:
        # Reference is a full file path
        bt2_base = os.path.basename(reference).split(".")[0] + "_index"
        if ix:
            script += "{}/bowtie2-build {} {}\n".format(align_bin,reference, bt2_base)
        ge_indices.append(bt2_base)

    log.info("Done... Aligning files")
    for f_name in files:

        if 'fastq' not in os.path.basename(f_name):
            print("Skipping {}, not a fastq file".format(f_name))
            continue
        if f_name.endswith(".gz"):
            script += "gunzip {}\n".format(f_name)
            name = f_name.split(".gz")[0]
        else:
            name = f_name

        sam_name = os.path.join(out_dir, os.path.basename(name).split('.')[0] + ".sam")

        if len(ge_indices) == 1:
            bt2_base = ge_indices[0]
        else:
            prefix = os.path.basename(name).split("_")[0]
            bt2_base = [g for g in ge_indices if prefix in g][0]

        script += "{}/bowtie2 -x {} -U {} -S {}\n".format(align_bin, bt2_base, name, sam_name)
    return script

