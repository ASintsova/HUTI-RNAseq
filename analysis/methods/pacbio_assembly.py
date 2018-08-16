#import generateJobs
import argparse
import configparser
import subprocess
import os


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", nargs='+', help="List of files", required=True)
    parser.add_argument("-o", "--out_dir", help="Location of the output directory", required=True)
    parser.add_argument('-gs', '--genome_size', help="Reference genome for alignments", required=False)

    return parser

def call_canu(files, config, out_dir, genome_size = 5.5):
    """
    files: full paths to pacbio reads
    :return: list of tuples (pbs_name, script) to be submitted for pacbio assembly
    """
    scripts = []

    canu = config.get("pacbio", "bin")
    for fi in files:
        prefix ="canu_" + os.path.basename(fi).split('.fastq')[0]
        script = "{} -p {} -d {} genomeSize = {}m "\
            "-pacbio-raw {}\n".format(canu, prefix, out_dir, str(genome_size), fi)
        scripts.append((prefix, script))
    return scripts


if __name__ == "__main__":

    config = configparser.ConfigParser()
    config.read("config")
    args = get_args().parse_args()
    if args.genome_size:
        scripts = call_canu(args.input, config, args.out_dir, args.genome_size)
    else:
        scripts = call_canu(args.input, config, args.out_dir)

    for script in scripts:
        pbs_name = os.path.join(args.out_dir, script[0] +".pbs")
        pbs = generateJobs.generatePBSScript(pbs_name, script[1],
                                pmem='47gb', walltime="76:00:00")
        subprocess.call(["qsub", pbs_name])