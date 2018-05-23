import datetime
import os
import argparse
import configparser

import subprocess
from methods import generateJobs as gJ
from methods import pacbio_assembly

from argparse import RawTextHelpFormatter


def getArgs():
    parser = argparse.ArgumentParser("RNASeq Pipeline\n", formatter_class=RawTextHelpFormatter)
    parser.add_argument("-a", "--analysis", help="Analysis options", required=True)
    parser.add_argument("-i", "--input", nargs='+', help="List of files", required=True)
    parser.add_argument("-o", "--out_dir", help="Location of the output directory", required=True)
    parser.add_argument('-ref', '--reference_genome', help="Reference genome for alignments", required=False)
    parser.add_argument('--no_index', help="Don't build bowtie2 index again", action="store_true", required=False)#explore metavar

    return parser


def pipeline():

    args = getArgs().parse_args()
    analysis = args.analysis

    if os.path.isdir(args.input[0]):
        files = [os.path.join(os.path.abspath(args.input[0]),fi) for fi in os.listdir(args.input[0])]
    elif args.input:
        files = [os.path.abspath(fi) for fi in args.input]
    else:
        raise IOError

    out_dir = os.path.abspath(args.out_dir)
    subprocess.call(["mkdir", "-p", out_dir])

    if os.path.exists(out_dir):
        print("Output directory found...")
    else:
        raise IOError

    config = configparser.ConfigParser()
    config.read("code/methods/config")

    # Need to pass the config to the various modules
#####################################################################################################
    if analysis == 'pacbio':

        scripts = pacbio_assembly.call_canu(files, config, out_dir)

        for script in scripts:
            pbs_name = os.path.join(out_dir, datetime.datetime.now().strftime("%Y-%m-%d") + '_' + script[0] +'.pbs')
            pbs = gJ.generatePBSScript(pbs_name, script[1],
                                                 pmem='47gb', walltime="76:00:00")
            #subprocess.call(["qsub", pbs_name])
            print(pbs_name)
            print(script[1])


####################################################################################################
"""


    if analysis == "fastqc":

        pbs_name = os.path.join(out_dir, datetime.datetime.now().strftime("%Y-%m-%d") + '_fastqc.pbs')
        gJ.generatePBSScript(pbs_name, fq.fastqc(files, out_dir, config))
        analysis_name = "FastQC Job"

        # should i delete these pbs scripts after?
    elif analysis == 'trim':
        pbs_name = os.path.join(out_dir, datetime.datetime.now().strftime("%Y-%m-%d") + '_trimmomatic.pbs')
        gJ.generatePBSScript(pbs_name, fq.Trimmomatic(files, out_dir, config))
        analysis_name = "Trimmomatic Job"

    elif analysis == 'align':
        log.info("Calling align module, fastqAligner function")
        ref = os.path.abspath(args.reference_genome)
        ix = False if args.no_index else True
        assert ref
        if os.path.isdir(ref):
            multi_ref = True
        else:
            multi_ref = False
        pbs_name = os.path.join(out_dir, datetime.datetime.now().strftime("%Y-%m-%d") + '_bowtie2.pbs')
        gJ.generatePBSScript(pbs_name, align.fasqAligner(files, ref, out_dir, multi_ref, ix, config))
        analysis_name = "Alignment Job"
    elif analysis == 'sam2bam':
        log.info("Calling samtools module")
        pbs_name = os.path.join(out_dir, datetime.datetime.now().strftime("%Y-%m-%d") + '_sam2bam.pbs')
        gJ.generatePBSScript(pbs_name, samtools.sam2bam(files, out_dir, config))
        analysis_name = "Sam to Bam Conversion Job"
    elif analysis == 'count':
        ref = os.path.abspath(args.reference_genome)
        pbs_name = os.path.join(out_dir, datetime.datetime.now().strftime("%Y-%m-%d") + '_htseq_count.pbs')
        gJ.generatePBSScript(pbs_name, samtools.countReads(files, out_dir, ref, config))
        analysis_name = "HTseq Count Job"
    else:
        return "Wrong Answer"
    #subprocess.call(["qsub", pbs_name])
    return "{} submitted!".format(analysis_name)


"""

if __name__ == '__main__':
    pipeline()


