import datetime
import os
import argparse
# import configparser
import subprocess
from modules import generateJobs as gJ
from modules import fastqc as fq
from argparse import RawTextHelpFormatter


def getArgs():
    parser = argparse.ArgumentParser("RNASeq Pipeline\n", formatter_class=RawTextHelpFormatter)
    parser.add_argument("-a", "--analysis", help="Analysis options", required=True)
    parser.add_argument("-i", "--input", nargs='+', help="List of files", required=True)
    parser.add_argument("-o", "--out_dir", help="Location of the output directory", required=True)
    parser.add_argument('-ref', '--reference_genome', help="Reference genome for alignments", required=False)
    parser.add_argument('--no_index', help="Don't build bowtie2 index again", action="store_true", required=False)

    return parser


# def ConfigSectionMap(section, Config):
#     dict1 = {}
#     if not Config.has_section(section):
#         # keep_logging('ERROR: Please Check the section name: \'{}\' in config file'.format(section), 'Please Check the section name: \'{}\' in config file'.format(section), logger, 'exception')
#         print ("ERROR: Please Check the section name: %s in config file" % section)
#         exit()
#     options = Config.options(section)
#     for option in options:
#         try:
#             dict1[option] = Config.get(section, option)
#             if dict1[option] == -1:
#                 DebugPrint("skip: %s" % option)
#         except:
#             print("exception on %s!" % option)
#             dict1[option] = None
#     return dict1


def pipeline():  # will eventually add a logger

    args = getArgs().parse_args()
    analysis = args.analysis
    if os.path.isdir(args.input[0]):
        files = [os.path.join(args.input[0],fi) for fi in os.listdir(args.input[0])]
    elif args.input:
        files = args.input
    else:
        raise IOError
    out_dir = args.out_dir
    subprocess.call(["mkdir", "-p", out_dir])
    assert os.path.exists(out_dir)

    # Config = configparser.ConfigParser()
    # Config.read(os.path.dirname(os.path.abspath(__file__)) + "/config")
    if analysis == "fastqc":
        pbs_name = os.path.join(out_dir, datetime.datetime.now().strftime("%Y-%m-%d") + '_fastqc.pbs')
        gJ.generatePBSScript(pbs_name, fq.fastqc(files, out_dir))
        analysis_name = "FastQC Job"
        # need to add logging to this
        # should i delete these pbs scripts after?
    elif analysis == 'trim':
        pbs_name = os.path.join(out_dir, datetime.datetime.now().strftime("%Y-%m-%d") + '_trimmomatic.pbs')
        gJ.generatePBSScript(pbs_name, fq.Trimmomatic(files, out_dir))
        analysis_name = "Trimmomatic Job"
    elif analysis == 'align':
        ref = args.reference_genome
        ix = False if args.no_index else True
        assert ref
        if os.path.isdir(ref):
            multi_ref = True
        else:
            multi_ref = False
        pbs_name = os.path.join(out_dir, datetime.datetime.now().strftime("%Y-%m-%d") + '_bowtie2.pbs')
        gJ.generatePBSScript(pbs_name, fq.BowTieAlign(files, ref, out_dir, multi_ref, ix))
        analysis_name = "BowTie2 Alignment Job"
    else:
        return "Wrong Answer"
    #subprocess.call(["qsub", pbs_name])
    return "{} submitted!".format(analysis_name)


print(pipeline())
