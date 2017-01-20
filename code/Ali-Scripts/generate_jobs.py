__author__ = 'alipirani'

import os
import argparse
import re

parser = argparse.ArgumentParser(description='Generate Assembly PBS scripts for Assembly Pipeline')
parser.add_argument('-dir', action='store', dest="dir", help='directory of fastq files')
parser.add_argument('-filenames', action='store', dest="filenames", help='These file should contain name of forward fastq files. \
One file per line. \
These can be obtained by running <ls *R1*.gz > filenames>')
parser.add_argument('-out_dir', action='store', dest="out_dir", help='Provide a path where you want to save the assembly output')
parser.add_argument('-pipeline', action='store', dest="pipeline", help='Generating Jobs for which pipeline? varcall/assembly/new_assembly/rna/ptr')
parser.add_argument('-reference', action='store', dest="reference", help='Reference Genome to be used for pipeline')
parser.add_argument('-reference_key', action='store', dest="reference_key", help='Reference Genome to be used for each sample')
args = parser.parse_args()

sample_files_directory = args.dir
filenames = args.filenames
main_out_directory = args.out_dir
pipeline = args.pipeline

# Save filenames in an array
filenames_array = []
with open(filenames) as fp:
    for line in fp:
        line = line.strip()
        line = sample_files_directory + "/" + line
        filenames_array.append(line)



# These is the Pbs initial resources assignment
# Spades requires high memory for base error correction
# Please change the email address, flux account and memory requirement accordingly


def create_RNA_seq_jobs():
    reference_key_dict = {}
    #with open(args.reference_key) as fp:
    fp = open(args.reference_key, 'r')
    for line in fp:
        line = line.strip()
        line_split = line.split('\t')
        line_split_file = sample_files_directory + "/" + line_split[0]
        reference_key_dict[line_split_file] = line_split[1]
    Pbs_model_lines = "#PBS -M apirani@med.umich.edu\n" \
    "#PBS -m abe\n" \
    "#PBS -V\n" \
    "#PBS -l nodes=1:ppn=4,pmem=4000mb,walltime=12:00:00\n" \
    "#PBS -q fluxod\n" \
    "#PBS -A esnitkin_fluxod\n" \
    "#PBS -l qos=flux\n"
    for file in reference_key_dict.keys():
        filename_base = os.path.basename(file)
        #first_part_split = filename_base.split('.')
        first_part_split = filename_base.split('_')
	first_part = str(first_part_split[0]) + "_"
        #first_part = str(first_part_split[0]) + "_"
	sample_output_dir = first_part
        #first_file = file
        #second_part = filename_base.replace("_1_", "_2_")
        #second_part = filename_base.replace("_R1", "_R2")
        #second_part = filename_base.replace("forward", "reverse")
        #second_part = filename_base.replace("1_sequence", "2_sequence")
        #second_file = sample_files_directory + "/" + second_part

        job_name = "./" + first_part + ".pbs"
        # Change these directory path to where your
        cd_command = "cd /nfs/esnitkin/bin_group/pipeline/RNA_seq_pipeline/"
        command = "/home/apirani/anaconda/bin/python pipeline.py -PE1 %s -o %s/%s -analysis %s -index %s -type SE -s yes" % (file, main_out_directory, first_part, first_part, reference_key_dict[file])
        with open(job_name, 'w') as out:
            job_title = "#PBS -N %s" % first_part
            out.write(job_title+'\n')
            out.write(Pbs_model_lines+'\n')
            out.write(cd_command+'\n')
            out.write(command+'\n')
    








# Function to create assembly pbs scripts for each fastq files
def create_assembly_jobs():
    Pbs_model_lines = "#PBS -M apirani@med.umich.edu\n" \
        "#PBS -m abe\n" \
        "#PBS -V\n" \
        "#PBS -l nodes=1:ppn=4,mem=47000mb,walltime=12:00:00\n" \
        "#PBS -q fluxod\n" \
        "#PBS -A esnitkin_fluxod\n" \
        "#PBS -l qos=flux\n"	
    for file in filenames_array:
        filename_base = os.path.basename(file)
        first_part_split = filename_base.split('.')
        first_part = first_part_split[0]
        sample_output_dir = first_part
        first_file = file
        #second_part = filename_base.replace("_1", "_2")
        #sample_output_dir = first_part.replace("_1", "")
	second_part = filename_base.replace("_R1_", "_R2_")
	#second_part = filename_base.replace("forward", "reverse")
        #second_part = filename_base.replace("1_sequence", "2_sequence")
	second_file = sample_files_directory + "/" + second_part

        job_name = "./" + first_part + ".pbs"
        # Change these directory path to where your
        cd_command = "cd /nfs/esnitkin/bin_group/pipeline/assembly_umich/"
	#cd_command = "cd /home2/apirani/bin/assembly_umich/"
        command = "~/anaconda/bin/python pipeline.py -f1 %s -f2 %s -o %s/%s -start_step 1 -end_step 4 -A spades -reference %s" % (first_file, second_file, main_out_directory, sample_output_dir, args.reference)
        with open(job_name, 'w') as out:
            job_title = "#PBS -N %s" % first_part
            out.write(job_title+'\n')
            out.write(Pbs_model_lines+'\n')
            out.write(cd_command+'\n')
            out.write(command+'\n')


def create_varcall_jobs():
    Pbs_model_lines = "#PBS -M apirani@med.umich.edu\n" \
    "#PBS -m abe\n" \
    "#PBS -V\n" \
    "#PBS -l nodes=1:ppn=4,pmem=4000mb,walltime=12:00:00\n" \
    "#PBS -q fluxod\n" \
    "#PBS -A esnitkin_fluxod\n" \
    "#PBS -l qos=flux\n"
    for file in filenames_array:
        filename_base = os.path.basename(file)
	print filename_base
        first_part_split = filename_base.split('_R1_001.fastq.gz')
        first_part = str(first_part_split[0])
	#print first_part
	#first_part_split = filename_base.split('_')
        #first_part = str(first_part_split[0]) + "_"
        sample_output_dir = first_part
        first_file = file
        #second_part = first.replace("_1", "_2")
        second_part = filename_base.replace("_R1_", "_R2_")
        #second_part = filename_base.replace("_1.", "_2.")
        #second_part = filename_base.replace("forward", "reverse")
	#second_part = filename_base.replace("1_sequence", "2_sequence")
	second_file = sample_files_directory + "/" + second_part
	#print job_name
        job_name = "./" + first_part + ".pbs"
	print job_name
        # Change these directory path to where your
        cd_command = "cd /nfs/esnitkin/bin_group/pipeline/Github/varcall_umich"
        config = "/nfs/esnitkin/bin_group/pipeline/Github/varcall_umich/config"
        command = "/home/apirani/anaconda/bin/python pipeline.py -PE1 %s -PE2 %s -o %s/%s -analysis %s -index %s -type PE -config %s" % (first_file, second_file, main_out_directory, sample_output_dir, first_part, args.reference, config)
        with open(job_name, 'w') as out:
            job_title = "#PBS -N %s" % first_part
            out.write(job_title+'\n')
            out.write(Pbs_model_lines+'\n')
            out.write(cd_command+'\n')
            out.write(command+'\n')

def create_PTR_jobs():
    Pbs_model_lines = "#PBS -M apirani@med.umich.edu\n" \
    "#PBS -m abe\n" \
    "#PBS -V\n" \
    "#PBS -l nodes=1:ppn=4,pmem=4000mb,walltime=12:00:00\n" \
    "#PBS -q fluxod\n" \
    "#PBS -A esnitkin_fluxod\n" \
    "#PBS -l qos=flux\n"
    for file in filenames_array:
        filename_base = os.path.basename(file)
        first_part_split = filename_base.split('.')
        first_part = first_part_split[0]
        sample_output_dir = first_part
        first_file = file
        #second_part = first.replace("_1", "_2")
        #second_part = filename_base.replace("R1", "R2")
        #second_part = filename_base.replace("_1.", "_2.")
        #second_part = filename_base.replace("forward", "reverse")
	#second_part = filename_base.replace("1_sequence", "2_sequence")
	#second_file = sample_files_directory + "/" + second_part

        job_name = "./" + first_part + ".pbs"
        # Change these directory path to where your
        cd_command = "cd /nfs/esnitkin/bin_group/pipeline/PTR_analysis/"
        command = "/home/apirani/anaconda/bin/python pipeline.py -SE %s -o %s/%s -analysis %s -index CFT073" % (first_file, main_out_directory, sample_output_dir, first_part)
        with open(job_name, 'w') as out:
            job_title = "#PBS -N %s" % first_part
            out.write(job_title+'\n')
            out.write(Pbs_model_lines+'\n')
            out.write(cd_command+'\n')
            out.write(command+'\n')

# Function to create assembly pbs scripts for each fastq files
def create_new_assembly_jobs():
    Pbs_model_lines = "#PBS -M apirani@med.umich.edu\n" \
        "#PBS -m abe\n" \
        "#PBS -V\n" \
        "#PBS -l nodes=1:ppn=4,mem=47000mb,walltime=12:00:00\n" \
        "#PBS -q fluxod\n" \
        "#PBS -A esnitkin_fluxod\n" \
        "#PBS -l qos=flux\n"
    for file in filenames_array:
        filename_base = os.path.basename(file)
        #first_part_split = filename_base.split('.')
        first_part_split = filename_base.split('_L001_R1_001.fastq.gz')
	first_part = first_part_split[0]
        sample_output_dir = first_part
        first_file = file
        second_part = filename_base.replace("_R1_", "_R2_")
        second_file = sample_files_directory + "/" + second_part
        job_name = "./" + first_part + ".pbs"
        # Change these directory path to where your
        cd_command = "cd /nfs/esnitkin/bin_group/pipeline/Github/assembly_umich/"
        config = "/nfs/esnitkin/bin_group/pipeline/Github/assembly_umich/config"
	    #cd_command = "cd /home2/apirani/bin/assembly_umich/"
        if args.reference:
	    command = "~/anaconda/bin/python pipeline.py -f1 %s -f2 %s -o %s/%s -start_step 1 -end_step 4 -A spades -reference %s -type PE -analysis %s -config %s" % (first_file, second_file, main_out_directory, sample_output_dir, args.reference, first_part, config)
        else:
	    command = "~/anaconda/bin/python pipeline.py -f1 %s -f2 %s -o %s/%s -start_step 1 -end_step 4 -A spades -type PE -analysis %s -config %s" % (first_file, second_file, main_out_directory, sample_output_dir, first_part, config)	    
	with open(job_name, 'w') as out:
            job_title = "#PBS -N %s" % first_part
            out.write(job_title+'\n')
            out.write(Pbs_model_lines+'\n')
            out.write(cd_command+'\n')
            out.write(command+'\n')









if pipeline == "varcall":
    create_varcall_jobs()
elif pipeline == "assembly":
    create_assembly_jobs()
elif pipeline == "new_assembly":
    create_new_assembly_jobs()
elif pipeline == "rna":
    create_RNA_seq_jobs()
elif pipeline == "ptr":
    create_PTR_jobs()
else:
    print "Please procide pipeline argument"

#second_part = filename_base.replace("_1", "_2")
#sample_output_dir = first_part.replace("_1", "")

#second_part = filename_base.replace("forward", "reverse")
#second_part = filename_base.replace("1_sequence", "2_sequence")
