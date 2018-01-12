import datetime
import os

### FASTQC, TRIMMOMATIC, and ALIGNMENT with BOWTIE2

### All of these need to be adjusted to accomodate PE situation

def fastqc(files, out_dir = (datetime.datetime.now().strftime("%Y-%m-%d") + "_fastqc_results"), multiqc = True):
    #list of full paths to fastq files and output output directory
    script = ''
    for f_name in files:
        if "fastq" not in os.path.basename(f_name):
            print("{} not a fastq file, skipping".format(f_name))
            continue
        elif ".gz" in f_name:
            script += "/home/annasint/bin/FastQC/fastqc -o {} {} --extract\n".format(out_dir, f_name)
        else:
            script += "/home/annasint/bin/FastQC/fastqc -o {} {}\n".format(out_dir, f_name)
    if multiqc:
        report_name = os.path.join(out_dir, (datetime.datetime.now().strftime("%Y-%m-%d") + "_multiqc_report"))
        script += "multiqc {} --force --filename {}\n".format(out_dir, report_name)
    return script

def Trimmomatic(files, out_dir = (datetime.datetime.now().strftime("%Y-%m-%d") + "_trimmed_results")):
    # List of full paths to fastq files and output output directory, can handle .gz
    script = ''
    for f_name in files:
        file_name = os.path.basename(f_name).split(".")[0] + "_trimmed.fastq"
        f_out = os.path.join(out_dir, file_name)
        script += "java -jar /home/annasint/bin/Trimmomatic-0.36/trimmomatic-0.36.jar " \
                  "SE {} {} ILLUMINACLIP:/home/annasint/bin/Trimmomatic-0.36/adapters/TruSeq3-SE.fa:2:30:10:8:true" \
                  " SLIDINGWINDOW:4:15 MINLEN:40 HEADCROP:0\n".format(f_name, f_out)
    return script

def BowTieAlign(files, reference, out_dir, multi_ref=False, ix=True):

    # Generate an index.

    script = "cd {}\n".format(out_dir)
    ge_indices = []
    if multi_ref:  # build all the indices
        genomes = [os.path.join(reference, r) for r in os.listdir(reference)]  # reference needs to be path to folder

        for ge in genomes:
            bt2_base = os.path.basename(ge).split(".")[0] + "_index"
            if ix:
                script += "bowtie2-build {} {}\n".format(ge, bt2_base)
            ge_indices.append(bt2_base)
    else:
        # Reference is a full file path
        bt2_base = os.path.basename(reference).split(".")[0] + "_index"
        if ix:
            script += "bowtie2-build {} {}\n".format(reference, bt2_base)
        ge_indices.append(bt2_base)


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

        script += "bowtie2 -x {} -U {} -S {}\n".format(bt2_base, name, sam_name)
    return script

