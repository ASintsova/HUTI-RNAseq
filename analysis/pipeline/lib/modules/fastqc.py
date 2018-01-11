import datetime
import os

### FASTQC and ALIGNMENT with BOWTIE2

### All of these need to be adjusted to accomodate PE situation

def fastqc(files, out_dir = (datetime.datetime.now().strftime("%Y-%m-%d") + "_fastqc_results"), multiqc = True):
    #list of fastq files and output output directory
    script = ''
    for f_name in files:
        if ".gz" in f_name:
            script += "/home/annasint/bin/FastQC/fastqc -o {} {} --extract\n".format(out_dir, f_name)
        else:
            script += "/home/annasint/bin/FastQC/fastqc -o {} {}\n".format(out_dir, f_name)
    if multiqc == True:
    #The name below needs to include outdir in it
        report_name = datetime.datetime.now().strftime("%Y-%m-%d") + "_multiqc_report"
        script += "multiqc {} --force --filename {}\n".format(out_dir, report_name)
    return script

def Trimmomatic(files, out_dir = (datetime.datetime.now().strftime("%Y-%m-%d") + "_trimmed_results")):
    #list of fastq files and output output directory
    script = ''
    suffix = "_trimmed.fastq"
    for f_name in files:
        file_name = f_name.split("/")[-1].split(".")[0]
        f_out = out_dir + "/" + file_name +suffix
        script += "java -jar /home/annasint/bin/Trimmomatic-0.36/trimmomatic-0.36.jar SE {} {} ILLUMINACLIP:/home/annasint/bin/Trimmomatic-0.36/adapters/TruSeq3-SE.fa:2:30:10:8:true SLIDINGWINDOW:4:15 MINLEN:40 HEADCROP:0\n".format(f_name, f_out)
    return script

def BowTieAlign(files, reference, out_dir, multi_ref=False):
    #generate an index. In the future should check if index already exists(maybe?)
    script = f"cd {out_dir}\n"
    ge_indices = []

    if multi_ref:#build all the indicids
        genomes = os.listdir(reference)#reference needs to be path to folder

        for ge in genomes:
            bt2_base = os.path.basename(ge).split(".")[0] + "_index"
            script += "bowtie2-build {} {}\n".format(os.path.abspath(ge), bt2_base)
            ge_indices.append(bt2_base)
    else:
        bt2_base = reference.split("/")[-1].split(".")[0] + "_index"
        script += "bowtie2-build {} {}\n".format(os.path.abspath(reference), bt2_base)
        ge_indices.append(bt2_base)
    for f_name in files:

        if 'fastq' not in f_name:
            print(f"Skipping {f_name}, not a fastq file")
            continue
        if f_name.endswith(".gz"):
            script+= "gunzip {}\n".format(os.path.abspath(f_name))
            name = f_name.split(".gz")[0]
        else:
            name = f_name

        sam_name = os.path.join(out_dir, os.path.basename(name).split('.')[0] + ".sam")

        if len(ge_indices) == 1:
            bt2_base = ge_indices[0]
        else:
            prefix = name.split("_")[0]
            bt2_base = [g for g in ge_indices if prefix in g][0]

        script += "bowtie2 -x {} -U {} -S {}\n".format(bt2_base, os.path.abspath(name), os.path.abspath(sam_name))
    return script

