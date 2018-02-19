import datetime
import os
import logging

#logger = logging.getLogger("module1")

### FASTQC, TRIMMOMATIC

### All of these need to be adjusted to accomodate PE situation


def fastqc(files, out_dir,  analysis_name, config, multiqc = True):
    #list of full paths to fastq files and output output directory
    log = logging.getLogger("{}.fp.fastqc".format(analysis_name))

    script = ''
    log.info("Checking the fastq files, writting fastqc commands")
    for f_name in files:
        if "fastq" not in os.path.basename(f_name):
            log.info("{} not a fastq file, skipping".format(f_name))
            continue
        elif ".gz" in f_name:
            script += "/home/annasint/bin/FastQC/fastqc -o {} {} --extract\n".format(out_dir, f_name)
        else:
            script += "/home/annasint/bin/FastQC/fastqc -o {} {}\n".format(out_dir, f_name)
    log.info("Writting multiqc command")
    if multiqc:
        report_name = os.path.join(out_dir, (datetime.datetime.now().strftime("%Y-%m-%d") + "_multiqc_report"))
        script += "multiqc {} --force --filename {}\n".format(out_dir, report_name)
    log.info("Fastqc script generated")
    return script

def Trimmomatic(files, out_dir = (datetime.datetime.now().strftime("%Y-%m-%d") + "_trimmed_results")):
    # List of full paths to fastq files and output output directory, can handle .gz
    log = logging.getLogger("{}.fp.trimmomatic".format(analysis_name))
    script = ''
    log.info("Writting trimmomatic script")
    for f_name in files:
        file_name = os.path.basename(f_name).split(".")[0] + "_trimmed.fastq"
        f_out = os.path.join(out_dir, file_name)
        script += "java -jar /home/annasint/bin/Trimmomatic-0.36/trimmomatic-0.36.jar " \
                  "SE {} {} ILLUMINACLIP:/home/annasint/bin/Trimmomatic-0.36/adapters/TruSeq3-SE.fa:2:30:10:8:true" \
                  " SLIDINGWINDOW:4:15 MINLEN:40 HEADCROP:0\n".format(f_name, f_out)
    log.info("Trimmomatic script generated")
    return script