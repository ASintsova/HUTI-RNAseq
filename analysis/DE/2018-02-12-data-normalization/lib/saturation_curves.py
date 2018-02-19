
import pysam
import random
import subprocess
import pandas as pd
import configparser
import sampleNames
import os


# How would I pick random reads with replacement
# Can use samtools for the same purpose (-s)

def subSampleBam(bam_file, subsample_bam, fraction): # without replacement

    bam = pysam.AlignmentFile(bam_file)
    output = pysam.AlignmentFile(subsample_bam, "wb", template=bam) #subsampled file
    for read in bam.fetch():
        if random.random() < fraction:
            output.write(read)
    bam.close()
    print("{} subsampled".format(bam_file))
    output.close()


def subSampleSamtools(bam_file, subsample_bam, fraction): # for some reason generates same subsample everytime
    ###
    cmd =  "/Users/annasintsova/tools/samtools-1.7/samtools view -b -s {} -o {} {}".format(fraction, subsample_bam, bam_file)
    subprocess.call(cmd, shell=True)
    print("{} subsampled".format(bam_file))


def coverageBedTools(bam_file, gff_file, coverage_file):

    cmd = "/Users/annasintsova/tools/bedtools2/bin/bedtools coverage -a {} -b {} > {}".format(gff_file, bam_file, coverage_file)
    subprocess.call(cmd, shell=True) # Read about security concerns
    print("Coverage calculated for {}".format(bam_file))
    # Produces gff file with an extra column that shows coverage

def highCoveredGenes(coverage_file, cutoff, reads = True):

    df = pd.read_csv(coverage_file,  sep="\t", header=None, names=["chr", "soft", "feat", "start", "end", "score", "str", "frame", "attr", "num_reads", "num_bases_covered", "len_gene", "coverage" ])
    if reads:
        return df.loc[df["num_reads"] >= cutoff].shape[0]
    else:
        return df.loc[df["coverage"] >= cutoff].shape[0]




def calculateSaturationCurveData(bam_file, gff_file, subsample_range, cutoff, reads, output_dir, iter):

    results = []
    for fraction in subsample_range:

        subsample_bam = os.path.basename(bam_file).split(".")[0] + "_{}_subsample.bam".format(str(fraction))
        subsample_path = os.path.join(output_dir, subsample_bam)
        coverage_file = os.path.basename(bam_file).split(".")[0] + "_{}_coverage.txt".format(str(fraction))
        coverage_path = os.path.join(output_dir, coverage_file)
        gene_num_in_iter = []
        for i in range(iter):
            #subSampleSamtools(bam_file, subsample_path, fraction)
            subSampleBam(bam_file, subsample_path, fraction)
            coverageBedTools(subsample_path, gff_file, coverage_path)
            num_high_genes = highCoveredGenes(coverage_path, cutoff, reads)
            #num_high_genes = i
            gene_num_in_iter.append(num_high_genes)
            print(gene_num_in_iter)
        results.append((fraction, sum(gene_num_in_iter)/len(gene_num_in_iter)))
        subprocess.call(["rm", subsample_path])
        subprocess.call(["rm", coverage_path])
    return results



def getDataforSaturationCurves(samples, config, subsample_range, cutoff, reads, output_file, i):
    out_dir = os.path.dirname(output_file)
    df = pd.DataFrame(index=[str(c * 100) + "%" for c in subsample_range])
    for sample in samples:
        genome = sample.split("_", 1)[0]
        bam_path = config.get("BAMS", "path")
        bam_file = os.path.join(bam_path, "{}_trimmed_sorted.bam".format(sample))
        gff_path = config.get("GFF", "path")
        gff_file = os.path.join(gff_path, "{}.gff".format(genome))
        results = calculateSaturationCurveData(bam_file, gff_file, subsample_range, cutoff, reads, out_dir, i)
        df[sample] = [y[1] for y in results]
    df.to_csv(os.path.join(output_file))


if __name__ == "__main__":

        samples = sampleNames.generateSampleNames(['HM06', 'HM07'], ["UR", "UTI"])
        print(samples)
        config = configparser.ConfigParser()
        config.read(os.path.dirname(os.path.abspath(__file__))+"/config")
        output_file = config.get("output", "file")
        getDataforSaturationCurves(samples, config, [0.001, 0.005, 0.01, 0.1, 0.5, 0.75, 1], 10, True, output_file, 5)

#Same problem as with gh program - need it to spit out intermediate results

#Also very slow :'(

#gff_file = "/Users/annasintsova/git_repos/HUTI-RNAseq/analysis/DE/2018-01-16-sync-genomes/data/gff_files/HM01.gff"
#bam_file = "/Users/annasintsova/git_repos/HUTI-RNAseq/data/alignments/HM01_UR_trimmed_sorted.bam"
#output_dir = "/Users/annasintsova/git_repos/HUTI-RNAseq/test"

#coverage_file = "test/test.1.coverage.txt"
#plotSaturationCurve(bam_file, gff_file, [0.1, 0.2, 0.3], 10, True, output_dir, iter=3)

#Before I forget output:

"""
10: number of reads overalapped with gene,
11: number of bases in the gene covered by reads
12: length of gene
13: % of gene covered

to the best of my ability to understand this

will try to look at this 2 different ways: 10 reads per gene / 10x coverage



"""

# ~/tools/bedtools2/bin/bedtools coverage -a $gff -b test/test.1.bam >test/test.1.coverage.txt
#After generating subsample call bedtools coverage on gff and subsampled bam. Need to first make sure bedtools is working

#Process bedtools output to see how many genes > 10x

# Add this number to pandas dataframe

# Plot with matplotlib



#################
"""
subsampling from iterator

def iter_sample_fast(iterable, samplesize):
    results = []
    iterator = iter(iterable)
    # Fill in the first samplesize elements:
    try:
        for _ in xrange(samplesize):
            results.append(iterator.next())
    except StopIteration:
        raise ValueError("Sample larger than population.")
    random.shuffle(results)  # Randomize their positions
    for i, v in enumerate(iterator, samplesize):
        r = random.randint(0, i)
        if r < samplesize:
            results[r] = v  # at a decreasing rate, replace random items
    return results
    
    
"""