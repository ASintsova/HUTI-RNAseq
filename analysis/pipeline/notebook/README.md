#### 2018-01-10
##### Summary: Organizing data, renaming fastq files, putting stuff on flux. Working on extending pipeline to work for this dataset

* Went back to hard drive and copied, renamed, concatenated fastq files with `fastq_search_and_name_edit.py`

* Concatenated with the following. Not very elegant but seems to work

```bash
for f in *001.fastq*; do pr=${f%%_*}; cat $pr*001* >> $pr".fastq.gz"; rm $pr*001*; done
```

* Put on flux with (started with urine samples, but did all)

```bash
scp UR* annasint@flux-xfer.arc-ts.umich.edu:/scratch/hmobley_fluxod/annasint/HUTI-RNAseq/data/reads/fastq
```

* Now have a dilemma - want to run these through the beginnings of my pipeline, but it's all in a different repo. What's the best way to deal with this?? Have a separate version here? I don't want them to diverge. But will in the end be different pipelines. So maybe good for them to diverge? 

* also have to get my repos synced somehow, this might be quite a challenge at this point

* going to copy files over from chipseq stuff

* going to add dealing with folders instead of list of files as well as .gz extension

* well this is awkward, but I want to rename these files yet again, to make it easier down the line to align them

```bash
for f in *.fastq*
do
pr=${f%%.*}
sf=${f#*.}
echo $pr
echo $sf
cond=${pr//[0-9]/}
strain=${pr//[A-Z]/}
mv $f "HM"$strain"_"$cond"."$sf
done

```

* Updating fastq module so that it could take a folder of matching references and map each fastq to a different genome
* Needs to be tested


#### 2018-01-11
##### Summary: Running trim, fastqc, and bowtie jobs on RNAseq datase. Needs debugging

* Everything seems to be working. Going to try to sync everything with flux via git. Wish me luck!

* Python f string formatting does not work on flux, so stop showing off. Only had to change in 2 spots.
* Path problems in trimmomatic to do with the fact that's passing a directory: Fixed.
* Following job submitted:
```bash
python analysis/pipeline/lib/pipeline.py -a trim -i data/reads/fastq/ -o data/reads/trim_results
```
* and
```bash
python analysis/pipeline/lib/pipeline.py -a trim -i data/reads/fastq/ -o data/reads/trim_results

```
* and

```bash
python analysis/pipeline/lib/pipeline.py -a align -i data/reads/trim_results/ -o data/human/ceacam3_align -ref data/ref/human/ceacam3_transcript_variant_1.fasta

```
* and

```bash
python analysis/pipeline/lib/pipeline.py -a align -i data/reads/trim_results -o data/bt2_align -ref data/ref/huti
```
This last one is not working, need to use path.join not path.abspath. Too tired right now


#### 2018-01-12
##### Summary: ASM abstract. Looked at multiqc report. Fixed a bunch of path related bugs. fasqc, trim and align modules working

* multiqc report still is appearing in the wrong place (root), instead of outdir. Might want to look into the parameters for multiqc

* dowload the report. Overall looks fine. I'm concerned about combining reads for re-sequence files. Want to see if better off analyzing seperately. If results will be differnt, so have to go back to uploading files. Yay. 
```bash
 scp annasint@flux-xfer.arc-ts.umich.edu:/scratch/hmobley_fluxod/annasint/HUTI-RNAseq/2018-01-11_multiqc_report.html ~/git_repos/HUTI-RNAseq/analysis/pipeline/
```

* Submitted all the alignment jobs. (Human completed yesterday - looks fine in terms of no human reads in urine samples, but need next module to actually look at the alignment files)

* Wrote ASM abstract


#### 2018-01-16
##### Summary: Wrote samtools and htseq modules. Need to be tested

* Using old version of samtools, need to make sure the commands actually work
* Need to make sure strandiness works


#### 2018-01-17
##### Summary: Testing samtools and htseq modules

* Uploaded original files that were sequenced 2x, going through the whole pipeline with them. First going to rename them. 

* Just as note seq1 and seq2 were used randomly, do not actually refer to first or second sequencing event. Easy enough to check in the future if need be. 

```bash
for f in *seq*
do
pr=${f%%_*}
sf=${f#*_}
echo $pr
echo $sf
cond=${pr//[0-9]/}
strain=${pr//[A-Z]/}
mv $f "HM"$strain"_"$cond"_"$sf
done

```

* So now going to run trim job, and fastqc job as before. Fingers crossed. 
```bash
reads=data/reads/fastq/*seq*
out_dir=out_dir=data/reads/trim_results/
python analysis/pipeline/lib/pipeline.py -a trim -i $reads -o $out_dir

```

* Then 

```bash
sam=data/bt2_align/*sam
out_dir=data/bt2_align/bam
python analysis/pipeline/lib/pipeline.py -a sam2bam -i $sam -o $out_dir
```

* Everything seems to be running so far. 


#### 2018-01-18
##### Summary: Running more flux jobs with double sequenced samples. Running htseq-count.

* Sam2bam module seems to be working. Slightly modified sorting command to be more uptodate.
* Going to test htseq count today.

* Jobs to be submitted:
    
    - Fastqc on reseq'ed samples that were trimmed yesterday
 
   ```bash
   files=data/reads/trim_results/*seq*
   python analysis/pipeline/lib/pipeline.py -a fastqc -i $files -o data/reads/trim_results/qc
  
  ```
* Running

    - **Align on the same samples** 
       
 ```bash
 python analysis/pipeline/lib/pipeline.py -a align -i $files -o data/bt2_align -ref data/ref/huti --no_index
```    
* Running

    - **Align same samples to human ceacam3**
```bash
python analysis/pipeline/lib/pipeline.py -a align -i $files -o data/human/ceacam3_align/ -ref data/ref/human/ceacam3_transcript_variant_1.fasta --no_index
```    
* Running
    
    - **htseq on aligned bam files**
* Had to fix a few minor bugs
* Running. 
* Using Ali-provided gff files, need to replace ID= with gene_id=:    

```bash
sed -i -- 's/ID=/gene_id=/g' *
python analysis/pipeline/lib/pipeline.py -a count -i data/bt2_align/bam/*.bam -o data/counts -ref data/annot/

```


* Realized using Ali's installation of bowtie - mine isn't working. Re-install mine to fix problem. Now need to get the proper path in the module file. 
* Everything is taking way too long to run. Need to find a way to speed up the process. 

     
#### 2018-01-19
##### Summary: All jobs finished, except for htseq-count: ran out of time.

* htseq job ran out of time, but other than that seems to have worked. Also using Ali's installation, need to fix that. Also using Ali's installation of python I think. Not sure why those things are permanently in my PATH. 


* want to extend htseq-count to also give me RPKMs as well as raw counts. 

#### 2018-01-21
##### Summary: In the process of restructuring. Added logging as well as config file so that can easily modify in the future. Have not tested restructured version. 


#### 2018-01-22
##### Summary: Setting up a bunch of other flux jobs (on the previous version of pipeline that's on flux)

* need to extend time for htseq-count

* download multiqc report for re'sequed files
  
  - **sam2bam on human alignments**
  
```bash
sam=data/human/ceacam3_align/*sam
out_dir=data/human/ceacam3_align/bam
python analysis/pipeline/lib/pipeline.py -a sam2bam -i $sam -o $out_dir
```

 -  **sam2bam on reseq'ed aligned samples**
```bash
sam=data/bt2_align/*seq*sam
out_dir=data/bt2_align/bam
python analysis/pipeline/lib/pipeline.py -a sam2bam -i $sam -o $out_dir
```
 

* Think I figured out a couple of mistakes with htseq: Took so long because I ran on both sorted and unsorted bam files. Extra files in the output because it automatically generates one, but I redirected the output? 
* Anyhow, goint to re-run htseq on all the **sorted** files, including re'sequed ones, extended walltime for all jobs to 24 h.
* Concurency issues must be solved at flux level. How to run commands in parallel. 


#### 2018-01-22
##### Summary: Re-organizing modules, adding logging and configpaser modules


#### 2018-01-24
##### Summary: Re-writting samtools module, using samtools 1.5. Re-running sam2bam job.

* Still need to merge local version of the pipeline with the one on flux!!!
```bash
sam=sam=data/bt2_align/*sam
# there were 40 sam samples
python analysis/pipeline/lib/pipeline.py -a sam2bam -i $sam -o data/bt2_align/bam
```
#### 2018-01-25
##### Summary: re-running the counting jobs

* Re-submitting human sam2bam

```bash
sam=data/human/ceacam3_align/*sam
python analysis/pipeline/lib/pipeline.py -a sam2bam -i $sam -o data/human/ceacam3_align/bam
```
 
* Re-submitting htseq job. Cannot find any way to do this in parallel, except for submitting a job for each file seperately. 

```bash
python analysis/pipeline/lib/pipeline.py -a count -i data/bt2_align/bam/sorted/*.bam -o data/counts -ref data/annot/
```

* Nest things that need to be done here:

- commit change to git
- run through everything on flux and make sure it's still working
- add more things to config file: i.e. htseq settings, so that could be easily changed
- learned more advanced flux options, i.e. how to schedule jobs to run after a certain other jobs completed, plus arraying jobs. 
- next step analyzing the count data: look into `rpy2`
- save all the `pbs` scripts to a log folder instead of output folder


#### 2018-01-26
##### Summary: running htseq on CEACAM3 alignments for fun. Htseq broke, had to fix it. 

* This is silly, but running htseq for CEACAM3 alignments, where ref is something I manually threw together just now

```bash
python analysis/pipeline/lib/pipeline.py -a count -i $bam -o data/human/ceacam3_align/counts -ref $ref
```


