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
##### Summary:

* multiqc report still is appearing in the wrong place (root), instead of outdir. Might want to look into the parameters for multiqc

* dowload the report. Overall looks fine. I'm concerned about combining reads for re-sequence files. Want to see if better off analyzing seperately. If results will be differnt, so have to go back to uploading files. Yay. 
```bash
 scp annasint@flux-xfer.arc-ts.umich.edu:/scratch/hmobley_fluxod/annasint/HUTI-RNAseq/2018-01-11_multiqc_report.html ~/git_repos/HUTI-RNAseq/analysis/pipeline/
```




