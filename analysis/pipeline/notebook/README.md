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
##### Summary:

* Everything seems to be working. Going to try to sync everything with flux via git. Wish me luck!

