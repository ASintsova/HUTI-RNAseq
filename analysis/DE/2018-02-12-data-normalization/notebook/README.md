#### 2018-02-12
##### Summary: Putting together `saturation_curves.py`

* Wrote function to subsample without replacement
* Calculate coverage on the the subsampled file
* Count how many genes are at or above certain coverage (or # of reads)


#### 2018-02-17
##### Summary: Testing and running `saturation_curves.py`

* Running for HM01 in UR and in UTI, looking at genes that have 10 or more reads at 0.1, 1, 10, 50 and 75%
* Running 5 iterations, seems to be pretty consistent


#### 2018-02-19
##### Summary: want to include 25% as well

* Running subsampling, but it's very slow. eh. 
* Want another point (25 % to see what that looks like)


#### 2018-02-20
##### Summary: Finished plotting the curves. What do they mean?

* Overall, I actually think they look better than expected
* Clear loser in HM66, how did not see that before?
* Have questions about HM60 genome
* Not hitting the saturation limits with about half of the samples. 


Notes from reading about RNAseq data normalization:

* Major assumptions that might not actually hold: symmetry: most genes are not changed, number of upregulated gene is about the same between the two conditions
* If there's a global shift in transcription, i.e. if there's just more RNA in one of the conditions will not notice this
* If anything our hypothesis should be there more mRNA in 
