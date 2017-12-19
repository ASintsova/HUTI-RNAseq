#### 2017-11-12
##### Summary: Reading about model seed and tutorial on FBA

**20802497**
Model building with SEED:

- Annotate with RAST
- Create a model based on annotations (and previously available models)
- Gap fill but selecting the smallest set of reactions needed to allow flux via biomass reactions
- FVA can be used to classify reactions as active/inactive/essential
- FBA to predict essential genes and nutrients essential for growth
- Optimization/manual curation is required to make it more accurate

- Need to go through their manual curation tutorial



**What is flux balance analysis?**

- have to know what is physiologically relevant in order to put proper constraints on the system
- going through the tutorial, trying to replicate results with python tool box


#### 2017-12-13
##### Summary: Reading up on metabolic model reconstruction methods

* Found a whole [**book**](https://link.springer.com/protocol/10.1007%2F978-1-4939-7528-0_11) on this!
  - Refs of interest:

  		- https:/ /doi.org/10.1073/pnas.1307797110
  		- https:/ /doi.org/10.1073/pnas.1523199113
  			- key idea with these two is to use multiple references for your reconstruciton
  		- 17
  		- 18 (key background piece)
  		- 31 for reference model (iJO1366)
  		
  		
  	- Amazing! This seems to provide a step by step process of model building. 
  	- Look at Note 1 script for how to download genomes
  	- Recommend prokka and Roary for pangenome analysis
  	- [Collection of models](https://github.com/opencobra/m_model_collection)
  	- For reconciliation use MNXref namespace from MetaNetX database
  	- Manual curation and merge the referecen GENRE into rGENRE
  	- [Orthology benchmark service]( http://orthology.benchmarkservice.org)
  	- They suggest Inparanoid
  	- Identify orthologs of genes present in rGENRE. 
  	- Get_homologs to identify orthologous genes
  	- Might need Gap-filling methods (GrowMatch and Smiley)
  	- Workflow:
  			- Unfortunately don't show how to reconcile multiple models for rGENRE, use only one
  			- Start with ortho matrix, import with pandas
  			- `delete_model_genes` module
  			- other than that they do the same thing I did! Good for me. 
  			- ok, their code is a little better than mine. 
  			- This is great, can't wait to play with the code
  			
 
#### 2017-12-14
#### Summary: Reading about environment specific models: found relevant CobraPy packages

**18483554**

* GIMME
* Binary: setting a threshold for absence/presence 
* Strains evolved to grow on certain carbon source - setting objective to maximize growth on this specific carbon source. 
* calculates 'consistency score', i.e. how consistent gene expression is with growth on specific substrate. Why would strains evolved on lactate have a higher consistency score compared to WT when grown on succinate or malate? What does this actually show?
* always comparing 2 expression sets and saying which one is more consistent with a certain growth condition
* I think implementation is only in MATLAB

**21081510**

* iMAT
* Web app - looks like it's no longer available? So sad


**[CobraPy](http://opencobra.github.io/cobrapy/packages/) Packages**

* From a brief look found 2 packages that can be used to integrate transcriptomics data into your model
* One is from Palsson group
* Exciting! Means I don't have to write my own :)


