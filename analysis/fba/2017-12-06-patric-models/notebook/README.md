## Metabolic Modeling with PATRIC

#### 2017-12-06

* Re-annotated HUTI genomes with RAST via PATRIC
* Ran metabolic model reconstruction with [Model SEED](http://modelseed.org/about/faq) via PATRIC ([details](https://docs.patricbrc.org//tutorial/metabolic_model_reconstruction/metabolic_model_reconstruction.html?highlight=metabolic%20model))

* Models and annotations can be found in `/data` folder




#### 2017-12-12
* Based on what learned from the tutorial yesterday went back to play with the draft PATRIC models. Again attempted to stimulate iron depletion as well as aerobic and anaerobic environments. Key take problem: no difference in biomass production at anaerobic or aerobic conditions. In fact there does not appear to be flux through oxygen uptake pathway under any condition tested. This is a problem. Need to go back and see if manual curation can help with that - but little hope, sounds like a ginormous task. Alternatively, back to my own model design. Although cannot say that it will be any better.

* Code saved in jupyter notebook in `/lib/Comparing_draft_PATRIC_models`
