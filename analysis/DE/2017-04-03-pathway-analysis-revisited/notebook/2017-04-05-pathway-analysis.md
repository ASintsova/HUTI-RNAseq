# Pathway Analysis

This was worked out over a week of 2017-03-29 and 2017-04-05

1. Received 13 files from Evan (one for each genome) that correlated genome specific PROKKA numbers with KO identifier (+ some extra info)
  - These can be found [here](/Users/annasintsova/git_repos/HUTI-RNAseq/data/annotations/kegg)

2. Made a file parser that identified all the prokka # for each of core gene_ids
[here](../lib/get_prokkas_for_gene_ids.py)

3. Made another parser that creates file linking each gene_id to KO +
other info [here](../lib/get_KOs_for_gene_ids.py) - poorly designed, and takes a long time to run

4. Played around a lot with what would be more useful for pathway analysis KOs or this other weird info provided in kegg annotations. After much back and forth stopped on using KO numbers.

5. Wrote a script that takes in DE expression matrix (from Jan-30 analysis), and adds KO numbers for each gene_id: [here](../lib/get_KO_for_DE_genes.py)

6. Now switching to R. Spent a lot of time here soul searching.
