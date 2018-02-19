import pandas as pd
import matplotlib as mpl
import os
import seaborn as sns; sns.set(color_codes=True)
import matplotlib.pyplot as plt
fimbriae_data = "/Users/annasintsova/git_repos/HUTI-RNAseq/analysis/comparative-genomics-huti/data/fimbria.csv"
OCG_data="/Users/annasintsova/git_repos/HUTI-RNAseq/analysis/DE/data/processed_counts/2018-01-29_combined_counts.csv"


ocg = pd.read_csv(OCG_data, index_col = 0)
ocg = ocg.dropna(subset=["CFT073"])
ocg.set_index("CFT073", inplace=True)
to_drop=[c for c in list(ocg.columns) if 'counts' in c or 'MG1655' in c or 'Unnamed' in c or 'seq' in c]
ocg = ocg.drop(to_drop, axis = 1)




fimb = pd.read_csv(fimbriae_data, index_col='genome_locus_tag')
fim_ocg = fimb.merge(ocg, how='inner', left_index=True, right_index=True)

# Presence/Absence
to_drop_for_pa = [c for c in list(fim_ocg.columns) if 'RPKM' in c] + ['gene_symbol', 'function', 'genome', 'gene_symbol']
fim_pa = fim_ocg.drop(to_drop_for_pa, axis = 1)
fim_type = fim_pa.pop("fimbrial_type")
fim_pa.fillna(0, inplace=True)
fim_pa.replace("PARALOGS", 2, inplace=True)
fim_pa.replace("PROKKA.*", 1, inplace=True, regex=True)
num_colors = len(fim_type.unique())
#print(num_colors)
colors = [mpl.cm.Set2(1.*i/9) for i in range(1, num_colors +1)]
lut = dict(zip(fim_type.unique(), colors))
row_colors = fim_type.map(lut)
cmap = mpl.colors.ListedColormap(["#f7f7f7", "#8b9dc3","#3b5998"])

g = sns.clustermap(fim_pa, row_cluster=False, row_colors=row_colors,cmap=cmap)
g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation = 0, fontsize = 8)


legend = [mpl.patches.Patch(color=c, label=l) for l,c in lut.items()]
l2=g.ax_heatmap.legend(handles=legend, bbox_to_anchor=(-0.1,0.5))
#l2=g.ax_heatmap.legend(loc='lower left',handles=legend,frameon=True)
l2.set_title(title='Fimbrial Type',prop={'size':10})

g.savefig("output.png")