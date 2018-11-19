from Bio.KEGG import REST



paths = ["path:eco00010", "path:eco00020"]
aa_biosynthesis = "eco00290,eco00300,eco00220,eco00400".split(",")
glycan_biosynthesis = "eco00550 eco00540".split()
lipid_biosynthesis = "eco00061 eco00121 eco01040".split()
def get_genes_for_pathway(pathway):
    genes = []
    pathway_file = REST.kegg_get(pathway).read()
    current_section = None
    for line in pathway_file.rstrip().split("\n"):
        section = line[:12].strip()  # section names are within 12 columns
        if not section == "":
            current_section = section

        if current_section == "GENE":
            print(line)
            gene_identifiers = line[12:].split("; ")[0]
            gene_id, gene_symbol = gene_identifiers.split()

            if not gene_symbol in genes:
                genes.append(gene_id)
    return genes


print(get_genes_for_pathway("path:eco00220"))

# # Filter all human pathways for repair pathways
# repair_pathways = []
# for line in human_pathways.rstrip().split("\n"):
#     entry, description = line.split("\t")
#     if "repair" in description:
#         repair_pathways.append(entry)
#
# # Get the genes for pathways and add them to a list
# repair_genes = []
# for pathway in repair_pathways:
#     pathway_file = REST.kegg_get(pathway).read()  # query and read each pathway
#
#     # iterate through each KEGG pathway file, keeping track of which section
#     # of the file we're in, only read the gene in each pathway
#     current_section = None

#
# print("There are %d repair pathways and %d repair genes. The genes are:" % \
#       (len(repair_pathways), len(repair_genes)))
# print(", ".join(repair_genes))