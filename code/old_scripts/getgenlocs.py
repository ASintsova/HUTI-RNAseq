import re
import sys, ast

"""
input: list of gene names, and gbk file

output: bed file with locations in the genome
"""


def getAllGenesFromGbk (filename):
	"""reads through gbk file, makes a list of tuples of
	all gene names present, and their locations in the genome"""
	gene_list = []
	z = ''
	y = ''
	with open(filename) as fh:
		while True:
			z = fh.readline()
			if len (z) == 0:#if reached the end of the file want to break
				break
			if ' gene' in z:
				y = fh.readline().strip()
				gene_list.append((z.strip(), y))
			else:
				pass
	return gene_list


def getYourGenesLocs(favgeneList, allgenelocs, fn):
	gen_locs = open(fn, 'w+')
	for gene in favgeneList:
		for ag in allgenelocs:
			if gene in ag[1]:
				if 'complement' in ag[0]:
					loc = ag[0].split('(')					
					print "gi|26111730|gb|AE014075.1|" +  "\t" + loc[1].split('..')[0] + '\t' + loc[1].split('..')[1][:-1]+ '\t' + str(gene)
					gen_locs.write("gi|26111730|gb|AE014075.1|"+ "\t" + loc[1].split('..')[0] + '\t' + loc[1].split('..')[1][:-1]+ '\t' + str(gene)+"\n")
				else:	
					loc = ag[0].split('            ')
					print "gi|26111730|gb|AE014075.1|" + "\t" + loc[1].split('..')[0] + '\t' + loc[1].split('..')[1]+ '\t' + str(gene)
					gen_locs.write("gi|26111730|gb|AE014075.1|" + "\t" + loc[1].split('..')[0] + '\t' + loc[1].split('..')[1]+ '\t' + str(gene)+"\n")
				
				
	
	
inputList = ast.literal_eval( sys.argv[1] )
gbk_file = str(sys.argv[2])
out_file = str(sys.argv[3])

print (getYourGenesLocs(inputList, getAllGenesFromGbk(gbk_file), out_file))

		# gene = ''
# 		gene = re.findall('iucA', line)
# 		if len(gene) != 0:
# 			fhand.readline()
# 			locs = fhand.readline()
# 			break
			
			