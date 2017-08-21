import re
def getAllGenesLocs (filename):
	
	x = []
	z = ''
	y = ''
	with open(filename) as fh:
		while True:
			z = fh.readline()
			if len (z) == 0:#if reached the end of the file want to break
				break
			if ' gene' in z:
				y = fh.readline().strip()
				x.append((z.strip(), y))
			else:
				pass
	return x


def getYourGenesLocs(geneList, allgenelocs, fn):
	gen_locs = open(fn, 'w+')
	for gene in geneList:
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
				
				
	
	
		
#print (getAllGenesLocs("test.gbk"))

print (getYourGenesLocs(['iucA', 'iutC', 'iucC', 'iucD', 'iutA', 'fitA'],getAllGenesLocs("CFT073.gbk"), "test2.bed"))

		# gene = ''
# 		gene = re.findall('iucA', line)
# 		if len(gene) != 0:
# 			fhand.readline()
# 			locs = fhand.readline()
# 			break
			
			