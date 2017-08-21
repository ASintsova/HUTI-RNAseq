
eutRbed = open("eutR.test.bed", 'w+')
strand = ''
for line in open("eutR.test"):
	line = line.rstrip()
	word = line.split("\t")
	if int(word[6]) < int(word[5]):
		strand = '-'
	else:
		strand = '+'
	eutRbed.write(word[1]+"\t"+ word[5]+"\t"+word[6]+"\t"+ strand +"\n")
    
    
