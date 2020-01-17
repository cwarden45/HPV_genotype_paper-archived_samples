import sys
import re
import os

alignmentFolder = ""
summaryFile = ""
hpv_cutoff = 5
geno_cutoff = 5
unclear_cutoff = 15
hpv_human_fc = 1.2
geno_human_fc = 1.0

#15% used instead of 5% for some results in paper, but this may not be necessary for some sample types

parameterFile = "parameters.txt"

inHandle = open(parameterFile)
lines = inHandle.readlines()
			
for line in lines:
	line = re.sub("\n","",line)
	line = re.sub("\r","",line)
	
	lineInfo = line.split("\t")
	param = lineInfo[0]
	value = lineInfo[1]
	
	if param == "Alignment_Folder":
		alignmentFolder = value
	
	if param == "Summary_File":
		summaryFile = value

	if param == "HPV_Freq_Cutoff":
		hpv_cutoff = int(value)

	if param == "Type_Freq_Cutoff":
		geno_cutoff = int(value)

	if param == "Unclear_Cutoff":
		unclear_cutoff = int(value)

	if param == "HPV_Human_Ref_FC":
		hpv_human_fc = float(value)
		
	if param == "Type_Human_Ref_FC":
		geno_human_fc = float(value)
		
if (alignmentFolder== "") or (alignmentFolder == "[required]"):
	print "Need to enter a value for 'Alignment_Folder'!"
	sys.exit()
	
if (summaryFile== "") or (summaryFile == "[required]"):
	print "Need to enter a value for 'Summary_File'!"
	sys.exit()

if (hpv_cutoff== "") or (hpv_cutoff == "[required]"):
	print "Need to enter a value for 'HPV_Freq_Cutoff'!"
	sys.exit()

if (geno_cutoff== "") or (geno_cutoff == "[required]"):
	print "Need to enter a value for 'Type_Freq_Cutoff'!"
	sys.exit()
	
if (hpv_human_fc== "") or (hpv_human_fc == "[required]"):
	print "Need to enter a value for 'HPV_Human_Ref_FC'!"
	sys.exit()
	
if (geno_human_fc== "") or (geno_human_fc == "[required]"):
	print "Need to enter a value for 'Type_Human_Ref_FC'!"
	sys.exit()
	
outHandle = open(summaryFile, 'w')
text = "Sample\tHPV.status\tHPV.percent\thuman.percent\tgenotype\tgenotype.percent\n"
outHandle.write(text)

fileResults = os.listdir(alignmentFolder)

for file in fileResults:
	result = re.search("(.*).bam$",file)
	
	if result:
		sample = result.group(1)
		
		alignmentStatFile = os.path.join(alignmentFolder, sample,"idxstats.txt")
		
		if os.path.isfile(alignmentStatFile):
			print sample

			inHandle = open(alignmentStatFile)
			line = inHandle.readline()
		
			genoHash={}
			hpvCount = 0
			humanCount = 0
			readCount = 0
		
			while line:
				line = re.sub("\n","",line)
				line = re.sub("\r","",line)
			
				lineInfo = line.split("\t")
				geno = lineInfo[0]
				alignedCount = int(lineInfo[2])
				unalignedCount = int(lineInfo[3])
				
				hpvHit = re.search("^HPV",geno)
				
				if geno == "*":
					readCount += unalignedCount
				elif hpvHit:
					genoCount = alignedCount - unalignedCount
					if(genoCount < 0):
						genoCount = 0
						
					hpvCount += genoCount
					readCount += genoCount
						
					genoHash[geno] = genoCount
				else:
					humanCount += alignedCount
					readCount += alignedCount
				line = inHandle.readline()
			
			inHandle.close()
			
			percentHPV = 100*float(hpvCount)/float(readCount)
			percentHuman = 100*float(humanCount)/float(readCount)
			
			hpvStatus = "neg"
			if (percentHPV > hpv_cutoff) and (hpvCount > hpv_human_fc * humanCount):
				hpvStatus = "pos"
			
			genoStatus = "NA"
			genoPercent= "NA"
			
			if hpvStatus == "pos":
				genoStatus = "unclear"
				
				revGenoHash = {}
				for geno in genoHash:
					genoCount = genoHash[geno]
					percentGeno = 100* float(genoCount)/float(readCount)
					if (percentGeno > geno_cutoff) and (genoCount > geno_human_fc * humanCount):
						revGenoHash[percentGeno]=geno
			
				#order by frequency
				typeAb = revGenoHash.keys()
				typeAb.sort(reverse=True)
				
				for percentGeno in typeAb:
					geno = revGenoHash[percentGeno]
					if genoStatus == "unclear":
						genoStatus = geno
						genoPercent = '{0:.3f}'.format(percentGeno) + "%"
					else:
						genoStatus = genoStatus + "," + geno
						genoPercent = genoPercent + "," + '{0:.3f}'.format(percentGeno)	+ "%"					
			else:
				if (percentHPV > unclear_cutoff):
					hpvStatus = "unclear"

			text = sample + "\t" + hpvStatus + "\t" + '{0:.3f}'.format(percentHPV) + "%\t"+ '{0:.3f}'.format(percentHuman) + "%\t" + genoStatus + "\t" + genoPercent + "\n"
			outHandle.write(text)
		
