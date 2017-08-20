import os
import sys
import re

finishedSamples = []

readsFolder = ""
hpvRef = ""
hpvAlignmentFolder = ""
threads = ""

parameterFile = "parameters.txt"

inHandle = open(parameterFile)
lines = inHandle.readlines()
			
for line in lines:
	line = re.sub("\n","",line)
	line = re.sub("\r","",line)
	
	lineInfo = line.split("\t")
	param = lineInfo[0]
	value = lineInfo[1]
	
	if param == "Reads_Folder":
		fullReadFolder = value
	
	if param == "BWA_Ref":
		hpvRef = value

	if param == "Alignment_Folder":
		hpvAlignmentFolder = value

	if param == "Threads":
		threads = value
		
if (fullReadFolder== "") or (fullReadFolder == "[required]"):
	print "Need to enter a value for 'Reads_Folder'!"
	sys.exit()

if (hpvRef== "") or (hpvRef == "[required]"):
	print "Need to enter a value for 'BWA_Ref'!"
	sys.exit()

if (hpvAlignmentFolder== "") or (hpvAlignmentFolder == "[required]"):
	print "Need to enter a value for 'Alignment_Folder'!"
	sys.exit()

if (threads== "") or (threads == "[required]"):
	print "Need to enter a value for 'Threads'!"
	sys.exit()
	
readsFolder = fullReadFolder + "/Cutadapt_Downsample_Reads"

command = "mkdir " + hpvAlignmentFolder
os.system(command)

fileResults = os.listdir(readsFolder)

for file in fileResults:
	result = re.search("(.*)_R1.fastq$",file)
	
	if result:
		sample = result.group(1)
		if sample not in finishedSamples:
			print sample

			read1 = readsFolder + "/" + file
			read2 = re.sub("_R1.fastq$","_R2.fastq",read1)
			
			#HPV alignment
			hpvSubfolder = hpvAlignmentFolder + "/" + sample
			command = "mkdir " + hpvSubfolder
			os.system(command)

			alnSam = hpvSubfolder + "/aligned.sam"
			command = "/opt/bwa/bwa mem -t "+ threads + " " + hpvRef + " " + read1 + " " + read2  + " > " + alnSam
			os.system(command)

			alnBam = hpvSubfolder + "/aligned.bam"
			command = "/opt/samtools-1.3/samtools view -b " + alnSam + " > " + alnBam
			os.system(command)

			command = "rm " + alnSam
			os.system(command)

			sortBam= hpvAlignmentFolder + "/" + sample + ".bam"
			command = "/opt/samtools-1.3/samtools sort " + alnBam + " -o " + sortBam
			os.system(command)

			command = "rm " + alnBam
			os.system(command)

			command = "/opt/samtools-1.3/samtools index " + sortBam
			os.system(command)

			statsFile = hpvSubfolder + "/alignment_stats.txt"
			command = "/opt/samtools-1.3/samtools flagstat " + sortBam + " > " + statsFile
			os.system(command)
			
			statFile = hpvSubfolder + "/idxstats.txt"
			command = "/opt/samtools-1.3/samtools idxstats " + sortBam + " > " + statFile
			os.system(command)
