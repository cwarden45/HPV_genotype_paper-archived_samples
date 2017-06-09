import os
import sys
import re
from Bio.Seq import Seq

fullReadFolder = ""
Fadapter = ""
Radapter = ""

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
	
	if param == "Forward_Primer":
		Fadapter = value

	if param == "Reverse_Primer":
		Radapter = value
		
if (fullReadFolder== "") or (fullReadFolder == "[required]"):
	print "Need to enter a value for 'Reads_Folder'!"
	sys.exit()

if (Fadapter== "") or (Fadapter == "[required]"):
	print "Need to enter a value for 'Forward_Primer'!"
	sys.exit()

if (Radapter== "") or (Radapter == "[required]"):
	print "Need to enter a value for 'Reverse_Primer'!"
	sys.exit()
	
inputFolder = fullReadFolder + "/Downsample_Reads"
outputFolder = fullReadFolder + "/Cutadapt_Downsample_Reads"

command = "mkdir " + outputFolder
os.system(command)

finishedSamples = []

fileResults = os.listdir(inputFolder)


for file in fileResults:
	result = re.search("(.*)_\w{6}_L\d{3}_R1_001.fastq$",file)
	
	if result:
		sample = result.group(1)
		if sample not in finishedSamples:
			print sample

			read1 = inputFolder + "/" + file
			read2 = re.sub("_R1_001.fastq$","_R2_001.fastq",read1)

			trim1 = outputFolder + "/" + sample + "_R1.fastq"
			trim2 = outputFolder + "/" + sample + "_R2.fastq"
	
			seqObj = Seq(Radapter)
			revcomR =  seqObj.reverse_complement()
			seqObj = Seq(Fadapter)
			revcomF =  seqObj.reverse_complement()
	
			command = "cutadapt --max-n 0 -a " + str(revcomR) + " -g " + Fadapter + " -A " + str(revcomF) + " -G " + Radapter + " -m 20 -o " + trim1 + " -p " + trim2 + " " + read1 + " " + read2
			os.system(command)
			