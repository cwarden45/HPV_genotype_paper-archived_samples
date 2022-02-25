import os
import sys
import re
from Bio.Seq import Seq

sampleFile = ""
fullReadFolder = ""
Fadapter = ""
Radapter = ""

finishedSamples = []
parameterFile = "parameters.txt"

inHandle = open(parameterFile)
lines = inHandle.readlines()
			
for line in lines:
	line = re.sub("\n","",line)
	line = re.sub("\r","",line)
	
	lineInfo = line.split("\t")
	param = lineInfo[0]
	value = lineInfo[1]

	if param == "Sample_Description_File":
		sampleFile = value
	
	if param == "Reads_Folder":
		fullReadFolder = value
	
	if param == "Forward_Primer":
		Fadapter = value

	if param == "Reverse_Primer":
		Radapter = value

if (sampleFile== "") or (sampleFile == "[required]"):
	print "Need to enter a value for 'Sample_Description_File'!"
	sys.exit()
		
if (fullReadFolder== "") or (fullReadFolder == "[required]"):
	print "Need to enter a value for 'Reads_Folder'!"
	sys.exit()

if (Fadapter== "") or (Fadapter == "[required]"):
	print "Need to enter a value for 'Forward_Primer'!"
	sys.exit()

if (Radapter== "") or (Radapter == "[required]"):
	print "Need to enter a value for 'Reverse_Primer'!"
	sys.exit()
	
inputFolder = fullReadFolder
outputFolder = fullReadFolder + "/Cutadapt_Full_Reads"

command = "mkdir " + outputFolder
os.system(command)

lineCount = 0
sampleIndex = -1
r1Index = -1
r2Index = -1

inHandle = open(sampleFile)
lines = inHandle.readlines()
			
for line in lines:
	line = re.sub("\n","",line)
	line = re.sub("\r","",line)
	
	lineInfo = line.split("\t")
	
	lineCount += 1
	
	if lineCount == 1:
		for i in range(0,len(lineInfo)):
			if lineInfo[i] == "SampleID":
				sampleIndex = i
			elif lineInfo[i] == "Forward_Read":
				r1Index = i
			elif lineInfo[i] == "Reverse_Read":
				r2Index = i

		if sampleIndex == -1:
			print "Need to have column 'SampleID' in " + sampleFile + " to map sampleID"
			sys.exit()
		if r1Index == -1:
			print "Need to have column 'Forward_Read' in " + sampleFile + " to map R1 read"
			sys.exit()
		if r2Index == -1:
			print "Need to have column 'Reverse_Read' in " + sampleFile + " to map R2 read"
			sys.exit()			
	else:
		sample = lineInfo[sampleIndex]
		if sample not in finishedSamples:
			print sample

			read1 = inputFolder + "/" + lineInfo[r1Index]
			read2 = inputFolder + "/" + lineInfo[r2Index]

			trim1 = outputFolder + "/" + sample + "_R1.fastq"
			trim2 = outputFolder + "/" + sample + "_R2.fastq"
	
			seqObj = Seq(Radapter)
			revcomR =  seqObj.reverse_complement()
			seqObj = Seq(Fadapter)
			revcomF =  seqObj.reverse_complement()
	
			command = "cutadapt --max-n 0 -a " + str(revcomR) + " -g " + Fadapter + " -A " + str(revcomF) + " -G " + Radapter + " -m 20 -o " + trim1 + " -p " + trim2 + " " + read1 + " " + read2
			os.system(command)
			
