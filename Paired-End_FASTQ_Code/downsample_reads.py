import os
import re
import sys

fullReadFolder = ""
downsampleCount = 10000

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
	
	if param == "Max_Reads":
		downsampleCount = int(value)

if (fullReadFolder== "") or (fullReadFolder == "[required]"):
	print "Need to enter a value for 'Reads_Folder'!"
	sys.exit()

if (downsampleCount== "") or (downsampleCount == "[required]"):
	print "Need to enter a value for 'Max_Reads'!"
	sys.exit()
	
downsampleReadFolder = fullReadFolder + "/Downsample_Reads"

command = "mkdir " + downsampleReadFolder
os.system(command)

fileResults = os.listdir(fullReadFolder)

lineLimit = 4 * downsampleCount

for file in fileResults:
	result = re.search("(.*).fastq$",file)
	
	if result:
		print file
		inFQ = fullReadFolder + "/" + file
		outFQ = downsampleReadFolder + "/" + file
		
		inHandle = open(inFQ)
		line = inHandle.readline()

		outHandle = open(outFQ, "w")

		lineCount = 0
		
		while line:	
			lineCount += 1
			
			if lineCount <= lineLimit:
				text = line
				outHandle.write(text)
			else:	
				break
		

			line = inHandle.readline()