import sys
import re
import os
import subprocess

#NOTE: this requires use of a PEAR-Merged alignment (rather than the paired-end alignment, used for most of the analysis)

#NOTE: this is used to create "HPV16_seq_stats_freq5.txt", "HPV18_seq_stats_freq5.txt", and "HPV58_seq_stats_freq5.txt" in "Downstream_R_Code/"

HPVtype = "HPV16"
#HPVtype = "HPV18"
#HPVtype = "HPV58"

#output_suffix = "_freq5.txt"
output_suffix = "_freq20.txt"

minFreq = 15

alignmentFolder = "../../../../../Paper_Draft_PLOS_Pathogens/Pipeline_Code/HPV_typing/PEAR_hg38_plus_35_HPV_Alignment"#these files are large and not on GitHub (put the output of this script is provided)
summaryFile = "../../Public_Input_Files/hg38_plus_35HPV_genotype_calls_freq5.txt"

seqStats = HPVtype + "_seq_stats" + output_suffix
outHandle = open(seqStats,"w")
text = "Sample\tRep.Seq\tRep.Percent\n"
outHandle.write(text)
	
inHandle = open(summaryFile)
line = inHandle.readline()

lineCount = 0

while line:
	line = re.sub("\r","",line)
	line = re.sub("\n","",line)
	
	lineCount += 1
	
	if lineCount > 1:
		lineInfo = line.split("\t")
		HPVresult = re.search(HPVtype,lineInfo[4])
		sampleID = lineInfo[0]
		
		if HPVresult:
			print HPVtype + " detected in " + sampleID
			text = sampleID
			
			bamFile = alignmentFolder + "/" + sampleID + ".bam"
			outputfolder = alignmentFolder + "/" + sampleID
			
			tempSam = outputfolder + "/temp.sam"
			
			outFA = outputfolder + "/"+HPVtype+"_"+sampleID + ".fa"
			
			command = "/opt/samtools-1.3/samtools view " + bamFile + " > " + tempSam
			os.system(command)
			
			inHandle2 = open(tempSam)
			line2 = inHandle2.readline()
			
			faHandle = open(outFA,"w")
			
			readHash = {}
			
			while line2:
				line2 = re.sub("\r","",line2)
				line2 = re.sub("\n","",line2)
				lineInfo2 = line2.split("\t")
				chr = lineInfo2[2]
				if chr == HPVtype:
					read = lineInfo2[0]
					seq = lineInfo2[9]
					
					if seq in readHash:
						readHash[seq] = readHash[seq] + 1
					else:
						readHash[seq] = 1
						
					faHandle.write(">"+read + "\n" + seq + "\n")
				line2 = inHandle2.readline()
			
			inHandle2.close()
			faHandle.close()
			
			revGenoHash = {}
			typeCount = 0
				
			for seq in readHash:
				seqCount = readHash[seq]
				typeCount += seqCount

			revSeqHash = {}
			for seq in readHash:
				seqCount = readHash[seq]
				percentSeq = 100* float(seqCount)/float(typeCount)
				if percentSeq >= minFreq:
					revSeqHash[percentSeq]=seq

			#order by frequency
			seqAb = revSeqHash.keys()
			seqAb.sort(reverse=True)
					
			seqText = ""
			freqText = ""
			
			for i in xrange(0,len(seqAb)):
				freq = seqAb[i]
				seq = revSeqHash[seqAb[i]]
				
				if i == 0:
					seqText=seq
					freqText='{0:.2f}'.format(freq)
				else:
					seqText=seqText + ", " + seq
					freqText=freqText + ", " + '{0:.2f}'.format(freq)					
					
			text = text + "\t" + seqText + "\t" + freqText + "\n"
			outHandle.write(text)
			
			command = "rm " + tempSam
			os.system(command)
	line = inHandle.readline()

inHandle.close()