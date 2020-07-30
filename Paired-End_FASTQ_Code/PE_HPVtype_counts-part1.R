#NOTE: This code assumes there is an extra variable in the parameter file called "Count_File" (for the output file of counts)
#NOTE: I set "Count_File" equal to "PE_HPVtype_counts.txt"

#NOTE: Next, use PE_HPVtype_counts-part2.R to change sample names

param.table = read.table("parameters_freq5.txt", header=T, sep="\t")
alignmentFolder=as.character(param.table$Value[param.table$Parameter == "Alignment_Folder"])
summaryFile=as.character(param.table$Value[param.table$Parameter == "Count_File"])

sample.folders = list.dirs(alignmentFolder, full.names = FALSE)
sample.folders = sample.folders[sample.folders!=""]

sampleID = sample.folders
input.files = paste(alignmentFolder, sample.folders,"idxstats.txt", sep="/")

total_counts = c()
human_counts = c()

for (i in 1:length(sampleID)){
	print(sampleID[i])
	input.table = read.table(input.files[i],head=F,sep="\t")
	
	#use slightly more complicated code for total counts to better match percentages from genotype table
	no_alignments = input.table$V4[input.table$V1 == "*"]
	human_counts[i] = sum(input.table$V3[grep("^chr",input.table$V1)])

	input.table=input.table[grep("^HPV",input.table$V1),]
	hpv.type = as.character(input.table$V1)
	aligned.PE.counts = input.table$V3
	unaligned.PE.counts = input.table$V4
	fragment.counts = aligned.PE.counts-unaligned.PE.counts
	fragment.counts[fragment.counts < 0]=0
	
	total_counts[i] = round((sum(fragment.counts) + no_alignments + human_counts[i])/2)
	human_counts[i] = round(human_counts[i] / 2)
	fragment.counts = round(fragment.counts/2)
	
	temp.mat = data.frame(fragment.counts=fragment.counts)
	rownames(temp.mat) = hpv.type
	
	if(i == 1){
		output.table = t(temp.mat)
	}else{
		output.table = rbind(output.table, t(temp.mat[match(colnames(output.table),hpv.type),]))
	}
}#end for (i in 1:length(sampleID))

rownames(output.table)=sampleID
output.table = t(output.table)
output.table=data.frame(HPV.type = rownames(output.table),output.table)
output.table$HPV.type = as.character(output.table$HPV.type)

numericHPV = gsub("b","",output.table$HPV.type)
numericHPV = as.numeric(gsub("HPV","",numericHPV))
output.table = output.table[order(numericHPV),]

output.table = rbind(output.table, c("HUMAN_REF",human_counts))
output.table = rbind(output.table, c("TOTAL",total_counts))
write.table(output.table, summaryFile, quote=F, sep="\t", row.names=F)
