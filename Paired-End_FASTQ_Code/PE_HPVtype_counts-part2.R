
#NOTE: the goal of this file is to re-label the counts, but it requires the output of the "Table1_Supplemental_TableS1_S2-LIMITED-INPUT.R" file in the "Downstream_R_Code" folder

input.file = "PE_HPVtype_counts.txt"
meta.file = "../Downstream_R_Code/Selected_Output_Files/combined_genotype_with_year_and_ethnicity_freq5.txt"
output.file = "PE_HPVtype_counts_final_names.txt"

count.table = read.table(input.file, head=T, sep="\t")
count.type = count.table$HPV.type
count.mat = count.table[,2:ncol(count.table)]
mappable.count = as.character(names(count.mat))
mappable.count = gsub("X","",mappable.count)

meta.table = read.table(meta.file, head=T, sep = "\t")
meta.table = meta.table[order(meta.table$batch, meta.table$Seq.IGC.ID),]
meta.table$Sample = as.character(meta.table$Sample)
meta.table$Sample = gsub("-",".",meta.table$Sample)
ethnicity.matchedIDs = meta.table$SAMPLEID

matchedIDs = mappable.count[match(meta.table$Sample, mappable.count, nomatch=0)]
unmatched.count = mappable.count[-match(matchedIDs,mappable.count)]
print(unmatched.count)
unmatched.meta = meta.table$Sample[-match(matchedIDs,meta.table$Sample)]
print(unmatched.meta)

count.mat = count.mat[,match(meta.table$Sample, mappable.count)]
colnames(count.mat)=ethnicity.matchedIDs

output.table = data.frame(HPV.type=count.type, count.mat)
write.table(output.table, output.file, row.names=F, sep="\t", quote=F)

#I may have previously calculated this somewhere else, but calculate 5% HPV types again (mixed sample was HPV-, but there was ~7% HPV reads)
rownames(count.mat)=count.type

convert.percent = function(arr){
	total = arr[length(arr)]
	arr = arr[-length(arr)]
	arr = arr[-length(arr)]#we don't want to count human samples
	return(arr/total)
}#end def convert.percent

percent.mat = apply(count.mat, 2, convert.percent)

for(i in 1:nrow(percent.mat)){
	HPV = rownames(percent.mat)[i]
	percent.arr = as.numeric(percent.mat[i,])
	print(paste(HPV," : ",length(percent.arr[percent.arr > 0.05])," samples",sep=""))
}#end for(i in 1:nrow(percent.mat))
