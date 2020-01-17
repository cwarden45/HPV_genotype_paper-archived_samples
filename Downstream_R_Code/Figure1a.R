count.file = "Public_Input_Files/PE_HPVtype_counts_final_names.txt"
meta.file = "Selected_Output_Files/combined_genotype_with_year_and_ethnicity_freq5.txt"

#I don't use the "genotype" column, so I can use either file with extended meta data (5% or 15%)

count.to.ab = function(counts, total){
	return(100*(counts/total))
}#end def count.to.ab

count.table = read.table(count.file, head=T, sep="\t")
total.counts = as.numeric(count.table[count.table$HPV.type == "TOTAL",2:ncol(count.table)])
count.mat = count.table[1:(nrow(count.table)-1),2:ncol(count.table)]
rownames(count.mat) = count.table$HPV.type[1:(nrow(count.table)-1)]
names(total.counts) = names(count.mat)
mappable.count = as.character(names(count.mat))


ab.mat = t(apply(count.mat, 1, count.to.ab, total=total.counts))

meta.table = read.table(meta.file, head=T, sep = "\t")
meta.table$batch = as.character(meta.table$batch)
meta.table$batch[meta.table$batch == "161007"] = "DNA"
meta.table$batch[meta.table$batch == "161206"] = "Frozen"
meta.table$batch[meta.table$batch == "170118"] = "FFPE"
meta.table$batch = factor(as.character(meta.table$batch), levels=c("DNA","Frozen","FFPE"))
mappable.meta = as.character(meta.table$SAMPLEID)
mappable.meta = gsub("-",".",mappable.meta)
mappable.meta = gsub("_",".",mappable.meta)

matchedIDs = mappable.count[match(mappable.meta, mappable.count, nomatch=0)]
unmatched.countID = mappable.count[-match(matchedIDs,mappable.count)]
unmatched.meta = mappable.meta[-match(matchedIDs,mappable.meta)]

ab.mat = ab.mat[,match(mappable.meta, mappable.count)]

color.palette = c("chartreuse4","orange","cyan")
group.levels = levels(meta.table$batch)
labelColors = rep("black",times=length(meta.table$batch))
for (i in 1:length(group.levels)){
	labelColors[meta.table$batch == as.character(group.levels[i])] = color.palette[i]
}#end for (i in 1:length(group.levels))

selected.HPV = c("HPV16","HPV18","HPV58","HPV45")

pdf("to_AI/Figure1a.pdf", width=30, height=10, useDingbats=FALSE)
par(mfcol=c(1,4))
for (i in 1:length(selected.HPV)){
	plot.type = selected.HPV[i]
	print(plot.type)
	subtype.freq = as.numeric(ab.mat[rownames(ab.mat) == plot.type,])
	par(mar=c(7,10,7,7))
	#plot(meta.table$batch, subtype.freq, cex.main=3, cex.axis=3,cex.lab=3,
	#	ylab="", xaxt="n",
	#	col=color.palette, pch=16, main=plot.type)
	boxplot(subtype.freq~meta.table$batch,cex.main=3, cex.axis=3,cex.lab=3,
		ylab="", xaxt="n", outline=FALSE, ylim=c(0,100),
		col=color.palette, pch=16, main=plot.type)
	mtext(levels(meta.table$batch), side=1, at =1:length(levels(meta.table$batch)), las=1, cex=2, line=3)
	mtext(paste("Percent ",plot.type," Reads",sep=""),2, cex=2, padj=-4)
	points(jitter(as.numeric(meta.table$batch), 1),subtype.freq, pch=16, cex=2)
}#end for (i in 1:nrow(ab.table))
dev.off()

