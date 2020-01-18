count.file = "../../Public_Input_Files/PE_HPVtype_counts_final_names.txt"
#NOTE: genotype column not used, so it doesn't really matter which extended meta table is used
meta.file = "../../Selected_Output_Files/combined_genotype_with_year_and_ethnicity_freq15.txt"

count.to.ab = function(counts, total){
	return(100*(counts/total))
}#end def count.to.ab

count.table = read.table(count.file, head=T, sep="\t")
total.counts = as.numeric(count.table[count.table$HPV.type == "TOTAL",2:ncol(count.table)])
count.mat = count.table[1:(nrow(count.table)-2),2:ncol(count.table)]
rownames(count.mat) = count.table$HPV.type[1:(nrow(count.table)-2)]

ab.mat = t(apply(count.mat, 1, count.to.ab, total=total.counts))

meta.table = read.table(meta.file, head=T, sep = "\t")
meta.table = meta.table[match(names(count.mat), meta.table$SAMPLEID),]
meta.table$batch = as.character(meta.table$batch)
meta.table$batch[meta.table$batch == "161007"] = "DNA"
meta.table$batch[meta.table$batch == "161206"] = "Frozen"
meta.table$batch[meta.table$batch == "170118"] = "FFPE"
meta.table$batch = factor(as.character(meta.table$batch), levels=c("DNA","Frozen","FFPE"))

meta.table$Age[meta.table$Age == "Unknown"]=NA
meta.table$Age = as.numeric(as.character(meta.table$Age))

color.palette = c("chartreuse4","orange","cyan")
group.levels = levels(meta.table$batch)
labelColors = rep("black",times=length(meta.table$batch))
for (i in 1:length(group.levels)){
	labelColors[meta.table$batch == as.character(group.levels[i])] = color.palette[i]
}#end for (i in 1:length(group.levels))

selected.HPV = c("HPV16","HPV18","HPV58")

min.year = min(meta.table$Age, na.rm=T)
max.year = max(meta.table$Age, na.rm=T)
pdf("HPV_genotype_by_Age.pdf", width=30, height=10, useDingbats=FALSE)
par(mfcol=c(1,3))
for (i in 1:length(selected.HPV)){
	plot.type = selected.HPV[i]
	print(plot.type)
	subtype.freq = as.numeric(ab.mat[rownames(ab.mat) == plot.type,])
		cor.coef = cor(meta.table$Age, subtype.freq, use="complete.obs")
		fit=lm(subtype.freq~ meta.table$Age)
		result = summary(fit)
		pvalue = result$coefficients[2,4]
		print(paste("r= ",round(cor.coef,digits=2),", lm p-value = ",signif(pvalue, digits=2),sep=""))

	par(mar=c(15,10,7,7))
	plot(meta.table$Age,subtype.freq,cex.main=3, cex.axis=3,cex.lab=3,
		ylab="", ylim=c(0,100), cex = 2, xlab="", bg=labelColors,
		col="black", pch=21, las=2, xlim=c(min.year,max.year))
	frozen.fit = lm(subtype.freq[meta.table$batch == "Frozen"]~ meta.table$Age[meta.table$batch == "Frozen"])
	abline(frozen.fit, col="orange")
	frozen.fit = lm(subtype.freq[meta.table$batch == "FFPE"]~ meta.table$Age[meta.table$batch == "FFPE"])
	abline(frozen.fit, col="cyan")
	
	text(28,112, labels=plot.type, xpd=T, cex=3, font=2)
	mtext("Age",1, cex=2, padj=6)
	mtext(paste("Percent ",plot.type," Reads",sep=""),2, cex=2, padj=-4)
	legend(62,115, legend = c("Frozen","FFPE"), col="black", pt.bg =c("orange","cyan"),
			ncol=2, pch=21,xpd=T, cex=2, inset=-0.1)
}#end for (i in 1:nrow(ab.table))
dev.off()