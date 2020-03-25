input.file = "../../Selected_Output_Files/combined_genotype_with_year_and_ethnicity_freq5.txt"
output.file = "human_reads_vs_HPV-genotype_freq5.pdf"

#input.file = "../../Selected_Output_Files/combined_genotype_with_year_and_ethnicity_freq5-FLAGGED.txt"
#output.file = "human_reads_vs_HPV-genotype_freq5-FLAGGED.pdf"

count.to.ab = function(counts, total){
	return(100*(counts/total))
}#end def count.to.ab

count.table = read.table("../../Public_Input_Files/PE_HPVtype_counts_final_names.txt",head=T, sep="\t")
total.counts = as.numeric(count.table[count.table$HPV.type == "TOTAL",2:ncol(count.table)])
count.mat = count.table[1:(nrow(count.table)-2),2:ncol(count.table)]
rownames(count.mat) = count.table$HPV.type[1:(nrow(count.table)-2)]

#I don't use the "genotype" column, so I can use either file with extended meta data (5% or 15%)
meta.table = read.table(input.file, head=T, sep="\t")
print(dim(meta.table))
meta.table = meta.table[meta.table$HPV.status != "qPCR Flag", ]
print(dim(meta.table))
percent.human = as.character(meta.table$human.percent)
percent.human = gsub("\\%","",percent.human)
percent.human = as.numeric(percent.human)

print(dim(count.mat))
total.counts = total.counts[match(meta.table$SAMPLEID, names(count.mat), nomatch=0)]
count.mat = count.mat[,match(meta.table$SAMPLEID, names(count.mat), nomatch=0)]
print(dim(count.mat))

batch = as.character(meta.table$batch)
batch[batch == "161007"] = "DNA"
batch[batch == "161206"] = "Frozen"
batch[batch == "170118"] = "FFPE"
batch = factor(as.character(batch), levels=c("DNA","Frozen","FFPE"))
batch = batch

color.palette = c("chartreuse4","orange","cyan")
group.levels = levels(batch)
labelColors = rep("black",times=length(batch))
for (i in 1:length(group.levels)){
	labelColors[batch == as.character(group.levels[i])] = color.palette[i]
}#end for (i in 1:length(group.levels))


ab.mat = t(apply(count.mat, 1, count.to.ab, total=total.counts))

selected.HPV = c("HPV16","HPV18","HPV58","HPV45")

#useDingbats=TRUE by default...testing for AI compatibility
pdf(output.file, width=7, height=2.2, pointsize=1, useDingbats=FALSE)
par(mfcol=c(1,4))
for (i in 1:length(selected.HPV)){
	plot.type = selected.HPV[i]
	print(plot.type)
	subtype.freq = as.numeric(ab.mat[rownames(ab.mat) == plot.type,])
	par(mar=c(5,5,5,3))
	plot(percent.human, subtype.freq,
		pch=21, col="black", bg=labelColors,
		xlim=c(0,100), ylim=c(0,100), las=2, lwd=1,
		main = "", cex=1.5, cex.main=1, cex.axis=1,cex.lab=1,
		xlab = "", ylab="")
	mtext("Percentage of Human Reads",1, cex=1, padj=3, font=2)
	mtext(paste("",plot.type," Read Fraction",sep=""),2, cex=1, padj=-4, font=2)
	DNA.freq = subtype.freq[batch == "DNA"]
	DNA.percent.human = percent.human[batch == "DNA"]
	frozen.freq = subtype.freq[batch == "Frozen"]
	frozen.percent.human = percent.human[batch == "Frozen"]
	FFPE.freq = subtype.freq[batch == "FFPE"]
	FFPE.percent.human = percent.human[batch == "FFPE"]
	
	cor.coef = cor(DNA.freq, DNA.percent.human)
	fit=lm(DNA.freq~ DNA.percent.human)
	result = summary(fit)
	pvalue = result$coefficients[2,4]
	print("DNA")
	print(paste("r= ",round(cor.coef,digits=2),", lm p-value = ",signif(pvalue, digits=2),sep=""))
	abline(fit,col=color.palette[1])
	
	cor.coef = cor(frozen.freq,frozen.percent.human)
	fit=lm(frozen.freq~ frozen.percent.human)
	result = summary(fit)
	pvalue = result$coefficients[2,4]
	print("Frozen")
	print(paste("r= ",round(cor.coef,digits=2),", lm p-value = ",signif(pvalue, digits=2),sep=""))	
	abline(lm(fit),col=color.palette[2])
	
	cor.coef = cor(FFPE.freq,FFPE.percent.human)
	fit=lm(FFPE.freq~ FFPE.percent.human)
	result = summary(fit)
	pvalue = result$coefficients[2,4]
	print("FFPE")
	print(paste("r= ",round(cor.coef,digits=2),", lm p-value = ",signif(pvalue, digits=2),sep=""))		
	abline(lm(fit),col=color.palette[3])
	text(10,115, labels=plot.type, xpd=T, cex=2, font=2)
	legend(40,120, legend = group.levels, col ="black", pt.bg=color.palette,
			ncol=3, pch=21,xpd=T, cex=1, inset=-0.4)
}#end for (i in 1:nrow(ab.table))
dev.off()