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

percent.human = as.character(meta.table$human.percent[match(names(count.mat), meta.table$SAMPLEID)])
percent.human = gsub("\\%","",percent.human)
percent.human = as.numeric(percent.human)

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

pdf("Figure2.pdf", width=30, height=20, useDingbats=FALSE)
par(mfrow=c(2,4))
#Code for Figure 2A
for (i in 1:length(selected.HPV)){
	plot.type = selected.HPV[i]
	print(plot.type)
	subtype.freq = as.numeric(ab.mat[rownames(ab.mat) == plot.type,])
	
	#Overall ANOVA
	fit = aov(subtype.freq ~ meta.table$batch)
	result = summary(fit)
	aov.pvalue = result[[1]][['Pr(>F)']][1]
	print(paste("Overall ANOVA p-value: ",aov.pvalue,sep=""))
	#FFPE vs DNA
	temp.freq  = subtype.freq[meta.table$batch != "Frozen"]
	temp.batch = meta.table$batch[meta.table$batch != "Frozen"]
	fit = aov(temp.freq ~ temp.batch)
	result = summary(fit)
	aov.pvalue = result[[1]][['Pr(>F)']][1]
	print(paste("-->FFPE vs DNA ANOVA p-value: ",aov.pvalue,sep=""))	
	#FFPE vs Frozen
	temp.freq  = subtype.freq[meta.table$batch != "DNA"]
	temp.batch = meta.table$batch[meta.table$batch != "DNA"]
	fit = aov(temp.freq ~ temp.batch)
	result = summary(fit)
	aov.pvalue = result[[1]][['Pr(>F)']][1]
	print(paste("-->FFPE vs Frozen ANOVA p-value: ",aov.pvalue,sep=""))	
	#Frozen vs DNA
	temp.freq  = subtype.freq[meta.table$batch != "FFPE"]
	temp.batch = meta.table$batch[meta.table$batch != "FFPE"]
	fit = aov(temp.freq ~ temp.batch)
	result = summary(fit)
	aov.pvalue = result[[1]][['Pr(>F)']][1]
	print(paste("-->Frozen vs DNA ANOVA p-value: ",aov.pvalue,sep=""))	
	
	par(mar=c(7,10,7,7))
	#plot(meta.table$batch, subtype.freq, cex.main=3, cex.axis=3,cex.lab=3,
	#	ylab="", xaxt="n",
	#	col=color.palette, pch=16, main=plot.type)
	boxplot(subtype.freq~meta.table$batch,cex.main=3, cex.axis=3,cex.lab=3,
		ylab="", xaxt="n", outline=FALSE, ylim=c(0,100),
		col=color.palette, pch=16, main=plot.type)
	mtext(levels(meta.table$batch), side=1, at =1:length(levels(meta.table$batch)), las=1, cex=2, line=3)
	mtext(paste("",plot.type," Read Fraction",sep=""),2, cex=2, padj=-3.5)
	points(jitter(as.numeric(meta.table$batch), 1),subtype.freq, pch=16, cex=2)

	if(i == 1){
		text(-0.1,110, labels="A.", xpd=T, cex=6, font=2)
	}#end if(i == 1)
}#end for (i in 1:nrow(ab.table))

##ANOVA for continuous values (Figure 2B)
fit = aov(percent.human ~ meta.table$batch)
result = summary(fit)
aov.pvalue = result[[1]][['Pr(>F)']][1]
print(paste("Human Read ANOVA p-value: ",aov.pvalue,sep=""))

##Fisher's Exact test for binary values
freq.cat = rep(NA, length(percent.human))
freq.cat[percent.human < 20]="low human"
freq.cat[(20 <= percent.human) & (percent.human <= 80)]="middle"
freq.cat[percent.human > 80]="low HPV"

fisher.mat = table(freq.cat, meta.table$batch)
print(fisher.test(fisher.mat))

for (i in 1:length(selected.HPV)){
	plot.type = selected.HPV[i]
	print(plot.type)
	subtype.freq = as.numeric(ab.mat[rownames(ab.mat) == plot.type,])
	#par(mar=c(5,5,5,3))
	par(mar=c(10,10,10,7))
	plot(percent.human, subtype.freq,
		pch=21, col="black", bg=labelColors,
		xlim=c(0,100), ylim=c(0,100), las=2, lwd=1,
		main = "", cex=3, cex.main=3, cex.axis=3,cex.lab=3,
		xlab = "", ylab="")
	mtext("Percentage of Human Reads",1, cex=2, padj=4)
	mtext(paste("",plot.type," Read Fraction",sep=""),2, cex=2, padj=-3.5)
	DNA.freq = subtype.freq[meta.table$batch == "DNA"]
	DNA.percent.human = percent.human[meta.table$batch == "DNA"]
	frozen.freq = subtype.freq[meta.table$batch == "Frozen"]
	frozen.percent.human = percent.human[meta.table$batch == "Frozen"]
	FFPE.freq = subtype.freq[meta.table$batch == "FFPE"]
	FFPE.percent.human = percent.human[meta.table$batch == "FFPE"]
	
	cor.coef = cor(DNA.freq, DNA.percent.human)
	fit=lm(DNA.freq~ DNA.percent.human)
	result = summary(fit)
	pvalue = result$coefficients[2,4]
	print("DNA")
	print(paste("r= ",round(cor.coef,digits=2),", lm p-value = ",signif(pvalue, digits=2),sep=""))
	abline(fit,col=color.palette[1], lwd=3)
	
	cor.coef = cor(frozen.freq,frozen.percent.human)
	fit=lm(frozen.freq~ frozen.percent.human)
	result = summary(fit)
	pvalue = result$coefficients[2,4]
	print("Frozen")
	print(paste("r= ",round(cor.coef,digits=2),", lm p-value = ",signif(pvalue, digits=2),sep=""))	
	abline(lm(fit),col=color.palette[2], lwd=3)
	
	cor.coef = cor(FFPE.freq,FFPE.percent.human)
	fit=lm(FFPE.freq~ FFPE.percent.human)
	result = summary(fit)
	pvalue = result$coefficients[2,4]
	print("FFPE")
	print(paste("r= ",round(cor.coef,digits=2),", lm p-value = ",signif(pvalue, digits=2),sep=""))		
	abline(lm(fit),col=color.palette[3], lwd=3)
	text(11,112, labels=plot.type, xpd=T, cex=3, font=2)
	legend(35,115, legend = group.levels, col ="black", pt.bg=color.palette,
			ncol=3, pch=21,xpd=T, cex=2, inset=-0.1)
			
	if(i == 1){
		text(-20,119, labels="B.", xpd=T, cex=6, font=2)
	}#end if(i == 1)
}#end for (i in 1:nrow(ab.table))
dev.off()