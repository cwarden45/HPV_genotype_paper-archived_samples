count.file = "../../Public_Input_Files/PE_HPVtype_counts_final_names.txt"
meta.file = "../../Selected_Output_Files/combined_genotype_with_year_and_ethnicity_freq15.txt"

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

median.insert.size = meta.table$L1.median.human.insert.size
max.insert.size = meta.table$L1.max.human.insert.size

color.palette = color.palette = c("green","orange","cyan")
group.levels = levels(meta.table$batch)
labelColors = rep("black",times=length(meta.table$batch))
for (i in 1:length(group.levels)){
	labelColors[meta.table$batch == as.character(group.levels[i])] = color.palette[i]
}#end for (i in 1:length(group.levels))

		fit = aov(median.insert.size ~ meta.table$batch)
		result = summary(fit)
		aov.pvalue = result[[1]][['Pr(>F)']][1]

print(paste("ANOVA Insert Size co-infection p-value: ",aov.pvalue,sep=""))

#, useDingbats=FALSE
png("Insert_Size_Figure.png", width=600, height=200)
par(mfcol=c(1,3))

	boxplot(median.insert.size~meta.table$batch,cex.main=0.8, cex.axis=0.8,cex.lab=0.8, cex=0.6,
		ylab="Median Insert Size", outline=FALSE, ylim=c(0,250),
		col=color.palette, pch=16, main="Median Insert Size by Sample Type", na.rm=T)
	points(jitter(as.numeric(meta.table$batch), 1),median.insert.size, pch=16, cex=0.6)

	#boxplot(max.insert.size~meta.table$batch,cex.main=1, cex.axis=1,cex.lab=1,
	#	ylab="", xaxt="n", outline=FALSE, ylim=c(0,1000),
	#	col=color.palette, pch=16, main="Max Insert Size by Batch", na.rm=T)
	#mtext("Max Insert Size",2, cex=1, padj=-3)
	#mtext(levels(meta.table$batch), side=1, at =1:length(levels(meta.table$batch)), las=1, cex=1, line=1)
	#points(jitter(as.numeric(meta.table$batch), 1),max.insert.size, pch=16, cex=1)

#percent human
percent.human = as.character(meta.table$human.percent)
percent.human = gsub("%","",percent.human)
percent.human = as.numeric(percent.human)

plot(median.insert.size, percent.human, cex=0.6,
		xlab="Median Insert Size",ylab="Percent Human",
		col=labelColors, pch=16)
#inset=-0.35
legend("top",group.levels, pch=16, col=color.palette,
			ncol=3, inset=-0.20, xpd=T, cex=0.8)

#percent HPV58 (if genotyped)
percent.HPV58 = rep(0, nrow(meta.table))
for(i in 1:nrow(meta.table)){
	temp.geno = as.character(meta.table$genotype[i])
	temp.percent = as.character(meta.table$genotype.percent[i])
	temp.percent = gsub("%","",temp.percent)
	
	temp.geno = temp.geno[grep("HPV58",temp.geno)]
	if(length(temp.geno)==1){
		geno.arr = unlist(strsplit(temp.geno,split=","))
		percent.arr = unlist(strsplit(temp.percent,split=","))
		
		percent.HPV58[i]=as.numeric(percent.arr[geno.arr == "HPV58"])
	}#end if(length(temp.geno)==1)
}#end for(i in 1:nrow(meta.table))

plot(median.insert.size, percent.HPV58, cex=0.6,
		xlab="Median Insert Size",ylab="Percent HPV58",
		col=labelColors, pch=16)
#inset=-0.35
legend("top",group.levels, pch=16, col=color.palette,
			ncol=3, inset=-0.20, xpd=T, cex=0.8)
dev.off()