library("RColorBrewer")

#I don't use the "genotype" column, so I can use either file with extended meta data (5% or 15%)
count.file = "Public_Input_Files/PE_HPVtype_counts_final_names.txt"
meta.file = "Selected_Output_Files/combined_genotype_with_year_and_ethnicity_freq5.txt"
meta.file.FLAGGED = "Selected_Output_Files/combined_genotype_with_year_and_ethnicity_freq5-FLAGGED.txt"

count.to.ab = function(counts, total){
	return(100*(counts/total))
}#end def count.to.ab

count.table = read.table(count.file, head=T, sep="\t")
human.counts = as.numeric(count.table[count.table$HPV.type == "HUMAN_REF",2:ncol(count.table)])
total.counts = as.numeric(count.table[count.table$HPV.type == "TOTAL",2:ncol(count.table)])
count.mat = count.table[1:(nrow(count.table)-2),2:ncol(count.table)]
rownames(count.mat) = count.table$HPV.type[1:(nrow(count.table)-2)]

normalize.to.human = function(arr, human.ref){
	print(arr)
	print(human.ref)
	ratio = (arr / human.ref)
	ratio[ratio < 1] = 0
	ratio[ratio > 100] = 100
	return(ratio)
}#end def normalize.to.human

#ab.mat = t(apply(count.mat, 1, normalize.to.human, human.ref = human.counts))
ab.mat = t(apply(count.mat, 1, count.to.ab, total=total.counts))

meta.table = read.table(meta.file, head=T, sep = "\t")
meta.table = meta.table[match(names(count.mat), meta.table$SAMPLEID),]
meta.table$batch = as.character(meta.table$batch)
meta.table$batch[meta.table$batch == "161007"] = "DNA"
meta.table$batch[meta.table$batch == "161206"] = "Frozen"
meta.table$batch[meta.table$batch == "170118"] = "FFPE"
meta.table$batch = factor(as.character(meta.table$batch), levels=c("DNA","Frozen","FFPE"))

meta.table.FLAGGED = read.table(meta.file.FLAGGED, head=T, sep = "\t")
meta.table.FLAGGED = meta.table.FLAGGED[match(names(count.mat), meta.table.FLAGGED$SAMPLEID),]
qPCR.flag = rep("OK",nrow(meta.table))
qPCR.flag[meta.table.FLAGGED$HPV.status == "qPCR Flag"]="qPCR Flag"
print(table(qPCR.flag))

color.palette = c("chartreuse4","orange","cyan")
group.levels = levels(meta.table$batch)
labelColors = rep("black",times=length(meta.table$batch))
for (i in 1:length(group.levels)){
	labelColors[meta.table$batch == as.character(group.levels[i])] = color.palette[i]
}#end for (i in 1:length(group.levels))

selected.HPV = c("HPV16","HPV18","HPV58")

#compare tumor-normal pairs with multiple sample types
#assumes only cause for >2 paried samples is tumor-normal pair in same batch, so I have to manually remove 3-tumor pairing (partialQCarray.pair1)

#manually specify pair order (think of transposed matrix)
legend.order = c(6,5,
				1,2,
				7,8,
				3,4)

pair.counts = table(meta.table$extra.pairID)
pair.counts = pair.counts[pair.counts > 2]
pair.counts = pair.counts[names(pair.counts) != "partialQCarray.pair1"]

frozen.only.tumor.normal.pairID = meta.table$extra.pairID[grep(".N",meta.table$SAMPLEID)]
frozen.only.tumor.normal.pairID = frozen.only.tumor.normal.pairID[-match(names(pair.counts),frozen.only.tumor.normal.pairID, nomatch=0)]
frozen.only.tumor.normal.pairID = sort(frozen.only.tumor.normal.pairID)

#compare tumor samples from multiple-batches
qPCR.flag = qPCR.flag[-grep(".N",meta.table$SAMPLEID)]
meta.table = meta.table[-grep(".N",meta.table$SAMPLEID),]
print(dim(meta.table))

pair.counts = table(meta.table$extra.pairID)
pair.counts = pair.counts[pair.counts != 1]

tumor.sample.x = c()
tumor.sample.y = c()
pair.color = c()

for (i in 1:length(pair.counts)){
	pair.table = meta.table[!is.na(meta.table$extra.pairID)&(meta.table$extra.pairID == names(pair.counts)[i]),]
	
	if(nrow(pair.table) == 2){
		tumor.sample.x = c(tumor.sample.x, as.character(pair.table$SAMPLEID[1]))
		tumor.sample.y = c(tumor.sample.y, as.character(pair.table$SAMPLEID[2]))
		
		if((pair.table$batch[1] == "Frozen")&(pair.table$batch[2] == "Frozen")){
			pair.color =c(pair.color,"orange")
		}else if((pair.table$batch[1] == "FFPE")&(pair.table$batch[2] == "FFPE")){
			pair.color =c(pair.color,"cyan")
		}else if((pair.table$batch[1] != "DNA")&(pair.table$batch[2] != "DNA")){
			pair.color =c(pair.color,"black")
		}else{
			#should be QCarray pair
			pair.color =c(pair.color,"brown")
		}
	}else if(names(pair.counts)[i] == "partialQCarray.pair1"){
		#manually add points to suspicious looking pair
		#1st add QCarray pairs (extracted DNA and Frozen)
		#2nd add reported pairs (Frozen and FFPE)
		tumor.sample.x = c(tumor.sample.x, "E373.65","S51458.T")
		tumor.sample.y = c(tumor.sample.y, "S51458.T","S16142.01.12")
		pair.color =c(pair.color,"brown","black")
	}else{
		stop("Revise code to handle sample pairs with more than 2 tumor samples")
	}
}#end for (i in 1:length(pair.counts)){

#legend.names = c("Frozen:Both","FFPE:Both","Mixed:Reported","Mixed:QCarray")
#legend.colors = c("orange","cyan","black","brown")
legend.names = c("FFPE:Both","Mixed:Reported","Mixed:QCarray")
legend.colors = c("cyan","black","brown")

pdf("FigureS1.pdf", width=7, height=6)
par(mfrow=c(2,3))
#################################################################################
### plot for all samples, with 20% read fraction genotypes, WITHOUT qPCR flag ###
#################################################################################
for (i in 1:length(selected.HPV)){
	plot.type = selected.HPV[i]
	print(plot.type)
	subtype.freq = as.numeric(ab.mat[rownames(ab.mat) == plot.type,])
	
	tumor.x.freq = subtype.freq[match(tumor.sample.x, colnames(ab.mat))]
	tumor.y.freq = subtype.freq[match(tumor.sample.y, colnames(ab.mat))]

	par(mar=c(8,5,6,3))
	plot(tumor.x.freq, tumor.y.freq, col=pair.color,
		cex.main=1.5, cex.axis=1,cex.lab=1, cex=1.5,
		ylim=c(0,100), xlim=c(0,100), las=2,
		xlab = "", ylab = "",main = "", pch=16)
	abline(lm(tumor.y.freq~tumor.x.freq),col="gray")
	#text(10,110,plot.type, cex=3, font=2, xpd=T)
	legend("bottom", legend=legend.names, col=legend.colors,
			xpd=T, inset=-0.6, ncol=3, pch=16, cex=0.6)
	mtext(paste(plot.type,", r = ",round(cor(tumor.x.freq, tumor.y.freq),digits=2), sep=""),
			side=3, cex=1, font = 2, padj=-0.6)
	mtext(paste("Tumor #1: Percent ",plot.type," Reads",sep=""),1, cex=0.6, padj=4)
	mtext(paste("Tumor #2: Percent ",plot.type," Reads",sep=""),2, cex=0.6, padj=-4.5)

	if(i == 1){
		text(-20,150, labels="A.", xpd=T, cex=3, font=2)
	}#end if(i == 1)
}#end for (i in 1:nrow(ab.table))

##############################################################################
### plot for all samples, with 20% read fraction genotypes, WITH qPCR flag ###
##############################################################################
meta.table = meta.table[qPCR.flag != "qPCR Flag",]
pair.counts = table(as.character(meta.table$extra.pairID))
pair.counts = pair.counts[pair.counts != 1]

tumor.sample.x = c()
tumor.sample.y = c()
pair.color = c()

for (i in 1:length(pair.counts)){
	pair.table = meta.table[!is.na(meta.table$extra.pairID)&(meta.table$extra.pairID == names(pair.counts)[i]),]
	
	if(nrow(pair.table) == 2){
		tumor.sample.x = c(tumor.sample.x, as.character(pair.table$SAMPLEID[1]))
		tumor.sample.y = c(tumor.sample.y, as.character(pair.table$SAMPLEID[2]))
		
		if((pair.table$batch[1] == "Frozen")&(pair.table$batch[2] == "Frozen")){
			pair.color =c(pair.color,"orange")
		}else if((pair.table$batch[1] == "FFPE")&(pair.table$batch[2] == "FFPE")){
			pair.color =c(pair.color,"cyan")
		}else if((pair.table$batch[1] != "DNA")&(pair.table$batch[2] != "DNA")){
			pair.color =c(pair.color,"black")
		}else{
			#should be QCarray pair
			pair.color =c(pair.color,"brown")
		}
	}else if(names(pair.counts)[i] == "partialQCarray.pair1"){
		#manually add points to suspicious looking pair
		#1st add QCarray pairs (extracted DNA and Frozen)
		#2nd add reported pairs (Frozen and FFPE)
		tumor.sample.x = c(tumor.sample.x, "E373.65","S51458.T")
		tumor.sample.y = c(tumor.sample.y, "S51458.T","S16142.01.12")
		pair.color =c(pair.color,"brown","black")
	}else{
		stop("Revise code to handle sample pairs with more than 2 tumor samples")
	}
}#end for (i in 1:length(pair.counts)){

for (i in 1:length(selected.HPV)){
	plot.type = selected.HPV[i]
	print(plot.type)
	subtype.freq = as.numeric(ab.mat[rownames(ab.mat) == plot.type,])
	
	tumor.x.freq = subtype.freq[match(tumor.sample.x, colnames(ab.mat))]
	tumor.y.freq = subtype.freq[match(tumor.sample.y, colnames(ab.mat))]

	par(mar=c(7,5,7,3))
	plot(tumor.x.freq, tumor.y.freq, col=pair.color,
		cex.main=1.3, cex.axis=1,cex.lab=1, cex=1.5,
		ylim=c(0,100), xlim=c(0,100), las=2,
		xlab = "", ylab = "",main = "", pch=16)
	abline(lm(tumor.y.freq~tumor.x.freq),col="gray")
	#text(10,110,plot.type, cex=3, font=2, xpd=T)
	legend("bottom", legend=legend.names, col=legend.colors,
			xpd=T, inset=-0.6, ncol=3, pch=16, cex=0.6)
	mtext(paste(plot.type,", r = ",round(cor(tumor.x.freq, tumor.y.freq),digits=2), sep=""),
			side=3, cex=1, font = 2, padj=-0.6)
	mtext(paste("Tumor #1: Percent ",plot.type," Reads",sep=""),1, cex=0.6, padj=4)
	mtext(paste("Tumor #2: Percent ",plot.type," Reads",sep=""),2, cex=0.6, padj=-4.5)

	if(i == 1){
		text(-20,150, labels="B.", xpd=T, cex=3, font=2)
	}#end if(i == 1)
}#end for (i in 1:nrow(ab.table))

dev.off()