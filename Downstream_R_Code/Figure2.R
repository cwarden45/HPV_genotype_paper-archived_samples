library("RColorBrewer")

count.file = "Public_Input_Files/PE_HPVtype_counts_final_names.txt"
meta.file = "Selected_Output_Files/combined_genotype_with_year_and_ethnicity.txt"

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

pdf("to_AI/Figure2b.pdf", width=7, height=2.2, pointsize=1, useDingbats=FALSE)
par(mfcol=c(1,3))
for (i in 1:length(selected.HPV)){
	plot.type = selected.HPV[i]
	print(plot.type)
	subtype.freq = as.numeric(ab.mat[rownames(ab.mat) == plot.type,])
	
	tumor.freq = c()
	normal.freq = c()
	plot.color = c()
	
	frozen.colors = colorRampPalette(c("orange","pink"))(n = length(pair.counts))
	ffpe.colors = colorRampPalette(c("cyan","blue"))(n = length(pair.counts))
	
	legend.names = c()
	legend.color = c()
	for (j in 1:length(pair.counts)){
		pair.samples = meta.table$SAMPLEID[(meta.table$extra.pairID == names(pair.counts)[j])&!is.na(meta.table$extra.pairID)]
		pair.type = meta.table$batch[(meta.table$extra.pairID == names(pair.counts)[j])&!is.na(meta.table$extra.pairID)]
		normal.frozen.sample = as.character(pair.samples[grep(".N$",pair.samples)])
		tumor.frozen.sample = as.character(pair.samples[(pair.samples != normal.frozen.sample)&(pair.type == "Frozen")])
		tumor.ffpe.sample = as.character(pair.samples[(pair.samples != normal.frozen.sample)&(pair.type == "FFPE")])
		
		#frozen
		tumor.freq = c(tumor.freq, subtype.freq[colnames(ab.mat) == tumor.frozen.sample])
		normal.freq = c(normal.freq, subtype.freq[colnames(ab.mat) == normal.frozen.sample])
		plot.color = c(plot.color, frozen.colors[j])
		legend.color = c(legend.color, frozen.colors[j])
		legend.names = c(legend.names,paste(names(pair.counts)[j],"Frozen",sep=":"))
		
		#FFPE
		tumor.freq = c(tumor.freq, subtype.freq[colnames(ab.mat) == tumor.ffpe.sample])
		normal.freq = c(normal.freq, subtype.freq[colnames(ab.mat) == normal.frozen.sample])
		plot.color = c(plot.color, ffpe.colors[j])
		legend.color = c(legend.color, ffpe.colors[j])
		legend.names = c(legend.names,paste(names(pair.counts)[j],"FFPE",sep=":"))
	}#end for (j in 1:length(pair.counts))
	
	legend.names = c(legend.names,"")
	legend.color = c(legend.color,"white")
	other.colors = c("purple","maroon")

	for (j in 1:length(frozen.only.tumor.normal.pairID)){
		pair.samples = meta.table$SAMPLEID[(meta.table$extra.pairID == frozen.only.tumor.normal.pairID[j])&!is.na(meta.table$extra.pairID)]
		pair.type = meta.table$batch[(meta.table$extra.pairID == frozen.only.tumor.normal.pairID[j])&!is.na(meta.table$extra.pairID)]
		normal.frozen.sample = as.character(pair.samples[grep(".N$",pair.samples)])
		tumor.frozen.sample = as.character(pair.samples[(pair.samples != normal.frozen.sample)&(pair.type == "Frozen")])
		tumor.ffpe.sample = as.character(pair.samples[(pair.samples != normal.frozen.sample)&(pair.type == "FFPE")])
		
		#frozen
		tumor.freq = c(tumor.freq, subtype.freq[colnames(ab.mat) == tumor.frozen.sample])
		normal.freq = c(normal.freq, subtype.freq[colnames(ab.mat) == normal.frozen.sample])
		plot.color = c(plot.color, other.colors[j])
		legend.color = c(legend.color, other.colors[j])
		legend.names = c(legend.names,paste(frozen.only.tumor.normal.pairID[j],"Frozen",sep=":"))
		
	}#end for (j in 1:length(pair.counts))
	
	legend.names = c(legend.names,"")
	legend.color = c(legend.color,"white")


	par(mar=c(10,5,5,3))
	plot(normal.freq, tumor.freq, col=plot.color,
		cex.main=2, cex.axis=1.5,cex.lab=1, cex=2,
		ylim=c(0,100), xlim=c(0,100),
		xlab = "", ylab = "",xaxt="n", main=plot.type,
		pch=16)
		
	#assumes 2 tumors, starting with frozen
	for(j in 1:length(pair.counts)){
		start.index = 2*j-1
		stop.index = 2*j
		arrows(normal.freq[start.index],tumor.freq[start.index],
				normal.freq[stop.index],tumor.freq[stop.index],
				lwd=1)
	}#end for(j in 1:(length(tumor.freq)/2))
	#text(10,110,plot.type, cex=3, font=2, xpd=T)
	#legend.order = as.vector(matrix(1:length(legend.names), ncol = 4, byrow = T))
	if(length(legend.names) != length(legend.order)){
		stop("Need to modify manual specification of legend order")
	}
	print(legend.names)
	legend("bottom", legend.names[legend.order], col=legend.color[legend.order],
			xpd=T, inset=-0.5, ncol=4, pch=16, cex=1)
	mtext(20*0:5,1,at=20*0:5, cex=1, padj=1)
	mtext(paste("Normal: Percent ",plot.type," Reads",sep=""),1, cex=1, padj=3.5)
	mtext(paste("Tumor: Percent ",plot.type," Reads",sep=""),2, cex=1, padj=-4)
}#end for (i in 1:nrow(ab.table))
dev.off()

#compare tumor samples from multiple-batches
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

pdf("to_AI/Figure2a.pdf", width=7, height=2.2, pointsize=1, useDingbats=FALSE)
par(mfcol=c(1,3))
for (i in 1:length(selected.HPV)){
	plot.type = selected.HPV[i]
	print(plot.type)
	subtype.freq = as.numeric(ab.mat[rownames(ab.mat) == plot.type,])
	
	tumor.x.freq = subtype.freq[match(tumor.sample.x, colnames(ab.mat))]
	tumor.y.freq = subtype.freq[match(tumor.sample.y, colnames(ab.mat))]

	par(mar=c(10,5,5,3))
	plot(tumor.x.freq, tumor.y.freq, col=pair.color,
		cex.main=2, cex.axis=1.5,cex.lab=1, cex=2,
		ylim=c(0,100), xlim=c(0,100),
		xlab = "", ylab = "",xaxt="n",
		main = paste(plot.type,", r = ",round(cor(tumor.x.freq, tumor.y.freq),digits=2), sep=""),
		pch=16)
	abline(lm(tumor.y.freq~tumor.x.freq),col="gray")
	#text(10,110,plot.type, cex=3, font=2, xpd=T)
	legend("bottom", legend=legend.names, col=legend.colors,
			xpd=T, inset=-0.5, ncol=3, pch=16, cex=1.2)
	mtext(20*0:5,1, at=20*0:5, cex=1, padj=1)
	mtext(paste("Tumor #1: Percent ",plot.type," Reads",sep=""),1, cex=1, padj=3.5)
	mtext(paste("Tumor #2: Percent ",plot.type," Reads",sep=""),2, cex=1, padj=-4)
}#end for (i in 1:nrow(ab.table))
dev.off()