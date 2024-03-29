#human.cutoff = 1.2
#human.cutoff = 1.0
#human.cutoff = 0.8
human.cutoff = 0.6

library("RColorBrewer")

meta.file = "../../Selected_Output_Files/combined_genotype_with_year_and_ethnicity_freq5.txt"#genotype assignments not directly used, so it doesn't matter which extended annotation file is used
count.file = "../../Public_Input_Files/PE_HPVtype_counts_final_names.txt"

count.to.ab = function(counts, total){
	ab = 100*(counts/total)
	ab[is.na(ab)]=0
	return(ab)
}#end def count.to.ab

count.table = read.table(count.file, head=T, sep="\t")
human.counts = as.numeric(count.table[count.table$HPV.type == "HUMAN_REF",2:ncol(count.table)])
total.counts = as.numeric(count.table[count.table$HPV.type == "TOTAL",2:ncol(count.table)])
count.mat = count.table[1:(nrow(count.table)-2),2:ncol(count.table)]
rownames(count.mat) = count.table$HPV.type[1:(nrow(count.table)-2)]

#using a ratio can create multiple HPV types with values greater than 100%
#so, censor counts less than percent of human reads
censor.human = function(arr, human.ref,cutoff){
	arr[arr < cutoff * human.ref]=0
	return(arr)
}#end def censor.human

print(dim(count.mat))
count.mat = t(apply(count.mat, 1, censor.human, human.ref = human.counts, cutoff=human.cutoff))
print(dim(count.mat))

total.counts = apply(count.mat, 2, sum)
ab.mat = t(apply(count.mat, 1, count.to.ab, total=total.counts))
print(dim(ab.mat))

meta.table = read.table(meta.file, head=T, sep = "\t")
meta.table = meta.table[match(colnames(count.mat), meta.table$SAMPLEID),]
meta.table$batch = as.character(meta.table$batch)
meta.table$batch[meta.table$batch == "161007"] = "DNA"
meta.table$batch[meta.table$batch == "161206"] = "Frozen"
meta.table$batch[meta.table$batch == "170118"] = "FFPE"
meta.table$batch = factor(as.character(meta.table$batch), levels=c("DNA","Frozen","FFPE"))

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

tumor.normal.plot = paste("Type_Human_Ref_FC",human.cutoff,"-Tumor_Normal_Divergence.pdf",sep="")
pdf(tumor.normal.plot, width=30, height=10, useDingbats=FALSE)
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


	par(mar=c(18,10,7,7))
	plot(normal.freq, tumor.freq, col=plot.color,
		cex.main=3, cex.axis=3,cex.lab=3, cex=3,
		ylim=c(0,100), xlim=c(0,100),
		xlab = "", ylab = "",xaxt="n",
		pch=16, main=plot.type)
		
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
			xpd=T, inset=-0.3, ncol=4, pch=16, cex=2)
	mtext(20*0:5,1,at=20*0:5, cex=2, padj=1)
	mtext(paste("Normal: Percent ",plot.type," Reads",sep=""),1, cex=2, padj=3.5)
	mtext(paste("Tumor: Percent ",plot.type," Reads",sep=""),2, cex=2, padj=-4)
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

#legend.names = c("Frozen:Both","FFPE#2:Both","Mixed:Reported","Mixed:QCarray")
#legend.colors = c("orange","cyan","black","brown")
legend.names = c("FFPE#2:Both","Mixed:Reported","Mixed:QCarray")
legend.colors = c("cyan","black","brown")

tumor.tumor.plot = paste("Type_Human_Ref_FC",human.cutoff,"-Tumor_Tumor_Divergence.pdf",sep="")
pdf(tumor.tumor.plot, width=30, height=10, useDingbats=FALSE)
par(mfcol=c(1,3))
for (i in 1:length(selected.HPV)){
	plot.type = selected.HPV[i]
	print(plot.type)
	subtype.freq = as.numeric(ab.mat[rownames(ab.mat) == plot.type,])
	
	tumor.x.freq = subtype.freq[match(tumor.sample.x, colnames(ab.mat))]
	tumor.y.freq = subtype.freq[match(tumor.sample.y, colnames(ab.mat))]

	par(mar=c(15,10,7,7))
	plot(tumor.x.freq, tumor.y.freq, col=pair.color,
		cex.main=3, cex.axis=3,cex.lab=3, cex=3,
		ylim=c(0,100), xlim=c(0,100),
		xlab = "", ylab = "",xaxt="n",
		main = paste(plot.type,", r = ",round(cor(tumor.x.freq, tumor.y.freq),digits=2), sep=""),
		pch=16)
	abline(lm(tumor.y.freq~tumor.x.freq),col="gray")
	#text(10,110,plot.type, cex=3, font=2, xpd=T)
	legend("bottom", legend=legend.names, col=legend.colors,
			xpd=T, inset=-0.25, ncol=3, pch=16, cex=2)
	mtext(20*0:5,1, at=20*0:5, cex=2, padj=1)
	mtext(paste("Tumor #1: Percent ",plot.type," Reads",sep=""),1, cex=2, padj=3.5)
	mtext(paste("Tumor #2: Percent ",plot.type," Reads",sep=""),2, cex=2, padj=-4)
}#end for (i in 1:nrow(ab.table))
dev.off()