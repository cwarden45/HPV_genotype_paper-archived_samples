meta.qPCR20 = read.table("../../Selected_Output_Files/combined_genotype_with_year_and_ethnicity_freq20-FLAGGED.txt",head=T, sep="\t")
meta.ALL5 = read.table("../../Selected_Output_Files/combined_genotype_with_year_and_ethnicity_freq5.txt",head=T, sep="\t")

print(dim(meta.qPCR20))
meta.qPCR20 = meta.qPCR20[meta.qPCR20$HPV.status != "qPCR Flag",]
print(dim(meta.qPCR20))

print(dim(meta.ALL5))
meta.ALL5 = meta.ALL5[grep("HPV58",meta.ALL5$genotype),]
print(dim(meta.ALL5))
meta.ALL5 = meta.ALL5[meta.ALL5$batch == 170118,]
print(dim(meta.ALL5))
meta.ALL5$hist.subtype = factor(meta.ALL5$hist.subtype,
								levels = c("Adenocarcinoma (Adeno)","Adenosquamous","Squamous Cell Carcinoma (SCC)"))
meta.ALL5 = meta.ALL5[!is.na(meta.ALL5$hist.subtype),]
print(dim(meta.ALL5))

pointCol = rep("red",nrow(meta.ALL5))
pointCol[match(meta.qPCR20$SAMPLEID, meta.ALL5$SAMPLEID, nomatch=0)]="black"

human.reads = as.character(meta.ALL5$human.percent)
human.reads = gsub("\\%","",human.reads)
human.reads=as.numeric(human.reads)

png("FFPE_HPV58_HumanReads_by_HistologicalSubtype.png")
par(mar=c(15,5,3,2))
plot(jitter(as.numeric(meta.ALL5$hist.subtype)), human.reads,
	main= "HPV58+ FFPE Samples",
	ylab="Percent Off-Target Human Reads",xaxt="n",xlab="",
	pch=16, col=pointCol)
mtext(c("Adenocarcinoma","Adenosquamous","Squamous Cell Carcinoma"), side=1, at =1:3, las=2, line=2)
legend("bottom",legend=c("Pass qPCR Filter","Fail qPCR Filter"),
		pch=16, col=c("black","red"), inset = -0.9, xpd=T, ncol=2)
dev.off()

fit = aov(human.reads ~ meta.ALL5$hist.subtype)
result = summary(fit)
aov.pvalue = result[[1]][['Pr(>F)']][1]
print(aov.pvalue)