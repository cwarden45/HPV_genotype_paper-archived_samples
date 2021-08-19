meta.qPCR20 = read.table("../../Selected_Output_Files/combined_genotype_with_year_and_ethnicity_freq20-FLAGGED.txt",head=T, sep="\t")
meta.ALL5 = read.table("../../Selected_Output_Files/combined_genotype_with_year_and_ethnicity_freq5.txt",head=T, sep="\t")

print(dim(meta.qPCR20))
meta.qPCR20 = meta.qPCR20[meta.qPCR20$HPV.status != "qPCR Flag",]
print(dim(meta.qPCR20))
meta.qPCR20$hist.subtype = factor(meta.qPCR20$hist.subtype,
								levels = c("Adenocarcinoma (Adeno)","Adenosquamous","Squamous Cell Carcinoma (SCC)"))
meta.qPCR20 = meta.qPCR20[!is.na(meta.qPCR20$hist.subtype),]
print(dim(meta.ALL5))

print(dim(meta.ALL5))
meta.ALL5$hist.subtype = factor(meta.ALL5$hist.subtype,
								levels = c("Adenocarcinoma (Adeno)","Adenosquamous","Squamous Cell Carcinoma (SCC)"))
meta.ALL5 = meta.ALL5[!is.na(meta.ALL5$hist.subtype),]
print(dim(meta.ALL5))

meta.ALL5$batch = as.character(meta.ALL5$batch)
meta.ALL5$batch[meta.ALL5$batch == 161007] = "DNA"
meta.ALL5$batch[meta.ALL5$batch == 161206] = "Frozen"
meta.ALL5$batch[meta.ALL5$batch == 170118] = "FFPE"
meta.ALL5$batch = factor(meta.ALL5$batch, levels=c("DNA","Frozen","FFPE"))

png("HistSubtype-all5.png")
plot.table = table(meta.ALL5$hist.subtype, meta.ALL5$batch)
barplot(plot.table, beside=T, col=c("red","magenta", "blue"), main = "All Samples, No qPCR Filter")
legend("bottom", legend = levels(meta.ALL5$hist.subtype), col=c("red","magenta", "blue"),
		pch=15, ncol=2, xpd=T, inset = -0.2, cex=0.7)
dev.off()

meta.qPCR20$batch = as.character(meta.qPCR20$batch)
meta.qPCR20$batch[meta.qPCR20$batch == 161007] = "DNA"
meta.qPCR20$batch[meta.qPCR20$batch == 161206] = "Frozen"
meta.qPCR20$batch[meta.qPCR20$batch == 170118] = "FFPE"
meta.qPCR20$batch = factor(meta.qPCR20$batch, levels=c("DNA","Frozen","FFPE"))

png("HistSubtype-qPCR20.png")
plot.table = table(meta.qPCR20$hist.subtype, meta.qPCR20$batch)
barplot(plot.table, beside=T, col=c("red","magenta", "blue"), main = "Samples Passing qPCR Filter")
legend("bottom", legend = levels(meta.qPCR20$hist.subtype), col=c("red","magenta", "blue"),
		pch=15, ncol=2, xpd=T, inset = -0.2, cex=0.7)
dev.off()