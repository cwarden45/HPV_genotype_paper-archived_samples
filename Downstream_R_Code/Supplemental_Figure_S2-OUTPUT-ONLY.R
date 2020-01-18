files = list.files("Y:\\cwarden_templates\\PhiX\\HPV_Cross-Contamination")
files = files[grep("^\\d+_",files)]#look for exact matches (not "fuzzy" matches, or other files)
files = files[-grep("Ogembo",files)]#remove samples that were the source of the HPV sequences, if re-analyzed
files = files[-grep("Unassigned",files)]#there was a considerable difference for "Undetermined" reads, but I wanted to be more fair and focus on the user samples
previous.exact.match = files[grep("-previous_",files)]

#earlier base calling

sampleID = c()
runID = c()
total.assembled = c()
absolute.HPV16 = c()
absolute.HPV18 = c()
absolute.HPV58 = c()
percent.HPV16 = c()
percent.HPV18 = c()
percent.HPV58 = c()

for (i in 1:length(previous.exact.match)){
	temp.runID = previous.exact.match[i]
	temp.runID = gsub("-previous_.*","",temp.runID)
	print(temp.runID)
	input.file = read.table(paste("Y:\\cwarden_templates\\PhiX\\HPV_Cross-Contamination\\",previous.exact.match[i],sep=""), head=T, sep="\t")
	temp.percent.HPV16 = 100 * input.file$HPV16 / input.file$Total.Assembled
	temp.percent.HPV18 = 100 * input.file$HPV18 / input.file$Total.Assembled
	temp.percent.HPV58 = 100 * input.file$HPV58 / input.file$Total.Assembled
	
	sampleID = c(sampleID, as.character(input.file$Sample))
	runID = c(runID, rep(temp.runID, nrow(input.file)))
	
	total.assembled = c(total.assembled, input.file$Total.Assembled)
	absolute.HPV16 = c(absolute.HPV16, input.file$HPV16)
	absolute.HPV18 = c(absolute.HPV18, input.file$HPV18)
	absolute.HPV58 = c(absolute.HPV58, input.file$HPV58)
	
	percent.HPV16 = c(percent.HPV16, temp.percent.HPV16)
	percent.HPV18 = c(percent.HPV18, temp.percent.HPV18)
	percent.HPV58 = c(percent.HPV58, temp.percent.HPV58)
}#end for (i in 1:length(previous.exact.match))

previous.sample.stats = data.frame(sampleID, runID,
								total.assembled, absolute.HPV16, absolute.HPV18, absolute.HPV58,
								percent.HPV16, percent.HPV18, percent.HPV58)

plot.mat = data.frame(Type=c(rep("HPV16",nrow(previous.sample.stats)),rep("HPV18",nrow(previous.sample.stats)),rep("HPV58",nrow(previous.sample.stats))),
						percent.HPV.leak = c(previous.sample.stats$percent.HPV16, previous.sample.stats$percent.HPV18, previous.sample.stats$percent.HPV58),
						absolute.HPV.leak = c(previous.sample.stats$absolute.HPV16, previous.sample.stats$absolute.HPV18, previous.sample.stats$absolute.HPV58))

pdf("Selected_Output_Files/Additional_File_07_Supplemental_Figure_S2.pdf")
par(mfcol=c(1,3))
plot(plot.mat$Type, plot.mat$percent.HPV.leak, pch=16,
	xlab="HPV L1 Cross-Contamination Genotype", ylab="Percent Cross-Contamination (%)")
plot(plot.mat$Type, plot.mat$absolute.HPV.leak, pch=16,
	xlab="HPV L1 Cross-Contamination Genotype", ylab="Absolute Cross-Contamination")
plot(plot.mat$Type, log10(plot.mat$absolute.HPV.leak+1), pch=16, col="gray",
	xlab="HPV L1 Cross-Contamination Genotype", ylab="log10(Absolute Cross-Contamination+1)")
dev.off()
