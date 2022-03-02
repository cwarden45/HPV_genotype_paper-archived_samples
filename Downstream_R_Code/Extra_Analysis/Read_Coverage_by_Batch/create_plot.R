meta.file = "../../Selected_Output_Files/combined_genotype_with_year_and_ethnicity_freq20.txt"

meta.table = read.table(meta.file, head=T, sep="\t")
qPCR = meta.table$qPCR.nM

sample.col = rep(NA,nrow(meta.table))
sample.col[meta.table$batch == "161007"]="darkgreen"
sample.col[meta.table$batch == "161206"]="orange"
sample.col[meta.table$batch == "170118"]="cyan"

#plot all
png("Total_Reads_vs_qPCR_Concentration.png")
#color by sample type
plot(meta.table$total.reads, meta.table$qPCR.nM,
	pch=16, col=sample.col,
	xlab="Total Reads", ylab="DNA concentration via qPCR (nM)")
legend("top", legend=c("DNA","Frozen","FFPE"), col = c("darkgreen","orange","cyan"),
	xpd=T, inset = -0.1, ncol=3, pch=16)
dev.off()
