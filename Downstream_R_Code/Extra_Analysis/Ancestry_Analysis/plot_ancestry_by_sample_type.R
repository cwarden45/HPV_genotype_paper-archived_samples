meta.file = "../../Selected_Output_Files/combined_genotype_with_year_and_ethnicity_freq20.txt"

meta.table = read.table(meta.file, head=T, sep = "\t")
meta.table$batch = as.character(meta.table$batch)
meta.table$batch[meta.table$batch == "161007"] = "DNA"
meta.table$batch[meta.table$batch == "161206"] = "Frozen"
meta.table$batch[meta.table$batch == "170118"] = "FFPE"
meta.table$batch = factor(as.character(meta.table$batch), levels=c("DNA","Frozen","FFPE"))

#HPV+ tumor
print(dim(meta.table))
meta.table = meta.table[meta.table$HPV.status=="pos",]
print(dim(meta.table))
meta.table = meta.table[-grep(".N",meta.table$SAMPLEID),]
print(dim(meta.table))

merged.ethnicity = as.character(meta.table$reported.race)
merged.ethnicity[merged.ethnicity == "White/Caucasian"]="Caucasian-Reported"
merged.ethnicity[merged.ethnicity == "Missing"]=NA
merged.ethnicity[merged.ethnicity == "Missing "]=NA
merged.ethnicity[merged.ethnicity == "Other"]=NA
merged.ethnicity[merged.ethnicity == "Unknown"]=NA
merged.ethnicity[merged.ethnicity == "Black"]="African"
#print(table(merged.ethnicity))
merged.ethnicity[is.na(merged.ethnicity) & (meta.table$ADMIXTURE.mixed=="AFR")]="African"
merged.ethnicity[is.na(merged.ethnicity) & (meta.table$ADMIXTURE.mixed=="EAS")]="Asian"
merged.ethnicity[(merged.ethnicity == "Caucasian-Reported") & (meta.table$ADMIXTURE.mixed=="EUR,AMR")]="Caucasian-AMR/EUR"
merged.ethnicity[(meta.table$ADMIXTURE.mixed=="EUR")]="Caucasian-EUR"
merged.ethnicity[(meta.table$ADMIXTURE.mixed=="AMR")]="Caucasian-AMR"
print(table(merged.ethnicity,meta.table$batch))

plot.ethnicity = merged.ethnicity
plot.ethnicity[grep("Caucasian",plot.ethnicity)]="Caucasian-Reported"
plot.table = table(plot.ethnicity, meta.table$batch)
png("predicted_ancestry_by_sample_type.png")
par(mar=c(5,5,5,12))
barplot(plot.table, beside=T, col=c("orange","green", "plum4"))
legend(13,30,legend=c("African / AFR", "Asian / EAS", "Caucasian-Reported","Caucasian-EUR","Caucasian-AMR","Caucasian-EUR/AMR"),
			col=c("orange", "green", "plum4","blue","red","black"), xpd=T, pch=15)

plot.table = table(merged.ethnicity, meta.table$batch)
plot.table = plot.table[match(c("Caucasian-EUR","Caucasian-AMR","Caucasian-AMR/EUR"),rownames(plot.table)),]

extra.colors = c("blue","red","black")
plot.index = 0
for (i in 1:ncol(plot.table)){
	plot.index = plot.index + 3
	type.count = 0
	for (j in 1:nrow(plot.table)){
		plot.count = plot.table[j,i]
		if (plot.count != 0){
			rect(plot.index, type.count, plot.index+1, type.count+plot.count, col=extra.colors[j])
			type.count = type.count + plot.count
		}
	}#end for (j in 1:ncol(plot.table))
	plot.index = plot.index + 1
}#end for (i in 1:ncol(plot.table))
dev.off()