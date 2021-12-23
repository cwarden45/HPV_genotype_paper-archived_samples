meta.file = "Selected_Output_Files/combined_genotype_with_year_and_ethnicity_freq20.txt"


#NOTE: These results were removed and back into the paper during iterations of drafts.  So, this content is also available within https://github.com/cwarden45/HPV_genotype_paper-archived_samples/tree/master/Downstream_R_Code/Extra_Analysis/Ancestry_Analysis/plot_ancestry_by_sample_type.R

meta.table = read.table(meta.file, head=T, sep = "\t")
meta.table$batch = as.character(meta.table$batch)
meta.table$batch[meta.table$batch == "161007"] = "DNA"
meta.table$batch[meta.table$batch == "161206"] = "Frozen"
meta.table$batch[meta.table$batch == "170118"] = "FFPE"
meta.table$batch = factor(as.character(meta.table$batch), levels=c("DNA","Frozen","FFPE"))

##This is presented as an introduction, so revise to **not* filter for HPV+ tumor
##Even though some samples are present in multiple sample types, filter for 1 sample per patient per group.
print(dim(meta.table))
meta.table = meta.table[meta.table$HPV.status == "pos",]
print(dim(meta.table))
meta.table = meta.table[(meta.table$sample.type == "Invasive Cervical Cancer")|(meta.table$sample.type == "Vulvar Cancer")|(meta.table$sample.type == "Endometrial + Cervical Cancer"),]
print(dim(meta.table))
keep.flag = rep(TRUE,nrow(meta.table))
for (batch_value in levels(meta.table$batch)){
	batch.table = meta.table[meta.table$batch == batch_value,]
	batch.table = batch.table[!is.na(batch.table$extra.pairID),]
	
	batch.pair.counts = table(batch.table$extra.pairID)
	rep.pairs = names(batch.pair.counts[batch.pair.counts > 1])
	for (rep.pair in rep.pairs){
		print(rep.pair)
		#arbitrarily remove 2nd sample from plot
		pair.table = batch.table[batch.table$extra.pairID == rep.pair,]
		#print(pair.table$SAMPLEID)
		
		dup.sample = as.character(pair.table$SAMPLEID[2])
		print(dup.sample)
		keep.flag[meta.table$SAMPLEID == dup.sample]=FALSE
	}#end for (rep.pair in rep.pairs)
	
}#end for (batch_value in meta.table$batch)

print(table(keep.flag))
meta.table = meta.table[keep.flag,]
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
pdf("to_AI/Figure1c.pdf", useDingbats=FALSE)
par(mar=c(5,5,12,5))
barplot(plot.table, beside=T, col=c("orange","green", "plum4"), main="",
			cex.axis=2, cex.names=1.5)
mtext("HPV+ Tumor",3, cex=1.5, padj=-2, font=2)
legend("top",legend=c("African / AFR", "Asian / EAS", "Caucasian-Reported","Caucasian-EUR","Caucasian-AMR","Caucasian-EUR/AMR"),
			col=c("orange", "green", "plum4","blue","red","black"), xpd=T, pch=15,
			cex=1.5, ncol=2, inset=-0.55)

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