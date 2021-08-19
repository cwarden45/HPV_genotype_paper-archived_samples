summary.file = "Pvalue_Table_Histological_Subtype.txt"
comp.type_2var = c("Overall-SCC_vs_AdenoTypes-2var-qPCR20","Overall-SCC_vs_AdenoTypes-2var-all5")
comp.type_1var = c("DNA-SCC_vs_AdenoTypes-1var-qPCR20","Frozen-SCC_vs_AdenoTypes-1var-qPCR20","FFPE-SCC_vs_AdenoTypes-1var-qPCR20",
					"DNA-SCC_vs_AdenoTypes-1var-all5","Frozen-SCC_vs_AdenoTypes-1var-all5","FFPE-SCC_vs_AdenoTypes-1var-all5")

for(i in 1:length(comp.type_1var)){
	limma.file= paste(comp.type_1var[i],"_differential_abundance_limma-voom.txt",sep="")
	fe.file= paste(comp.type_1var[i],"_differential_status_FE_test.txt",sep="")
	binomialGLM.file= paste(comp.type_1var[i],"_differential_abundance_binomialGLM.txt",sep="")
	
	limma.table = read.table(limma.file, head=T, sep="\t")
	fe.table = read.table(fe.file, head=T, sep="\t")
	fe.table = fe.table[match(limma.table$HPV.subtype,fe.table$HPV.subtype),]
	binomialGLM.table = read.table(binomialGLM.file, head=T, sep="\t")
	binomialGLM.table = binomialGLM.table[match(limma.table$HPV.subtype,binomialGLM.table$HPV.subtype),]
	
	genotypes = as.character(limma.table$HPV.subtype)

	FE.detection.diff = paste(as.character(fe.table$detection.diff),"%",sep="")
	FE.pvalue = fe.table$fe.pvalue
	FE.fdr = fe.table$fe.fdr

	limma.voom.pvalue = limma.table$limma.pvalue
	limma.voom.fdr = limma.table$limma.fdr

	binomialGLM.pvalue = binomialGLM.table$bGLM.pvalue
	binomialGLM.fdr = binomialGLM.table$bGLM.fdr
	
	if(i == 1){
		output.table = data.frame(HPV.type=genotypes, comparison = rep(comp.type_1var[i], length(genotypes)),
									limma.voom.pvalue, limma.voom.fdr,
									FE.detection.diff, FE.pvalue, FE.fdr,
									binomialGLM.pvalue, binomialGLM.fdr)
		output.table$comparison = as.character(output.table$comparison)
	}else{
		temp.table = data.frame(HPV.type=genotypes, comparison = rep(comp.type_1var[i], length(genotypes)),
									limma.voom.pvalue, limma.voom.fdr,
									FE.detection.diff, FE.pvalue, FE.fdr,
									binomialGLM.pvalue, binomialGLM.fdr)
		temp.table$comparison = as.character(temp.table$comparison)
		output.table = rbind(output.table, temp.table)
	}#end else
}#endfor(i in 1:length(comp.type_1var))

for(i in 1:length(comp.type_2var)){
	limma.file= paste(comp.type_2var[i],"_differential_abundance_limma-voom.txt",sep="")
	binomialGLM.file= paste(comp.type_2var[i],"_differential_abundance_binomialGLM.txt",sep="")
	
	limma.table = read.table(limma.file, head=T, sep="\t")
	binomialGLM.table = read.table(binomialGLM.file, head=T, sep="\t")
	binomialGLM.table = binomialGLM.table[match(limma.table$HPV.subtype,binomialGLM.table$HPV.subtype),]
	
	genotypes = as.character(limma.table$HPV.subtype)

	FE.detection.diff = rep("skipped",nrow(limma.table))
	FE.pvalue = rep("skipped",nrow(limma.table))
	FE.fdr = rep("skipped",nrow(limma.table))

	limma.voom.pvalue = limma.table$limma.pvalue
	limma.voom.fdr = limma.table$limma.fdr

	binomialGLM.pvalue = binomialGLM.table$bGLM.pvalue
	binomialGLM.fdr = binomialGLM.table$bGLM.fdr
	
		temp.table = data.frame(HPV.type=genotypes, comparison = rep(comp.type_2var[i], length(genotypes)),
									limma.voom.pvalue, limma.voom.fdr,
									FE.detection.diff, FE.pvalue, FE.fdr,
									binomialGLM.pvalue, binomialGLM.fdr)
		temp.table$comparison = as.character(temp.table$comparison)
		output.table = rbind(output.table, temp.table)
}#endfor(i in 1:length(comp.type_2var))

write.table(output.table, summary.file, quote=F, sep="\t", row.names=F)