summary.file = "Pvalues-Sample_Type.txt"
stat.folder = "."
comp.type = c("FFPE_vs_DNA","FFPE_vs_Frozen","Frozen_vs_DNA")

for(i in 1:length(comp.type)){
	min5.FE.file = paste(stat.folder,"/",comp.type[i],"/differential_status_FE_test.txt",sep="")
	min20.FE.file = paste(stat.folder,"/",comp.type[i],"_freq20/differential_status_FE_test.txt",sep="")
	limma.file= paste(stat.folder,"/",comp.type[i],"/differential_abundance_limma-voom_2var.txt",sep="")
	binomialGLM.file= paste(stat.folder,"/",comp.type[i],"/differential_status_binomial_GLM.txt",sep="")
	linear.lm.1var.file= paste(stat.folder,"/",comp.type[i],"/differential_status_linear-linear_regression_1var.txt",sep="")
	linear.lm.2var.file= paste(stat.folder,"/",comp.type[i],"/differential_status_linear-linear_regression_2var.txt",sep="")
	log.lm.1var.file= paste(stat.folder,"/",comp.type[i],"/differential_status_log2-linear_regression_1var.txt",sep="")
	log.lm.2var.file= paste(stat.folder,"/",comp.type[i],"/differential_status_log2-linear_regression_2var.txt",sep="")
	
	min5.FE.table = read.table(min5.FE.file, head=T, sep="\t")
	min20.FE.table = read.table(min20.FE.file, head=T, sep="\t")
	min20.FE.table = min20.FE.table[match(min5.FE.table$HPV.subtype,min20.FE.table$HPV.subtype),]
	limma.table = read.table(limma.file, head=T, sep="\t")
	limma.table = limma.table[match(min5.FE.table$HPV.subtype,limma.table$HPV.subtype),]
	binomialGLM.table = read.table(binomialGLM.file, head=T, sep="\t")
	binomialGLM.table = binomialGLM.table[match(min5.FE.table$HPV.subtype,binomialGLM.table$HPV.subtype),]
	linear.lm.1var.table = read.table(linear.lm.1var.file, head=T, sep="\t")
	linear.lm.1var.table = linear.lm.1var.table[match(min5.FE.table$HPV.subtype,linear.lm.1var.table$HPV.subtype),]
	linear.lm.2var.table = read.table(linear.lm.2var.file, head=T, sep="\t")
	linear.lm.2var.table = linear.lm.2var.table[match(min5.FE.table$HPV.subtype,linear.lm.2var.table$HPV.subtype),]
	log.lm.1var.table = read.table(log.lm.1var.file, head=T, sep="\t")
	log.lm.1var.table = log.lm.1var.table[match(min5.FE.table$HPV.subtype,log.lm.1var.table$HPV.subtype),]
	log.lm.2var.table = read.table(log.lm.2var.file, head=T, sep="\t")
	log.lm.2var.table = log.lm.2var.table[match(min5.FE.table$HPV.subtype,log.lm.2var.table$HPV.subtype),]
		
	genotypes = as.character(min5.FE.table$HPV.subtype)
	min5.detection.diff = paste(as.character(min5.FE.table$detection.diff),"%",sep="")
	min5.FE.pvalue = min5.FE.table$fe.pvalue
	min5.FE.fdr = min5.FE.table$fe.fdr
		
	min20.detection.diff = paste(as.character(min20.FE.table$detection.diff),"%",sep="")
	min20.FE.pvalue = min20.FE.table$fe.pvalue
	min20.FE.fdr = min20.FE.table$fe.fdr

	limma.voom.pvalue = limma.table$limma.pvalue
	limma.voom.fdr = limma.table$limma.fdr
	
	binomialGLM.pvalue = binomialGLM.table$bGLM.pvalue
	binomialGLM.fdr = binomialGLM.table$bGLM.fdr
	
	linear.lm.1var.pvalue = linear.lm.1var.table$lm.pvalue
	linear.lm.1var.fdr = linear.lm.1var.table$lm.fdr
	
	linear.lm.2var.pvalue = linear.lm.2var.table$lm.pvalue
	linear.lm.2var.fdr = linear.lm.2var.table$lm.fdr

	log.lm.1var.pvalue = log.lm.1var.table$lm.pvalue
	log.lm.1var.fdr = log.lm.1var.table$lm.fdr
	
	log.lm.2var.pvalue = log.lm.2var.table$lm.pvalue
	log.lm.2var.fdr = log.lm.2var.table$lm.fdr
	
	if(i == 1){
		output.table = data.frame(HPV.type=genotypes, comparison = rep(comp.type[i], length(genotypes)),
									min5.detection.diff, min5.FE.pvalue, min5.FE.fdr,
									min20.detection.diff, min20.FE.pvalue, min20.FE.fdr,
									limma.voom.pvalue, limma.voom.fdr,
									binomialGLM.pvalue, binomialGLM.fdr,
									linear.lm.1var.pvalue, linear.lm.1var.fdr,
									linear.lm.2var.pvalue, linear.lm.2var.fdr,
									log.lm.1var.pvalue, log.lm.1var.fdr,
									log.lm.2var.pvalue, log.lm.2var.fdr)
		output.table$comparison = as.character(output.table$comparison)
	}else{
		temp.table = data.frame(HPV.type=genotypes, comparison = rep(comp.type[i], length(genotypes)),
									min5.detection.diff, min5.FE.pvalue, min5.FE.fdr,
									min20.detection.diff, min20.FE.pvalue, min20.FE.fdr,
									limma.voom.pvalue, limma.voom.fdr,
									binomialGLM.pvalue, binomialGLM.fdr,
									linear.lm.1var.pvalue, linear.lm.1var.fdr,
									linear.lm.2var.pvalue, linear.lm.2var.fdr,
									log.lm.1var.pvalue, log.lm.1var.fdr,
									log.lm.2var.pvalue, log.lm.2var.fdr)
		temp.table$comparison = as.character(temp.table$comparison)
		output.table = rbind(output.table, temp.table)
	}#end else
}#endfor(i in 1:length(comp.type))

write.table(output.table, summary.file, quote=F, sep="\t", row.names=F)