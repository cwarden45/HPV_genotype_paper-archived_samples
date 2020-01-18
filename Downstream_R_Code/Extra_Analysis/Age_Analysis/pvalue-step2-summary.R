summary.file = "Pvalue_Table_Age.txt"
comp.type = c("Overall","Frozen","FFPE")

for(i in 1:length(comp.type)){
	bGLM.file = paste(comp.type[i],"_differential_status_binomial_GLM.txt",sep="")
	lm.file = paste(comp.type[i],"_differential_status_linear_regression.txt",sep="")
	limma.file= paste(comp.type[i],"_differential_abundance_limma-voom.txt",sep="")
	fe.file= paste(comp.type[i],"_differential_status_FE_test.txt",sep="")
	ks.file= paste(comp.type[i],"_differential_status_KS_test.txt",sep="")
	
	bGLM.table = read.table(bGLM.file, head=T, sep="\t")
	lm.table = read.table(lm.file, head=T, sep="\t")
	lm.table = lm.table[match(bGLM.table$HPV.subtype,lm.table$HPV.subtype),]
	limma.table = read.table(limma.file, head=T, sep="\t")
	limma.table = limma.table[match(bGLM.table$HPV.subtype,limma.table$HPV.subtype),]
	fe.table = read.table(fe.file, head=T, sep="\t")
	fe.table = fe.table[match(bGLM.table$HPV.subtype,fe.table$HPV.subtype),]
	ks.table = read.table(ks.file, head=T, sep="\t")
	ks.table = ks.table[match(bGLM.table$HPV.subtype,ks.table$HPV.subtype),]
	
	genotypes = as.character(bGLM.table$HPV.subtype)
	logistic.regression.pvalue = bGLM.table$bGLM.pvalue
	logistic.regression.fdr = bGLM.table$bGLM.fdr
		
	hpv.type.cor = limma.table$hpv.cor
	linear.regression.pvalue = lm.table$lm.pvalue
	linear.regression.fdr = lm.table$lm.fdr

	limma.voom.pvalue = limma.table$limma.pvalue
	limma.voom.fdr = limma.table$limma.fdr
	
	if(i == 1){
		output.table = data.frame(HPV.type=genotypes, comparison = rep(comp.type[i], length(genotypes)),
									logistic.regression.pvalue, logistic.regression.fdr,
									hpv.type.cor, linear.regression.pvalue, linear.regression.fdr,
									limma.voom.pvalue, limma.voom.fdr,
									fe.percent.detection.gt50=fe.table$gt50, fe.percent.detection.lt50=fe.table$lt50,
									fe.detection.diff=fe.table$detection.diff, fe.pvalue=fe.table$fe.pvalue, fe.fdr=fe.table$fe.fdr,
									ks.avg.pos=ks.table$avg.pos, ks.avg.neg=ks.table$avg.neg,
									ks.detection.diff=ks.table$avg.diff, ks.pvalue=ks.table$ks.pvalue, ks.fdr=ks.table$ks.fdr)
		output.table$comparison = as.character(output.table$comparison)
	}else{
		temp.table = data.frame(HPV.type=genotypes, comparison = rep(comp.type[i], length(genotypes)),
									logistic.regression.pvalue, logistic.regression.fdr,
									hpv.type.cor, linear.regression.pvalue, linear.regression.fdr,
									limma.voom.pvalue, limma.voom.fdr,
									fe.percent.detection.gt50=fe.table$gt50, fe.percent.detection.lt50=fe.table$lt50,
									fe.detection.diff=fe.table$detection.diff, fe.pvalue=fe.table$fe.pvalue, fe.fdr=fe.table$fe.fdr,
									ks.avg.pos=ks.table$avg.pos, ks.avg.neg=ks.table$avg.neg,
									ks.detection.diff=ks.table$avg.diff, ks.pvalue=ks.table$ks.pvalue, ks.fdr=ks.table$ks.fdr)
		temp.table$comparison = as.character(temp.table$comparison)
		output.table = rbind(output.table, temp.table)
	}#end else
}#endfor(i in 1:length(comp.type))

write.table(output.table, summary.file, quote=F, sep="\t", row.names=F)