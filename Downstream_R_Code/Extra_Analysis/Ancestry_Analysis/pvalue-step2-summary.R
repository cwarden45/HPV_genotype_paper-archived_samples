summary.file = "Pvalue_Table_Ancestry.txt"
comp.type = c("AMR_vs_EUR")

for(i in 1:length(comp.type)){
	limma.1var.file= paste(comp.type[i],"_1var","_differential_abundance_limma-voom.txt",sep="")
	fe.file= paste(comp.type[i],"_1var","_differential_status_FE_test.txt",sep="")
	limma.2var.file= paste(comp.type[i],"_2var","_differential_abundance_limma-voom.txt",sep="")
	
	limma.1var.table = read.table(limma.1var.file, head=T, sep="\t")
	fe.table = read.table(fe.file, head=T, sep="\t")
	fe.table = fe.table[match(limma.1var.table$HPV.subtype,fe.table$HPV.subtype),]
	limma.2var.table = read.table(limma.2var.file, head=T, sep="\t")
	limma.2var.table = limma.2var.table[match(limma.1var.table$HPV.subtype,limma.2var.table$HPV.subtype),]
	
	genotypes = as.character(limma.1var.table$HPV.subtype)

	FE.detection.diff = paste(as.character(fe.table$detection.diff),"%",sep="")
	FE.pvalue = fe.table$fe.pvalue
	FE.fdr = fe.table$fe.fdr

	limma.voom.1var.pvalue = limma.1var.table$limma.pvalue
	limma.voom.1var.fdr = limma.1var.table$limma.fdr

	limma.voom.2var.pvalue = limma.2var.table$limma.pvalue
	limma.voom.2var.fdr = limma.2var.table$limma.fdr
	
	if(i == 1){
		output.table = data.frame(HPV.type=genotypes, comparison = rep(comp.type[i], length(genotypes)),
									limma.voom.1var.pvalue, limma.voom.1var.fdr,
									FE.detection.diff, FE.pvalue, FE.fdr,
									limma.voom.2var.pvalue, limma.voom.2var.fdr)
		output.table$comparison = as.character(output.table$comparison)
	}else{
		temp.table = data.frame(HPV.type=genotypes, comparison = rep(comp.type[i], length(genotypes)),
									limma.voom.1var.pvalue, limma.voom.1var.fdr,
									FE.detection.diff, FE.pvalue, FE.fdr,
									limma.voom.2var.pvalue, limma.voom.2var.fdr)
		temp.table$comparison = as.character(temp.table$comparison)
		output.table = rbind(output.table, temp.table)
	}#end else
}#endfor(i in 1:length(comp.type))

write.table(output.table, summary.file, quote=F, sep="\t", row.names=F)