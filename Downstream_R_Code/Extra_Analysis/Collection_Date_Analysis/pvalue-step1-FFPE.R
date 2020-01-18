compID = "FFPE"

avgGroupExpression = function (geneExpr, groups) {
	avg.expr = tapply(geneExpr, groups, mean)
	return(avg.expr)
}#end def avgGroupExpression

ab2fc = function(trt.ab, ref.ab, min.ab)
{
	ratio = (trt.ab+min.ab)/(ref.ab+min.ab)
	if(ratio >= 1){
		return(ratio)
	} else {
		return (-1/ratio)
	}
}#end def ab2fc

count.to.ab = function(arr, total.reads)
{
	return(round(100*(arr/total.reads), digits=1))
}#end def ratio2fc

count.defined.values = function(arr, expr.cutoff)
{
	return(length(arr[arr > expr.cutoff]))
}#end def count.values

max.defined.values = function(arr, groups, expr.cutoff){
	temp.defined.values = tapply(arr,groups, count.defined.values, expr.cutoff=expr.cutoff)
	total.counts = table(groups)
	grp.freq = temp.defined.values / table(groups)
	return(max(grp.freq))
}#end def max.defined.values

gene.lm = function(arr, var1, var2=c(), var3=c())
{	
	if (length(var2) == 0){
		fit = lm(as.numeric(arr) ~ var1)
		result = summary(fit)
		pvalue = result$coefficients[2,4]
	} else if (length(var3) == 0){
		fit = lm(as.numeric(arr) ~ var1 + var2)
		result = summary(fit)
		pvalue = result$coefficients[2,4]
	} else {
		fit = lm(as.numeric(arr) ~ var2*var3 + var2 + var3)
		result = summary(fit)
		pvalue = result$coefficients[4,4]
	}
	return(pvalue)
}#end def gene.lm

calc.fe.pvalue = function(binary.arr, grp){
	fe.mat = table(binary.arr, grp)
	result = fisher.test(fe.mat)
	return(result$p.value)
}#end def calc.fe.pvalue

#use 15% read fraction genotypes
meta.table = read.table("../../Selected_Output_Files/combined_genotype_with_year_and_ethnicity_freq15.txt", head=T, sep = "\t")
print(dim(meta.table))
meta.table = meta.table[meta.table$HPV.status == "pos",]
print(dim(meta.table))
meta.table = meta.table[meta.table$genotype != "unclear",]
print(dim(meta.table))
meta.table = meta.table[-grep(".N",meta.table$SAMPLEID),]
print(dim(meta.table))
meta.table$batch = as.character(meta.table$batch)
meta.table$batch[meta.table$batch == "161007"] = "DNA"
meta.table$batch[meta.table$batch == "161206"] = "Frozen"
meta.table$batch[meta.table$batch == "170118"] = "FFPE"
meta.table$batch = factor(as.character(meta.table$batch), levels=c("DNA","Frozen","FFPE"))

### start comparison-specific code ###
meta.table = meta.table[!is.na(meta.table$collection.year)&(meta.table$batch == "FFPE"),]
print(dim(meta.table))
### end  comparison-specific code (until p-value calculation)###

count.table = read.table("../../Public_Input_Files/PE_HPVtype_counts_final_names.txt",head=T, sep="\t")
total.counts = as.numeric(count.table[count.table$HPV.type == "TOTAL",2:ncol(count.table)])
count.mat = count.table[1:(nrow(count.table)-2),2:ncol(count.table)]
rownames(count.mat) = count.table$HPV.type[1:(nrow(count.table)-2)]

ab.mat = t(apply(count.mat, 1, count.to.ab, total=total.counts))

ab.mat = ab.mat[,match(meta.table$SAMPLEID,colnames(ab.mat))]
total.counts = total.counts[match(meta.table$SAMPLEID,names(count.mat))]
count.mat = count.mat[,match(meta.table$SAMPLEID,names(count.mat))]

fdr.cutoff = 0.05

group = meta.table$collection.year
percent.human = meta.table$human.percent
percent.human = as.character(percent.human)
percent.human = gsub("%","",percent.human)
percent.human = as.numeric(percent.human)

#FE-code
meta.table$genotype = as.factor(as.character(meta.table$genotype))
mixed.geno = levels(meta.table$genotype)[grep(",",levels(meta.table$genotype))]

genotypes = levels(meta.table$genotype)[-grep(",",levels(meta.table$genotype))]
print(length(genotypes))
for (geno in mixed.geno){
	geno.arr = unlist(strsplit(geno,split=","))
	for (j in 1:length(geno.arr)){
		if(!(geno.arr[j] %in% genotypes)){
			genotypes = c(genotypes, geno.arr[j])
		}
	}#end for (j in 1:length(geno.arr))
}#end for (geno in levels(meta.table$genotype))
print(length(genotypes))
genotypes = sort(genotypes)

status.mat = matrix(0,ncol=length(genotypes),nrow=nrow(meta.table))
colnames(status.mat)=genotypes
rownames(status.mat)=rep("",nrow(meta.table))

for (i in 1:nrow(meta.table)){
	geno.text = as.character(meta.table$genotype[i])
	geno.arr = unlist(strsplit(geno.text,split=","))
	for (j in 1:length(geno.arr)){
		status.mat[i,colnames(status.mat)==geno.arr[j]]=1
	}#end for (j in 1:length(geno.arr))
}#end for (i in 1:nrow(meta.table))

FE.group = rep(NA,nrow(meta.table))
FE.group[!is.na(meta.table$collection.year) & (meta.table$collection.year >= 2000)]="gt2000"
FE.group[!is.na(meta.table$collection.year) & (meta.table$collection.year < 2000)]="lt2000"

percent.detected.by.group = function(status, group){
	counts = tapply(status, group, sum)
	total = tapply(status, group, length)
	return(counts / total)
}#end def percent.detected.by.group

grp.percent.detected = round(100*data.frame(t(apply(status.mat, 2, percent.detected.by.group, group=FE.group))), digits=1)
detection.trt = grp.percent.detected$gt2000
detection.ref = grp.percent.detected$lt2000
detection.diff = detection.trt - detection.ref

fe.pvalue = apply(status.mat, 2, calc.fe.pvalue, FE.group)
fe.fdr = p.adjust(fe.pvalue, "fdr")
fe.status = rep("No Change",length(fe.fdr))
fe.status[(fe.fdr < fdr.cutoff)&(detection.diff > 10)]=">2000 Up"
fe.status[(fe.fdr < fdr.cutoff)&(detection.diff < -10)]="<2000 Down"

print(genotypes[fe.fdr<0.05])

deg.table = data.frame(HPV.subtype = genotypes, grp.percent.detected, detection.diff,
						fe.pvalue=fe.pvalue,fe.fdr=fe.fdr, fe.status=fe.status)
write.table(deg.table,paste(compID,"_differential_status_FE_test.txt",sep=""), row.names=F, sep="\t", quote=F)

#KS-Test

avg.date.by.status = function(status, dates){
	return.arr = c(NA,NA)
	names(return.arr)=c("avg.pos","avg.neg")
	return.arr[1]=round(mean(dates[status == 1], na.rm=T))
	return.arr[2]=round(mean(dates[status == 0], na.rm=T))
	return(return.arr)
}#end def avg.date.by.status

avg.date = data.frame(t(apply(status.mat, 2, avg.date.by.status, dates=meta.table$collection.year)))
avg.trt = avg.date$avg.pos
avg.ref = avg.date$avg.neg
avg.diff = avg.trt - avg.ref

ks.date = function(status, dates){
	#compare to overall distribution
	pos.dates = dates[status == 1]
	result = ks.test(pos.dates, dates)
	return(result$p.value)
}#end def ks.date


ks.pvalue = apply(status.mat, 2, ks.date, meta.table$collection.year)
ks.fdr = p.adjust(ks.pvalue, "fdr")
ks.status = rep("No Change",length(ks.fdr))
ks.status[(ks.fdr < fdr.cutoff)&(avg.diff > 0)]=">2000 Up"
ks.status[(ks.fdr < fdr.cutoff)&(avg.diff < 0)]="<2000 Down"

print(genotypes[ks.fdr<0.05])

deg.table = data.frame(HPV.subtype = genotypes, avg.date, avg.diff,
						ks.pvalue=ks.pvalue,ks.fdr=ks.fdr, ks.status=ks.status)
write.table(deg.table,paste(compID,"_differential_status_KS_test.txt",sep=""), row.names=F, sep="\t", quote=F)

#binomial GLM
calc.bGLM.pvalue = function(arr, var1, var2=c(), var3=c()){
	if (length(var2) == 0){
		fit = glm(as.numeric(arr) ~ var1, family="binomial")
		result = summary(fit)
		pvalue = result$coefficients[2,4]
	} else if (length(var3) == 0){
		fit = glm(as.numeric(arr) ~ var1 + var2, family="binomial")
		result = summary(fit)
		pvalue = result$coefficients[2,4]
	} else {
		fit = glm(as.numeric(arr) ~ var1 + var2 + var3, family="binomial")
		result = summary(fit)
		pvalue = result$coefficients[2,4]
	}
	return(pvalue)
}#end def percent.detected.by.group

bGLM.pvalue = apply(status.mat, 2, calc.bGLM.pvalue, var1=meta.table$collection.year)
bGLM.fdr = p.adjust(bGLM.pvalue, "fdr")
bGLM.status = rep("No Change",length(bGLM.fdr))
bGLM.status[(bGLM.fdr < fdr.cutoff)]="Variable"

print(genotypes[bGLM.fdr<0.05])

deg.table = data.frame(HPV.subtype = genotypes, 
						bGLM.pvalue=bGLM.pvalue, bGLM.fdr=bGLM.fdr, bGLM.status=bGLM.status)
write.table(deg.table,paste(compID,"_differential_status_binomial_GLM.txt",sep=""), row.names=F, sep="\t", quote=F)

#limma-voom code
ab.mat = ab.mat[match(genotypes, rownames(ab.mat)),]
count.mat = count.mat[match(genotypes, rownames(count.mat)),]

calc.gene.cor = function(arr, indep.var)
{	
	na.count = length(arr[!is.na(arr)])
	if((na.count >= 3) & (sd(arr) != 0)){
		gene.cor.coef = cor(arr,indep.var)
	} else {
		gene.cor.coef = NA
	}
	return(gene.cor.coef)
}#end def calc.gene.cor

hpv.cor = apply(ab.mat, 1, calc.gene.cor, indep.var=group)

library("edgeR")
design = model.matrix(~group+percent.human)
y <- DGEList(counts=count.mat, genes=genotypes, lib.size=total.counts)
png(paste(compID,"_voom_plot.png",sep=""))
v <- voom(y,design,plot=TRUE)
dev.off()
fit <- lmFit(v,design)
fit = eBayes(fit)
pvalue.mat = data.frame(fit$p.value)
limma.pvalue = pvalue.mat[,2]
limma.fdr = p.adjust(limma.pvalue, "fdr")
limma.status = rep("No Change",length(limma.fdr))
limma.status[(limma.fdr < fdr.cutoff)&(hpv.cor > 0.2)]="Increased Abundance"
limma.status[(limma.fdr < fdr.cutoff)&(hpv.cor < -0.2)]="Decreased Abundance"

print(genotypes[limma.fdr<0.05])

deg.table = data.frame(HPV.subtype = genotypes, hpv.cor,
						limma.pvalue=limma.pvalue,limma.fdr=limma.fdr, limma.status=limma.status)
write.table(deg.table,paste(compID,"_differential_abundance_limma-voom.txt",sep=""), row.names=F, sep="\t", quote=F)

#lm code

calc.lm.pvalue = function(arr, var1, var2=c(), var3=c())
{	
	if (length(var2) == 0){
		fit = lm(as.numeric(arr) ~ var1)
		result = summary(fit)
		pvalue = result$coefficients[2,4]
	} else if (length(var3) == 0){
		fit = lm(as.numeric(arr) ~ var1 + var2)
		result = summary(fit)
		pvalue = result$coefficients[2,4]
	} else {
		fit = lm(as.numeric(arr) ~ var1 + var2+var3)
		result = summary(fit)
		pvalue = result$coefficients[2,4]
	}
	return(pvalue)
}#end def calc.lm.pvalue

lm.pvalue = apply(ab.mat, 1, calc.lm.pvalue, var1=meta.table$collection.year, var2=percent.human)
lm.pvalue[is.na(lm.pvalue)]=1
lm.fdr = p.adjust(lm.pvalue, "fdr")
lm.status = rep("No Change",length(lm.fdr))
lm.status[(lm.fdr < fdr.cutoff)]="Variable"

print(genotypes[lm.fdr<0.05])

deg.table = data.frame(HPV.subtype = genotypes, 
						lm.pvalue=lm.pvalue, lm.fdr=lm.fdr, lm.status=lm.status)
write.table(deg.table,paste(compID,"_differential_status_linear_regression.txt",sep=""), row.names=F, sep="\t", quote=F)