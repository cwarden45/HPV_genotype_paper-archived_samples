compID = "FFPE_vs_Frozen"

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

count2detectedStatus = function(arr, min.ab){
	subtype.status = rep(0,length(arr))
	subtype.status[arr > min.ab] =1
	return(subtype.status)
}#end def count2detectedStatus

calc.fe.pvalue = function(binary.arr, grp){
	fe.mat = table(binary.arr, grp)
	result = fisher.test(fe.mat)
	return(result$p.value)
}#end def calc.fe.pvalue

#library(metagenomeSeq)
library(gplots)

pvalue.dir = compID
dir.create(pvalue.dir)

#run analysis with lower read fraction and **without** qPCR flag to filter samples
meta.table = read.table("../../Selected_Output_Files/combined_genotype_with_year_and_ethnicity_freq5.txt", head=T, sep = "\t")
print(dim(meta.table))
meta.table = meta.table[meta.table$HPV.status == "pos",]
print(dim(meta.table))
meta.table = meta.table[meta.table$genotype != "unclear",]
print(dim(meta.table))
meta.table = meta.table[(meta.table$sample.type == "Invasive Cervical Cancer")|(meta.table$sample.type == "Vulvar Cancer")|(meta.table$sample.type == "Endometrial + Cervical Cancer"),]
print(dim(meta.table))
meta.table$batch = as.character(meta.table$batch)
meta.table$batch[meta.table$batch == "161007"] = "DNA"
meta.table$batch[meta.table$batch == "161206"] = "Frozen"
meta.table$batch[meta.table$batch == "170118"] = "FFPE"
meta.table$batch = factor(as.character(meta.table$batch), levels=c("DNA","Frozen","FFPE"))

### start comparison-specific code ###
meta.table = meta.table[(meta.table$batch == "FFPE") | (meta.table$batch == "Frozen"),]
meta.table$batch = as.character(meta.table$batch)
### end  comparison-specific code (until p-value calculation)###

count.table = read.table("../../Public_Input_Files/PE_HPVtype_counts_final_names.txt",head=T, sep="\t")
total.counts = as.numeric(count.table[count.table$HPV.type == "TOTAL",2:ncol(count.table)])
count.mat = count.table[1:(nrow(count.table)-2),2:ncol(count.table)]
rownames(count.mat) = count.table$HPV.type[1:(nrow(count.table)-2)]

ab.mat = t(apply(count.mat, 1, count.to.ab, total=total.counts))

ab.mat = ab.mat[,match(meta.table$SAMPLEID,colnames(ab.mat))]
total.counts = total.counts[match(meta.table$SAMPLEID,names(count.mat))]
count.mat = count.mat[,match(meta.table$SAMPLEID,names(count.mat))]

#portion of code that was fixed on 3/14/2018 --> in this case, can really just use v2 script
#count.mat = count.mat[,match(meta.table$SAMPLEID,names(count.mat))]
#ab.mat = ab.mat[,match(meta.table$SAMPLEID,names(count.mat))]
#total.counts = total.counts[match(meta.table$SAMPLEID,names(count.mat))]

min.ab.sample = 5
min.ab.group = 1
min.freq = 0.01
trt.group = "FFPE"

fc.cutoff = 1.2
fdr.cutoff = 0.05

hpv.subtypes = rownames(count.mat)
group = meta.table$batch
percent.human = meta.table$human.percent
percent.human = as.character(percent.human)
percent.human = gsub("%","",percent.human)
percent.human = as.numeric(percent.human)

print(dim(count.mat))
max.sample.freq = apply(ab.mat, 2, max)
count.mat = count.mat[,max.sample.freq > min.ab.sample]
ab.mat = ab.mat[,max.sample.freq > min.ab.sample]
group = group[max.sample.freq > min.ab.sample]
total.counts = total.counts[max.sample.freq > min.ab.sample]
percent.human = percent.human[max.sample.freq > min.ab.sample]
max.sample.freq  = max.sample.freq[max.sample.freq > min.ab.sample]

print(dim(count.mat))
percent.detected = apply(ab.mat, 1, max.defined.values, expr.cutoff=min.ab.sample, groups=group)
count.mat = count.mat[percent.detected > min.freq,]
hpv.subtypes = hpv.subtypes[percent.detected > min.freq]
ab.mat = ab.mat[percent.detected > min.freq,]
percent.detected=percent.detected[percent.detected > min.freq]
print(dim(count.mat))

grp.ab = data.frame(t(apply(ab.mat, 1, avgGroupExpression, groups=group)))
trt.ab = grp.ab[,"FFPE"]
ref.ab = grp.ab[,"Frozen"]
grp.fc = round(mapply(ab2fc, trt.ab = trt.ab, ref.ab=ref.ab, min.ab=min.ab.group), digits=2)
colnames(grp.ab)=paste(colnames(grp.ab),"grp.avg.ab",sep=".")

#metagenomeSeq code
#pheno = data.frame(group)
#rownames(pheno) = colnames(count.mat)
#features = data.frame(subtype = hpv.subtypes)
#rownames(features) = rownames(count.mat)
#mrObj = newMRexperiment(count.mat, phenoData=AnnotatedDataFrame(pheno), featureData=AnnotatedDataFrame(features))
#scalePercentile = cumNormStatFast(mrObj)
#mrObj = cumNorm(mrObj, p = scalePercentile)
#design = model.matrix(~group)
#fit = fitFeatureModel(mrObj,design)
#test.pvalue = as.numeric(fit$pvalues)

#limma code - percent
#design = model.matrix(~group)
#fit = lmFit(ab.mat, design)
#fit = eBayes(fit)
#pvalue.mat = data.frame(fit$p.value)
#limma.pvalue = pvalue.mat[,2]
#limma.fdr = p.adjust(limma.pvalue, "fdr")

#limma code - log2(count +1)
#pseudocount = log2(count.mat+1)
#design = model.matrix(~group)
#fit = lmFit(pseudocount, design)
#fit = eBayes(fit)
#pvalue.mat = data.frame(fit$p.value)
#limma.pvalue = pvalue.mat[,2]
#limma.fdr = p.adjust(limma.pvalue, "fdr")

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

percent.detected.by.group = function(status, group){
	counts = tapply(status, group, sum)
	total = tapply(status, group, length)
	return(counts / total)
}#end def percent.detected.by.group

grp.percent.detected = round(100*data.frame(t(apply(status.mat, 2, percent.detected.by.group, group=meta.table$batch))), digits=1)
detection.trt = grp.percent.detected$FFPE 
detection.ref = grp.percent.detected$Frozen
detection.diff = detection.trt - detection.ref

fe.pvalue = apply(status.mat, 2, calc.fe.pvalue, grp=meta.table$batch)
fe.fdr = p.adjust(fe.pvalue, "fdr")
fe.status = rep("No Change",length(fe.fdr))
fe.status[(fe.fdr < fdr.cutoff)&(detection.diff > 10)]=paste(trt.group,"Up")
fe.status[(fe.fdr < fdr.cutoff)&(detection.diff < -10)]=paste(trt.group,"Down")

print(genotypes[fe.fdr<0.05])

deg.table = data.frame(HPV.subtype = genotypes, detection.diff,
						fe.pvalue=fe.pvalue,fe.fdr=fe.fdr, fe.status=fe.status)
write.table(deg.table,paste(pvalue.dir,"/differential_status_FE_test.txt",sep=""), row.names=F, sep="\t", quote=F)