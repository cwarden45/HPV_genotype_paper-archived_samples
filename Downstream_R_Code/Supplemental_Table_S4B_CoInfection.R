#without qPCR flags
#input.file = "Selected_Output_Files/combined_genotype_with_year_and_ethnicity_freq5.txt"

#with qPCR flags
input.file = "Selected_Output_Files/combined_genotype_with_year_and_ethnicity_freq5-FLAGGED.txt"

parse.coinfection.count = function(string){
	if(is.na(string)){
		return(0)
	}else if(string=="unclear"){
		return(0)
	}else{
		HPV.arr = unlist(strsplit(string, split=","))
		return(length(HPV.arr))
	}
}#end def parse.coinfection.count

parse.coinfection.count2 = function(string, threshold){
	if(is.na(string)){
		return(0)
	}else if(string=="unclear"){
		return(0)
	}else{
		HPV.arr = unlist(strsplit(string, split=","))
		HPV.arr = gsub("\\%","",HPV.arr)
		HPV.arr = as.numeric(HPV.arr)
		HPV.arr = HPV.arr[HPV.arr > threshold]
		return(length(HPV.arr))
	}
}#end def parse.coinfection.count

input.table = read.table(input.file, head=T, sep="\t")
print(dim(input.table))
input.table = input.table[-grep("prostate", input.table$sample.type),]
print(dim(input.table))
input.table = input.table[-grep("adjacent normal", input.table$sample.type),]
print(dim(input.table))
input.table = input.table[input.table$HPV.status != "qPCR Flag",]
print(dim(input.table))

batch = rep(NA,nrow(input.table))
batch[grep("^13",as.character(input.table$Sample))]="161007"
batch[grep("^14",as.character(input.table$Sample))]="161206"
batch[grep("^15",as.character(input.table$Sample))]="170118"

batch.count = tapply(batch,batch,length)
print(batch.count)

#5%
coinfection.count = sapply(as.character(input.table$genotype), parse.coinfection.count)

HPV.counts = table(coinfection.count, batch)
print(HPV.counts)
coinfection.percent = t(apply(HPV.counts, 1, function(type, counts){return(sprintf(type/counts,fmt="%#.3f"))},counts = batch.count))
print(coinfection.percent)

		fit = aov(coinfection.count ~ batch)
		result = summary(fit)
		aov.pvalue = result[[1]][['Pr(>F)']][1]

print(paste("ANOVA Sample Type co-infection p-value: ",aov.pvalue,sep=""))

#20%
coinfection.count = sapply(as.character(input.table$genotype.percent), parse.coinfection.count2, threshold=20)

HPV.counts = table(coinfection.count, batch)
print(HPV.counts)
coinfection.percent = t(apply(HPV.counts, 1, function(type, counts){return(sprintf(type/counts,fmt="%#.3f"))},counts = batch.count))
print(coinfection.percent)

#50%
coinfection.count = sapply(as.character(input.table$genotype.percent), parse.coinfection.count2, threshold=50)

HPV.counts = table(coinfection.count, batch)
print(HPV.counts)
coinfection.percent = t(apply(HPV.counts, 1, function(type, counts){return(sprintf(type/counts,fmt="%#.3f"))},counts = batch.count))
print(coinfection.percent)