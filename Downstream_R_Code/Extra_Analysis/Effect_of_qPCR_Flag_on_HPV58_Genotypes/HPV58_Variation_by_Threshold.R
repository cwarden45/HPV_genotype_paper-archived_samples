#without qPCR flags
#input.file = "../../Selected_Output_Files/combined_genotype_with_year_and_ethnicity_freq5.txt"

#with qPCR flags
input.file = "../../Selected_Output_Files/combined_genotype_with_year_and_ethnicity_freq5-FLAGGED.txt"


parse.HPV58.status = function(string){
	if(length(grep("^HPV58",string))!=0){
		return("Plurality")
	}else if(length(grep("HPV58",string))!=0){
		return("Positive")
	}else{
		return("Negative")
	}
}#end def parse.HPV58.status

parse.HPV58.percent = function(arr){
	HPV.text = arr[9]
	percent.text = arr[10]
	if(length(grep("HPV58",HPV.text))==0){
		return("Negative")
	}else{
		percent.text = gsub("\\%","",percent.text)
		HPV.info = unlist(strsplit(HPV.text, split=","))
		percent.info = as.double(unlist(strsplit(percent.text, split=",")))
		for (i in 1:length(HPV.info)){
			if(HPV.info[i] == "HPV58"){
				if(percent.info[i] >= 50){
					return("Majority")
				}else{
					return("Minority")
				}
			}
		}#end for (i in 1:length(HPV.info))
	}#end else
}#end def parse.HPV58.percent

parse.HPV58.percent.20 = function(arr){
	HPV.text = arr[9]
	percent.text = arr[10]
	if(length(grep("HPV58",HPV.text))==0){
		return("Negative")
	}else{
		percent.text = gsub("\\%","",percent.text)
		HPV.info = unlist(strsplit(HPV.text, split=","))
		percent.info = as.double(unlist(strsplit(percent.text, split=",")))
		for (i in 1:length(HPV.info)){
			if(HPV.info[i] == "HPV58"){
				if(percent.info[i] >= 20){
					return("Pos.2")
				}else{
					return("Neg.2")
				}
			}
		}#end for (i in 1:length(HPV.info))
	}#end else
}#end def parse.HPV58.percent.20

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

HPV.status = sapply(as.character(input.table$genotype), parse.HPV58.status)

batch.count = tapply(batch,batch,length)

#sum "plurality" + "positive" to get full positive count
status.counts = table(HPV.status,batch)
print(status.counts)
status.percent = apply(status.counts, 1, function(type, counts){return(sprintf(type/counts,fmt="%#.3f"))},counts = batch.count)
print(t(status.percent))

FE.mat.freq5 = matrix(c(status.counts[1,1],sum(status.counts[2:3,1]),
						status.counts[1,3],sum(status.counts[2:3,3])),
						ncol=2)
colnames(FE.mat.freq5) = c("DNA","FFPE")
rownames(FE.mat.freq5) = c("Negative","Positive")
print(FE.mat.freq5)

	result = fisher.test(FE.mat.freq5)
	print(paste("FFPE vs DNA, FE P-value (5%) :", result$p.value, sep=""))

FE.mat.plurality = matrix(c(status.counts[1,1],status.counts[2,1],
						status.counts[1,3],status.counts[2,3]),
						ncol=2)
colnames(FE.mat.plurality) = c("DNA","FFPE")
rownames(FE.mat.plurality) = c("Negative","Positive")
print(FE.mat.plurality)

	result = fisher.test(FE.mat.plurality)
	print(paste("FFPE vs DNA, FE P-value (Plurality) :", result$p.value, sep=""))

mod.detect = apply(input.table, 1, parse.HPV58.percent.20)
status.counts = table(mod.detect, batch)
print(status.counts)
status.percent = apply(status.counts, 1, function(type, counts){return(sprintf(type/counts,fmt="%#.3f"))},counts = batch.count)
print(t(status.percent))

FE.mat.freq15 = matrix(c(sum(status.counts[1:2,1]),status.counts[3,1],
						sum(status.counts[1:2,3]),status.counts[3,3]),
						ncol=2)
colnames(FE.mat.freq15) = c("DNA","FFPE")
rownames(FE.mat.freq15) = c("Negative","Positive")
print(FE.mat.freq15)

	result = fisher.test(FE.mat.freq15)
	print(paste("FFPE vs DNA, FE P-value (20%) :", result$p.value, sep=""))


#alternative strategy for calling plurality
HPV.plurality = input.table[grep("^HPV58",input.table$genotype),]
status.counts = table(HPV.plurality$batch)
#print(status.counts)
status.percent = status.counts / batch.count
#print(status.percent)

#calculate majority
majority.percent = apply(input.table, 1, parse.HPV58.percent)
status.counts = table(majority.percent,batch)
print(status.counts)
status.percent = apply(status.counts, 1, function(type, counts){return(sprintf(type/counts,fmt="%#.3f"))},counts = batch.count)
print(t(status.percent))

FE.mat.majority = matrix(c(status.counts[3,1],status.counts[1,1],
						status.counts[3,3],status.counts[1,3]),
						ncol=2)
colnames(FE.mat.majority) = c("DNA","FFPE")
rownames(FE.mat.majority) = c("Negative","Positive")
print(FE.mat.majority)

	result = fisher.test(FE.mat.majority)
	print(paste("FFPE vs DNA, FE P-value (Majority) :", result$p.value, sep=""))

#alternative counts from https://www.nature.com/articles/6601024/tables/2
FE.mat.majority = matrix(c((5646 + 2175)-(0.03 * 5646 + 0.069 * 2175),0.03 * 5646 + 0.069 * 2175,
						status.counts[3,1],status.counts[1,1]),
						ncol=2)
FE.mat.majority = round(FE.mat.majority)
colnames(FE.mat.majority) = c("Clifford","FFPE")
rownames(FE.mat.majority) = c("Negative","Positive")
#print(FE.mat.majority)

	result = fisher.test(FE.mat.majority)
	#print(paste("DNA vs Public SCC+HSIL, FE P-value (Majority) :", result$p.value, sep=""))
	
FE.mat.majority = matrix(c((5646 + 2175)-(0.03 * 5646 + 0.069 * 2175),0.03 * 5646 + 0.069 * 2175,
						status.counts[3,3],status.counts[1,3]),
						ncol=2)
FE.mat.majority = round(FE.mat.majority)
colnames(FE.mat.majority) = c("Clifford","FFPE")
rownames(FE.mat.majority) = c("Negative","Positive")
#print(FE.mat.majority)

	result = fisher.test(FE.mat.majority)
	#print(paste("FFPE vs Public SCC+HSIL, FE P-value (Majority) :", result$p.value, sep=""))
#the p-value for this comparison is so low that I wonder if a different strategy should be used for the p-value calculation (although the DNA p-value was ~1.00)
