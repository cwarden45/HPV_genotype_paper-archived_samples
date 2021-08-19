parse.HPV.status = function(string, type){
	if(length(grep(type,string))!=0){
		return("Positive")
	}else{
		return("Negative")
	}
}#end def parse.HPV.status

input.table = read.table("Selected_Output_Files/combined_genotype_with_year_and_ethnicity_freq20-FLAGGED.txt", head=T, sep="\t")
print(dim(input.table))
input.table = input.table[-grep("Prostate", input.table$sample.type),]
print(dim(input.table))
input.table = input.table[input.table$sample.type != "Adjacent Normal Cervix",]
print(dim(input.table))
input.table = input.table[input.table$sample.type != "Non-Malignant Vagina",]
print(dim(input.table))
input.table = input.table[input.table$HPV.status != "qPCR Flag",]
input.table$HPV.status = as.factor(as.character(input.table$HPV.status))
print(dim(input.table))

batch = rep(NA,nrow(input.table))
batch[grep("^13",as.character(input.table$Sample))]="161007"
batch[grep("^14",as.character(input.table$Sample))]="161206"
batch[grep("^15",as.character(input.table$Sample))]="170118"

#overall
total.samples = table(batch)
overall.table = table(input.table$HPV.status, batch)
print(overall.table)
print(overall.table[2,]/total.samples)

#per HPV-type
HPV.types = c("HPV16","HPV18","HPV31","HPV33","HPV45","HPV58","HPV59","HPV67","HPV73","unclear")

batches = levels(as.factor(batch))
count.table = matrix(0, ncol=3, nrow=length(HPV.types))
colnames(count.table) = batches
rownames(count.table) = HPV.types
percent.table = matrix(0, ncol=3, nrow=length(HPV.types))
colnames(percent.table) = batches
rownames(percent.table) = HPV.types

for (i in 1:length(batches)){
	batch.genos = input.table$genotype[input.table$batch ==batches[i]]
	for (j in 1:length(HPV.types)){
		HPV.status = sapply(as.character(batch.genos), parse.HPV.status, type=HPV.types[j])
		type.count = length(HPV.status[HPV.status == "Positive"])
		count.table[j,i]=type.count
		percent.table[j,i]=type.count / length(batch.genos)
	}#end for (j in 1:length(HPV.types))
}#end for (i in 1:length(batches))

print(count.table)
print(percent.table)

#co-infections
count.coinfections = function(char){
	if(is.na(char)){
		return(0)
	}else{
		types = unlist(strsplit(char,split=","))
		return(length(types))
	}
}#end def count.coinfections

coinfection.count = 0:4

batches = levels(as.factor(batch))
count.table = matrix(0, ncol=3, nrow=length(coinfection.count))
colnames(count.table) = batches
rownames(count.table) = coinfection.count
percent.table = matrix(0, ncol=3, nrow=length(coinfection.count))
colnames(percent.table) = batches
rownames(percent.table) = coinfection.count

for (i in 1:length(batches)){
	batch.genos = input.table$genotype[input.table$batch ==batches[i]]
	infection.count = unlist(sapply(as.character(batch.genos), count.coinfections))
	for (j in 1:length(coinfection.count)){
		level.count = length(infection.count[infection.count == coinfection.count[j]])
		count.table[j,i]=level.count
		percent.table[j,i]=level.count / length(batch.genos)
	}#end for (j in 1:length(HPV.types))
}#end for (i in 1:length(batches))

print(count.table)
print(percent.table)

#overall
count.arr = rep(0, length(HPV.types))
percent.arr = rep(0, length(HPV.types))
names(count.arr) = HPV.types
names(percent.arr) = HPV.types

	for (j in 1:length(HPV.types)){
		HPV.status = sapply(as.character(input.table$genotype), parse.HPV.status, type=HPV.types[j])
		type.count = length(HPV.status[HPV.status == "Positive"])
		count.arr[j]=type.count
		percent.arr[j]=type.count / nrow(input.table)
	}#end for (j in 1:length(HPV.types))

print(count.arr)
print(percent.arr)

count.arr = rep(0, length(coinfection.count))
percent.arr = rep(0, length(coinfection.count))
names(count.arr) = coinfection.count
names(percent.arr) = coinfection.count

	infection.count = unlist(sapply(as.character(input.table$genotype), count.coinfections))
	for (j in 1:length(coinfection.count)){
		level.count = length(infection.count[infection.count == coinfection.count[j]])
		count.arr[j]=level.count
		percent.arr[j]=level.count/ nrow(input.table)
	}#end for (j in 1:length(HPV.types))

print(count.arr)
print(percent.arr)