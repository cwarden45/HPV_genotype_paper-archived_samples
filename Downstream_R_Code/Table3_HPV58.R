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
	HPV.text = arr[8]
	percent.text = arr[9]
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

parse.HPV58.percent.15 = function(arr){
	HPV.text = arr[8]
	percent.text = arr[9]
	if(length(grep("HPV58",HPV.text))==0){
		return("Negative")
	}else{
		percent.text = gsub("\\%","",percent.text)
		HPV.info = unlist(strsplit(HPV.text, split=","))
		percent.info = as.double(unlist(strsplit(percent.text, split=",")))
		for (i in 1:length(HPV.info)){
			if(HPV.info[i] == "HPV58"){
				if(percent.info[i] >= 15){
					return("Pos.2")
				}else{
					return("Neg.2")
				}
			}
		}#end for (i in 1:length(HPV.info))
	}#end else
}#end def parse.HPV58.percent.15

input.table = read.table("Selected_Output_Files/combined_genotype_with_year_and_ethnicity_freq5.txt", head=T, sep="\t")
print(dim(input.table))
input.table = input.table[-grep("prostate", input.table$sample.type),]
print(dim(input.table))
input.table = input.table[-grep("adjacent normal", input.table$sample.type),]
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

mod.detect = apply(input.table, 1, parse.HPV58.percent.15)
status.counts = table(mod.detect, batch)
print(status.counts)
status.percent = apply(status.counts, 1, function(type, counts){return(sprintf(type/counts,fmt="%#.3f"))},counts = batch.count)
print(t(status.percent))

#alternative strategy for calling plurality
HPV.plurality = input.table[grep("^HPV58",input.table$genotype),]
status.counts = table(HPV.plurality$batch)
#print(status.counts)
status.percent = status.counts / batch.count
#print(status.percent)

majority.percent = apply(input.table, 1, parse.HPV58.percent)
status.counts = table(majority.percent,batch)
print(status.counts)
status.percent = apply(status.counts, 1, function(type, counts){return(sprintf(type/counts,fmt="%#.3f"))},counts = batch.count)
print(t(status.percent))