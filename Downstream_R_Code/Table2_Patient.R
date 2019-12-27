parse.HPV.status = function(string, type){
	if(length(grep(type,string))!=0){
		return("Positive")
	}else{
		return("Negative")
	}
}#end def parse.HPV.status

input.table = read.table("Selected_Output_Files/combined_genotype_with_year_and_ethnicity.txt", head=T, sep="\t")
print(dim(input.table))
input.table = input.table[-grep("prostate", input.table$sample.type),]
print(dim(input.table))
input.table = input.table[-grep("adjacent normal", input.table$sample.type),]
print(dim(input.table))

batch = rep(NA,nrow(input.table))
batch[grep("^13",as.character(input.table$Sample))]="161007"
batch[grep("^14",as.character(input.table$Sample))]="161206"
batch[grep("^15",as.character(input.table$Sample))]="170118"

input.table$extra.pairID

patientID = as.character(input.table$extra.pairID)
patientID[is.na(patientID)]= as.character(input.table$SAMPLEID[is.na(patientID)])

unique.patientID = unique(patientID)
print(length(unique.patientID))

HPV.status=c()
HPV.geno=c()

for (i in 1:length(unique.patientID)){
	test.patient = unique.patientID[i]
	patient.mat = input.table[patientID == test.patient,]
	
	if(nrow(patient.mat)==1){
		if(patient.mat$HPV.status == "pos"){
			HPV.status[i]="pos"
			HPV.geno[i]=as.character(patient.mat$genotype)
		}else if(patient.mat$HPV.status == "unclear"){
			HPV.status[i]="unclear"
			HPV.geno[i]=NA
		}else if(patient.mat$HPV.status == "neg"){
			HPV.status[i]="neg"
			HPV.geno[i]=NA
		}else{
			print("Parse Result:")
			print(patient.mat)
			stop()
		}
	}else{
		if(length(patient.mat$HPV.status[grep("pos",patient.mat$HPV.status)])>=1){
			HPV.status[i]="pos"
			patient.mat = patient.mat[patient.mat$HPV.status=="pos",]
			HPV.geno[i]=paste(patient.mat$genotype,collapse=",")
			HPV.geno.arr = unlist(strsplit(HPV.geno[i],split=","))
			HPV.geno.arr=unique(HPV.geno.arr)
			HPV.geno[i] = paste(HPV.geno.arr, collapse=",")
		}else if(patient.mat$HPV.status == rep("neg",nrow(patient.mat))){
			HPV.status[i]="neg"
			HPV.geno[i]=NA
			print("Negative Result?")
			print(patient.mat)
			stop()
		}else{
			HPV.status[i]="unclear"
			HPV.geno[i]=NA
		}#end else
	}#end else
}#end for (test.patient in unique.patientID)

print(table(HPV.status))

#per HPV-type
HPV.types = c("HPV16","HPV18","HPV31","HPV33","HPV45","HPV58","HPV59","HPV67","HPV73","unclear")

for (type in HPV.types){
	print(paste(type," : ",
			length(HPV.geno[grep(type,HPV.geno)]),
			" (",round(100 * length(HPV.geno[grep(type,HPV.geno)]) / length(HPV.geno),digits=1),"%)",sep=""))
}#end for (type in HPV.types)

