meta.ALL5 = read.table("../../Selected_Output_Files/combined_genotype_with_year_and_ethnicity_freq5.txt",head=T, sep="\t")
meta.qPCR20 = read.table("../../Selected_Output_Files/combined_genotype_with_year_and_ethnicity_freq20-FLAGGED.txt",head=T, sep="\t")

create_genotype_table = function(genotype_str){
	genotype_arr = as.factor(as.character(genotype_str))
	mixed.geno = levels(genotype_arr)[grep(",",levels(genotype_arr))]

	if(length(mixed.geno) > 0){
		genotypes = levels(genotype_arr)[-grep(",",levels(genotype_arr))]
	}else{
		genotypes = levels(genotype_arr)
	}#end else
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

	status.mat = matrix(0,ncol=length(genotypes),nrow=length(genotype_arr))
	colnames(status.mat)=genotypes
	rownames(status.mat)=names(genotype_str)

	for (i in 1:length(genotype_arr)){
		geno.text = as.character(genotype_arr[i])
		geno.arr = unlist(strsplit(geno.text,split=","))
		for (j in 1:length(geno.arr)){
			status.mat[i,colnames(status.mat)==geno.arr[j]]=1
		}#end for (j in 1:length(geno.arr))
	}#end for (i in 1:nrow(meta.table))
	return(status.mat)
}#end def create_genotype_table

########################################
###  5% assignments for all samples  ###
########################################

print(dim(meta.ALL5))
meta.ALL5$hist.subtype = factor(meta.ALL5$hist.subtype,
								levels = c("Squamous Cell Carcinoma (SCC)","Adenocarcinoma (Adeno)","Adenosquamous"))
meta.ALL5 = meta.ALL5[!is.na(meta.ALL5$hist.subtype),]
print(dim(meta.ALL5))

meta.ALL5$batch = as.character(meta.ALL5$batch)
meta.ALL5$batch[meta.ALL5$batch == 161007] = "DNA"
meta.ALL5$batch[meta.ALL5$batch == 161206] = "Frozen"
meta.ALL5$batch[meta.ALL5$batch == 170118] = "FFPE"
meta.ALL5$batch = factor(meta.ALL5$batch, levels=c("DNA","Frozen","FFPE"))

geno_obj=meta.ALL5$genotype
names(geno_obj)= meta.ALL5$SAMPLEID
status_table.ALL5 = data.frame(create_genotype_table(geno_obj))

geno_summary.ALL5 = rep("Other",nrow(status_table.ALL5))
geno_summary.ALL5[(status_table.ALL5$HPV16 == 1)]="HPV16"
geno_summary.ALL5[(status_table.ALL5$HPV18 == 1)]="HPV18"
geno_summary.ALL5[(status_table.ALL5$HPV16 == 1) & (status_table.ALL5$HPV18 == 1)]="HPV16+HPV18"
geno_summary.ALL5 = factor(geno_summary.ALL5, levels=c("HPV16","HPV18","HPV16+HPV18","Other"))

print(table(meta.ALL5$hist.subtype,meta.ALL5$batch))

print("DNA")
print(table(meta.ALL5$hist.subtype[meta.ALL5$batch == "DNA"],geno_summary.ALL5[meta.ALL5$batch == "DNA"]))
print("Frozen")
print(table(meta.ALL5$hist.subtype[meta.ALL5$batch == "Frozen"], geno_summary.ALL5[meta.ALL5$batch == "Frozen"]))
print("FFPE")
print(table(meta.ALL5$hist.subtype[meta.ALL5$batch == "FFPE"], geno_summary.ALL5[meta.ALL5$batch == "FFPE"]))

###################################################
###  20% assignments for qPCR-filtered samples  ###
###################################################

print(dim(meta.qPCR20))
meta.qPCR20 = meta.qPCR20[meta.qPCR20$HPV.status != "qPCR Flag",]
print(dim(meta.qPCR20))
meta.qPCR20$hist.subtype = factor(meta.qPCR20$hist.subtype,
								levels = c("Squamous Cell Carcinoma (SCC)","Adenocarcinoma (Adeno)","Adenosquamous"))
meta.qPCR20 = meta.qPCR20[!is.na(meta.qPCR20$hist.subtype),]
print(dim(meta.qPCR20))

meta.qPCR20$batch = as.character(meta.qPCR20$batch)
meta.qPCR20$batch[meta.qPCR20$batch == 161007] = "DNA"
meta.qPCR20$batch[meta.qPCR20$batch == 161206] = "Frozen"
meta.qPCR20$batch[meta.qPCR20$batch == 170118] = "FFPE"
meta.qPCR20$batch = factor(meta.qPCR20$batch, levels=c("DNA","Frozen","FFPE"))

geno_obj=meta.qPCR20$genotype
names(geno_obj)= meta.qPCR20$SAMPLEID
status_table.qPCR20 = data.frame(create_genotype_table(geno_obj))

geno_summary.qPCR20 = rep("Other",nrow(status_table.qPCR20))
geno_summary.qPCR20[(status_table.qPCR20$HPV16 == 1)]="HPV16"
geno_summary.qPCR20[(status_table.qPCR20$HPV18 == 1)]="HPV18"
geno_summary.qPCR20[(status_table.qPCR20$HPV16 == 1) & (status_table.qPCR20$HPV18 == 1)]="HPV16+HPV18"
geno_summary.qPCR20 = factor(geno_summary.qPCR20, levels=c("HPV16","HPV18","HPV16+HPV18","Other"))

print(table(meta.qPCR20$hist.subtype,meta.qPCR20$batch))

print("DNA")
print(table(meta.qPCR20$hist.subtype[meta.qPCR20$batch == "DNA"],geno_summary.qPCR20[meta.qPCR20$batch == "DNA"]))
print("Frozen")
print(table(meta.qPCR20$hist.subtype[meta.qPCR20$batch == "Frozen"], geno_summary.qPCR20[meta.qPCR20$batch == "Frozen"]))
print("FFPE")
print(table(meta.qPCR20$hist.subtype[meta.qPCR20$batch == "FFPE"], geno_summary.qPCR20[meta.qPCR20$batch == "FFPE"]))
