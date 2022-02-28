HPV35_genotypes = "../../Selected_Output_Files/combined_genotype_with_year_and_ethnicity_freq20-FLAGGED.txt"
PAVE220_genotypes = "hg38_plus_PaVE_220HPV_genotype_calls.txt"
output.file = "genotype_concordance.png"

library("gplots")

HPV35.table = read.table(HPV35_genotypes, head=T, sep="\t")
print(dim(HPV35.table))
HPV35.table = HPV35.table[!is.na(HPV35.table$genotype),]
print(dim(HPV35.table))
HPV35.table = HPV35.table[HPV35.table$genotype != "qPCR Flag",]
print(dim(HPV35.table))
HPV35.table = HPV35.table[order(HPV35.table$SAMPLEID),]
HPV35.table = HPV35.table[order(HPV35.table$batch),]

PAVE220.table = read.table(PAVE220_genotypes, head=T, sep="\t")

unique.genotypes = c()
for(i in 1:nrow(HPV35.table)){
	temp_genotype_arr = unlist(strsplit(HPV35.table$genotype[i], split=","))
	
	unique.genotypes = union(unique.genotypes, temp_genotype_arr)
}#end for(i in nrow(HPV35.table))
for(i in 1:nrow(HPV35.table)){
	temp_genotype_arr = unlist(strsplit(PAVE220.table$genotype[i], split=","))
	
	unique.genotypes = union(unique.genotypes, temp_genotype_arr)
}#end for(i in nrow(PAVE220.table))
unique.genotypes = sort(unique.genotypes)
unique.genotypes = unique.genotypes[-grep("unclear",unique.genotypes)]

archive_type = rep(NA,nrow(HPV35.table))
archive_type[HPV35.table$batch == "161007"]="Archived DNA"
archive_type[HPV35.table$batch == "161206"]="Frozen Tissue"
archive_type[HPV35.table$batch == "170118"]="FFPE Tissue"

archiveCol = rep("gray",nrow(HPV35.table))
archiveCol[archive_type == "Archived DNA"]="darkgreen"
archiveCol[archive_type == "Frozen Tissue"]="orange"
archiveCol[archive_type == "FFPE Tissue"]="cyan"

concordance.mat = matrix(0, ncol=length(unique.genotypes), nrow=nrow(HPV35.table))
rownames(concordance.mat)=rep("",length(archive_type))
colnames(concordance.mat)=unique.genotypes

for(i in 1:nrow(HPV35.table)){
	HPV35_genotype = unlist(strsplit(HPV35.table$genotype[i], split=","))
	PAVE220_genotype = unlist(strsplit(PAVE220.table$genotype[PAVE220.table$Sample == HPV35.table$Sample[i]], split=","))
	
	if(all.equal(HPV35_genotype,PAVE220_genotype)){
		concordance.mat[i,match(HPV35_genotype, unique.genotypes)]=1
	}else{
		stop("Write code for display of discordant genotypes!")
	}
}#end for(i in nrow(HPV35.table))

#If there were any discordant genotypes, use this:
#colorPanel = colorpanel(3, low="red", mid="black", high="green")
#because there are no discordant genotypes, use this:
colorPanel = colorpanel(2, low="black", high="green")

png(file = output.file)
heatmap.2(concordance.mat,  dendrogram="none", Rowv=FALSE, Colv=FALSE,
				  col=colorPanel, density.info="none", key=FALSE,
				  trace="none", margins = c(5,12),
				  RowSideColors=archiveCol)
legend(0.2,0.95,legend=c("No Genotype","","Concordant Genotype","Discordant Genotype"), col=c("black","white","green","red"), pch=15, ncol=2)
legend("bottomright",legend=c("DNA","Frozen","FFPE"), col=c("darkgreen","orange","cyan"), pch=15)
dev.off()