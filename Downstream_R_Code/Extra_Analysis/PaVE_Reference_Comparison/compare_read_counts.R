HPV35_genotypes = "../../Public_Input_Files/PE_HPVtype_counts_final_names.txt"
PAVE220_genotypes = "PaVE_220HPV_PE_HPVtype_counts-final_names.txt"
output.file = "read_count_concordance.png"

HPV35.table = read.table(HPV35_genotypes, head=T, sep="\t")
HPV35.table = HPV35.table[HPV35.table$HPV.type != "HUMAN_REF",]
HPV35.table = HPV35.table[HPV35.table$HPV.type != "TOTAL",]

PAVE220.table = read.table(PAVE220_genotypes, head=T, sep="\t")
PAVE220.table = PAVE220.table[PAVE220.table$HPV.type != "HUMAN_REF",]
PAVE220.table = PAVE220.table[PAVE220.table$HPV.type != "TOTAL",]

pairID = c()
HPV35.count = c()
PAVE220.count = c()
plotCol = c()

for (i in 2:ncol(PAVE220.table)){
	for (j in 1:nrow(PAVE220.table)){
		pairID = c(pairID, paste(colnames(PAVE220.table)[i], PAVE220.table$HPV.type[j],sep=" : "))
		
		temp_PAVE220_count = PAVE220.table[j,i]
		if(PAVE220.table$HPV.type[j] %in% HPV35.table$HPV.type){
			temp_HPV35_count = HPV35.table[HPV35.table$HPV.type==PAVE220.table$HPV.type[j], colnames(HPV35.table) == colnames(PAVE220.table)[i]]
		
			HPV35.count = c(HPV35.count, temp_HPV35_count)
			PAVE220.count = c(PAVE220.count, temp_PAVE220_count)
			plotCol = c(plotCol, "black")
		}else{
			temp_HPV35_count = 0
			
			HPV35.count = c(HPV35.count, temp_HPV35_count)
			PAVE220.count = c(PAVE220.count, temp_PAVE220_count)
			plotCol = c(plotCol, "blue")
		}#end else
	}#end for (j in 1:nrow(PAVE220.table))
}#end for (i in 2:ncol(PAVE220.table))

png(output.file)
plot(PAVE220.count, HPV35.count, col=plotCol, pch=16,
		xlab="PaVE 220 HPV + hg38, Adjusted Count",
		ylab="35 HPV + hg38, Adjusted Count",
		main = paste("r = ",signif(cor(PAVE220.count, HPV35.count), digits=6),sep=""))
legend("topleft", legend=c("35 Common HPV","Additional PaVE HPV"), col=c("black","blue"), pch=16)
dev.off()
