meta.file = "Selected_Output_Files/combined_genotype_with_year_and_ethnicity_freq20.txt"

meta.table = read.table(meta.file, head=T, sep="\t")
qPCR = meta.table$qPCR.nM

avg.archivedDNA.neg = mean(qPCR[(meta.table$batch == "161007")&(meta.table$HPV.status == "neg")])

percent.human = as.character(meta.table$human.percent)
percent.human = gsub("%","",percent.human)
percent.human = as.numeric(percent.human)

sample.col = rep(NA,nrow(meta.table))
sample.col[meta.table$batch == "161007"]="darkgreen"
sample.col[meta.table$batch == "161206"]="orange"
sample.col[meta.table$batch == "170118"]="cyan"

calculate.num.coinfections = function(value){
	if(is.na(value)){
		return(0)
	}else{
		genotypes = unlist(strsplit(as.character(value),split=","))
		return(length(genotypes))
	}
}#end def calculate.num.coinfections

num.coinfections = sapply(meta.table$genotype, calculate.num.coinfections)

coinfection.col = rep("black",nrow(meta.table))
coinfection.col[num.coinfections == 0]="gray"
coinfection.col[num.coinfections == 2]="orange"
coinfection.col[num.coinfections == 3]="red"

HPV58.col = rep("black",nrow(meta.table))
HPV58.col[grep("HPV58",meta.table$genotype)]="red"

#confirm that amount of 
	fit = aov(qPCR ~ sample.col)
	result = summary(fit)
	aov.pvalue = result[[1]][['Pr(>F)']][1]

print(paste("qPCR ANOVA p-value = ",aov.pvalue,sep=""))

#plot all
pdf("Figure1.pdf")
#color by sample type
plot(percent.human, qPCR,
	pch=16, col=sample.col,
	xlab="Percent Human", ylab="qPCR (nM)")
abline(h=avg.archivedDNA.neg, col="gray")
rect(-20, -10, 120, 2, col=rgb(red=1, green=0, blue=0, alpha=0.2), border=NA)
legend("top", legend=c("DNA","Frozen","FFPE"), col = c("darkgreen","orange","cyan"),
	xpd=T, inset = -0.1, ncol=3, pch=16)
dev.off()

pdf("Additional_File_05_Supplemental_Figure_S1.pdf")
par(mfrow=c(2,2))	
#color by number of co-infections (all)
plot(percent.human, qPCR,
	pch=16, col=coinfection.col)
abline(h=avg.archivedDNA.neg, col="gray")
rect(-20, -10, 120, 2, col=rgb(red=1, green=0, blue=0, alpha=0.1), border=NA)
legend("top", legend=0:3, col = c("gray","black","orange","red"),
	xpd=T, inset = -0.2, ncol=4, pch=16, cex=0.8)

#color by HPV58 status (all)
plot(percent.human, qPCR,
	pch=16, col=HPV58.col)
abline(h=avg.archivedDNA.neg, col="gray")
rect(-20, -10, 120, 2, col=rgb(red=1, green=0, blue=0, alpha=0.1), border=NA)
legend("top", legend=c("HPV58-","HPV58+"), col = c("black","red"),
	xpd=T, inset = -0.2, ncol=2, pch=16, cex=0.8)

#plot FFPE
FFPE.table = meta.table[meta.table$batch == "170118",]
qPCR = qPCR[meta.table$batch == "170118"]
percent.human = percent.human[meta.table$batch == "170118"]
coinfection.col = coinfection.col[meta.table$batch == "170118"]
HPV58.col = HPV58.col[meta.table$batch == "170118"]

#color by number of co-infections (FFPE)
plot(percent.human, qPCR,
	pch=16, col=coinfection.col)
abline(h=avg.archivedDNA.neg, col="gray")
rect(-20, -10, 120, 2, col=rgb(red=1, green=0, blue=0, alpha=0.1), border=NA)
legend("top", legend=0:3, col = c("gray","black","orange","red"),
	xpd=T, inset = -0.2, ncol=4, pch=16, cex=0.8)
	
#color by HPV58 status (FFPE)
plot(percent.human, qPCR,
	pch=16, col=HPV58.col)
abline(h=avg.archivedDNA.neg, col="gray")
rect(-20, -10, 120, 2, col=rgb(red=1, green=0, blue=0, alpha=0.1), border=NA)
legend("top", legend=c("HPV58-","HPV58+"), col = c("black","red"),
	xpd=T, inset = -0.2, ncol=2, pch=16, cex=0.8)
dev.off()