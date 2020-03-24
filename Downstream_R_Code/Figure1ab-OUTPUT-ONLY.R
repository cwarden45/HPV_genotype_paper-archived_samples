min.call.rate = 85

combined.meta.file="../../Copied_Files/QC_Array_combined_sample_description.txt"
allele.count.table="../../Copied_Files/QC_Array_combined_allele_counts.txt"

#NOTE: These results were removed and back into the paper during iterations of drafts.  So, this is also available on the following page: https://github.com/cwarden45/QCarray_SuperPop/tree/master/1000_Genomes_Benchmarks/pca_plots.R

calc.percent.called = function(arr){
	return(length(arr[!is.na(arr)]))
}#end def calc.percent.called

input.file = paste(allele.count.table,sep="")

allele.mat = read.table(input.file, head=T, sep="\t")
snpIDs = allele.mat[,1]
allele.mat = allele.mat[,2:ncol(allele.mat)]
rownames(allele.mat) = snpIDs
print(dim(allele.mat))

sample.call.rate = 100 * apply(allele.mat, 2, calc.percent.called) / nrow(allele.mat)

allele.mat = allele.mat[,sample.call.rate > min.call.rate]
kept.samples = names(allele.mat)
print(dim(allele.mat))

meta.table = read.table(combined.meta.file, head=T, sep="\t")
meta.table = meta.table[match(kept.samples, meta.table$subject),]

super.pop = c("AFR","AMR","EAS","EUR","SAS")
color.pallete = c("orange","red","green","blue","purple")
labelColors =rep("black",nrow(meta.table))
for (i in 1:length(super.pop)){
	labelColors[meta.table$super.pop == super.pop[i]] = color.pallete[i]
}

allele.mat = allele.mat[,match(meta.table$subject,names(allele.mat),nomatch=0)]
rownames(allele.mat) = snpIDs

print(dim(allele.mat))
temp.mat = na.omit(data.matrix(allele.mat))
print(dim(temp.mat))

pca.values = prcomp(temp.mat)
pc.values = data.frame(pca.values$rotation)
variance.explained = (pca.values $sdev)^2 / sum(pca.values $sdev^2)
pca.table = data.frame(PC = 1:length(variance.explained), percent.variation = variance.explained, t(pc.values))

pdf(paste("to_AI/Figure1a.pdf",sep=""), useDingbats=FALSE)
par(mar=c(5,5,3,3), mfrow=c(2,2))
plot(pc.values$PC1, pc.values$PC2, col = labelColors, xlab = paste("PC1 (",round(100* variance.explained[1] , digits = 2),"%)", sep = ""),
	ylab = paste("PC2 (",round(100* variance.explained[2] , digits = 2),"%)", sep = ""),
	pch=19, cex=0.5)
box(which="figure", xpd=T)

plot(pc.values$PC3, pc.values$PC2, col = labelColors,
	xlab = paste("PC3 (",round(100* variance.explained[3] , digits = 2),"%)", sep = ""),
	ylab = paste("PC2 (",round(100* variance.explained[2] , digits = 2),"%)", sep = ""),
	pch=19, cex=0.5)
box(which="figure", xpd=T)

plot(pc.values$PC1, pc.values$PC3, col = labelColors,
	xlab = paste("PC1 (",round(100* variance.explained[1] , digits = 2),"%)", sep = ""),
	ylab = paste("PC3 (",round(100* variance.explained[3] , digits = 2),"%)", sep = ""),
	pch=19, cex=0.5)
box(which="figure", xpd=T)
		
par(mar=c(0,0,0,0))
plot(NA, xlim=0:1, ylim=0:1, xaxt="n",yaxt="n")
legend("center",legend=c(super.pop,"QC Array"),col=c(color.pallete,"black"),
		pch=19, xpd=T, inset=-0.25)
dev.off()

pdf(paste("to_AI/Figure1b.pdf",sep=""), useDingbats=FALSE)
plot(pc.values$PC4, pc.values$PC3, col = labelColors,
	xlab = paste("PC4 (",round(100* variance.explained[4] , digits = 2),"%)", sep = ""),
	ylab = paste("PC3 (",round(100* variance.explained[3] , digits = 2),"%)", sep = ""),
	pch=19, cex=0.7)
legend("top",legend=c(super.pop,"QC Array"),col=c(color.pallete,"black"),ncol=6, xpd=T, inset=-0.1, pch=19)
dev.off()