meta.file = "Selected_Output_Files/combined_genotype_with_year_and_ethnicity.txt"

meta.table = read.table(meta.file, head=T, sep="\t")
meta.table = meta.table[!is.na(meta.table$QCarray.call.rate),]

call.rate = meta.table$QCarray.call.rate
LRR.SD = meta.table$QCarray.LRR.SD

color.palette = c("chartreuse4","orange","cyan")
meta.table$batch = as.character(meta.table$batch)
meta.table$batch[meta.table$batch == "161007"] = "DNA"
meta.table$batch[meta.table$batch == "161206"] = "Frozen"
#meta.table$batch[meta.table$batch == "170118"] = "FFPE"
#meta.table$batch = factor(as.character(meta.table$batch), levels=c("DNA","Frozen","FFPE"))
meta.table$batch = factor(as.character(meta.table$batch), levels=c("DNA","Frozen"))

group.levels = levels(meta.table$batch)
labelColors = rep("black",times=length(meta.table$batch))
for (i in 1:length(group.levels)){
	labelColors[meta.table$batch == as.character(group.levels[i])] = color.palette[i]
}#end for (i in 1:length(group.levels))

#Call Rate vs batch
fit = aov(call.rate ~ meta.table$batch)
result = summary(fit)
aov.pvalue = result[[1]][['Pr(>F)']][1]
print(paste("ANOVA p-value = ",signif(aov.pvalue,digits=2),sep=""))

pdf("to_AI/Figure4a.pdf", width=2, height=2, useDingbats=FALSE)
par(mar=c(2.2,2,1.5,1))
plot(meta.table$batch, call.rate,
	cex=0.6, cex.axis=0.5,cex.lab=0.5,cex.main=0.37,
	yaxt="n",ylab="", ylim=c(0,100),
	main="QCarray Call Rate vs. Sample Type",
	col =c("chartreuse4","orange"))
axis(2, mgp=c(1, .5, 0), cex.axis=0.4)
mtext("QCarray Call Rate",2, cex=0.4, padj=-5)
dev.off()


percent.human = as.character(meta.table$human.percent)
percent.human = as.numeric(gsub("%","",percent.human))

#plot LRR-SD and percent human
print("plot LRR-SD and percent human")
		fit = lm(percent.human ~ LRR.SD)
		result = summary(fit)
		pvalue = result$coefficients[2,4]
		print(paste("lm p-value = ",pvalue,sep=""))
		
pdf("to_AI/Figure4c.pdf", width=2, height=2, useDingbats=FALSE)
par(mar=c(2.2,2,1.5,1))
plot(LRR.SD, percent.human,
	col="black", bg=labelColors, xlim=c(0,max(LRR.SD)),
	cex=0.5, cex.axis=0.4,cex.lab=0.4,
	mgp=c(1.2,0.5,0),
	pch=21, xlab="QCArray LRR-SD", ylab="L1 Percent Human")
abline(lm(percent.human~LRR.SD),col="black")
abline(lm(percent.human[meta.table$batch=="DNA"]~LRR.SD[meta.table$batch=="DNA"]),col="chartreuse4", lwd=1)
abline(lm(percent.human[meta.table$batch=="Frozen"]~LRR.SD[meta.table$batch=="Frozen"]),col="orange", lwd=1)
#text(0.08,110, labels=paste("r = ",round(cor(LRR.SD, percent.human, use="pairwise.complete.obs"),digits=2),sep=""), xpd=T, cex=0.6, font=2)
legend(0.4,120, legend = c("DNA","Frozen"), col="black", pt.bg =c("chartreuse4","orange"),
			ncol=2, pch=21,xpd=T, cex=0.4, inset=-0.1)
dev.off()

#plot call rate and percent human
print("plot call rate and percent human")
		fit = lm(percent.human ~ call.rate)
		result = summary(fit)
		pvalue = result$coefficients[2,4]
		print(paste("lm p-value = ",pvalue,sep=""))
		
pdf("to_AI/Figure4b.pdf", width=2, height=2, useDingbats=FALSE)
par(mar=c(2.2,2,1.5,1))
plot(call.rate, percent.human,
	col="black", bg=labelColors, xlim=c(0,100),
	cex=0.5, cex.axis=0.4,cex.lab=0.4,
	mgp=c(1.2,0.5,0),
	pch=21, xlab="QCArray Call Rate", ylab="L1 Percent Human")
abline(lm(percent.human~call.rate),col="black")
abline(lm(percent.human[meta.table$batch=="DNA"]~call.rate[meta.table$batch=="DNA"]),col="chartreuse4", lwd=1)
abline(lm(percent.human[meta.table$batch=="Frozen"]~call.rate[meta.table$batch=="Frozen"]),col="orange", lwd=1)
#text(8,110, labels=paste("r = ",round(cor(call.rate, percent.human, use="pairwise.complete.obs"),digits=2),sep=""), xpd=T, cex=0.6, font=2)
legend(50,120, legend = c("DNA","Frozen"), col="black", pt.bg =c("chartreuse4","orange"),
			ncol=2, pch=21,xpd=T, cex=0.4, inset=-0.1)
dev.off()