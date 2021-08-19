allele.mat = read.table("../../../../Copied_Files/QC_Array_combined_allele_counts-IBD.txt", head=T, sep="\t")
plink.table = read.table("../../../../Copied_Files/pruned.genome", head=T)
meta.1KG.table = read.table("../../../../Copied_Files/QC_Array_1000_Genomes_sample_description.txt", head=T, sep="\t")

meta.QCarray.table = read.table("../../../../Copied_Files/QC_Array_combined_sample_description.txt", head=T, sep="\t")
colnames(meta.QCarray.table) = c("subject","population","super.population")
meta.QCarray.table = meta.QCarray.table[is.na(meta.QCarray.table$population),]
meta.QCarray.table$subject = as.character(meta.QCarray.table$subject)
meta.QCarray.table$subject[meta.QCarray.table$subject == "S51540.T"] = "S51450.T" #fix typo

plink.table$IID1 = as.character(plink.table$IID1)
plink.table$IID1[plink.table$IID1 == "S51540.T"] = "S51450.T" #fix typo
plink.table$IID2 = as.character(plink.table$IID2)
plink.table$IID2[plink.table$IID2 == "S51540.T"] = "S51450.T" #fix typo
plink.pairID = paste(plink.table$IID1,plink.table$IID2,sep="-")

sampleIDs = colnames(allele.mat)
sampleIDs[sampleIDs == "S51540.T"] = "S51450.T" #fix typo
colnames(allele.mat) = sampleIDs

k1g.parent.PI_HAT = c()
k1g.child.PI_HAT = c()
QCarray.PI_HAT = c()
k1g.parent.IBD.dist = c()
k1g.child.IBD.dist = c()
QCarray.IBD.dist = c()

#S51423.T
combinedIDs = c(as.character(meta.1KG.table$subject),as.character(meta.QCarray.table$subject))

ref.alleles = allele.mat[,match(combinedIDs,names(allele.mat))]

familyID = levels(meta.1KG.table$family)

expected.father = c()
expected.mother = c()

k1g.parent.dist = c()
k1g.child.dist = c()

for (i in 1:length(familyID)){
	family.mat = meta.1KG.table[meta.1KG.table$family == familyID[i],]
	
	motherID = as.character(family.mat$subject[family.mat$relationship.1kg == "mother"])
	fatherID = as.character(family.mat$subject[family.mat$relationship.1kg == "father"])
	childID = as.character(family.mat$subject[family.mat$relationship.1kg == "child"])
	
	mother.alleles = as.numeric(ref.alleles[,motherID])
	father.alleles = as.numeric(ref.alleles[,fatherID])
	child.alleles = as.numeric(ref.alleles[,childID])

	parent.dist = dist(t(data.frame(mother.alleles, father.alleles)))
	child.dist1 = dist(t(data.frame(child.alleles, father.alleles)))
	child.dist2 = dist(t(data.frame(child.alleles, mother.alleles)))
	
	k1g.parent.dist[i] = parent.dist
	k1g.child.dist[i] = (child.dist1 + child.dist2)/2
	
	if(paste(motherID,fatherID,sep="-") %in% plink.pairID){
		parent.plink.stats = plink.table[plink.pairID == paste(motherID,fatherID,sep="-"),]
	}else if(paste(fatherID,motherID,sep="-") %in% plink.pairID){
		parent.plink.stats = plink.table[plink.pairID == paste(fatherID,motherID,sep="-"),]
	}else{
		stop(paste("Problem parsing ",motherID,"/",fatherID,sep=""))
	}

	if(paste(motherID,childID,sep="-") %in% plink.pairID){
		MtoC.plink.stats = plink.table[plink.pairID == paste(motherID,childID,sep="-"),]
	}else if(paste(childID,motherID,sep="-") %in% plink.pairID){
		MtoC.plink.stats = plink.table[plink.pairID == paste(childID,motherID,sep="-"),]
	}else{
		stop(paste("Problem parsing ",motherID,"/",childID,sep=""))
	}

	if(paste(fatherID,childID,sep="-") %in% plink.pairID){
		PtoC.plink.stats = plink.table[plink.pairID == paste(fatherID,childID,sep="-"),]
	}else if(paste(childID,fatherID,sep="-") %in% plink.pairID){
		PtoC.plink.stats = plink.table[plink.pairID == paste(childID,fatherID,sep="-"),]
	}else{
		stop(paste("Problem parsing ",fatherID,"/",childID,sep=""))
	}
	
	k1g.parent.PI_HAT[i] = parent.plink.stats$PI_HAT
	k1g.parent.IBD.dist[i] = parent.plink.stats$DST

	k1g.child.PI_HAT[i] = (MtoC.plink.stats$PI_HAT + PtoC.plink.stats$PI_HAT)/2
	k1g.child.IBD.dist[i] = (MtoC.plink.stats$DST + PtoC.plink.stats$DST)/2
	
}#end for (i in 1:length(familyID))

QCarraySamples = as.character(meta.QCarray.table$subject)

QCarrayPair = c()
QCarrayDist = c()

for (i in 1:(length(QCarraySamples)-1)){
	nextIndex = i + 1
	for (j in nextIndex:length(QCarraySamples)){
		sampleID1 = QCarraySamples[i]
		sampleID2 = QCarraySamples[j]
		
		alleles1 = ref.alleles[,sampleID1]
		alleles2 = ref.alleles[,sampleID2]
		QCarrayDist[length(QCarrayPair)+1] = dist(t(data.frame(alleles1, alleles2)))

		if(paste(sampleID1,sampleID2,sep="-") %in% plink.pairID){
			QCarray.plink.stats = plink.table[plink.pairID == paste(sampleID1,sampleID2,sep="-"),]
		}else if(paste(sampleID2,sampleID1,sep="-") %in% plink.pairID){
			QCarray.plink.stats = plink.table[plink.pairID == paste(sampleID2,sampleID1,sep="-"),]
		}else{
			stop(paste("Problem parsing ",sampleID1,"/",sampleID2,sep=""))
		}

		
		QCarray.PI_HAT[length(QCarrayPair)+1] = QCarray.plink.stats$PI_HAT
		QCarray.IBD.dist[length(QCarrayPair)+1] = QCarray.plink.stats$DST

		QCarrayPair[length(QCarrayPair)+1]=paste(sampleID1,sampleID2,sep="-")		
	}#end for (j in nextIndex:length(cohSamples))
}#end for (i in 1:length(cohSamples))

QCarray.status = rep("QCarray",length(QCarrayPair))
QCarray.status[QCarrayPair == "S51468.T-S51468.N"]="QCarray:tumor-normal pair"
QCarray.status[QCarrayPair == "S51468.N-S51468.T"]="QCarray:tumor-normal pair"

QCarray.status[QCarrayPair == "S51490.T-S51491.N"]="QCarray:tumor-normal pair"
QCarray.status[QCarrayPair == "S51491.N-S51490.T"]="QCarray:tumor-normal pair"

QCarray.status[QCarrayPair == "S51466.T-S51466.N"]="QCarray:tumor-normal pair"
QCarray.status[QCarrayPair == "S51466.N-S51466.T"]="QCarray:tumor-normal pair"

QCarray.status[QCarrayPair == "S51480.T-S51480.N"]="QCarray:tumor-normal pair"
QCarray.status[QCarrayPair == "S51480.N-S51480.T"]="QCarray:tumor-normal pair"

print(table(QCarray.status))
combinedID = c(paste(familyID,"parents",sep="."),
				paste(familyID,"child",sep="."),
				QCarrayPair)
combined.dist = c(k1g.parent.dist, k1g.child.dist, QCarrayDist)
plink.IBD = c(k1g.parent.PI_HAT, k1g.child.PI_HAT, QCarray.PI_HAT)
plink.dist = c(k1g.parent.IBD.dist, k1g.child.IBD.dist, QCarray.IBD.dist)
relationship.type = c(rep("1KG.parent.to.parent",length(familyID)),
						rep("1KG.child.to.parent.avg",length(familyID)),
						QCarray.status)

#need to find bug in code, but plink statistics match earlier allele distance
#output.table = data.frame(PairID = combinedID, type=relationship.type,
#							plink.IBD.PI_HAT=plink.IBD, plink.IBD.dist=plink.dist, allele.count.dist = combined.dist)
						
output.table = data.frame(PairID = combinedID, type=relationship.type,
							plink.IBD.PI_HAT=plink.IBD)
write.table(output.table,"PI_HAT-pairwise_values.txt", sep="\t", quote=F, row.names=F)

png("PI_HAT-distributions.png")
par(mar=c(13,5,2,5))
plot(output.table$type, output.table$plink.IBD.PI_HAT,
	ylab="PI_HAT",col="gray", las=2)
dev.off()

#pdf("FigureS3b.pdf")
#par(mar=c(12,5,2,5))
#plot(output.table$type, output.table$plink.IBD.dist,
#	ylab="IBD Distance", col="gray", las=2)
#dev.off()

#pdf("FigureS3c.pdf")
#par(mar=c(12,5,2,5))
#plot(output.table$type, output.table$alleleDist,
#	ylab="Allele Count Distance", col="gray", las=2)
#dev.off()

