set.seed(0)

#NOTE: Using the 15% frequency causes us to lose the 2 lower freqeuncy HPV18 sequences (while keeping the overall trend for HPV16)

HPV16.seqs = c("ACGCAGTACAAATATGTCATTATGTGCTGCCATATCTACTTCAGAAACTACATATAAAAATACTAACTTTAAGGAGTACCTACGACATGGGGAG",
				"ACGCAGTACAAATATGTCATTATGTGCTGCCATATCTACTTCAGAACCTACATATAAAAATACTAACTTTAAAGAGTACCTACGACATGGGGAG",
				"ACGCAGTACAAATATGTCATTATGTGCTGCCATATCTACTTCAGAAACTACATATAAAAATACTAACTTTAAGGAGTACCTACGACATGGGGAA",
				"ACGCAGTACAAATATGTCATTATGTGCTGCCATATCTACTTCAGAACCTACATATAAAAATACTAACTTTAAAGAGTACCTATGACATGGGGAG",
				"ACGCAGTACAAATATGTCATTATGTGCTGCCATATCTACTTCAGAAACTACATATAAAAATACTAACTTTAAAGAGTACCTACGACATGGGGAG")
HPV18.seqs = c("TCGCAGTACCAATTTAACAATATGTGCTTCTACACAGTCTCCTGTACCTGGGCAATATGATGCTACCAAATTTAAGCAGTATAGCAGACATGTTGAG",
				"TCGTAGTACCAATTTAACAATATGTGCTTCTACACAGTCTCCTGTACCTGGGCAATATGATGCTACCAAATTTAAGCAGTATAGCAGACATGTTGAA")
HPV58.seqs = c("TCGTAGCACTAATATGACATTATGCACTGAAGTAACTAAGGAAGGTACATATAAAAATGATAATTTTAAGGAATATGTACGTCATGTTGAA")

HPV16.table = read.table("Public_Input_Files/HPV16_seq_stats_freq15.txt", head=T, sep="\t")
HPV18.table = read.table("Public_Input_Files/HPV18_seq_stats_freq15.txt", head=T, sep="\t")
HPV58.table = read.table("Public_Input_Files/HPV58_seq_stats_freq15.txt", head=T, sep="\t")

meta.table = read.table("Selected_Output_Files/combined_genotype_with_year_and_ethnicity_freq15.txt", head=T, sep="\t")


HPV16.sample = c()
HPV16.rep.seqs = c()
HPV16.freq = c()

for (i in 1:length(HPV16.seqs)){
	temp.table = HPV16.table[grep(HPV16.seqs[i],HPV16.table$Rep.Seq),]
	
	for (j in 1:nrow(temp.table)){
		sampleID = as.character(temp.table$Sample[j])
		seq.text = as.character(temp.table$Rep.Seq[j])
		seq.percent = as.character(temp.table$Rep.Percent[j])
		
		HPV16.sample = c(HPV16.sample, sampleID)
		
		seq.arr = unlist(strsplit(seq.text, split=", "))
		if(length(seq.arr) == 1){
			HPV16.rep.seqs = c(HPV16.rep.seqs, HPV16.seqs[i])
			HPV16.freq = c(HPV16.freq, as.numeric(seq.percent))
		}else{
			percent.arr = unlist(strsplit(seq.percent, split=", "))
			matched.percent = percent.arr[seq.arr == HPV16.seqs[i]]
			if(length(matched.percent) > 1){
				stop("Issue with Matched Percent")
			}else if(length(matched.percent) == 0){
				stop("No matching sequence?")
			}
			
			HPV16.rep.seqs = c(HPV16.rep.seqs, HPV16.seqs[i])
			HPV16.freq = c(HPV16.freq, as.numeric(matched.percent))
		}
	}#end for (j in 1:nrow(temp.table))
}#end for (i in 1:length(HPV16.seqs))

HPV16.rep.seqs = factor(HPV16.rep.seqs, HPV16.seqs)
print(table(HPV16.rep.seqs))
print(tapply(HPV16.freq, HPV16.rep.seqs, mean))
print(tapply(HPV16.freq, HPV16.rep.seqs, sd))

HPV18.sample = c()
HPV18.rep.seqs = c()
HPV18.freq = c()

for (i in 1:length(HPV18.seqs)){
	temp.table = HPV18.table[grep(HPV18.seqs[i],HPV18.table$Rep.Seq),]
	
	for (j in 1:nrow(temp.table)){
		sampleID = as.character(temp.table$Sample[j])
		seq.text = as.character(temp.table$Rep.Seq[j])
		seq.percent = as.character(temp.table$Rep.Percent[j])
		
		HPV18.sample = c(HPV18.sample, sampleID)
		
		seq.arr = unlist(strsplit(seq.text, split=", "))
		if(length(seq.arr) == 1){
			HPV18.rep.seqs = c(HPV18.rep.seqs, HPV18.seqs[i])
			HPV18.freq = c(HPV18.freq, as.numeric(seq.percent))
		}else{
			percent.arr = unlist(strsplit(seq.percent, split=", "))
			matched.percent = percent.arr[seq.arr == HPV18.seqs[i]]
			if(length(matched.percent) > 1){
				stop("Issue with Matched Percent")
			}else if(length(matched.percent) == 0){
				stop("No matching sequence?")
			}
			
			HPV18.rep.seqs = c(HPV18.rep.seqs, HPV18.seqs[i])
			HPV18.freq = c(HPV18.freq, as.numeric(matched.percent))
		}
	}#end for (j in 1:nrow(temp.table))
}#end for (i in 1:length(HPV18.seqs))

HPV18.rep.seqs = factor(HPV18.rep.seqs, HPV18.seqs)
print(table(HPV18.rep.seqs))
print(tapply(HPV18.freq, HPV18.rep.seqs, mean))
print(tapply(HPV18.freq, HPV18.rep.seqs, sd))

pdf("to_AI/Additional_File_05_Supplemental_Figure_S1.pdf", height = 8, width = 12)
par(mfrow=c(2,3))
###########################
### HPV16 - Sample Type ###
###########################

seqCol=rep("gray",length(HPV16.sample))
seqCol[grep("^13",HPV16.sample)]="green"
seqCol[grep("^14",HPV16.sample)]="orange"
seqCol[grep("^15",HPV16.sample)]="cyan"

repSeq= as.character(HPV16.rep.seqs)
for (i in 1:length(HPV16.seqs)){
	repSeq[HPV16.rep.seqs == HPV16.seqs[i]]=paste("Seq",i,sep="")
}#end for (i in 1:length(HPV16.seqs))
repSeq = factor(repSeq, levels=paste("Seq",1:length(HPV16.seqs),sep=""))

jitter.HPV16 = jitter(as.numeric(repSeq), factor=1.5)

par(mar=c(8,5,3,2))
plot(jitter.HPV16, HPV16.freq,
		xaxt="n", xlim=c(0,length(levels(repSeq))+1), ylim=c(0,100), cex=0.8, cex.main=1.4,
		xlab="", ylab="Percent Identical Reads", pch=19, col=seqCol, main = "HPV16 Representative Sequences")
mtext(levels(repSeq), side=1, at =1:length(levels(repSeq)), las=2, line=2)
legend("bottom",legend=c("DNA","Frozen","FFPE"),col=c("green","orange","cyan"),
			cex=1, pch=16, inset=0.025, xpd=T, ncol=3)

###########################
### HPV18 - Sample Type ###
###########################

seqCol=rep("gray",length(HPV18.sample))
seqCol[grep("^13",HPV18.sample)]="green"
seqCol[grep("^14",HPV18.sample)]="orange"
seqCol[grep("^15",HPV18.sample)]="cyan"

repSeq= as.character(HPV18.rep.seqs)
for (i in 1:length(HPV18.seqs)){
	repSeq[HPV18.rep.seqs == HPV18.seqs[i]]=paste("Seq",i,sep="")
}#end for (i in 1:length(HPV18.seqs))
repSeq = factor(repSeq, levels=paste("Seq",1:length(HPV18.seqs),sep=""))

jitter.HPV18 = jitter(as.numeric(repSeq), factor=1)

par(mar=c(8,5,3,2))
plot(jitter.HPV18, HPV18.freq,
		xaxt="n", xlim=c(0,length(levels(repSeq))+1), ylim=c(0,100), cex=0.8, cex.main=1.4,
		xlab="", ylab="Percent Identical Reads", pch=19, col=seqCol, main = "HPV18 Representative Sequences")
mtext(levels(repSeq), side=1, at =1:length(levels(repSeq)), las=2, line=2)
legend("bottom",legend=c("DNA","Frozen","FFPE"),col=c("green","orange","cyan"),
			cex=1, pch=16, inset=0.025, xpd=T, ncol=3)

###########################
### HPV58 - Sample Type ###
###########################

seqCol=rep("gray",nrow(HPV58.table))
seqCol[grep("^13",HPV58.table$Sample)]="green"
seqCol[grep("^14",HPV58.table$Sample)]="orange"
seqCol[grep("^15",HPV58.table$Sample)]="cyan"

repSeq= rep("Seq1",nrow(HPV58.table))
repSeq = factor(repSeq, levels=c("Seq1"))

print(nrow(HPV58.table))
print(mean(HPV58.table$Rep.Percent))
print(sd(HPV58.table$Rep.Percent))

jitter.HPV58 = jitter(as.numeric(repSeq), factor=3)

par(mar=c(8,5,3,2))
plot(jitter.HPV58, HPV58.table$Rep.Percent,
		xaxt="n", xlim=c(0,length(levels(repSeq))+1), ylim=c(0,100), cex=0.8, cex.main=1.4,
		xlab="", ylab="Percent Identical Reads", pch=19, col=seqCol, main = "HPV58 Representative Sequences")
mtext(levels(repSeq), side=1, at =1:length(levels(repSeq)), las=2, line=2)
legend("bottom",legend=c("DNA","Frozen","FFPE"),col=c("green","orange","cyan"),
			cex=1, pch=16, inset=0.025, xpd=T, ncol=3)

##################################
### HPV16 : HPV58 Co-Infection ###
##################################

seqCol=rep("black",length(HPV16.sample))
for(i in 1:length(HPV16.sample)){
	temp.subtype = as.character(meta.table$genotype[meta.table$Sample == HPV16.sample[i]])
	temp.subtype.arr = unlist(strsplit(temp.subtype, split=","))
	if(length(temp.subtype.arr) > 1){
		#print(temp.subtype.arr)
		if("HPV58" %in% temp.subtype.arr){
			seqCol[i]="red"
		}#end if("HPV58" %in% temp.subtype.arr)
	}#end if(length(temp.subtype.arr) > 1)
}#end for(i in 1:length(HPV16.sample))


repSeq= as.character(HPV16.rep.seqs)
for (i in 1:length(HPV16.seqs)){
	repSeq[HPV16.rep.seqs == HPV16.seqs[i]]=paste("Seq",i,sep="")
}#end for (i in 1:length(HPV16.seqs))
repSeq = factor(repSeq, levels=paste("Seq",1:length(HPV16.seqs),sep=""))

par(mar=c(8,5,3,2))
plot(jitter.HPV16, HPV16.freq,
		xaxt="n", xlim=c(0,length(levels(repSeq))+1), ylim=c(0,100), cex=0.8, cex.main=1.4,
		xlab="", ylab="Percent Identical Reads", pch=19, col=seqCol, main = "HPV16 Representative Sequences")
mtext(levels(repSeq), side=1, at =1:length(levels(repSeq)), las=2, line=2)
legend("bottom",legend=c("No HPV58 Co-Infection","HPV58 Co-Infection"),col=c("black","red"),
			cex=1, pch=16, inset=0.025, xpd=T, ncol=2)

##################################
### HPV18 : HPV58 Co-Infection ###
##################################

seqCol=rep("black",length(HPV18.sample))
for(i in 1:length(HPV16.sample)){
	temp.subtype = as.character(meta.table$genotype[meta.table$Sample == HPV16.sample[i]])
	temp.subtype.arr = unlist(strsplit(temp.subtype, split=","))
	if(length(temp.subtype.arr) > 1){
		print(temp.subtype.arr)
		if("HPV58" %in% temp.subtype.arr){
			seqCol[i]="red"
		}#end if("HPV58" %in% temp.subtype.arr)
	}#end if(length(temp.subtype.arr) > 1)
}#end for(i in 1:length(HPV16.sample))

repSeq= as.character(HPV18.rep.seqs)
for (i in 1:length(HPV18.seqs)){
	repSeq[HPV18.rep.seqs == HPV18.seqs[i]]=paste("Seq",i,sep="")
}#end for (i in 1:length(HPV18.seqs))
repSeq = factor(repSeq, levels=paste("Seq",1:length(HPV18.seqs),sep=""))

par(mar=c(8,5,3,2))
plot(jitter.HPV18, HPV18.freq,
		xaxt="n", xlim=c(0,length(levels(repSeq))+1), ylim=c(0,100), cex=0.8, cex.main=1.4,
		xlab="", ylab="Percent Identical Reads", pch=19, col=seqCol, main = "HPV18 Representative Sequences")
mtext(levels(repSeq), side=1, at =1:length(levels(repSeq)), las=2, line=2)
legend("bottom",legend=c("No HPV58 Co-Infection","HPV58 Co-Infection"),col=c("black","red"),
			cex=1, pch=16, inset=0.025, xpd=T, ncol=2)

##################################
### HPV58 : HPV58 Co-Infection ###
##################################
	
seqCol=rep("black",nrow(HPV58.table))
for(i in 1:length(HPV16.sample)){
	temp.subtype = as.character(meta.table$genotype[meta.table$Sample == HPV16.sample[i]])
	temp.subtype.arr = unlist(strsplit(temp.subtype, split=","))
	if(length(temp.subtype.arr) > 1){
		print(temp.subtype.arr)
		if("HPV58" %in% temp.subtype.arr){
			seqCol[i]="red"
		}#end if("HPV58" %in% temp.subtype.arr)
	}#end if(length(temp.subtype.arr) > 1)
}#end for(i in 1:length(HPV16.sample))

repSeq= rep("Seq1",nrow(HPV58.table))
repSeq = factor(repSeq, levels=c("Seq1"))

print(nrow(HPV58.table))
print(mean(HPV58.table$Rep.Percent))
print(sd(HPV58.table$Rep.Percent))

par(mar=c(8,5,3,2))
plot(jitter.HPV58, HPV58.table$Rep.Percent,
		xaxt="n", xlim=c(0,length(levels(repSeq))+1), ylim=c(0,100), cex=0.8, cex.main=1.4,
		xlab="", ylab="Percent Identical Reads", pch=19, col=seqCol, main = "HPV58 Representative Sequences")
mtext(levels(repSeq), side=1, at =1:length(levels(repSeq)), las=2, line=2)
legend("bottom",legend=c("No HPV58 Co-Infection","HPV58 Co-Infection"),col=c("black","red"),
			cex=1, pch=16, inset=0.025, xpd=T, ncol=2)
	
dev.off()
