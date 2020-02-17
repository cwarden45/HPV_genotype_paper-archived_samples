HPV.type = "HPV16"
#HPV.type = "HPV18"
#HPV.type = "HPV58"

#suffix = "_freq5_conservative.txt"#not on GitHub
#suffix = "_freq5.txt"#not on GitHub
#suffix = "_freq5_liberal.txt"#not on GitHub
#suffix = "_freq5_liberal2.txt"#not on GitHub

#suffix = "_freq15_conservative.txt"#not on GitHub
suffix = "_freq15.txt"
#suffix = "_freq15_liberal.txt"#not on GitHub
#suffix = "_freq15_liberal2.txt"#not on GitHub

input.file = paste("Public_Input_Files/",HPV.type,"_seq_stats",suffix,sep="")
output.file = paste("Selected_Output_Files/",HPV.type,"_seq_stats_update_names_MERGE_PAIRS",suffix,sep="")

#NOTE: doesn't use genotype column, so any extended annotation file will work
meta.file = "Selected_Output_Files/combined_genotype_with_year_and_ethnicity_freq15.txt"

seq.table = read.table(input.file, head=T, sep="\t")
mappableID= as.character(seq.table$Sample)
mappableID = gsub("X","",mappableID)

meta.table = read.table(meta.file, head=T, sep = "\t")
meta.table = meta.table[order(meta.table$batch, meta.table$Seq.IGC.ID),]
meta.table$Sample = as.character(meta.table$Sample)
#meta.table$Sample = gsub("-",".",meta.table$Sample)

matchedIDs = mappableID[match(meta.table$Sample, mappableID, nomatch=0)]
unmatched.count = mappableID[-match(matchedIDs,mappableID)]
print(unmatched.count)
ethnicity.matchedIDs = as.character(meta.table$SAMPLEID[match(mappableID,meta.table$Sample)])
pairID = as.character(meta.table$extra.pairID[match(mappableID,meta.table$Sample)])

seq.table$Sample = as.character(seq.table$Sample)
seq.table$Rep.Seq = as.character(seq.table$Rep.Seq)
seq.table$Sample = ethnicity.matchedIDs

write.table(seq.table, output.file, row.names=F, sep="\t", quote=F)

tmp.table =seq.table[is.na(pairID),]
print(dim(tmp.table))

pair.counts = table(pairID)

seq.table = seq.table[!is.na(pairID),]
pairID = pairID[!is.na(pairID)]
tmp.table = rbind(tmp.table, seq.table[match(names(pair.counts[pair.counts == 1]),pairID),])
print(dim(tmp.table))

pair.counts = pair.counts[pair.counts > 1]
for (pair.level in names(pair.counts)){
	pair.table = seq.table[pairID == pair.level,]
	
	pair.total = list()
	pair.count = list()
	
	for (i in 1:nrow(pair.table)){
		rep.seqs = as.character(unlist(strsplit(as.character(pair.table$Rep.Seq[i]),split=",")))
		rep.seqs = gsub(" ","",rep.seqs)
		rep.percent  = as.numeric(as.character(unlist(strsplit(as.character(pair.table$Rep.Percent[i]),split=","))))
		for (j in 1:length(rep.seqs)){
			if (rep.seqs[j] %in% names(pair.total)){
				pair.total[rep.seqs[j]]=as.numeric(pair.total[rep.seqs[j]])+rep.percent[j]
				pair.count[rep.seqs[j]]=as.numeric(pair.count[rep.seqs[j]])+1
			}else{
				pair.total[rep.seqs[j]]=rep.percent[j]
				pair.count[rep.seqs[j]]=1
			}			
		}#end for (j in 1:length(rep.seqs))
	}#end for (i in 1:nrow(pair.table))
	
	pair.seqs = ""
	pair.percent = ""
	
	for (rep.seq in names(pair.total)){
		temp.count = as.numeric(pair.count[rep.seq])
		temp.total = as.numeric(pair.total[rep.seq])
		temp.avg = temp.total / temp.count
		if(pair.seqs == ""){
			pair.seqs = rep.seq
			pair.percent = as.character(temp.avg)
		}else{
			pair.seqs = paste(pair.seqs, rep.seq, sep=", ")
			pair.percent = paste(pair.percent, as.character(temp.avg), sep=", ")
		}
	}#end for (rep.seq in names(pair.total))
	
	pair.table = data.frame(Sample=pair.level, Rep.Seq=pair.seqs,Rep.Percent=pair.percent)
	tmp.table = rbind(tmp.table, pair.table)
}#end for (pair.level in names(pair.counts))

print(dim(tmp.table))

seq.table = tmp.table

seq.count.all = list()
seq.count.majority = list()
seq.count.plurality = list()

for (i in 1:nrow(seq.table)){
	seq.char = as.character(seq.table$Rep.Seq[i])
	count.char = as.character(seq.table$Rep.Percent[i])
	seq.arr = unlist(strsplit(seq.char,split=", "))
	count.arr = unlist(strsplit(count.char,split=", "))
	
	majority.flag = 0
	for (j in 1:length(seq.arr)){
		temp.seq = seq.arr[j]
		temp.count = as.numeric(count.arr[j])
		
		if (temp.seq %in% names(seq.count.majority)){
			seq.count.all[temp.seq]=as.numeric(seq.count.all[temp.seq])+1
		}else{
			seq.count.all[temp.seq]=1
		}
		
		if (j == 1){
			if (temp.seq %in% names(seq.count.majority)){
				seq.count.plurality[temp.seq]=as.numeric(seq.count.plurality[temp.seq])+1
			}else{
				seq.count.plurality[temp.seq]=1
			}			
		}#end if (j == 1)
		
		if (temp.count >= 50){
			majority.flag = 1
			if (temp.seq %in% names(seq.count.majority)){
				seq.count.majority[temp.seq]=as.numeric(seq.count.majority[temp.seq])+1
			}else{
				seq.count.majority[temp.seq]=1
			}		
		}#end if (temp.count >= 50)
	}#end for (j in 1:length(seq.arr))
	
	if(majority.flag == 0){
		print(paste(seq.table$Sample[i]," did not have a representative majority sequence",sep=""))
	}
}#end for (i in 1:nrow(seq.table))

print("Majority")
print(seq.count.majority)

print("Plurality")
print(seq.count.plurality)

print("All")
print(seq.count.all)

print(nrow(seq.table))
