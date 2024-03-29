###Bowtie1 Alignment Similar to Table S4 (> 1.2x  Human Overall, > 1.0x Human Specific Genotype) --> 20% frequency
qPCR_flag = FALSE
input.file = "Bowtie1-hg38_plus_35HPV_genotype_calls_freq20.txt"
output.file = "Bowtie1-combined_genotype_with_year_and_ethnicity_freq20.txt"

###Bowtie1 Alignment Similar to Table S3 + S4 (> 1.2x  Human Overall, > 1.0x Human Specific Genotype) --> 5% frequency
#qPCR_flag = FALSE
#input.file = "Bowtie1-hg38_plus_35HPV_genotype_calls_freq5.txt"
#output.file = "Bowtie1-combined_genotype_with_year_and_ethnicity_freq5.txt"

get.igcid = function(char.name){
	sample.info = unlist(strsplit(as.character(char.name), split="_"))
	return(sample.info[1])
}#end def get.igcid

input.table = read.table(input.file,head=T, sep="\t")

batchID = rep(NA,nrow(input.table))
batchID[grep("^13",as.character(input.table$Sample))]="161007"
batchID[grep("^14",as.character(input.table$Sample))]="161206"
batchID[grep("^15",as.character(input.table$Sample))]="170118"

igcID = sapply(input.table$Sample, get.igcid)
pathID = as.character(input.table$Sample)
pathID = gsub("^\\d+_","",pathID, perl=T)
pathID[batchID == "161007"]=gsub("^\\d+_","",pathID[batchID == "161007"])

#add qPCR values
qPCR.nM = rep(NA, length(igcID))

FileMaker_files = c("FileMaker_ArchivedDNA_qPCR.tab","FileMaker_Frozen_qPCR.tab","FileMaker_FFPE_qPCR.tab")

for (temp.inputfile in FileMaker_files){
	temp.table = read.table(paste("../../../../Copied_Files/",temp.inputfile,sep=""), head=F, sep="\t")
	print(dim(temp.table))
	if("15969" %in% temp.table$V2){
		#remove sample that we withdrew
		temp.table = temp.table[temp.table$V2 != "15969",]
		print(dim(temp.table))
	}#end if("15969" %in% temp.table$V2)
	
	qPCR.nM[match(temp.table$V2, igcID)]=temp.table$V12
}#end for (temp.inputfile in FileMaker_files)

if(qPCR_flag){
	#flag samples with qPCR amplification value less than 2 nM, instead of providing genotype
	input.table$HPV.status = as.character(input.table$HPV.status)
	input.table$HPV.status[(qPCR.nM < 2)&(input.table$HPV.status != "neg")]="qPCR Flag"
	input.table$HPV.percent = as.character(input.table$HPV.percent)
	input.table$HPV.percent[(qPCR.nM < 2)&(input.table$HPV.status != "neg")]="qPCR Flag"
	
	input.table$genotype = as.character(input.table$genotype)
	input.table$genotype[(qPCR.nM < 2)&(input.table$HPV.status != "neg")]="qPCR Flag"
	input.table$genotype.percent = as.character(input.table$genotype.percent)
	input.table$genotype.percent[(qPCR.nM < 2)&(input.table$HPV.status != "neg")]="qPCR Flag"
}#end if(qPCR_flag)


meta.table = data.frame(Seq.IGC.ID=igcID, seq.pathID=pathID,
						Sample = as.character(input.table$Sample), qPCR.nM, human.percent = as.character(input.table$human.percent),
						input.table[,c(2:3,5:6)],
						batch=batchID)

meta.table$seq.pathID = as.character(meta.table$seq.pathID)
meta.table$seq.pathID[meta.table$seq.pathID == "16142-01-1"]="16142-01-01"
meta.table$seq.pathID[meta.table$seq.pathID == "16142-01-2"]="16142-01-02"
meta.table$seq.pathID[meta.table$seq.pathID == "16142-01-3"]="16142-01-03"
meta.table$seq.pathID[meta.table$seq.pathID == "16142-01-4"]="16142-01-04"
meta.table$seq.pathID[meta.table$seq.pathID == "16142-01-5"]="16142-01-05"
meta.table$seq.pathID[meta.table$seq.pathID == "16142-01-6"]="16142-01-06"
meta.table$seq.pathID[meta.table$seq.pathID == "16142-01-7"]="16142-01-07"

SAMPLEID=meta.table$seq.pathID
SAMPLEID[grep("^\\d",SAMPLEID)]=paste("S",SAMPLEID[grep("^\\d",SAMPLEID)],sep="")
SAMPLEID = gsub("-",".",SAMPLEID)
meta.table = data.frame(SAMPLEID, meta.table)

geno.table = table(meta.table$genotype, meta.table$batch)
print(geno.table)

sum.types = function(arr, coinfections, type){
	arr = arr[grep(type,coinfections)]
	return(sum(arr))
}

HPV16.counts = apply(geno.table, 2, sum.types, coinfections = rownames(geno.table), type = "HPV16")
HPV18.counts = apply(geno.table, 2, sum.types, coinfections = rownames(geno.table), type = "HPV18")
HPV58.counts = apply(geno.table, 2, sum.types, coinfections = rownames(geno.table), type = "HPV58")
print(HPV18.counts)
print(HPV58.counts)
print(HPV58.counts/HPV18.counts)

print(paste("Total HPV16+: ",sum(HPV16.counts),sep=""))
print(paste("Total HPV18+: ",sum(HPV18.counts),sep=""))
print(paste("Total HPV58+: ",sum(HPV58.counts),sep=""))

count.coinfections = function(char){
	char.arr = unlist(strsplit(char,split=","))
	return(length(char.arr))
}
print(max(sapply(as.character(rownames(geno.table)), count.coinfections)))
coinfections.per.sample = sapply(as.character(meta.table$genotype), count.coinfections)

#Pathology Department Archived DNA
temp.table= read.table("../../../../Copied_Files/HPV_cases_Ogembo_Sharon.txt", head=T,sep="\t")
temp.table$Sample.name = as.character(temp.table$Sample.name)
temp.table$Sample.name = gsub("-",".",temp.table$Sample.name)

sample.type = as.character(temp.table$Tumor[match(SAMPLEID,temp.table$Sample.name)])
sample.type = gsub("\\s","",sample.type)
sample.type[grep(".N$",meta.table$SAMPLEID)]="adjacent normal"

prev.HPV.status = as.character(temp.table$HPV.PCR[match(SAMPLEID,temp.table$Sample.name)])
prev.HPV.status = gsub("\\s","",prev.HPV.status)

meta.table = data.frame(meta.table, sample.type, prev.HPV.status)

prev.HPV.table = table(paste(meta.table$sample.type,meta.table$HPV.status,meta.table$genotype,sep=" "), prev.HPV.status)
print(prev.HPV.table)
print(sum(prev.HPV.table))

#Pathology Core Frozen
temp.table= read.table("../../../../Copied_Files/16142LOI01frozenrelease.txt", head=T)
pairID = temp.table[,1]
pathID = as.character(temp.table[,2])
collection.year = temp.table[,3]

path.table = data.frame(pathID=pathID, collection.year=collection.year, pairID=pairID)

#Pathology Core FFPE
temp.table = read.table("../../../../Copied_Files/16142LOI01paraffin.txt", head=T)
pathID = as.character(rownames(temp.table))
collection.year = temp.table[,1]
pairID = rep(NA, length(pathID))

temp.table = data.frame(pathID=pathID, collection.year=collection.year, pairID=pairID)
path.table = rbind(path.table, temp.table)

meta.table = cbind(meta.table, path.table[match(meta.table$seq.pathID, path.table$pathID),])

#add histological subtype information
temp.table = read.table("../../../../Copied_Files/Histology_Cnext_pathologyreports-REFORMAT.txt", head=T, sep="\t")

temp.table$DEID = as.character(temp.table$DEID)
temp.table$DEID = gsub("-\\(JOG-\\d+\\)","",temp.table$DEID)

normal.rows = temp.table[grep("-N-T",temp.table$DEID),]
temp.table$DEID = gsub("-N-T","-T",temp.table$DEID)
normal.rows$DEID = gsub("-N-T","-N",normal.rows$DEID)
temp.table = rbind(temp.table, normal.rows)

hist.table = temp.table[match(meta.table$seq.pathID, temp.table$DEID),]
print(meta.table$seq.pathID[is.na(hist.table$Histological.Subtype)])

hist.table$Histological.Subtype = as.character(hist.table$Histological.Subtype)
hist.table$Histological.Subtype[meta.table$sample.type == "adjacent normal"] = "Adjacent Normal Cervix"

meta.table$sample.type = as.character(meta.table$sample.type)
meta.table$sample.type[hist.table$Histological.Subtype == "Adenocarcinoma (Adeno)"] = "Invasive Cervical Cancer"
meta.table$sample.type[hist.table$Histological.Subtype == "Adenosquamous"] = "Invasive Cervical Cancer"
meta.table$sample.type[hist.table$Histological.Subtype == "Cervical Small Cell Carcinoma"] = "Invasive Cervical Cancer"
meta.table$sample.type[hist.table$Histological.Subtype == "Squamous Cell Carcinoma (SCC)"] = "Invasive Cervical Cancer"
meta.table$sample.type[hist.table$Histological.Subtype == "Unavailable"] = "Invasive Cervical Cancer"


meta.table$sample.type[meta.table$seq.pathID == "16142-01-55"] = "Non-Malignant Vagina"
meta.table$sample.type[meta.table$sample.type == "nlprostate"] = "Normal Prostate"
meta.table$sample.type[meta.table$sample.type == "adjacent normal"] = "Adjacent Normal Cervix"

meta.table$sample.type[meta.table$sample.type == "cxca"] = "Invasive Cervical Cancer"
meta.table$sample.type[meta.table$sample.type == "vulvaca"] = "Vulvar Cancer"
meta.table$sample.type[hist.table$Histological.Subtype == "Endometroid + Serous Papillary Carcinoma"] = "Endometrial + Cervical Cancer"

meta.table$sample.type[meta.table$sample.type == "prostate"] = "Prostate Cancer"
meta.table$sample.type = as.factor(meta.table$sample.type)

meta.table = data.frame(meta.table, hist.subtype=hist.table$Histological.Subtype)

print(table(meta.table$hist.subtype))
overall.hist.subtype = table(meta.table$batch, meta.table$hist.subtype)
print(overall.hist.subtype)

hist.FE.mat = matrix(ncol=2,nrow=3)
colnames(hist.FE.mat) = c("Adeno","SSC")
rownames(hist.FE.mat) = c("DNA","Frozen","FFPE")

#hist.FE.mat[1,1]=overall.hist.subtype[1,1]
#hist.FE.mat[2,1]=overall.hist.subtype[2,1]
#hist.FE.mat[3,1]=overall.hist.subtype[3,1]

hist.FE.mat[1,1]=overall.hist.subtype[1,1] + overall.hist.subtype[1,2]
hist.FE.mat[2,1]=overall.hist.subtype[2,1] + overall.hist.subtype[2,2]
hist.FE.mat[3,1]=overall.hist.subtype[3,1] + overall.hist.subtype[3,2]

hist.FE.mat[1,2]=overall.hist.subtype[1,8]
hist.FE.mat[2,2]=overall.hist.subtype[2,8]
hist.FE.mat[3,2]=overall.hist.subtype[3,8]

print(hist.FE.mat)
print(fisher.test(hist.FE.mat))

print(hist.FE.mat[c(1,3),])
print(fisher.test(hist.FE.mat[c(1,3),]))

#print(hist.FE.mat[c(2,3),])
#print(fisher.test(hist.FE.mat[c(2,3),]))

#print(hist.FE.mat[c(1,2),])
#print(fisher.test(hist.FE.mat[c(1,2),]))

temp.HPV58.table = meta.table[grep("HPV58",meta.table$genotype),]
HPV58.hist.subtype = table(temp.HPV58.table$batch,temp.HPV58.table$hist.subtype)

FFPE.HPV58.hist.FE.mat = matrix(ncol=2,nrow=2)
colnames(FFPE.HPV58.hist.FE.mat) = c("Adeno","SSC")
rownames(FFPE.HPV58.hist.FE.mat) = c("Positive","Negative")

FFPE.HPV58.hist.FE.mat[1,1]=HPV58.hist.subtype[3,1]
FFPE.HPV58.hist.FE.mat[1,2]=HPV58.hist.subtype[3,8]

FFPE.HPV58.hist.FE.mat[2,1]=overall.hist.subtype[3,1] - HPV58.hist.subtype[3,1]
FFPE.HPV58.hist.FE.mat[2,2]=overall.hist.subtype[3,8] - HPV58.hist.subtype[3,8]

print(FFPE.HPV58.hist.FE.mat)
print(fisher.test(FFPE.HPV58.hist.FE.mat))

#parse collection year information
print(mean(meta.table$collection.year, na.rm=T))
print(sd(meta.table$collection.year, na.rm=T))

print(tapply(meta.table$collection.year, meta.table$batch, mean, na.rm=T))
print(tapply(meta.table$collection.year, meta.table$batch, sd, na.rm=T))

#collection date density plot
Frozen.years = meta.table$collection.year[meta.table$batch == "161206"]
Frozen.years = Frozen.years[!is.na(Frozen.years)]
print(paste("Frozen (1990-2010): ",length(Frozen.years[(Frozen.years >= 1990)&(Frozen.years <=2010)]),
						" / ",length(Frozen.years),sep=""))
FFPE.years = meta.table$collection.year[meta.table$batch == "170118"]
FFPE.years = FFPE.years[!is.na(FFPE.years)]
print(paste("FFPE (1990-2010): ",length(FFPE.years[(FFPE.years >= 1990)&(FFPE.years <=2010)]),
						" / ",length(FFPE.years),sep=""))
print(ks.test(Frozen.years, FFPE.years))
print(wilcox.test(Frozen.years, FFPE.years))
						
collection.start = 1980
collection.stop = 2020

temp.years = meta.table$collection.year
names(temp.years) = meta.table$genotype
temp.batch = meta.table$batch
temp.batch = temp.batch[!is.na(temp.years)]
temp.years = temp.years[!is.na(temp.years)]

HPV58.status = rep("Negative", length(temp.years))
HPV58.status[grep("HPV58",names(temp.years))]="Positive"
HPV58.status = factor(HPV58.status, levels=c("Positive","Negative"))

year.bin = rep("gt2000",length(temp.years))
year.bin[temp.years < 2000]="lt2000"
year.bin = factor(year.bin, levels = c("lt2000","gt2000"))
print(table(HPV58.status,year.bin))
print(table(HPV58.status[temp.batch == 161206],year.bin[temp.batch == 161206]))
print(table(HPV58.status[temp.batch == 170118],year.bin[temp.batch == 170118]))

		fit = lm(coinfections.per.sample ~ meta.table$collection.year)
		result = summary(fit)
		pvalue = result$coefficients[2,4]
print(paste("Co-Infection vs Year (1-var) p-value: ",pvalue,sep=""))

		fit = lm(coinfections.per.sample ~ meta.table$collection.year + meta.table$batch)
		result = summary(fit)
		pvalue = result$coefficients[2,4]
print(paste("Co-Infection vs Year (2-var) p-value: ",pvalue,sep=""))

#extra meta.data (Pathology Core)
extra.pairs = read.table("../../../../Copied_Files/Yafan_same_samples.txt", head=T, sep="\t")
extra.pairs = extra.pairs[extra.pairs$Pair1 != "16142-01-22",]
extra.pairs$Pair2 = as.character(extra.pairs$Pair2)
#correct label is 51450-T --> do not change back to match minor typo (51540-T)

extra.pairID = rep(NA, nrow(meta.table))
batch2.TNcounts = table(meta.table$pairID)
batch2.TNcounts=batch2.TNcounts[batch2.TNcounts > 1]

pair.count = 0

for(pairID in names(batch2.TNcounts)){
	pair.count = pair.count + 1
	extra.pairID[meta.table$pairID == pairID]=paste("pair",pair.count,sep="")
}

for (i in 1:nrow(extra.pairs)){
	paired.samples = extra.pairs[i,]
	paired.samples = as.character(paired.samples[!is.na(paired.samples)])
	paired.samples[grep("^\\d",paired.samples)]=paste("S",paired.samples,sep="")
	paired.samples = gsub("-",".",paired.samples)
	
	prev.label = extra.pairID[match(paired.samples,meta.table$SAMPLEID)]
	prev.label = prev.label[!is.na(prev.label)]
	if(length(prev.label) > 0){
		newID=prev.label[1]
	}else{
		pair.count = pair.count + 1
		newID = paste("pair",pair.count,sep="")
	}
	for (j in 1:length(paired.samples)){
		if(paired.samples[j] %in% meta.table$SAMPLEID){
			extra.pairID[meta.table$SAMPLEID == paired.samples[j]]=newID
		}else{
			stop(paste("Need to map ",paired.samples[j],sep=""))
		}
	}#end for (j in 1:length(paired.samples))
}#end for (i in 1:nrow(extra.pairs))

QCarray.table = read.table("../../../../Copied_Files/IBD_Table-Revised_Names.txt", head=T, sep="\t")
QCarray.table = QCarray.table[(QCarray.table$type == "QCarray")&(QCarray.table$plink.IBD.PI_HAT > 0.9),]

QCarray.count = 0
double.pair.count = 0

for(i in 1:nrow(QCarray.table)){
	QCarray.pairs = unlist(strsplit(as.character(QCarray.table$PairID[i]),split="-"))
	if(QCarray.pairs[1] %in% meta.table$SAMPLEID){
		if(QCarray.pairs[2] %in% meta.table$SAMPLEID){
			if(is.na(extra.pairID[meta.table$SAMPLEID == QCarray.pairs[1]])&is.na(extra.pairID[meta.table$SAMPLEID == QCarray.pairs[2]])){
				QCarray.count = QCarray.count + 1
				extra.pairID[match(QCarray.pairs,meta.table$SAMPLEID)]=paste("QCarray.pair",QCarray.count,sep="")
			}else{
				double.pair.count = double.pair.count + 1
				prevID1 = extra.pairID[meta.table$SAMPLEID == QCarray.pairs[1]]
				prevID2 = extra.pairID[meta.table$SAMPLEID == QCarray.pairs[2]]
				extra.pairID[match(QCarray.pairs,meta.table$SAMPLEID)]=paste("partialQCarray.pair",double.pair.count,sep="")
				if(!is.na(prevID1)){
					extra.pairID[extra.pairID == prevID1]=paste("partialQCarray.pair",double.pair.count,sep="")
				}
				if(!is.na(prevID2)){
					extra.pairID[extra.pairID == prevID2]=paste("partialQCarray.pair",double.pair.count,sep="")
				}
			}
		}else{
			stop(paste("Need to revise code to map QCarray pair sample: ",QCarray.pairs[2]))
		}
	}else{
		stop(paste("Need to revise code to map QCarray pair sample: ",QCarray.pairs[1]))
	}
}#end for(i in 1:nrow(QCarray.table$PairID))

#add QC Array call rate
call.rate.table = read.table("../../../../Copied_Files/R_call_rate.txt", head=T, sep="\t")
call.rate.table$SampleID = as.character(call.rate.table$SampleID)
call.rate.table$SampleID[call.rate.table$SampleID == "S51540.T"] = "S51450.T" #minor typo

QCarray.call.rate = as.character(call.rate.table$overall.call.rate[match(meta.table$SAMPLEID,call.rate.table$SampleID)])

meta.table = data.frame(meta.table, extra.pairID, QCarray.call.rate)

#add age
age.table = read.table("../../../../Copied_Files/OgemboStudy_DemographicV2_DEID.txt", head=T, sep="\t")
age.table$PathID = as.character(age.table$PathID)
age.table$PathID = paste("S",age.table$PathID,sep="")
#correct label is 51450-T --> do not change back to match minor typo (51540-T)
#S16142.01.22 was not processed
print(dim(age.table))
#print(age.table$DEID[-match(meta.table$seq.pathID, age.table$DEID, nomatch=0)])
#print(age.table$PathID[-match(as.character(meta.table$SAMPLEID), as.character(age.table$PathID), nomatch=0)])
mapped.age.table = age.table[match(meta.table$SAMPLEID, age.table$PathID),3:5]

#table also includes collection date and race
#when present, collection dates always matched (and this table less complete than earlier table)
#there were some discrepancies between race, typically "Other" or "Missing" (but I'll check that in later part of code)

mapped.age.table$Age[mapped.age.table$Age == "Unknown"]=NA
mapped.age.table$Age = as.numeric(as.character(mapped.age.table$Age))

meta.table = data.frame(meta.table, Age = mapped.age.table$Age)

print(mean(meta.table$Age, na.rm=T))
print(sd(meta.table$Age, na.rm=T))
print(tapply(meta.table$Age, meta.table$batch, mean, na.rm=T))
print(tapply(meta.table$Age, meta.table$batch, sd, na.rm=T))

		fit = lm(coinfections.per.sample ~ meta.table$Age)
		result = summary(fit)
		pvalue = result$coefficients[2,4]
print(paste("Co-Infection vs Age (1-var) p-value: ",pvalue,sep=""))

		fit = lm(coinfections.per.sample ~ meta.table$Age + meta.table$batch)
		result = summary(fit)
		pvalue = result$coefficients[2,4]
print(paste("Co-Infection vs Age (2-var) p-value: ",pvalue,sep=""))

#L1 Amplicon-Seq Insert Size

L1.median.human.insert.size = c()
L1.max.human.insert.size = c()

for (i in 1:nrow(meta.table)){
	#print(meta.table$Sample[i])
	
	input.file = paste("../../../../Copied_Files/Human_Only_Alignment/",meta.table$Sample[i],"/insert_size_metrics.txt",sep="")
	if (file.exists(input.file)){
		#based upon example from https://stackoverflow.com/questions/39110755/skip-specific-rows-using-read-csv-in-r
		input.table = read.delim(input.file,head=T, sep="\t", skip=6, nrows = 1, as.is = T)
		
		L1.median.human.insert.size[i]=input.table$MEDIAN_INSERT_SIZE
		L1.max.human.insert.size[i]=input.table$MAX_INSERT_SIZE
	}else{
		L1.median.human.insert.size[i]=NA
		L1.max.human.insert.size[i]=NA
	}
}#end for (i in 1:nrow(meta.table))

meta.table = data.frame(meta.table,
						L1.median.human.insert.size, L1.max.human.insert.size)

#ethnicity information - ADMIXTURE
array.table = read.table("../../../../Copied_Files/ADMIXTURE_super_pop_assignments.txt", head=T, sep="\t")					
array.table$sample = as.character(array.table$sample)
array.table$sample[array.table$sample == "S51540.T"] = "S51450.T" #minor typo

ADMIXTURE.top = array.table$max.assignment[match(meta.table$SAMPLEID,array.table$sample)]
ADMIXTURE.mixed = array.table$mixed.assignment[match(meta.table$SAMPLEID,array.table$sample)]
meta.table = data.frame(meta.table,ADMIXTURE.top,ADMIXTURE.mixed)

print(table(ADMIXTURE.mixed, meta.table$batch))
print(table(ADMIXTURE.mixed))

#ethnicity information - bootstrap
array.table = read.table("../../../../Copied_Files/overall_distance_bootstrap_assignments.txt", head=T, sep="\t")					
array.table$sample = as.character(array.table$sample)
array.table$sample[array.table$sample == "S51540.T"] = "S51450.T" #minor typo

bootstrap.ethnicity = array.table$min.overall.dist[match(meta.table$SAMPLEID,array.table$sample)]
meta.table = data.frame(meta.table,bootstrap.ethnicity)

print(table(bootstrap.ethnicity, ADMIXTURE.top))
print(table(bootstrap.ethnicity, ADMIXTURE.mixed))

ADMIXTURE.mixed.AMR.EUR = rep(NA,length(ADMIXTURE.mixed))
ADMIXTURE.mixed.AMR.EUR[ADMIXTURE.mixed=="AMR"]="AMR"
ADMIXTURE.mixed.AMR.EUR[ADMIXTURE.mixed=="EUR"]="EUR"

AMR.EUR.table = table(ADMIXTURE.mixed.AMR.EUR, meta.table$batch)
fisher.mat =AMR.EUR.table[1:2,1:2]
result = fisher.test(fisher.mat)
print(result)

#reported ethnicity
reported.race.table = read.table("../../../../Copied_Files/Ogembo_Cohort_Demographics.txt", head=T, sep="\t")					
reported.race.table$PathID = paste("S",as.character(reported.race.table$PathID),sep="")
#correct label is 51450-T --> do not change back to match minor typo (51540-T)

#matchedIDs = reported.race.table$PathID[match(meta.table$SAMPLEID,reported.race.table$PathID, nomatch=0)]
#print(meta.table$SAMPLEID[-match(matchedIDs,meta.table$SAMPLEID)])
#print(reported.race.table$PathID[-match(matchedIDs,reported.race.table$PathID)])

reported.race = reported.race.table$Race[match(meta.table$SAMPLEID,reported.race.table$PathID)]
reported.race = as.character(reported.race)

#compare to earlier metadata table, and update

#all discrepanices change from "Missing" or "Other" to "White/Caucasian"
#print(table(reported.race, mapped.age.table$lookup_Display))
mapped.age.table$lookup_Display = as.character(mapped.age.table$lookup_Display)
mapped.age.table$lookup_Display[mapped.age.table$lookup_Display== "Native American"]="Other"
mapped.age.table$lookup_Display[mapped.age.table$lookup_Display== "Pacific Islander"]="Other"

print(table(reported.race))
reported.race[(reported.race == "Missing ")&!is.na(mapped.age.table$lookup_Display)&!is.na(reported.race)]=as.character(mapped.age.table$lookup_Display[(reported.race == "Missing ")&!is.na(mapped.age.table$lookup_Display)&!is.na(reported.race)])
print(table(reported.race))
reported.race[(reported.race == "Other")&!is.na(mapped.age.table$lookup_Display)&!is.na(reported.race)]=as.character(mapped.age.table$lookup_Display[(reported.race == "Other")&!is.na(mapped.age.table$lookup_Display)&!is.na(reported.race)])
print(table(reported.race))

reported.race.table = table(reported.race, meta.table$batch)
print(reported.race.table)

meta.table = data.frame(meta.table, reported.race)

reported.race.table = table(reported.race, meta.table$batch)
asian.white.table = reported.race.table[c(1,5),2:3]
result = fisher.test(asian.white.table)
print(result)

ADMIXTURE.table = table(meta.table$ADMIXTURE.mixed, meta.table$batch)
print(ADMIXTURE.table)
AMR.EUR.table = ADMIXTURE.table[c(3,5),1:2]
result = fisher.test(AMR.EUR.table)
print(result)

#barcode and read counts
barcode1 = read.table("../../../../Copied_Files/barcode_and_read_count_161007.txt", head=T, sep="\t")
barcode1$SampleID = as.character(barcode1$SampleID)
barcode1$barcode = as.character(barcode1$barcode)
print(paste("DNA: ratio=",max(barcode1$TotalReads)/min(barcode1$TotalReads),", diff=",max(barcode1$TotalReads)-min(barcode1$TotalReads),sep=""))

barcode2 = read.table("../../../../Copied_Files/barcode_and_read_count_161206.txt", head=T, sep="\t")
barcode2$SampleID = as.character(barcode2$SampleID)
barcode2$barcode = as.character(barcode2$barcode)
print(paste("Frozen: ratio=",max(barcode2$TotalReads)/min(barcode2$TotalReads),", diff=",max(barcode2$TotalReads)-min(barcode2$TotalReads),sep=""))

barcode3 = read.table("../../../../Copied_Files/barcode_and_read_count_170118.txt", head=T, sep="\t")
barcode3$SampleID = as.character(barcode3$SampleID)
barcode3$barcode = as.character(barcode3$barcode)
print(paste("FFPE: ratio=",max(barcode3$TotalReads)/min(barcode3$TotalReads),", diff=",max(barcode3$TotalReads)-min(barcode3$TotalReads),sep=""))

barcode.table = rbind(barcode1, barcode2)
barcode.table = rbind(barcode.table, barcode3)
print(dim(barcode.table))
barcode.table = barcode.table[match(meta.table$Seq.IGC.ID, barcode.table$SeqID),]
print(dim(barcode.table))

#add lane information
lane.files = c("../../../../Copied_Files/read_stats_161007_by_sample.txt",
				"../../../../Copied_Files/read_stats_161206_by_sample.txt",
				"../../../../Copied_Files/read_stats_170118_by_sample.txt")
				
for (i in 1:length(lane.files)){
	if(i == 1){
		lane.table = read.table(lane.files[i], head=T, sep="\t")
	}else{
		temp.table = read.table(lane.files[i], head=T, sep="\t")
		temp.table = apply(temp.table, 2, as.character)
		lane.table = apply(lane.table, 2, as.character)
		lane.table = rbind(lane.table, temp.table)
	}
}#end for (i in 1:length(lane.files))

lane.table = data.frame(lane.table)
runInfo = lane.table$sample.runInfo[match(meta.table$Sample, lane.table$sampleID)]

meta.table = data.frame(meta.table, barcode=barcode.table$barcode, runInfo,
						total.reads = barcode.table$TotalReads, checkID=barcode.table$SampleID)

write.table(meta.table, output.file, quote=F, sep="\t", row.names=F)