mean.pop = function(arr, pop){
	return(tapply(arr,pop,mean))
}

max.assignment = function(arr, pop.IDs, cutoff){
	assignment = pop.IDs[arr > cutoff]
	if(length(assignment) == 1){
		return(assignment)
	}else{
		return(NA)
	}
}#end def max.assignment

min.conf = 95
K=5

meta.table = read.table("../../Copied_Files/MXL_test_unrelated_sample_description.txt", head=T, sep="\t")
sample.order = read.table("../../Copied_Files/MXL_test_1KG_unrelated.fam", head=F, sep="\t")

censored.pop.code = read.table("../../Copied_Files/MXL_test_1KG_unrelated.pop", head=F)
censored.pop.code[censored.pop.code == "-"]=NA
censored.pop.code = as.character(censored.pop.code[[1]])

pop.code = meta.table$ethnicity[match(sample.order[,2],meta.table$subject)]
super.pop.code = meta.table$super.pop[match(sample.order[,2],meta.table$subject)]

boot.assignment.table = read.table("../../Copied_Files/MXL_test_overall_bootstrap_assignments.txt", head=T, sep="\t")
bootstrap.assignment = boot.assignment.table$min.overall.dist

bootstrap.count.mat = boot.assignment.table[,2:ncol(boot.assignment.table)]

extract.bootstrap.confidence = function(arr){
	temp.assignment = arr[1]
	bootstrap.counts = as.numeric(arr[2:length(arr)])
	level.names = names(arr)[2:length(arr)]
	level.names = gsub(".bootstrap","",level.names)
	return(bootstrap.counts[level.names==temp.assignment])
}#end def extract.bootstrap.confidence

bootstrap.confidence = apply(bootstrap.count.mat, 1, extract.bootstrap.confidence)

output.table = data.frame(Sample=boot.assignment.table$sample,
							assigned.population = pop.code[match(boot.assignment.table$sample,sample.order[,2])],
							assigned.super.pop = super.pop.code[match(boot.assignment.table$sample,sample.order[,2])],
							used.pop.code = censored.pop.code[match(boot.assignment.table$sample,sample.order[,2])],
							bootstrap.assignment = bootstrap.assignment,
							bootstrap.confidence = bootstrap.confidence,
							boot.assignment.table)
#write.table(output.table,
#			paste("unrelated_1KG_bootstrap_concordance.txt",sep=""),
#			sep="\t", quote=F, row.names=F)

#validation set
final.assignments = output.table$bootstrap.assignment
final.assignments[bootstrap.confidence < min.conf]=NA
known.super.pop = output.table$assigned.super.pop
known.pop = output.table$assigned.population

print("Overall Super-Pop")
print(length(final.assignments))
super.pop.assignments = table(known.super.pop, final.assignments)
print(sum(super.pop.assignments))

print(t(super.pop.assignments))

accuracy = (super.pop.assignments[1,1]
			+super.pop.assignments[2,2]
			+super.pop.assignments[3,3]
			+super.pop.assignments[4,4]
			+super.pop.assignments[5,5])/sum(super.pop.assignments)
print(accuracy)

print("AMR")
temp.known = known.super.pop[known.super.pop == "AMR"]
print(length(temp.known))
temp.assigned = final.assignments[known.super.pop == "AMR"]
temp.table = table(temp.known, temp.assigned)
print(t(temp.table))
print(sum(temp.table))

AMR.counts = temp.table[,"AMR"]
accuracy = (AMR.counts[2])/sum(temp.table)
print(accuracy)

print("CLM")
temp.known = known.super.pop[known.pop == "CLM"]
print(length(temp.known))
temp.assigned = final.assignments[known.pop == "CLM"]
temp.table = table(temp.known, temp.assigned)
print(t(temp.table))
print(sum(temp.table))

AMR.counts = temp.table[,"AMR"]
accuracy = (AMR.counts[2])/sum(temp.table)
print(accuracy)

print("MXL")
temp.known = known.super.pop[known.pop == "MXL"]
print(length(temp.known))
temp.assigned = final.assignments[known.pop == "MXL"]
temp.table = table(temp.known, temp.assigned)
print(t(temp.table))
print(sum(temp.table))

AMR.counts = temp.table[,"AMR"]
accuracy = (AMR.counts[2])/sum(temp.table)
print(accuracy)

print("PEL")
temp.known = known.super.pop[known.pop == "PEL"]
print(length(temp.known))
temp.assigned = final.assignments[known.pop == "PEL"]
temp.table = table(temp.known, temp.assigned)
print(t(temp.table))
print(sum(temp.table))

AMR.counts = temp.table[,"AMR"]
accuracy = (AMR.counts[2])/sum(temp.table)
print(accuracy)

print("PUR")
temp.known = known.super.pop[known.pop == "PUR"]
print(length(temp.known))
temp.assigned = final.assignments[known.pop == "PUR"]
temp.table = table(temp.known, temp.assigned)
print(t(temp.table))
print(sum(temp.table))

AMR.counts = temp.table[,"AMR"]
accuracy = (AMR.counts[2])/sum(temp.table)
print(accuracy)