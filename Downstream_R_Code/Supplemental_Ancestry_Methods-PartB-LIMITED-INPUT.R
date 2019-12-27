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

##Vary the minimum fraction for an ADMIXTURE assignment
#component.threshold = 0.8
component.threshold = 0.5

K=5

meta.table = read.table("../../Copied_Files/MXL_test_unrelated_sample_description.txt", head=T, sep="\t")
sample.order = read.table("../../Copied_Files/MXL_test_1KG_unrelated.fam", head=F, sep="\t")

censored.pop.code = read.table("../../Copied_Files/MXL_test_1KG_unrelated.pop", head=F)
censored.pop.code[censored.pop.code == "-"]=NA
censored.pop.code = as.character(censored.pop.code[[1]])

pop.code = meta.table$ethnicity[match(sample.order[,2],meta.table$subject)]
super.pop.code = meta.table$super.pop[match(sample.order[,2],meta.table$subject)]

admixture.assignment.table = read.table(paste("../../Copied_Files/MXL_test_1KG_unrelated.",K,".Q",sep=""), head=F, sep=" ")
super.pops = c("AFR","EAS","EUR","AMR","SAS")
colnames(admixture.assignment.table) = super.pops

ref.pop.freq = apply(admixture.assignment.table, 2, mean.pop, pop = censored.pop.code)
print("Average Frequency per Column:")
print(ref.pop.freq)

admixture.max.assignment = apply(admixture.assignment.table, 1, max.assignment, pop.IDs=super.pops, cutoff=component.threshold)
output.table = data.frame(Sample=sample.order[,2],
							assigned.population = pop.code, assigned.super.pop = super.pop.code,
							used.pop.code = censored.pop.code,
							admixture.assignment = admixture.max.assignment, admixture.assignment.table)
#write.table(output.table,
#			paste("unrelated_1KG_ADMIXTURE_assignments_",component.threshold,".txt",sep=""),
#			sep="\t", quote=F, row.names=F)

admixture.assignment.table = admixture.assignment.table[!is.na(pop.code),]
pop.code = meta.table$super.pop[!is.na(pop.code)]

pop.freq = apply(admixture.assignment.table, 2, mean.pop, pop = pop.code)
#write.table(pop.freq, "ADMIXTURE_super_pop_means.txt",sep="\t", quote=F)

#validation set
output.table = output.table[is.na(output.table$used.pop.code),]

final.assignments = output.table$admixture.assignment
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