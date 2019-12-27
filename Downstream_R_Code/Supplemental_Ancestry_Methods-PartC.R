meta.file = "Selected_Output_Files/combined_genotype_with_year_and_ethnicity.txt"

meta.table = read.table(meta.file, head=T, sep = "\t")
meta.table$reported.race[meta.table$reported.race == "Unknown"]=NA

ADMIXTURE.max.confusion.table = table(meta.table$reported.race, meta.table$ADMIXTURE.top)
print(t(ADMIXTURE.max.confusion.table))

bootstrap.confusion.table = table(meta.table$reported.race, meta.table$bootstrap.ethnicity)
print(t(bootstrap.confusion.table))

ADMIXTURE.mixed.confusion.table = table(meta.table$reported.race, meta.table$ADMIXTURE.mixed)
print(t(ADMIXTURE.mixed.confusion.table))

#check for discrepancies
print(table(meta.table$ADMIXTURE.top[!is.na(meta.table$reported.race)], meta.table$bootstrap.ethnicity[!is.na(meta.table$reported.race)]))