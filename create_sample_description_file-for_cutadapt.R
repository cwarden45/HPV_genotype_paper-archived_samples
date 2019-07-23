#NOTE: this code works with the format of bcl2fastq version 2.17 (or later).  You'll need to modify code (or manually create file) if reads are formatted differently.

input.folder = "/path/to/Downsample_Reads"
output.file = "sample_description_file.txt"

files = list.files(input.folder)

print(length(files))
Forward_Read = files[grep("_R1_001.fastq",files)]
print(length(Forward_Read))

Reverse_Read = gsub("_R1_001.fastq","_R2_001.fastq",Forward_Read)

##if you want to compare samples from multiple lanes
#SampleID = gsub("_R1_001.fastq","",Forward_Read)
SampleID = gsub("_S\\d+_L\\d{3}_R1_001.fastq","",Forward_Read)

output.table = data.frame(SampleID, Forward_Read, Reverse_Read)

write.table(output.table, output.file, quote=F, sep="\t", row.names=F)