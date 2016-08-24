# Dalliance
files associated with Dalliance


###How to liftOver from Assembly 2.4 to 2.5

1. Obtain chain file for the 2 assemblies and chromsizes for new assmebly
2. `liftOver` only takes bed or gff files (bed is recommended)
3. Convert bigBED to bed using `bigBedToBed old.bb old.bed`
4. Run `liftOver old.bed over.chainfile new.bed unMapped`
5. Sort new bed file `sort -k1,1 -k2,2n new.bed > sorted.new.bed`
6. Convert bed to bigBED `bedToBigBed sorted.new.bed 2.5_chrom_sizes new.bb`
7. The new.bb file can now be used to replace old.bb file

###How to liftOver SNPs from Assembly 2.4 to 2.5

SNPs cannot have same start and end location when lifting. Running liftOver on a SNP bed file will result in all SNPs being excluded. To fix, subtract 1 from every start position.

```
test <- read.table("old.bed", header = F, as.is = T) ##Read in bed file, may need to remove any comments at top of file and readd later
options("scipen" = 100, "digits" = 4) #Prevent scientific notation
test[,2] <- test[,2]-1 #Bed uses base 0, this will allow for remapping
write.table(test,"minus1bed",col.names=F,row.names=F,sep="\t",quote=F)
liftOver minus1bed over.chainfile liftedminus1bed unMapped
test <- read.table("liftedminus1bed", header = F, as.is = T)
test[,2] <- test[,2]+1 #Shifting position back
write.table(test,"new.bed",col.names=F,row.names=F,sep="\t",quote=F)
```
