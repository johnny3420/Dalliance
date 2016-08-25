# Dalliance
files associated with Dalliance

###How to liftOver from Assembly 2.4 to 2.5

##Lifting over bigBED files

1. Obtain chain file for the 2 assemblies and chromsizes for new assembly
2. `liftOver` only takes bed or gff files (bed is recommended)
3. Convert bigBED to bed using `bigBedToBed old.bb old.bed`
4. Run `liftOver old.bed over.chainfile new.bed unMapped`
5. Sort new bed file `sort -k1,1 -k2,2n new.bed > sorted.new.bed`
6. Convert bed to bigBED `bedToBigBed sorted.new.bed 2.5_chrom_sizes new.bb`
7. The new.bb file can now be used to replace old.bb file

##How to liftOver SNP bedfiles from Assembly 2.4 to 2.5

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

##To liftover bigBED files within Dalliance (Only bigBED files can be directly remapped)

Have to convert over.chain to bigBED.

Use instructions at https://github.com/dasmoth/hubbub

Example setup for chain
```
chains: {
        itag24Toitag25: {uri: '/ITAG2.4ToITAG2.5.bb',
                        type: 'bigbed',
                        coords: {
                        speciesName: 'Solanum lycopersicum',
                              taxon: 4081,
                              auth: ' ',
                              version: "ITAG2.3"
                         }        
               }      
      },
      
sources: [  {name: 'Genes',
                      desc: 'on-the-fly mapping of genes from ITAG2.3',
                      bwgURI: '/itag23_gene_modelsnew.bb',
                      mapping: 'itag24Toitag25'
                     },
                ]
                
 ``` 
  
## Lifting over bigwig files

```
bigWigToBedGraph youfile.bw yourfile.bed
liftOver yourfile.bed overchain.file yourliftedfile.bed unMapped
sort -k1,1 -k2,2n yourliftedfile.bed > yourliftedfile.sorted.bed
bedGraphToBigWig yourliftedfile.sorted.bed 2.5_chromo_sizes yourliftedfile.bw
```
