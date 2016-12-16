# Creating Chain.over file for (OLD) ITAG2.4 to (NEW) ITAG2.5 (Can be modified to work with other organisms)

## Download and install necessary tools

Downloading all tools at once isn't working so will install one by one

```
cd /usr/local/bin
sudo -s
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/blat/blat
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/faSplit
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/chainMergeSort
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/chainNet
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/twoBitInfo
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/liftUp
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/axtChain
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/netChainSubset
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/faToTwoBit
###Set all to executable
chmod 755 blat faSplit chainMergeSort chainNet twoBitInfo liftUp axtChain netChainSubset faToTwoBit
```

## Setup workspace and download genome files from solgenomics.net


```
mkdir LiftOver
cd LiftOver
mkdir -p ITAG2.4 ITAG2.5 lift raw psl scratch/split
cd ITAG2.4
wget ftp://ftp.solgenomics.net/tomato_genome/assembly/build_2.40/S_lycopersicum_chromosomes.2.40.fa.gz
gunzip S_lycopersicum_chromosomes.2.40.fa.gz
cd ../ITAG2.5
wget ftp://ftp.solgenomics.net/tomato_genome/assembly/build_2.50/S_lycopersicum_chromosomes.2.50.fa.gz
gunzip S_lycopersicum_chromosomes.2.50.fa.gz
```

## Get chromosome sizes for later

```
cd LiftOver/ITAG2.4
faToTwoBit S_lycopersicum_chromosomes.2.40.fa S_lycopersicum_chromosomes.2.40.2bit
twoBitInfo S_lycopersicum_chromosomes.2.40.2bit 2.40.chrom.sizes
cd LiftOver/ITAG2.5
faToTwoBit S_lycopersicum_chromosomes.2.50.fa S_lycopersicum_chromosomes.2.50.2bit
twoBitInfo S_lycopersicum_chromosomes.2.50.2bit 2.50.chrom.sizes
```

## Split the (NEW) ITAG2.5 genome by chromosome then into 3kb chunks

```
cd LiftOver
for i in 00 01 02 03 04 05 06 07 08 09 10 11 12; do faSplit -lift=lift/SL2.50ch${i}.lft size ITAG2.5/S_lycopersicum_chromosomes.2.50.fa -oneFile 3000 scratch/split/SL2.50ch${i}; done
cd ..
```

## BLAT query sequences from (NEW) ITAG2.5 against (OLD) ITAG2.4

```
cd psl
for i in 00 01 02 03 04 05 06 07 08 09 10 11 12; do blat ../ITAG2.4/S_lycopersicum_chromosomes.2.40.2bit ../scratch/split/SL2.50ch${i}.fa SL2.40ch${i}.psl -tileSize=12 -minScore=100 -minIdentity=98 -fastMap; done
```

## Use liftUp to change the coordinate system. Requires .lft files

```
cd raw
for i in 00 01 02 03 04 05 06 07 08 09 10 11 12; do liftUp -pslQ SL2.50ch$i.psl ../lift/SL2.50ch$i.lft warn ../psl/SL2.50ch$i.psl; done

```


# Dalliance
files associated with Dalliance

###How to liftOver from Assembly 2.4 to 2.5

##Lifting over bed files
```
liftOver yourfile.bed over.chainfile youfilelifted.bed unMapped
```

##Lifting over bigBED files

1. Obtain chain file for the 2 assemblies and chromsizes for new assembly
2. `liftOver` only takes bed or gff files (bed is recommended)
3. Convert bigBED to bed using `bigBedToBed old.bb old.bed`
4. Run `liftOver old.bed over.chainfile new.bed unMapped`
5. Sort new bed file `sort -k1,1 -k2,2n new.bed > sorted.new.bed`
6. Convert bed to bigBED `bedToBigBed sorted.new.bed 2.5_chrom_sizes new.bb`
7. The new.bb file can now be used to replace old.bb file

##How to liftOver SNP bedfiles from Assembly 2.4 to 2.5

SNPs cannot have same start and end location when lifting. Running liftOver on a SNP bed file will result in all SNPs being excluded. To fix, subtract 1 from every start position. The below example is being performed in R

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
## Lifting over bam files

Depending on what information you wish to keep, you may need to change the options on the utilities
This is simply to show how to liftOver the positions
```
bedtool bamtobed -i yourfile.bam > yourfile.bed
liftOver yourfile.bed over.chainfile yourfilelifted.bed unMapped
sort -k1,1 -k2,2n yourfilelifted.bed > sorted.yourfilelifted.bed
bedtools bedtobam -i sorted.yourfilelifted.bed -g 2.5_chromo_sizes > yourfilelifted.bam
```
