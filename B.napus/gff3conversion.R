### load needed library
library(dplyr)
options(scipen = 100)
### load gff3 converted to bed
bed <- read.delim("Brassica_napus.annotation_v5.bed", header = F)

### simple changes

modified.bed <- subset(bed, V8 == "mRNA")
modified.bed$V2 <- modified.bed$V2 - 1 #bed files are 0 indexed, gff3 files are 1 indexed
modified.bed$V9 <- "255,0,0"
modified.bed$V5 <- 1000

### harder changes
## block numbers

blocknumbers <- unlist(bed %>% group_by(V4) %>%
  filter(V8 != "mRNA") %>%
  mutate(counts = n()) %>%
  select(V4,counts) %>%
  unique() %>%
  ungroup() %>%
  select(counts))
modified.bed$V10 <- blocknumbers

## block sizes

blocksizes <- bed %>% group_by(V4) %>%
  filter(V8 != "mRNA") %>%
  mutate(sizes = V3-V2, allblocks = list(sizes)) %>%
  select(V4,allblocks) %>%
  unique() %>%
  ungroup() %>%
  select(allblocks)
blocksizes$allblocks <- gsub("[c\\(\\) ]","",blocksizes$allblocks)
modified.bed$V11 <- blocksizes$allblocks

## block starts

CDS <- bed %>% 
  filter(V8 != "mRNA") %>%
  select(V2,V4)

mRNA <- modified.bed$V2
blockstarts <- data.frame(CDS, rep(mRNA,blocknumbers))
blockstarts$new_start <- blockstarts$V2 - blockstarts$rep.mRNA..blocknumbers.

blockstarts <- blockstarts %>% group_by(V4) %>%
  mutate(allstarts = list(new_start)) %>%
  select(V4,allstarts) %>%
  unique() %>%
  ungroup() %>%
  select(allstarts)
blockstarts$allstarts <- gsub("[c\\(\\) ]","",blockstarts$allstarts)
blockstarts$allstarts <- sub("^[0-9]+","0",blockstarts$allstarts)
modified.bed$V12 <- blockstarts$allstarts

## Adjust block sizes for 1 block features

modified.bed <- within(modified.bed, V11[V12 == '0'] <- as.character(as.numeric(V11[V12 == '0']) + 1))

## Adjust for weird stuff

modified.bed$V11 <- sub(":",",",modified.bed$V11)

## Add thick start and end

#thickstart <- as.numeric(sub(",.*","",sub("^0,","",blockstarts$allstarts)))
#modified.bed$V7 <- thickstart + modified.bed$V2
modified.bed$V7 <- modified.bed$V2

#thickend <- as.numeric(gsub(".*,","",blockstarts$allstarts))
#modified.bed$V8 <- thickend + modified.bed$V2
modified.bed$V8 <- modified.bed$V3

## write new bed file

write.table(modified.bed,"B.napus_V4.1_gene_models.bed", row.names = F, quote = F, col.names = F, sep = "\t")

## Need to manually change BnaA06g20800D due to too many features. Output prints on multiple lines
