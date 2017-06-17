#### Script collection used in
#### "Reevaluation of quadruplex propensity with G4Hunter"
#### L. Lacroix, laurent.lacroix@inserm.fr
#### 2015-09-28

## developed for R version 3.2.2 on mac OS X 10.10 running on a macbook pro 2013

#### Functions are located in the file function_G4Hunter.r

workingpath <-  '~/Desktop/G4HunterPaper/'
setwd(workingpath)
source(paste(workingpath,'function_G4Hunter.r',sep=''))

#### Example to compute G4Hscore for a single sequence
A22g <- 'AGGGTTAGGGTTAGGGTTAGGG'
G4Hscore(A22g)

#### Example to compute G4Hscore for a list of sequcences

RefSet <- read.table('ReferenceDataset.txt',header=T)
RefSet[,4] <- sapply(RefSet[,2],function(x) signif(G4Hscore(x),3))
names(RefSet)[4] <- 'G4Hscore'
write.table(as.data.frame(RefSet),'RefSetCalc.txt',sep='\t',col.names=NA)

#### use of the function QPtest on the reference dataset to test the presence of Quadparser type sequences

RefSet[,5] <- sapply(RefSet[,2],QPtest)
RefSet[,6] <- sapply(RefSet[,2],QPtest,lmax=12)
RefSet[,7] <- sapply(RefSet[,2],QPtest,gmin=2)
colnames(RefSet)[5:7] <- c('QP37','QP312','QP27')




#### Example to extract G4FS form the mitochondrial genome
#### Working with the release hg19 for the human genome
library("BSgenome.Hsapiens.UCSC.hg19")
genome <- BSgenome.Hsapiens.UCSC.hg19

toto=G4hunt(i=25,k=25,hl=1,gen=genome,masked=5)
totolist=G4huntlist(i=25,k=25,hl=c(1,1.2,1.5),gen=genome,masked=5)
titi=G4huntrefined(toto,gen=genome,i=25)
write.table(as.data.frame(titi),'MitoRef.txt',sep='\t',col.names=NA)


#### Example to extract G4FS from SacCer3 with a threshold of 1.5
library("BSgenome.Scerevisiae.UCSC.sacCer3")
genome <- BSgenome.Scerevisiae.UCSC.sacCer3
chrlist <- 1:17

G4H_sc3_1.5 <- do.call(c,lapply(chrlist,G4hunt,hl=1.2))
seqinfo(G4H_sc3_1.2) <- seqinfo(genome)
# refining results
G4H_sc3_1.2_ref <- sort(do.call(c,lapply(chrlist,G4huntrefined,gen=genome,Res_hl=G4H_sc3_1.2)),ignore.strand=T)

# this result can be then export to bed or txt, use to calculate a coverage by 1kb bin to generate a bigwig
# this result can also be filtered for the width to eliminate the very long or very short G4FS



#### Example to extract G4FS from hg19 with a threshold of 1.75
library("BSgenome.Hsapiens.UCSC.hg19")
genome <- BSgenome.Hsapiens.UCSC.hg19

chrlist <- 1:25

G4H_hg19_1.75 <- do.call(c,lapply(chrlist,G4hunt,hl=1.75))
seqinfo(G4H_hg19_1.75) <- seqinfo(genome)
seqlevels(G4H_hg19_1.75, force=T) =seqlevels(genome) [1:25]

# refining results
G4H_hg19_1.75_ref <- sort(do.call(c,lapply(chrlist,G4huntrefined,gen=genome,Res_hl=G4H_hg19_1.75)),ignore.strand=T)

# this result can be then export to bed or txt, use to calculate a coverage by 1kb bin to generate a bigwig
# this result can also be filtered for the width to eliminate the very long or very short G4FS
save(G4H_hg19_1.75_ref,file='G4H_hg19_1.75_ref.RData')



#### exporting to bed, txt, bigwig with 1kb resolution


name2export <- deparse(substitute(G4H_sc3_1.2_ref))
var2export=get(name2export)
bed2save=paste(workingpath,name2export,'.bed',sep='')
export(var2export,con=bed2save)

bigwig2export <- G4GR2bigwig(var2export,gen=genome,mcore=mcore)
export(bigwig2export,paste(workingpath,name2export,'.bigwig',sep=''))

write.table(as.data.frame(var2export),paste(workingpath,name2export,'.txt',sep=''),sep='\t',col.names=NA)

# bed and bigwig file can be used with IGV

### Examples with Quadparser

library("BSgenome.Scerevisiae.UCSC.sacCer3")
genome <- BSgenome.Scerevisiae.UCSC.sacCer3
chrlist <- 1:17

QP37_sc3 <- do.call(c,lapply(chrlist,quadparser,gen=genome))

# this result can be then export to bed or txt, use to calculate a coverage by 1kb bin to generate a bigwig
# this result can also be filtered for the width to eliminate the very long or very short G4FS


# examples using other parameters for quadparser: gmin/gmax are the minimal/maximal sizes for the G-run
# lmin/max are the minimal/maximal sizes for loops
QP27_sc3 <- do.call(c,lapply(chrlist,quadparser,gen=genome,gmin=2,gmax='',lmin=1,lmax=7))
QP312_sc3 <- do.call(c,lapply(chrlist,quadparser,gen=genome,gmin=3,gmax='',lmin=1,lmax=12))




















