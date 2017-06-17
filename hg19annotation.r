#### Homo sapiens hg19 annotation
#### L. Lacroix, laurent.lacroix@inserm.fr
#### 2015-09-30

################################################################################
### Home made annotation of the human genome to avoid ambiguity ################
################################################################################

require(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
`%ni%`<- Negate(`%in%`)

### limitation to chr 1 to Y

seqnames <- c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9",
"chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19",
"chr20","chr21","chr22","chrX","chrY")
seqlevels(txdb,force=T) <- seqnames

gene <- genes(txdb)	#23274 
# only genes mapping to a unique location are kept 
transcripts <- transcripts(txdb)	#78807
names(transcripts) <- transcripts$tx_name
transBgen <- transcriptsBy(txdb, by="gene")	#23434 more than gene number!
transBgenClean <- transBgen[names(transBgen) %in% names(gene)]	#23274
transclean <- transcripts[transcripts$tx_name %in% unlist(transBgenClean)$tx_name]

txnames <- transclean$tx_name	#68842	tx_names uniques from single locus transc
genenames <- names(gene)		#23274	geneID from single locus gene

# promoters
prom <- promoters(txdb, upstream=1000, downstream=0)	#78807 (49052 uniques)
names(prom) <- prom$tx_name

promclean <- prom[prom$tx_name %in% txnames]	#68842

# exons
exonallnames <- exonsBy(txdb,by='tx', use.name=T)	#78807
exonallclean <- exonallnames[names(exonallnames) %in% txnames]	#68842
ulexonall <- unlist(exonallclean)
exon1 <- unlist(exonallclean)[cumsum(c(1,elementLengths(exonallclean)[1:length(elementLengths(exonallclean))-1]))]
length(exon1) 	#68842 (40839 uniques)
exonlast <- unlist(exonallclean)[cumsum(elementLengths(exonallclean))]
#68842 (35842 uniques)
uniqex <- exonlast[exonlast$exon_rank==1]	#3373 (3373 uniques)
Texonlast <- exonlast[exonlast$exon_rank!=1]	#65469 (32469 uniques)
Texon1 <- exon1[exon1 %ni% uniqex]	#65467 (37466 uniques)
Uuniqex=unique(uniqex)
UTexon1 <- unique(Texon1)
Uexonall <- unique(ulexonall)	#254964
# we have Uexon1=UTexon1+Uuniqex and Uexonlast=UTexonlast+Uuniqex but there are 
# exonlast that superimposed with first exon but with different rank!
TTexonlast <- exonlast[exonlast$exon_id %ni% exon1$exon_id]	#65160
UTTexonlast <- unique(TTexonlast)	#32286
exonNOTother <- c(Uuniqex,UTexon1,UTTexonlast)
Uexonother <- Uexonall[Uexonall$exon_id %ni% exonNOTother$exon_id]	#181839
#Uexonother2 <- Uexonall[!overlapsAny(Uexonall,exonNOTother, type='equal')]
AnnoExon <- GRangesList(Uuniqex,UTexon1,UTTexonlast,Uexonother)
names(AnnoExon) <- c('Exon_unique','Exon_1','Exon_last','Exon_other')
length(Uexonall)==length(Uuniqex)+length(UTexon1)+length(UTTexonlast)+length(Uexonother)
# TRUE

# introns
intronall <- intronsByTranscript(txdb, use.names=T)	
#78807
intronallclean <- intronall[names(intronall) %in% txnames]	
#68842
intronall2 <- intronallclean[elementLengths(intronallclean)!=0]	
#65469
length(intronallclean[elementLengths(intronallclean)==0])	
#3373, like the number of uniqex

ulintronall <- unlist(intronall2) #613262 (239396 uniques)
mtutu <- ulintronall[strand(ulintronall)=='-']
ptutu <- ulintronall[strand(ulintronall)=='+']
pnom <- unique(names(ptutu))	#33275
mnom <- unique(names(mtutu))	#32194
pintron <- intronall2[names(intronall2) %in% pnom]
mintron <- intronall2[names(intronall2) %in% mnom]
length(intronall2)==length(pintron)+length(mintron)	#TRUE
pintron1 <- unlist(pintron)[cumsum(c(1,elementLengths(pintron)[1:length(elementLengths(pintron))-1]))]	#33275 
mintron1 <- unlist(mintron)[cumsum(elementLengths(mintron))]	#32194
intron1 <- c(pintron1,mintron1)	#65469, same length as intronall2
pintronlast <- unlist(pintron)[cumsum(elementLengths(pintron))]
mintronlast <- unlist(mintron)[cumsum(c(1,elementLengths(mintron)[1:length(elementLengths(mintron))-1]))]
intronlast <- c(pintronlast,mintronlast)	#65469
uniqintr <- unlist(intronall2[elementLengths(intronall2)==1])	#3601
Uintron1 <- unique(intron1)	#43432
Uintronlast <- unique(intronlast)	#35942
Uuniqintr <- unique(uniqintr)	#3590
TUintron1 <- Uintron1[Uintron1 %ni% Uuniqintr]	#39842
TUintronlast <- Uintronlast[Uintronlast %ni% Uuniqintr]	#32352 
# but it remains cases at the intersection between TUintron1 and TUintronlast
TTUintronlast <- TUintronlast[TUintronlast %ni% TUintron1]	#32005
intronNOother <- c(Uuniqintr,TUintron1,TTUintronlast)	#75437 (uniques)
Uintronother <- unique(ulintronall)[unique(ulintronall) %ni% intronNOother]

AnnoIntron <- GRangesList(Uuniqintr,TUintron1,TTUintronlast,Uintronother)
names(AnnoIntron) <- c('Intron_unique','Intron_1','Intron_last','Intron_other')

# UTR and CDS
utr5p <- fiveUTRsByTranscript(txdb,use.name=T)	#58290
utr3p <- threeUTRsByTranscript(txdb,use.name=T)	#57673
cds <- cdsBy(txdb,by='tx',use.name=T)	#60396
utr5pclean <- utr5p[names(utr5p) %in% txnames]	#56859
utr3pclean <- utr3p[names(utr3p) %in% txnames]	#56301
cdsclean <- cds[names(cds) %in% txnames]	#58118
ulcdsclean <- unlist(cdsclean)	#558559 (224004 uniques)
Uutr5p <- unique(unlist(utr5pclean))	#66866
Uutr3p <- unique(unlist(utr3pclean))	#31745
CDS <- unique(ulcdsclean)	#224004
Uprom <- unique(promclean) #39692
Utrans <- unique(transclean)	#48207
TSS <- resize(promclean,fix='end',width=1)	#68842
UTSS <- unique(TSS)	#39692
TES <- resize(transclean,fix='end',width=1)	#68842
UTES <- unique(TES)	#34736

Trans <- GRangesList(Utrans,UTSS,UTES)
names(Trans)=c('Transcript','TSS','TES')
Prom <- Uprom

UTR <- GRangesList(Uutr5p,Uutr3p)
names(UTR) <- c('UTR5p','UTR3p')

# saving results
var2save <- c('Prom','Trans','AnnoIntron','AnnoExon','CDS','UTR')
save(list=var2save,file='hg19Anno.RData')

################################################################################

