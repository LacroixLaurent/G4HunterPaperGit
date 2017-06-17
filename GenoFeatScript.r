#### Script to look for in genomic features
#### L. Lacroix, laurent.lacroix@inserm.fr
#### 2015-09-30

workingpath <- '~/Desktop/G4HunterPaper/'
setwd(workingpath)
# loading HG19 annotation data
load(paste(workingpath,'hg19Anno.RData',sep=''))
######

# evaluation of the genome size without N
# chrM excluded because it lacks annotations

library(BSgenome.Hsapiens.UCSC.hg19)
genome <- BSgenome.Hsapiens.UCSC.hg19
N=0
for (i in 1:24)
	{ N=N+alphabetFrequency(genome[[i]])[15]}
	
genomesize <- sum(as.numeric(seqlengths(genome)[1:24]))
genomesizeNoN <- genomesize-N

### loading the list of G4FS as a GRanges

load('G4H_hg19_1.75_ref.RData')
# removing G4FS from chrM for which no annotation are present in our annotation file

G4data <- G4H_hg19_1.75_ref[seqnames(G4H_hg19_1.75_ref)!='chrM']


# fraction of genome in promoter region
# note that we use reduce(Prom,ignore.strand=T) to define 'promoter region'
fProm <- sum(width(reduce(Prom,ignore.strand=T)))/genomesizeNoN
# 1.2%
# fraction of G4FS found in promoter region
fG4inPromNotStr <- sum(countOverlaps(G4data,reduce(Prom),ignore.strand=T)!=0)/length(G4data)
# 5.6%
# fraction of promoter with at least one G4FS
fPromwithG4 <- sum(countOverlaps(reduce(Prom,ignore.strand=T),G4data,ignore.strand=T)!=0)/length(reduce(Prom,ignore.strand=T))
# 51.9%
# G4FS density global and in promoter (per kb)
gbG4dens <- length(G4data)*1000/genomesizeNoN
# 0.25/kb
promG4dens <- 1000*sum(countOverlaps(G4data,reduce(Prom,ignore.strand=T),ignore.strand=T))/sum(width(reduce(Prom,ignore.strand=T)))
# 1.12/kb

###########################
findNgaps <- function(x)
# x is a DNAString
	{ y=Rle(strsplit(as.character(x),NULL)[[1]])
	  y2=ranges(Views(y,y=='N'))
	  return(y2)	# y2 is a list of IRanges
	}


listGaps=list()
for (i in 1:24)
	{
		ra= findNgaps(genome[[i]])
		listGaps[[i]]=GRanges(seqnames=seqnames(genome)[i],ranges=ra,seqinfo=seqinfo(genome))
		print(i)
	}	
Ngaps=do.call(c,listGaps)
Ngaps2=Ngaps[width(Ngaps)>=20]

PromN=setdiff(Prom,Ngaps2,ignore.strand=T)
# but using PromN instead of reduce(Prom,ignore.strand=T) does not change much the result

sum(countOverlaps(G4data,PromN,ignore.strand=T)!=0)/length(G4data)
sum(countOverlaps(reduce(PromN,ignore.strand=T),G4data,ignore.strand=T)!=0)/length(reduce(PromN,ignore.strand=T))
1000*sum(countOverlaps(G4data,reduce(PromN,ignore.strand=T),ignore.strand=T))/sum(width(reduce(PromN,ignore.strand=T)))

# defining all the genomic features used in Table S6-8
junction_ei=unique(flank(AnnoIntron,start=T,both=T,width=50))
junction_ie=unique(flank(AnnoIntron,start=F,both=T,width=50))
names(junction_ei) <- paste(names(AnnoIntron),'_j_ei',sep='')
names(junction_ie) <- paste(names(AnnoIntron),'_j_ie',sep='')
FL <- list()
FL[[1]] <- Prom
FL <- c(FL,CDS, as.list(UTR),as.list(AnnoIntron),as.list(AnnoExon),as.list(junction_ei),as.list(junction_ie))
names(FL)[1] <- 'Prom'
names(FL)[2] <- 'CDS'

### defining intergenic region
FLgene <- FL[2:12]
for (i in 1:length(FLgene))
	{ mcols(FLgene[[i]])=NULL
	names(FLgene[[i]])=NULL}
names(FLgene)=NULL
allgene2 <- reduce(do.call(c,FLgene),ignore.strand=T)
intergen <- gaps(allgene2)
intergen <- intergen[strand(intergen)=='*']

## defining intergenic region that are also not found in promoters
FLgeneP <- FL[1:12]
for (i in 1:length(FLgeneP))
	{ mcols(FLgeneP[[i]])=NULL
	names(FLgeneP[[i]])=NULL}
names(FLgeneP)=NULL
allgeneP <- reduce(do.call(c,FLgeneP),ignore.strand=T)
intergennoP <- gaps(allgeneP)
intergennoP <- intergennoP[strand(intergennoP)=='*']


FL2=c(FL,intergennoP,intergen)
names(FL2)=c(names(FL),'intergennoP','intergen')

# FL2 is the final list of genomic features we are going to work with
# but there are some overlaps with Ngaps2 regions
sapply(FL2,function(x) sum(overlapsAny(x,Ngaps2)))
# so we will substract the Ngaps2 region from the feature list
FL <- FL2
FL2=sapply(FL,function(x) {setdiff(x,Ngaps2,ignore.strand=T)})
# note the setdiff also perform a reduce with ignore.strand=T

# for Table S6
fFeatwG4 <- sapply(FL2,function(x){sum(countOverlaps(x,G4data,ignore.strand=T)!=0)/length(x)})

# for Table S7
fFeatgb <- sapply(FL2,function(x) sum(as.numeric(width(x)))/genomesizeNoN)
names(fFeatgb) <- names(FL2)
fG4inFeat <- sapply(FL2,function(x){sum(countOverlaps(G4data,x,ignore.strand=T)!=0)/length(G4data)})

# for Table S8
G4d <- sapply(FL2,function(x) 1000*sum(countOverlaps(G4data,x,ignore.strand=T))/sum(as.numeric(width(x))))

##### resampled background

## shuffleGR, a function to resample GRanges in the genome keeping the chromosomal and strand distribution
## but avoiding long gaps of uncharacterized genome

# function to resample on a given chromosom keeping strand information
shuffleGR3=function(gen=genome,chrnb=24,inputGR=G4data,gap=Ngaps2)
	{	require(GenomicRanges)
		hit <- inputGR[seqnames(inputGR)==seqnames(gen)[chrnb]]
		gapchr=gap[seqnames(gap)==seqnames(gen)[chrnb]]
		rgap <- ranges(gapchr)
		ravail <- gaps(rgap)
		st_avail <- unlist(as.vector(ravail))
		st_rdgr <- sample(st_avail,length(hit))
		if (length(hit)==1)		
				{
				wi_rdgr <- width(hit)
				}else{
				wi_rdgr <- sample(width(hit))
				#necessary if only one range sample(width()) choose a number
				#betwen in 1:width() rather than one width
				}
		ra_rdgr <- sort(IRanges(start=st_rdgr,width=wi_rdgr))
		keep <- IRanges()
		ra_rdgr2 <- IRanges()
		while ((sum(overlapsAny(ra_rdgr,rgap))!=0) | (sum(overlapsAny(ra_rdgr2,keep))!=0))
			{
			keep <- ra_rdgr[overlapsAny(ra_rdgr,rgap)==0]
			hit2 <- ra_rdgr[overlapsAny(ra_rdgr,rgap)!=0]
			st_rdgr2 <- sample(st_avail,length(hit2))
			if (length(hit2)==1)
				{
				wi_rdgr2 <- width(hit2)
				}else{
				wi_rdgr2 <- sample(width(hit2))
				}
			ra_rdgr2 <- IRanges(start=st_rdgr2,width=wi_rdgr2)	
			ra_rdgr <- c(keep,ra_rdgr2)
			}
		rdgr <- sort(GRanges(seqnames=Rle(rep(seqnames(gen)[chrnb],length(hit))),ranges=ra_rdgr,strand=Rle(rep('*',length(hit))),seqinfo=seqinfo(inputGR)))
		return(rdgr)
	}

# function to resample on a genome

shuffleGRgen <- function(dummy=1,gen2=genome,inputGR2=G4data,gap2=Ngaps2,chrlist=1:24)
	{
	rdlist=GRangesList()
	for (i in chrlist) {rdlist[[i]] <- shuffleGR3(gen=gen2,chrnb=i,inputGR=inputGR2,gap=gap2)}
	y<- do.call(c,rdlist)
	return(y)
	}
	
# example to resample 4 times the G4data list.
# for the manuscript, this has been done 1000 times for each G4data list
# using mclapply is strongly advised


rdlength <- 1:4
# not the for the manuscript, rdlength <- 1:1000
rdlist2 <- lapply(rdlength,function(x) {shuffleGRgen(x,gen=genome,inputGR=G4data,gap=Ngaps2,chrlist=1:24)})

mcG4inFeatbg <- function(shufGR=rdlist2,FeatList=FL2,mc=4)
	{
	toto <- simplify2array(mclapply(shufGR,function(y) {sapply(FeatList,function(x) {sum(countOverlaps(y,x,ignore.strand=T)!=0)/length(y)})},mc.cores=mc))
	return(toto)
	}

mcFeatwG4bg <- function(shufGR=rdlist2,FeatList=FL2,mc=4)
	{
	toto <- simplify2array(mclapply(shufGR,function(y) {sapply(FeatList,function(x) {sum(countOverlaps(x,y,ignore.strand=T)!=0)/length(x)})},mc.cores=mc))
	return(toto)
	}

mcG4dbg <- function(shufGR=rdlist2,FeatList=FL2,mc=4)
	{
	x <- simplify2array(mclapply(shufGR,function(x) {sapply(FeatList,function(y) {1000*sum(countOverlaps(x,y,ignore.strand=T))/sum(as.numeric(width(y)))})},mc.cores=4))
	return(x)
	}
mc_used=4L	
bgG4inFeat <- mcG4inFeatbg(shufGR=rdlist2,FeatList=FL2,mc=mc_used)
bgFeatwG4 <- mcFeatwG4bg(shufGR=rdlist2,FeatList=FL2,mc=mc_used)
bgG4d <- mcG4dbg(shufGR=rdlist2,FeatList=FL2,mc=mc_used)


testbgG4inFeat=apply(apply(bgG4inFeat,2,function(x) (x>=fG4inFeat)),1,sum)
testbgG4d=apply(apply(bgG4d,2,function(x) (x>=G4d)),1,sum)
testbgFeatwG4=apply(apply(bgFeatwG4,2,function(x) (x>=fFeatwG4)),1,sum)	

final_bgG4inFeat <- as.data.frame(cbind(apply(bgG4inFeat,1,mean),apply(bgG4inFeat,1,sd),testbgG4inFeat/length(rdlength)))
names(final_bgG4inFeat) <- c('mean','sd',paste('p-value (n=',length(rdlength),')',sep=''))
final_bgG4d <- as.data.frame(cbind(apply(bgG4d,1,mean),apply(bgG4d,1,sd),testbgG4d/length(rdlength)))
names(final_bgG4d) <- c('mean','sd',paste('p-value (n=',length(rdlength),')',sep=''))
final_bgFeatwG4 <- as.data.frame(cbind(apply(bgFeatwG4,1,mean),apply(bgFeatwG4,1,sd),testbgFeatwG4/length(rdlength)))
names(final_bgFeatwG4) <- c('mean','sd',paste('p-value (n=',length(rdlength),')',sep=''))

### profiles	
## loading the list of G4FS as a GRanges
load('G4H_hg19_1.75_ref.RData')
# removing G4FS from chrM for which no annotation are present in our annotation file
G4data <- G4H_hg19_1.75_ref[seqnames(G4H_hg19_1.75_ref)!='chrM']

load('hg19Anno.RData')
TSS_ranges <-resize(resize(Trans$TSS,fix='end',width=1),fix='center',width=4001)
Feat <- TSS_ranges

pG4 <- G4data[strand(G4data)=='+']
pG4cov <- coverage(pG4)
mG4 <- G4data[strand(G4data)=='-']
mG4cov <- coverage(mG4)
pFeat <- Feat[strand(Feat)=='+']	
mFeat <- Feat[strand(Feat)=='-']	

################################################################################
# this function was developed by PGP Martin
# Pascal.Martin@toulouse.inra.fr
profcomp2_tx=function(covr,gr)
#Compute profiles from a genome-wide coverage and a set of windows
#covr is a genome-wide coverage (typically obtained with the coverage function)
#gr is a GRanges object containing the windows over which to compute the profiles
#Note that names of covr and seqnames of gr must match
{
prof=covr[gr]
prof[strand(gr)=='-']=lapply(prof[strand(gr)=='-'],rev)
return(prof)
}
################################################################################
profFeat_pp <- profcomp2_tx(pG4cov,pFeat)
profFeat_mm <- profcomp2_tx(mG4cov,mFeat)
profFeat_pm <- profcomp2_tx(mG4cov,pFeat)
profFeat_mp <- profcomp2_tx(pG4cov,mFeat)
profFeat_c <- c(profFeat_pp,profFeat_mm)
profFeat_nc <- c(profFeat_pm,profFeat_mp)
prof_G4nearFeat_c <- do.call("rbind",lapply(profFeat_c,function(x){as.vector(x)}))
prof_G4nearFeat_nc <- do.call("rbind",lapply(profFeat_nc,function(x){as.vector(x)}))

xlarg <- 1000L
sp <- 0.5
suff <- 'TSS'
ylableg='Fraction of G4FS'
plot(smooth.spline(apply(prof_G4nearFeat_c[,(2001L-xlarg):(2001L+xlarg)],2,mean),spar=sp,x=((-xlarg):xlarg)),lwd=2,type='l',xlab=paste('Dist to ',suff,' (bp)',sep=''),ylab=ylableg,col='red',ylim=c(0,0.07),cex.axis=1.4,cex.lab=1.4,font=2,font.lab=2)
lines(smooth.spline(apply(prof_G4nearFeat_nc[,(2001L-xlarg):(2001L+xlarg)],2,mean),spar=sp,x=((-xlarg):xlarg)),lwd=2,type='l',col='blue')
abline(v=0,col='purple',lty=2)
legend('topleft',c('G4H1.75','coding','non-coding'),text.col=c('black','red','blue'),bty='n',text.font=2,cex=1.4)