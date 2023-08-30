##################################################################################################################
########## Additional version of R analyses for other ways of processing 16s community sequencing data
######### Please see "duckweed GxGxE analysis.R" for fully annotated script, with all analyses and with 16s data reported in manuscript
##############################################################################################################


library(MCMCglmm)
library(vegan) #mantel test
library(tidyverse)#to allow qiime2R to work
#devtools::install_github("jbisanz/qiime2R")
library(qiime2R)

std.error  <- function(x) {sd(x[!is.na(x)])/sqrt(length(x[!is.na(x)]))} #calculate standard error

range01=function(x){
newnums=(x-min(x,na.rm=T))/(max(x,na.rm=T)-min(x,na.rm=T))
return(newnums)
} #rescales a vector to fall between 0 and 1, used for plotting

bufferX <- function(x,p) { 
	r<- range(x,na.rm=T)
	add <- c(-1,1)*p*(r[2]-r[1])
	return(r+add)
	} #finds the range of a vector and adds a proportion to either side, used for plotting

#read in data files
poptab <- read.csv("R inputs/popdata.csv",header=T)#data on sites at which duckweed lines and associated microbes were collected from
duckdat <- read.csv("R inputs/inputdata.csv",header=T) #experiment data, responses measured and treatments applied
duckdat$colint <- 100-duckdat$meangray #reverse the ImageJ "mean" metric so that it is color intensity rather then whiteness
duckdat$invOD <- 1/duckdat$OD600

#items from main set of scripts that are needed, but unchanged by different microbiome input files
tfitmat <- data.frame(sqmmE=duckdat$sqmmE, invOD=duckdat$invOD, colint=duckdat$colint, agg = duckdat$areapperE , rnd=duckdat$meanround)
mx <- mean(duckdat$x)
mx2 <- mean(duckdat$x^2)
my <- mean(duckdat$y)
my2 <- mean(duckdat$y^2)
mcol <- mean(duckdat$col)
mcol2 <- mean(duckdat$col^2)
mrow <- mean(duckdat$row)
mrow2 <- mean(duckdat$row^2) 
#calculate genotype means in each microbe (averaged across zinc levels) 
DmnsinM <- lapply(1:ncol(tfitmat), function(z) sapply(sort(unique(duckdat$micr)),function(m) 
		tapply(tfitmat[duckdat$micr== m,z],duckdat$duck[duckdat$micr==m],mean,na.rm=T)  ) )# duckweed are rows, microbes are columns each trait is own matrix
#calculate the correlation for every pairwise comparison of genotypic means across microbial community treatments
Bgb.Bu2 <- lapply(DmnsinM, function(z) cor(z)) #tables of duckweed genetic correlations across pairwise microbiome treatment, 1 matrix for each trait. rows and columns are microbe origins
#generally high, but some evidence of dissimilarity of expressed traits/fitness across duckweed genotypes when measured in different microbes -- esp for roundness.



#weighted unifrac distance matrices as output from qiime2 analyses
feat.wunfrc.M 	<- read_qza("R inputs/Lemna_GxGxE_coremetrics_deb_M/weighted_unifrac_distance_matrix.qza")$data #unifrac distances when analyses conducted on only cultured communities
feat.wunfrc.B 	<- read_qza("R inputs/Lemna_GxGxE_coremetrics_deb/weighted_unifrac_distance_matrix.qza")$data #unifrac distances when analyses conducted on cultured and field communities
#for overlapping distances, (e.g. within cultured community pairs), these unifrac distances are nearly perfectly correlated between the two analyses
# re-arrange into matrices
wunfrc.M <- diag(x=0,10)
wunfrc.M[lower.tri(wunfrc.M)] <- feat.wunfrc.M 
wunfrc.M <- Matrix::forceSymmetric(wunfrc.M,uplo="L")
wu.Mladpt<- as.matrix(wunfrc.M[order(attr(feat.wunfrc.M,"Labels")),order(attr(feat.wunfrc.M,"Labels"))])
#re-arrange all samples into matrix
wunfrc.Bx <- diag(x=0,20)
wunfrc.Bx[lower.tri(wunfrc.Bx)] <- feat.wunfrc.B 
wunfrc.Bx <- Matrix::forceSymmetric(wunfrc.Bx,uplo="L")
wunfrc.B <- as.matrix(wunfrc.Bx[order(attr(feat.wunfrc.B,"Labels")),order(attr(feat.wunfrc.B,"Labels"))])

#mantel test, using distances on the same 0,1 scale (wu.Mladpt above are stretched out wrt to msubB)
fsubB <- wunfrc.B[c(1,3,5,7,9,11,13,15,17,19),c(1,3,5,7,9,11,13,15,17,19)]
msubB <- wunfrc.B[c(2,4,6,8,10,12,14,16,18,20),c(2,4,6,8,10,12,14,16,18,20)]
mantel(fsubB,msubB)
#not significant

#make pairwise unifrac distance figure for all inocula and field communities
wbk <- colorRampPalette(c(rgb(1,1,1),rgb(0,0,0)))
pdf("pwisedistanceall_deb10.pdf",height=4.5,width=5.5)
layout(matrix(1:2,ncol=2),widths = c(1,0.125) )
par(mar=c(6,6,1,1))
orderdistm <- c(1,3,5,7,9,11,13,15,17,19,2,4,6,8,10,12,14,16,18,20)[c(order(poptab$Nimperv),order(poptab$Nimperv)+10)]
image(as.matrix(wunfrc.B[orderdistm,orderdistm]),xaxt="n",yaxt="n",xlab="",ylab="",col=wbk(100),zlim=c(0,1))
abline(h=.5); abline(v=.5)
axis(side=1, at = seq(from=0,to=1,length.out=20), labels=round(rep(sort(poptab$Nimperv), times=2), digits=1), las=2 )
axis(side=2, at = seq(from=0,to=1,length.out=20), labels=round(rep(sort(poptab$Nimperv), times=2), digits=1), las=2 )
mtext("Inocula",side=1,at=0.75,line=3.25); mtext("Field",side=1,at=0.25, line=3.25)
mtext("Inocula",side=2,at=0.75,line=3.25); mtext("Field",side=2,at=0.25, line=3.25)
mtext("% permeable surface area",side=1,line=4.5); mtext("% permeable surface area",side=2,line=4.5)
par(mar=c(7,0,4,2))
image(matrix(seq(from=0,to=1,length.out=20),ncol=20), xaxt="n",yaxt="n",xlab="",ylab="",col=wbk(100),zlim=c(0,1))
axis(side=4,at=c(0,1))
mtext("Unifrac distance",side=4)
dev.off()

#get pairwise distance from each field site community to each cultured community; add vector to data for distance between duckweed host's field community and the experimentally inoculated community
fseq <- seq(from=1,to=19,by=2)
mseq <- seq(from = 2, to = 20, by =2)
wunf.m2f <- wunfrc.B[fseq,mseq]
wuncdat <- data.frame(unf_fm = as.vector(wunf.m2f), popM = as.character(rep(poptab$pop,each=10)), popF = as.character(rep(poptab$pop,times=10)), matched = as.vector(matchedsmpl <- diag(1,nrow=10)))
duckdat$distItoF <- sapply(1:nrow(duckdat), function(z) wuncdat$unf_fm[ wuncdat$popM==duckdat$micr[z] &  wuncdat$popF ==duckdat$duck[z]   ])


######################################################
### do sympatric combinations result in greater fitness and does contamination disrupt this?
###only for models with unifrac distance between inoculated and field community
summary(MCMCglmm(sqmmE ~ x + I(x^2) + y + I(y^2) + col + row + I(row^2) +
 distItoF + numzinc + distItoF:numzinc, 
 family="gaussian", data=duckdat,verbose=FALSE,pr=TRUE,nitt=11000, thin=10, burnin=1000)) 
summary(MCMCglmm(sqmmE ~ x + I(x^2) + y + I(y^2) + col + row + I(row^2) +
 distItoF + numzinc , 
 family="gaussian", data=duckdat,verbose=FALSE,pr=TRUE,nitt=11000, thin=10, burnin=1000)) 
#distItoF  sig, and one other term must be removed. zinc & zincXdistItoF models apparently DIC equivalent. w/o zinc:distItoF and with main effect of zinc seems to have better behaviour
summary(MCMCglmm(invOD ~ x + I(x^2) + y + I(col^2)+ row + I(row^2) + 
 distItoF + numzinc + distItoF:numzinc, 
 family="gaussian", data=duckdat,verbose=FALSE,pr=TRUE,nitt=11000, thin=10, burnin=1000)) 
summary(MCMCglmm(invOD ~ x + I(x^2) + y + I(col^2)+ row + I(row^2) + 
 distItoF, 
 family="gaussian", data=duckdat,verbose=FALSE,pr=TRUE,nitt=11000, thin=10, burnin=1000)) 
#only distItoF, still pos.

#refit best models and extract predictions for reporting and figures
AE.if <-(MCMCglmm(sqmmE ~ x + I(x^2) + y + I(y^2) + col + row + I(row^2) +
 distItoF + numzinc , 
 family="gaussian", data=duckdat,verbose=FALSE,pr=TRUE,nitt=101000, thin=10, burnin=1000)) 
Aifsol <- AE.if$Sol
if.s <- seq(from=min(duckdat$distItoF),to=max(duckdat$distItoF),length.out=1000)
aif.ci <- sapply(if.s, function(z) HPDinterval(as.mcmc(
		Aifsol[,1] + Aifsol[,2]*mx + Aifsol[,3]*mx2 + Aifsol[,4]*my + Aifsol[,5]*my2 + Aifsol[,6]*mcol + Aifsol[,7]*mrow + Aifsol[,8]*mrow2 +
				 Aifsol[,9]*z + Aifsol[,10]*0.5  ),prob=0.95) )
aif.mn <-  sapply(if.s, function(z) 	mean(Aifsol[,1] + Aifsol[,2]*mx + Aifsol[,3]*mx2 + Aifsol[,4]*my + Aifsol[,5]*my2 + Aifsol[,6]*mcol + Aifsol[,7]*mrow + Aifsol[,8]*mrow2 +
				 Aifsol[,9]*z + Aifsol[,10]*0.5 ) )

invOD.if <- MCMCglmm(invOD ~ x + I(x^2) + y + I(col^2)+ row + I(row^2) + 
 distItoF ,
 family="gaussian", data=duckdat,verbose=FALSE,pr=TRUE,nitt=101000, thin=10, burnin=1000) 
ODifsol <- invOD.if$Sol
micr.if.ci <- sapply(if.s, function(z)    HPDinterval(as.mcmc(
	ODifsol[,1] + ODifsol[,2]*mx + ODifsol[,3]*mx2 + ODifsol[,4]*my + ODifsol[,5]*mcol2+ ODifsol[,6]*mrow + ODifsol[,7]*mrow2 + 
 		ODifsol[,8]*z ),prob=0.95) ) 
micr.if.mn <- sapply(if.s, function(z)    mean(
	ODifsol[,1] + ODifsol[,2]*mx + ODifsol[,3]*mx2 + ODifsol[,4]*my + ODifsol[,5]*mcol2+ ODifsol[,6]*mrow + ODifsol[,7]*mrow2 + 
 		ODifsol[,8]*z ) ) 

#calculate means for figure
mH2A.m  <- tapply(duckdat$sqmmE, paste(duckdat$duck,duckdat$micr,duckdat$numzinc),mean)
mH2A.se <- tapply(duckdat$sqmmE, paste(duckdat$duck,duckdat$micr,duckdat$numzinc),std.error)
mH2.m  <- tapply(duckdat$invOD, paste(duckdat$duck,duckdat$micr,duckdat$numzinc),mean)
mH2.se <- tapply(duckdat$invOD, paste(duckdat$duck,duckdat$micr, duckdat$numzinc),std.error)
mH2if <- tapply(duckdat$distItoF, paste(duckdat$duck,duckdat$micr, duckdat$numzinc),mean)

pdf("zincANDsymp_deb10.pdf",height=5,width=2.8)
par(mar=c(1.5,2.5,1.5,0))
par(oma=c(2.5,1.5,0,1.5))
layout(matrix(1:2,ncol=1))
plot(1~1,pch=NA,xlim=bufferX(if.s,0.05),ylim=bufferX(c(mH2A.m),0.1),ylab="",xlab="",main="" )
	polygon (x= c(if.s, rev(if.s) ), y= c( aif.ci[1,], rev(aif.ci[2,]) ) ,col=rgb(0,0,0,alpha=0.2),border=NA)
	points(mH2A.m ~ mH2if,pch=1,col=rgb(0,0,0),cex=0.5)
 	lines(aif.mn~if.s)
	mtext(expression("Duckweed frond area, mm"^2),side=2,line=2)	
plot(1~1,pch=NA,xlim=bufferX(if.s,0.05),ylim=bufferX(1/mH2.m,0.1),ylab="",xlab="",main="" )
	polygon (x= c(if.s, rev(if.s) ), y=1/ c( micr.if.ci[1,], rev(micr.if.ci[2,]) ) ,col=rgb(0,0,0,alpha=0.2),border=NA)
	points(1/mH2.m ~ mH2if,pch=1,col=rgb(0,0,0),cex=0.5)
 	lines(1/micr.if.mn~if.s)
	mtext("Microbial density, OD 600 nm",side=2,line=2)	
	mtext("Field-inocula Unifrac distance",side=1,line=2.5)
dev.off()


#########################################################
########## CONSISTENCY OF FITNESS, TRAIT EXPRESSSION
#######################################################

##use the matrix to test whether microbial community distance explains the variation in heritability (note for clones, heritability is covariance in genotype means)
UD.s <- seq(from = 0, to =max(wu.Mladpt), length.out=1000)
#rearrange unifrac distance and pairwise genotype mean correlations into a data frame
UnifracGCdat <- lapply(1:5, function(z) data.frame(GC = Bgb.Bu2[[z]][lower.tri(Bgb.Bu2[[z]])==T], Unifrac = wu.Mladpt[lower.tri(wu.Mladpt)==T]) )
#fit a model for each trait, correlation is dependent variable, independent is unifrac distance, extract predictions and HPDI
UnifracGCmods <-	lapply(1:5, function(z) MCMCglmm(GC ~ Unifrac, data=UnifracGCdat[[z]], burnin=1000 , nitt =101000 , thin=10, verbose=F ) )
UnifracGCpreds <-	lapply(1:5, function(z) sapply(UD.s, function(D) 		mean(UnifracGCmods[[z]]$Sol[,1] + UnifracGCmods[[z]]$Sol[,2]*D) ) )
UnifracGCHPDIs <-	lapply(1:5, function(z) sapply(UD.s, function(D) HPDinterval(UnifracGCmods[[z]]$Sol[,1] + UnifracGCmods[[z]]$Sol[,2]*D,prob=0.95 ) ) )

##plot results 
tnames<- c("Duckweed area","Microbial density","Color Intensity","Aggregation","Roundness")
##main figure
pdf("duckweed_pairwise_Vg_vs_microbeUnifrac_deb10.pdf",height=2.25,width=7.75)
par(mfrow=c(1,5))
par(mar=c(4,1,2.5,1.5))
par(oma=c(0,4,0,0))
for(z in 1:5){
 	plot(UnifracGCdat[[z]]$GC ~ UnifracGCdat[[z]]$Unifrac,ylab="",xlab="",ylim=c(ifelse(z!=5,-0.1,-0.5),1.1),xlim=c(0,0.8))
	if(z==3){mtext("Inoculated microbiome pairwise Unifrac distance",line=2.5,side=1)}
	if(z==1){mtext("Cross-microbiome",line=3,side=2)}
	if(z==1){mtext("trait similarity",line=2,side=2)}
	if(z!=5){lines(UnifracGCpreds[[z]]~UD.s)}
	abline(h=0,lty=3)
	if(z!=5){polygon(c(UD.s, rev(UD.s)),c(UnifracGCHPDIs[[z]][1,],rev(UnifracGCHPDIs[[z]][2,])),col=rgb(0,0,0,alpha=0.3))}
	mtext(tnames[z],side=3,line=0.5)
}
dev.off()


##############
###INVESTIGATE MICROBIOME COMPOSITION
############
##read in full taxonomic data
library(iNEXT)

feat.tab.qz 	<- read_qza("R inputs/Lemna_GxGxE_sepp_keptfeatures_table_debM.qza") #is filtered to rm streptophyta
feat.tabFM.qz 	<- read_qza("R inputs/Lemna_GxGxE_sepp_keptfeatures_table_deb.qza") #is filtered to rm streptophyta
feat.tax.qz 	<- read_qza("R inputs/Lemna_GxGxE_taxonomy_deb.qza")##NOT FILTERED TO kept features, or even to rm streptophyta (chloroplast, which are removed from the dataset, along with mitochondrial sequences)

Ltab.s	<- feat.tab.qz$data[,order(colnames(feat.tab.qz$data))] #inocula only
LtabFM.s	<- feat.tabFM.qz$data[,order(colnames(feat.tabFM.qz$data))] #both field and inocula
##sequences remaining after ALL filtering are just sums of LtabFM
 sum(LtabFM.s)# deb10:  508,642
 range(colSums(LtabFM.s)) #deb10: 	9421 49166	
 mean(colSums(LtabFM.s))# deb10:	25432.1
 mean(colSums(LtabFM.s)[c(1,3,5,7,9)])# deb10: 	15276  field
 mean(colSums(LtabFM.s)[c(2,4,6,8,10)]) #deb10: 36379.2	 inocula
 
#split tax info into columns
 feat.tax 	<- feat.tax.qz$data ##ALL taxa - this is not filtered to taxa that occur in samples, or even to remove streptophyta
 taxlevs <- sapply(1:nrow(feat.tax), function(z) strsplit(as.character(feat.tax$Taxon[z]),split="; "))
 tax <- matrix(NA,ncol=7,nrow=nrow(feat.tax))
 for(i in 1:nrow(tax)){ tax[i,1:length(taxlevs[[i]])]<- taxlevs[[i]]}
 rownames(tax) <- as.character(feat.tax$Feature.ID)
##subset
 Ltax.s <- 	 	tax[sapply(1:nrow(Ltab.s), function(z) which(rownames(tax)==rownames(Ltab.s)[z])),]
 LtaxFM.s <- 	 	tax[sapply(1:nrow(LtabFM.s), function(z) which(rownames(tax)==rownames(LtabFM.s)[z])),]#
 Ltax <- Ltax.s
 Ltab <- Ltab.s[rownames(Ltax.s),]
 rm(feat.tab.qz,feat.tax.qz, feat.tabFM.qz, feat.tax, taxlevs,tax)
#reset colsums to be 1
 Lprp <- Ltab
 for(i in 1:ncol(Ltab)){Lprp[,i] <- Ltab[,i]/sum(Ltab[,i])}
 LprpFM <- LtabFM.s
 for(i in 1:ncol(LtabFM.s)){LprpFM[,i] <- LtabFM.s[,i]/sum(LtabFM.s[,i])}
#organized the same
 sum(rownames(Ltax)==rownames(Lprp) )
 dim(Lprp)
#quick stats for ms
 colSums(sign(Lprp))
 #deb10:  Ars.M.101   BrBr.M.105    Cam.M.107    Cdv.M.109 HumBPE.M.119    KSR.M.121   MocT.M.129    Rat.M.133    UTM.M.145  WestD.M.148 
#           33           15           14           20           16           14           32           31           27           17 

  
##coverage estimation
DataInfo(LtabFM.s,datatype="abundance")
#there is a warning running on the full table, so this is a check to see that good's coverage is not altered
sapply(1:ncol(LtabFM.s), function(z) DataInfo(as.vector(LtabFM.s[ which(LtabFM.s[,z]>0) ,z]),datatype="abundance")$SC) 
#deb10:[1] 0.9990 1.0000 0.9997 1.0000 0.9982 1.0000 1.0000 1.0000 0.9996 1.0000 0.9991 1.0000 0.9993 1.0000 0.9993 1.0000 0.9990 1.0000 0.9990 1.0000


#field-inocula comparison
mtab <- LtabFM.s[,seq(from = 2, to = 20, by =2)]
ftab <- LtabFM.s[,seq(from = 1, to = 19, by =2)]
sum(sign(rowSums(mtab)))#deb10: 169
sum(sign(rowSums(ftab)))#deb10: 2284
sum( rowSums(ftab) >0 & rowSums(mtab)>0 ) #deb10: 86		overlap from field and master
FcapanyM <- colSums(ftab[ rowSums(ftab) >0 & rowSums(mtab)>0,])/colSums(ftab[ ,])#percent of field community reads belonging to taxa in any inocula community, 
#deb10 ranges from 0.3-18.9%
FcapbyM <- sapply(1:10, function(z) sum(ftab[ (mtab[,z]) >0,z])/sum(ftab[ ,z])) #percent of field community reads belonging to taxa in the matched inocula community, 
#deb10 ranges from ~0-6.9%
MpresinanyF <- colSums(mtab[ rowSums(ftab) >0 ,])/colSums(mtab[ ,])#percent of reads in inocula communities present in any field community
#deb 10: 2.1-97%
colSums(sign(mtab[ rowSums(ftab) >0 ,]))/colSums(sign(mtab[ ,]))# percent of ASV found in any field sample
#   Ars.M.101   BrBr.M.105    Cam.M.107    Cdv.M.109 HumBPE.M.119    KSR.M.121   MocT.M.129    Rat.M.133    UTM.M.145  WestD.M.148 
#    0.3939394    0.4000000    0.4285714    0.3500000    0.6875000    0.7142857    0.6250000    0.7741935    0.5555556    0.6470588 
MpresinF <- sapply(1:10, function(z) sum(mtab[ (ftab[,z]) >0,z])/sum(mtab[ ,z])) # same, but for matched field site, 
#deb10 ranges from very small to 18.9%; 
sapply(1:10, function(z) sum(sign(mtab[ (ftab[,z]) >0,z]))/sum(sign(mtab[ ,z])))#now number of ASV also in matched site / total ASV again
#  [1] 0.06060606 0.06666667 0.21428571 0.00000000 0.37500000 0.35714286 0.37500000 0.54838710 0.14814815 0.41176471

#same types of comparisons between field or inocula communities, but at genus level 
genbyS <- lapply(mseq, function(z) LtaxFM.s[LtabFM.s[,z]>0,6])
genbySrmNA <- lapply(genbyS, function(z) z[!z%in%c("g__",NA)])
geninf <- lapply(1:10, function(z) genbySrmNA[[z]] [which(genbySrmNA[[z]] %in% LtaxFM.s[ LtabFM.s[,fseq[z]]>0,6 ] )] )
MorGpresinF <- sapply(1:10, function(z)  #  (excluding taxa not identified to this level) percent of inocula community where the same genus is present in the matched field community
			sum( mtab[ftab[, z ] >0 | LtaxFM.s[,6]%in%geninf[[z]]    , z]) / 
			sum( mtab[,z ])  )
Fgens <-   (LtaxFM.s[ rowSums(ftab) >0 ,6 ]) [!(LtaxFM.s[ rowSums(ftab) >0 ,6 ]%in%c(NA,"g__"))]
MorGpresinanyF <- sapply(1:10, function(z)  #same but for the any field community
			sum( mtab[ rowSums(ftab) >0 | LtaxFM.s[,6]%in%Fgens  ,z]) /
			sum(mtab[,z] ) )
fieldug <- lapply(1:10, function(z) unique( LtaxFM.s[ ftab[,z]>0 ,6]))
sapply(1:10, function(z) sum(unlist(fieldug[z])%in%unlist(fieldug[-z])) )# how many taxa in one field site are found in any of the other 9 (note this still includes "g__" and "NA")
#deb10  88  48  84  13  73  70 105 112  62  71
 sapply(1:10, function(z) sum(unlist(fieldug[z])%in%unlist(fieldug[-z])) )/colSums(sign(ftab))
#deb10:    Ars.F.13   BrBr.F.20    Cam.F.22    Cdv.F.26 HumBPE.F.40    KSR.F.43   MocT.F.57    Rat.F.68    UTM.F.90  WestD.F.93 
#   0.1156373   0.2008368   0.1546961   0.3513514   0.2067989   0.1386139   0.1707317   0.1540578   0.1577608   0.1621005 

#continue with more in depth consideration of the taxa across communites
sum(rowSums(mtab)>0 ) # deb10:  total in inocula communities
sapply(0:10, function(z) length(which(rowSums(sign(mtab)) == z))) # how many entries in mtab occur in 0, 1, ...all 10 of the inocula?
#deb10: 2198  138   18    8    4    1    0    0    0    0    0
sum(rowSums(sign(mtab)) >= 2 & rowSums(ftab)>0 )  # how many taxa in two or more inocula also occur in a field community?
#deb10: 21 (of 31)
sum(rowSums(sign(mtab)) >= 3 & rowSums(ftab)>0 )  # how many taxa in three or more inocula also occur in a field community?
# deb10: 11 (of 13)
sum(rowSums(sign(mtab)) >= 4 & rowSums(ftab)>0 )  # how many taxa in four or more inocula also occur in a field community?
#deb10: 4 of 5.
sapply(1:10, function(z) sum(ftab[,z]>0 & mtab[,z] > 0) ) # found in mtab column and corresponding field site
#deb10: 2  1  3  0  6  5 12 17  4 7  
sapply(1:10, function(z) sum(mtab[,z] > 0 & rowSums(ftab)>0) ) # found in mtab column and any field sample
#deb10: 13  6  6  7 11 10 20 24 15 11    
sapply(1:10, function(z) sum(mtab[,z] > 0) ) # found in mtab column
#deb10: 33 15 14 20 16 14 32 31 27 17    
unique(LtaxFM.s[rowSums(mtab)>0,6] [which(!LtaxFM.s[rowSums(mtab)>0,6]%in%LtaxFM.s[rowSums(ftab)>0,6])]) #genus names in inocula but not field
#deb10:
# [1] "g__Tolumonas"         "g__Bacillus"          "g__Rahnella"          "g__Enterobacter"      "g__Gluconacetobacter"
# [6] "g__Kaistia"           "g__Streptomyces"     
unique(LtaxFM.s[rowSums(mtab)>0,6] [which(LtaxFM.s[rowSums(mtab)>0,6]%in%LtaxFM.s[rowSums(ftab)>0,6])]) # genus names in inocula and field
# deb10:
#  [1] "g__Pseudomonas"       "g__Rheinheimera"      "g__Agrobacterium"     "g__"                  "g__Shewanella"       
#  [6] "g__Elstera"           "g__Devosia"           NA                     "g__Methylibium"       "g__Acinetobacter"    
# [11] "g__Variovorax"        "g__Aeromonas"         "g__Rhodoferax"        "g__Azospirillum"      "g__Rhodobacter"      
# [16] "g__Chryseobacterium"  "g__Rhizobium"         "g__Mycoplana"         "g__Hydrogenophaga"    "g__Rhodococcus"      
# [21] "g__planctomycete"     "g__Sphingopyxis"      "g__Flavobacterium"    "g__Acidovorax"        "g__Sphingomonas"     
# [26] "g__Pseudoxanthomonas" "g__Caulobacter"       "g__Janthinobacterium" "g__Pedobacter"        "g__Stenotrophomonas" 
# [31] "g__Exiguobacterium"   "g__Pelomonas"         "g__Reyranella"        "g__Herbaspirillum"    "g__Dyadobacter"      
# [36] "g__Novosphingobium"   "g__Kineosporia"       "g__Flectobacillus"    "g__Runella"           "g__Paucibacter"      
# [41] "g__Phenylobacterium" 
length(unique(LtaxFM.s[rowSums(mtab)>0,6] [which(LtaxFM.s[rowSums(mtab)>0,6]%in% LtaxFM.s[rowSums(ftab)>0,6])]))-2 # same as above but count, here minus 2 for NA and "g__"
length(unique(LtaxFM.s[rowSums(mtab)>0,6] [which(!LtaxFM.s[rowSums(mtab)>0,6]%in% LtaxFM.s[rowSums(ftab)>0,6])])) 
#deb10: of 39 named genera in inocula communities, all but 7 (32) are in field communities
unique(LtaxFM.s[rowSums(mtab)>0,5] [which(!LtaxFM.s[rowSums(mtab)>0,5]%in% LtaxFM.s[rowSums(ftab)>0,5])])#same but for families
#deb10: [1] "f__Bacillaceae"       "f__Xanthobacteraceae" "f__Streptomycetaceae"
unique(LtaxFM.s[rowSums(mtab)>0,5] [which(LtaxFM.s[rowSums(mtab)>0,5]%in% LtaxFM.s[rowSums(ftab)>0,5])]) 
#deb10:  [1] "f__Pseudomonadaceae"    "f__[Chromatiaceae]"     "f__Rhizobiaceae"        "f__Sinobacteraceae"     "f__Aeromonadaceae"     
#  [6] "f__Shewanellaceae"      "f__Rhodospirillaceae"   "f__Hyphomicrobiaceae"   "f__Comamonadaceae"      "f__Moraxellaceae"      
# [11] "f__Saprospiraceae"      "f__Chitinophagaceae"    "f__Enterobacteriaceae"  "f__"                    "f__C111"               
# [16] "f__Rhodobacteraceae"    "f__[Weeksellaceae]"     "f__Caulobacteraceae"    "f__Nocardiaceae"        "f__Opitutaceae"        
# [21] "f__Pirellulaceae"       "f__Sphingomonadaceae"   "f__Flavobacteriaceae"   "f__Phyllobacteriaceae"  "f__Xanthomonadaceae"   
# [26] "f__Carnobacteriaceae"   "f__Oxalobacteraceae"    "f__Sphingobacteriaceae" "f__[Exiguobacteraceae]" "f__Cellulomonadaceae"  
# [31] NA                       "f__Cytophagaceae"       "f__Kineosporiaceae"    
GnotinF <- unique(LtaxFM.s[rowSums(mtab)>0,6] [which(!LtaxFM.s[rowSums(mtab)>0,6]%in% LtaxFM.s[rowSums(ftab)>0,6])]) #makes named object of above genus list in inocula but not field communites
colSums(mtab[LtaxFM.s[,6]%in%GnotinF,]) / colSums(mtab) # what is the abundance, by sample, of genera in inocula that are not in any field community
#deb10:   Ars.M.101   BrBr.M.105    Cam.M.107    Cdv.M.109 HumBPE.M.119    KSR.M.121   MocT.M.129    Rat.M.133    UTM.M.145  WestD.M.148 
#		 0.808648253  0.013325046  0.000000000  0.134719044  0.000000000  0.000000000  0.000000000  0.000000000  0.008566925  0.009037047 
GnotinFtab <- (LtabFM.s[LtaxFM.s[,6]%in%GnotinF,]) # inspect the abundance of taxa for the above
LtaxFM.s[rownames( GnotinFtab)[sapply(1:nrow(GnotinFtab), function(z) any(GnotinFtab[z,]/(colSums(mtab))>0.01) )], ]
#deb 10 : [1] "k__Bacteria"    "p__Firmicutes"  "c__Bacilli"     "o__Bacillales"  "f__Bacillaceae" "g__Bacillus"    NA              

#make summary table of abundance at family level, all samples. including taxa not in inocula
ufams <- sort(na.omit(c(unique(LtaxFM.s [,5]))))
famsum <- matrix(NA, ncol = ncol(LprpFM),nrow = length(ufams) )
rows <- list()
for(i in 1:length(ufams)) {
	rows[[i]] <- which(LtaxFM.s[,5] == ufams[i])
	subtax <- LprpFM[rows[[i]],]
	if(!is.null(nrow(subtax)) ){	
	famsum[i,] <- colSums(subtax)
	} else{
	famsum[i,] <- subtax	
	}
}
famsum[1,]<- famsum[1,] + colSums( LprpFM[ is.na(LtaxFM.s[,5]),] )  #leaving the NA in ufams fails to include NA, merging with "f__"
prettyfamnames <- c("Unidentified",sapply(2:length(ufams), function(z) strsplit(ufams[z],split="__")[[1]][[2]]))
lowabundfam <-  !sapply(1:nrow(famsum),function(z) any(famsum[z,]>0.05))
hiabundFsum <- famsum[!lowabundfam,]
simplerfamsum <- rbind(hiabundFsum[order(rowSums(hiabundFsum),decreasing=T),],colSums(famsum[lowabundfam,]) )
simplerfamsumnames <- c(prettyfamnames[!lowabundfam][order(rowSums(hiabundFsum),decreasing=T)], "Other")
simplerfamsumnames[which(simplerfamsumnames=="[Chromatiaceae]")] <- "Chromatiaceae"
simplerfamsumnames[which(simplerfamsumnames=="[Exiguobacteraceae]")] <- "Exiguobacteraceae"

 #set color
library(Polychrome)
load( file="R inputs/colorpal.RData")
siteorder <- c(order(poptab$Nimperv)*2-1,order(poptab$Nimperv)*2)

uniDdat <- data.frame(uniD=wuncdat[wuncdat$matched==1,1], nimp=poptab$Nimperv)
uniDm <- (MCMCglmm(uniD ~ nimp , family="gaussian", data=uniDdat,verbose=FALSE,pr=TRUE,nitt=101000, thin=10, burnin=1000)) #all kept
uniDs<-uniDm$Sol
n.s <- seq(from = min(uniDdat$nimp),to=max(uniDdat$nimp),length.out=1000)
uniD.p <- sapply(n.s, function(z) mean(uniDs[,1]+uniDs[,2]*z) )
uniD.ci <- sapply(n.s, function(z) HPDinterval(as.mcmc(uniDs[,1]+uniDs[,2]*z),prob=.95) )

##main microbiome composition figure
pdf("ladpt_microbiome_sim_deb10.pdf",height=6,width=9)
layout(matrix(c(1,1,2,3,4,2),nrow=2, byrow = TRUE),heights=c(3,3),widths=c(1.1,.75,.4))
par(oma=c(0,0,2,0))
par(mar=c(6,4,0,0))
bg <- barplot(simplerfamsum[,siteorder], col= newpal,ylab="",xlab="")#,names.arg=rep(round(poptab$Nimperv),times=2) )
axis(  at=c( (bg[5]+bg[6])/2, (bg[15]+bg[16])/2), side=3, labels = c("Field","Inocula") ,lty=0,line=-0.5,cex.axis=1.25)
axis( at = bg, side=1, labels = rep(sort(round(poptab$Nimperv)),times=2),lty=0, line=-1)
#axis( at = bg[c(1,10,11,20)], side=1, labels = rep(c("most urban","most rural"),times=2),lty=0, line=0,las=2)
abline(v = (bg[10]+bg[11])/2, lty=1, lwd=3)
mtext("% permeable surface area",side=1,line=4.5)
mtext("a.",side=3, adj = -0.075, line = 0.5) #,adj=-0.12,line=-0.5)
mtext("Proportion",side=2,line=2)
par(mar=c(0,0,0,0))
plot(rep(1,times=3)~c(1:3), pch=NA,bty="n",xaxt="n",yaxt="n",ylab="",xlab="")
legend(x=0.9,y=1.4,legend=simplerfamsumnames,fill=newpal,bty="n",ncol=1,cex=1.1)
legend(x=0.9,y=0.8,legend=c("ASV, matched site","Genus, matched","ASV, any site","Genus, any site"),fill=c("blue","gray","blue","gray"),density=c(100,100,50,50),bty="n",ncol=1,cex=1)
par(mar=c(6,4,1,2))
bg<- barplot(100*t(cbind(MpresinF,MorGpresinF,MpresinanyF,MorGpresinanyF))[,order(poptab$Nimperv)],beside=T,xlab="",
 		col=c("blue","gray"),density=c(100,100,50,50), ,names.arg=rep("",times=10),ylim=c(0,103))
axis( at = (bg[2,]+bg[3,])/2, side=1, labels = sort(round(poptab$Nimperv)),lty=0, line=-1)
mtext("% represented in field",side=2,line=2)
mtext("% permeable surface area",side=1,line=4.5)
mtext("b.",side=3,adj=-0.125,line = 0.25)
plot(wuncdat[wuncdat$matched==1,1]~poptab$Nimperv,pch=NA,ylab="",xlab="", ,cex.axis=1.25,ylim=c(0.25,0.9))#,
polygon(c(n.s,rev(n.s)),c(uniD.ci[1,],rev(uniD.ci[2,])),col=rgb(0,0,0,alpha=0.25),border=NA)
#text(,,"p<0.05")
lines(uniD.p~n.s)
points(uniDdat$uniD~uniDdat$nimp,pch=16,cex=2)
# mtext("Field vs inoc. Unifrac dist.",side = 2, line = 2.5)
 mtext("Field-inocula Unifrac distance",side = 2, line = 2.5,at=0.525)
 mtext("% permeable surface area", side =1, line = 4.5)
mtext("c.",side=3,adj=-0.2,line = 0.25)
 dev.off()

