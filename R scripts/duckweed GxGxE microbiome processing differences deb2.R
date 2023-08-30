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
feat.wunfrc.M 	<- read_qza("R inputs/Lemna_GxGxE_coremetrics_deb2_M/weighted_unifrac_distance_matrix.qza")$data #unifrac distances when analyses conducted on only cultured communities
feat.wunfrc.B 	<- read_qza("R inputs/Lemna_GxGxE_coremetrics_deb2/weighted_unifrac_distance_matrix.qza")$data #unifrac distances when analyses conducted on cultured and field communities
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
pdf("pwisedistanceall_deb2.pdf",height=4.5,width=5.5)
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
#get comparisons of cultured communities (culture columns) to the field communities (field rows)
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

pdf("zincANDsymp_deb2.pdf",height=5,width=2.8)
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
pdf("duckweed_pairwise_Vg_vs_microbeUnifrac_deb2.pdf",height=2.25,width=7.75)
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

feat.tab.qz 	<- read_qza("R inputs/Lemna_GxGxE_sepp_keptfeatures_table_deb2M.qza") #is filtered to rm streptophyta
feat.tabFM.qz 	<- read_qza("R inputs/Lemna_GxGxE_sepp_keptfeatures_table_deb2.qza") #is filtered to rm streptophyta
feat.tax.qz 	<- read_qza("R inputs/Lemna_GxGxE_taxonomy_deb2.qza")##NOT FILTERED TO kept features, or even to rm streptophyta (chloroplast, which are removed from the dataset, along with mitochondrial sequences)

Ltab.s	<- feat.tab.qz$data[,order(colnames(feat.tab.qz$data))] #inocula only
LtabFM.s	<- feat.tabFM.qz$data[,order(colnames(feat.tabFM.qz$data))] #both field and inocula
##sequences remaining after ALL filtering are just sums of LtabFM
 sum(LtabFM.s)# deb2:530,137     
 range(colSums(LtabFM.s)) #deb2: 9482 49338		
 mean(colSums(LtabFM.s))# deb2: 26506.85	
 mean(colSums(LtabFM.s)[c(1,3,5,7,9)])# deb2:17,079.6 	  field
 mean(colSums(LtabFM.s)[c(2,4,6,8,10)]) #deb2:36,477.6 	 inocula
 
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
 #deb2:    Ars.M.101   BrBr.M.105    Cam.M.107    Cdv.M.109 HumBPE.M.119    KSR.M.121   MocT.M.129    Rat.M.133    UTM.M.145  WestD.M.148
 #        105           27           35           93           35           52           89           84          121 		39
 
##coverage estimation
DataInfo(LtabFM.s,datatype="abundance")
#there is a warning running on the full table, so this is a check to see that good's coverage is not altered
sapply(1:ncol(LtabFM.s), function(z) DataInfo(as.vector(LtabFM.s[ which(LtabFM.s[,z]>0) ,z]),datatype="abundance")$SC) 
#deb2  0.9987 1.0000 0.9998 1.0000 0.9975 1.0000 1.0000 1.0000 0.9992 1.0000 0.9990 1.0000 0.9989 0.9998 0.9991 1.0000 0.9990 1.0000 0.9986 1.0000


#field-inocula comparison
mtab <- LtabFM.s[,seq(from = 2, to = 20, by =2)]
ftab <- LtabFM.s[,seq(from = 1, to = 19, by =2)]
sum(sign(rowSums(mtab)))#deb2:628
sum(sign(rowSums(ftab)))#deb2:8511
sum( rowSums(ftab) >0 & rowSums(mtab)>0 ) #deb2: 94 		overlap from field and master
FcapanyM <- colSums(ftab[ rowSums(ftab) >0 & rowSums(mtab)>0,])/colSums(ftab[ ,])#percent of field community reads belonging to taxa in any inocula community, 
#deb2 ranges from 0.2-16.6%
FcapbyM <- sapply(1:10, function(z) sum(ftab[ (mtab[,z]) >0,z])/sum(ftab[ ,z])) #percent of field community reads belonging to taxa in the matched inocula community, 
#deb2 ranges from ~0-6%
MpresinanyF <- colSums(mtab[ rowSums(ftab) >0 ,])/colSums(mtab[ ,])#percent of reads in inocula communities present in any field community
#deb2 ranges from 2-97%
colSums(sign(mtab[ rowSums(ftab) >0 ,]))/colSums(sign(mtab[ ,]))# percent of ASV found in any field sample
#    Ars.M.101   BrBr.M.105    Cam.M.107    Cdv.M.109 HumBPE.M.119    KSR.M.121   MocT.M.129    Rat.M.133    UTM.M.145  WestD.M.148 
#   0.14285714   0.22222222   0.17142857   0.07526882   0.34285714   0.21153846   0.24719101   0.29761905   0.12396694   0.30769231 
MpresinF <- sapply(1:10, function(z) sum(mtab[ (ftab[,z]) >0,z])/sum(mtab[ ,z])) # same, but for matched field site, 
#deb2 ranges from very small to 18.8%; 
sapply(1:10, function(z) sum(sign(mtab[ (ftab[,z]) >0,z]))/sum(sign(mtab[ ,z])))#now number of ASV also in matched site / total ASV again
#0.01904762 0.03703704 0.08571429 0.00000000 0.17142857 0.09615385 0.14606742 0.20238095 0.03305785 0.20512821

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
#deb2 140  62 125  19 106 120 144 159  83  98
 sapply(1:10, function(z) sum(unlist(fieldug[z])%in%unlist(fieldug[-z])) )/colSums(sign(ftab))
#   Ars.F.13   BrBr.F.20    Cam.F.22    Cdv.F.26 HumBPE.F.40    KSR.F.43   MocT.F.57    Rat.F.68    UTM.F.90  WestD.F.93 
#  0.05636071  0.18452381  0.09897070  0.33928571  0.15727003  0.09287926  0.09568106  0.09064994  0.12556732  0.09057301 

#continue with more in depth consideration of the taxa across communites
sum(rowSums(mtab)>0 ) # deb2: 198 total in inocula communities
sapply(0:10, function(z) length(which(rowSums(sign(mtab)) == z))) # how many entries in mtab occur in 0, 1, ...all 10 of the inocula?
#deb2: 8417  595   20    8    4    1    0    0    0    0    0
sum(rowSums(sign(mtab)) >= 2 & rowSums(ftab)>0 )  # how many taxa in two or more inocula also occur in a field community?
#deb2: 21 (of 33)
sum(rowSums(sign(mtab)) >= 3 & rowSums(ftab)>0 )  # how many taxa in three or more inocula also occur in a field community?
# deb2: 11 (of 13)
sum(rowSums(sign(mtab)) >= 4 & rowSums(ftab)>0 )  # how many taxa in four or more inocula also occur in a field community?
#deb2: 4 of 5.
sapply(1:10, function(z) sum(ftab[,z]>0 & mtab[,z] > 0) ) # found in mtab column and corresponding field site
#deb2: 2  1  3  0  6  5 13 17  4  8
sapply(1:10, function(z) sum(mtab[,z] > 0 & rowSums(ftab)>0) ) # found in mtab column and any field sample
#deb2: 15  6  6  7 12 11 22 25 15 12
sapply(1:10, function(z) sum(mtab[,z] > 0) ) # found in mtab column
#deb2: 105  27  35  93  35  52  89  84 121  39
unique(LtaxFM.s[rowSums(mtab)>0,6] [which(!LtaxFM.s[rowSums(mtab)>0,6]%in%LtaxFM.s[rowSums(ftab)>0,6])]) #genus names in inocula but not field
#deb2:
# [1] "g__Averyella"         "g__Streptomyces"      "g__Blastobacter"      "g__Achromobacter"     "g__Ancylobacter"     
# [6] "g__Serratia"          "g__Gluconacetobacter" "g__Kaistia"          
unique(LtaxFM.s[rowSums(mtab)>0,6] [which(LtaxFM.s[rowSums(mtab)>0,6]%in%LtaxFM.s[rowSums(ftab)>0,6])]) # genus names in inocula and field
# deb2:
#  [1] "g__"                  NA                     "g__Paucibacter"       "g__Mycobacterium"     "g__Rhodococcus"      
#  [6] "g__Shewanella"        "g__Aeromonas"         "g__Bacillus"          "g__Rahnella"          "g__Pseudomonas"      
# [11] "g__Rheinheimera"      "g__Stenotrophomonas"  "g__Rhodoferax"        "g__Roseomonas"        "g__Devosia"          
# [16] "g__Elstera"           "g__Sphingomonas"      "g__Methylibium"       "g__Flavobacterium"    "g__Dyadobacter"      
# [21] "g__Chryseobacterium"  "g__Agrobacterium"     "g__Variovorax"        "g__Pedobacter"        "g__Acidovorax"       
# [26] "g__Runella"           "g__Caulobacter"       "g__Rhodobacter"       "g__Cellulomonas"      "g__Leptothrix"       
# [31] "g__Brevundimonas"     "g__Pleomorphomonas"   "g__Phenylobacterium"  "g__Rhizobium"         "g__Opitutus"         
# [36] "g__Azospirillum"      "g__Planctomyces"      "g__Tolumonas"         "g__Rubrivivax"        "g__Dechloromonas"    
# [41] "g__Paenibacillus"     "g__Kineosporia"       "g__Enterobacter"      "g__Mycoplana"         "g__Exiguobacterium"  
# [46] "g__Janthinobacterium" "g__Novosphingobium"   "g__Pseudoxanthomonas" "g__Acinetobacter"     "g__Pelomonas"        
# [51] "g__Flectobacillus"    "g__Herbaspirillum"    "g__Yonghaparkia"      "g__Sphingopyxis"      "g__Reyranella"       
# [56] "g__planctomycete"     "g__Hydrogenophaga"   

length(unique(LtaxFM.s[rowSums(mtab)>0,6] [which(LtaxFM.s[rowSums(mtab)>0,6]%in% LtaxFM.s[rowSums(ftab)>0,6])]))-2 # same as above but count, here minus 2 for NA and "g__"
length(unique(LtaxFM.s[rowSums(mtab)>0,6] [which(!LtaxFM.s[rowSums(mtab)>0,6]%in% LtaxFM.s[rowSums(ftab)>0,6])])) 
#deb2: of 55 named genera in inocula communities, all but 8 (47) are in field communities
unique(LtaxFM.s[rowSums(mtab)>0,5] [which(!LtaxFM.s[rowSums(mtab)>0,5]%in% LtaxFM.s[rowSums(ftab)>0,5])])#same but for families
#deb2: none
unique(LtaxFM.s[rowSums(mtab)>0,5] [which(LtaxFM.s[rowSums(mtab)>0,5]%in% LtaxFM.s[rowSums(ftab)>0,5])]) 
#deb2: all 48 families in inocula are found in field communities
GnotinF <- unique(LtaxFM.s[rowSums(mtab)>0,6] [which(!LtaxFM.s[rowSums(mtab)>0,6]%in% LtaxFM.s[rowSums(ftab)>0,6])]) #makes named object of above genus list in inocula but not field communites
colSums(mtab[LtaxFM.s[,6]%in%GnotinF,]) / colSums(mtab) # what is the abundance, by sample, of genera in inocula that are not in any field community
#deb2:    Ars.M.101   BrBr.M.105    Cam.M.107    Cdv.M.109 HumBPE.M.119    KSR.M.121   MocT.M.129    Rat.M.133    UTM.M.145  WestD.M.148 
#		4.013134e-03 0.000000e+00 0.000000e+00 3.466205e-04 0.000000e+00 0.000000e+00 8.088325e-05 0.000000e+00 6.461706e-04 1.749902e-03 
GnotinFtab <- (LtabFM.s[LtaxFM.s[,6]%in%GnotinF,]) # inspect the abundance of taxa for the above
LtaxFM.s[rownames( GnotinFtab)[sapply(1:nrow(GnotinFtab), function(z) any(GnotinFtab[z,]/(colSums(mtab))>0.01) )], ]
 # There are 0 ASV over 1% where the genus is not present in the field

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
pdf("ladpt_microbiome_sim_deb2.pdf",height=6,width=9)
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

