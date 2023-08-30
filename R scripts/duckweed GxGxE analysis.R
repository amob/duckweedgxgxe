##################################################################################################################
########## R Analyses for manuscript
########## Microbiomes will affect host evolution: Experimental demonstration of microbiome effects on trait expression, heritability, and fitness of hosts
##############################################################################################################

#########################################################
########## LIBRARIES, FUNCTIONS, READ IN DATA, INSPECT UNIFRAC DATA
#######################################################

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

#weighted unifrac distance matrices as output from qiime2 analyses
feat.wunfrc.M 	<- read_qza("R inputs/Lemna_GxGxE_coremetrics_den_M/weighted_unifrac_distance_matrix.qza")$data #unifrac distances when analyses conducted on only cultured communities
feat.wunfrc.B 	<- read_qza("R inputs/Lemna_GxGxE_coremetrics_den/weighted_unifrac_distance_matrix.qza")$data #unifrac distances when analyses conducted on cultured and field communities
#for unifrac distances included in both objects/analyses, (e.g. within cultured community pairs), distances are nearly perfectly correlated

# re-arrange into matrices. First for subset with just cultured communities
wunfrc.M <- diag(x=0,10)
wunfrc.M[lower.tri(wunfrc.M)] <- feat.wunfrc.M 
wunfrc.M <- Matrix::forceSymmetric(wunfrc.M,uplo="L")
wu.Mladpt<- as.matrix(wunfrc.M[order(attr(feat.wunfrc.M,"Labels")),order(attr(feat.wunfrc.M,"Labels"))])
#re-arrange all samples into matrix
wunfrc.Bx <- diag(x=0,20)
wunfrc.Bx[lower.tri(wunfrc.Bx)] <- feat.wunfrc.B 
wunfrc.Bx <- Matrix::forceSymmetric(wunfrc.Bx,uplo="L")
wunfrc.B <- as.matrix(wunfrc.Bx[order(attr(feat.wunfrc.B,"Labels")),order(attr(feat.wunfrc.B,"Labels"))])

#mantel test for pairwise distance matrix of field communities and pairwise distance matrix of cultured communities
#using distances on the same 0,1 scale (wu.Mladpt above are stretched out wrt to msubB)
fsubB <- wunfrc.B[c(1,3,5,7,9,11,13,15,17,19),c(1,3,5,7,9,11,13,15,17,19)]
msubB <- wunfrc.B[c(2,4,6,8,10,12,14,16,18,20),c(2,4,6,8,10,12,14,16,18,20)]
mantel(fsubB,msubB)
#not significant

#make pairwise unifrac distance figure for all inocula and field communities
wbk <- colorRampPalette(c(rgb(1,1,1),rgb(0,0,0)))
pdf("pwisedistanceall_den.pdf",height=4.5,width=5.5)
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

#get pairwise distance from each field site community to each cultured community
#add vector to data for distance between duckweed host's field community and the experimentally inoculated community
fseq <- seq(from=1,to=19,by=2)
mseq <- seq(from = 2, to = 20, by =2)
cor(as.vector(wunfrc.B[mseq,mseq]), as.vector( wunfrc.M)) #note these are not identical: #but almost. p = 0.99985
#get comparisons of cultured communities (culture columns) to the field communities (field rows)
#add vector to experimental dataframe
wunf.m2f <- wunfrc.B[fseq,mseq]
wuncdat <- data.frame(unf_fm = as.vector(wunf.m2f), popM = as.character(rep(poptab$pop,each=10)), popF = as.character(rep(poptab$pop,times=10)), matched = as.vector(matchedsmpl <- diag(1,nrow=10)))
duckdat$distItoF <- sapply(1:nrow(duckdat), function(z) wuncdat$unf_fm[ wuncdat$popM==duckdat$micr[z] &  wuncdat$popF ==duckdat$duck[z]   ])

#########################################################
########## SPACE MODELS and INSPECTING NORMALITY, RESIDUALS
#######################################################

###normality test of response variables:
shapiro.test(duckdat$sqmmE) #W = 0.9889
shapiro.test(duckdat$lnOD600) #W = 0.96623
shapiro.test(duckdat$colint) # W = 0.97972
shapiro.test(duckdat$meanround) # W = 0.98363
shapiro.test(duckdat$areapperE) #W = 0.9971
shapiro.test(1/duckdat$OD600) #W = 0.95073
shapiro.test(duckdat$OD600) #W = 0.91474

####
# fit space effects first, remove n.s. terms one by one. 
#models shown are full space effects model and best fit space effects model, for fitness and each trait.
#only for the first two space models is the full stepwise procedure shown, however, the stepwise result is repeatable for all models,
	#it can be replicated by the code user if desired by removing least significant terms in turn
	#elsewhere in the code we similarly show only the full and best model, but in all places, the stepwise result is repeatable
summary(MCMCglmm(sqmmE ~ x + I(x^2) + y + I(y^2) + col + I(col^2)+ row + I(row^2), 
 family="gaussian", data=duckdat,verbose=FALSE,pr=TRUE,nitt=11000, thin=10, burnin=1000)) #first remove col^2
summary(MCMCglmm(sqmmE ~ x + I(x^2) + y + I(y^2) + col + row + I(row^2), 
 family="gaussian", data=duckdat,verbose=FALSE,pr=TRUE,nitt=11000, thin=10, burnin=1000)) #then keep the rest
#inverted OD 
summary(MCMCglmm(sqmmE ~ x + I(x^2) + y + I(y^2) + col + I(col^2)+ row + I(row^2), 
 family="gaussian", data=duckdat,verbose=FALSE,pr=TRUE,nitt=11000, thin=10, burnin=1000)) #first remove y^2
summary(MCMCglmm(sqmmE ~ x + I(x^2) + y + col + I(col^2)+ row + I(row^2), 
 family="gaussian", data=duckdat,verbose=FALSE,pr=TRUE,nitt=11000, thin=10, burnin=1000)) #then remove col
summary(MCMCglmm(sqmmE ~ x + I(x^2) + y + I(col^2)+ row + I(row^2), 
 family="gaussian", data=duckdat,verbose=FALSE,pr=TRUE,nitt=11000, thin=10, burnin=1000)) #then test remove y
summary(MCMCglmm(invOD ~ x + I(x^2) +I(col^2)+ row + I(row^2),
 family="gaussian", data=duckdat,verbose=FALSE,pr=TRUE,nitt=11000, thin=10, burnin=1000))#replace y, as DIC is lower with y
 #last model is worse than second last, second last used
#aggregation
summary(MCMCglmm(areapperE ~ x + I(x^2) + y + I(y^2) + col + I(col^2)+ row + I(row^2), 
 family="gaussian", data=duckdat,verbose=FALSE,pr=TRUE,nitt=11000, thin=10, burnin=1000)) #full spatial effects is best fit
#colint
summary(MCMCglmm(colint ~ x + I(x^2) + y + I(y^2) + col + I(col^2)+ row + I(row^2), 
 family="gaussian", data=duckdat,verbose=FALSE,pr=TRUE,nitt=11000, thin=10, burnin=1000)) #full spatial effects is best fit
#roundness
summary(MCMCglmm(meanround ~ x + I(x^2) + y + I(y^2) + col + I(col^2)+ row + I(row^2), 
 family="gaussian", data=duckdat,verbose=FALSE,pr=TRUE,nitt=11000, thin=10, burnin=1000)) 
summary(MCMCglmm(meanround ~ x + I(x^2) + y + I(y^2) + col + I(col^2) , 
 family="gaussian", data=duckdat,verbose=FALSE,pr=TRUE,nitt=11000, thin=10, burnin=1000))

###This section for inspecting residuals of space models, as a check of data distribution assumptions
#fit models with more iterations, compare predictions and residuals
set.seed(100)
SMODsqmmE <- (MCMCglmm(sqmmE ~ x + I(x^2) + y + I(y^2) + col + row + I(row^2) ,
  family="gaussian", data=duckdat,verbose=FALSE,pr=TRUE,nitt=101000, thin=10, burnin=1000))  
SMODareapperE <- (MCMCglmm(areapperE ~ x + I(x^2) + y + I(y^2) + col + I(col^2)+ row + I(row^2), 
 family="gaussian", data=duckdat,verbose=FALSE,pr=TRUE,nitt=101000, thin=10, burnin=1000))
SMODcolint <- (MCMCglmm(colint ~ x + I(x^2) + y + I(y^2) + col + I(col^2)+ row + I(row^2), 
 family="gaussian", data=duckdat,verbose=FALSE,pr=TRUE,nitt=101000, thin=10, burnin=1000))
SMODmeanround <- (MCMCglmm(meanround ~ x + I(x^2) + y + I(y^2) + col + I(col^2) , 
 family="gaussian", data=duckdat,verbose=FALSE,pr=TRUE,nitt=101000, thin=10, burnin=1000))
SMODinvOD <- (MCMCglmm(invOD ~ x + I(x^2) + y + I(col^2)+ row + I(row^2), 
 family="gaussian", data=duckdat,verbose=FALSE,pr=TRUE,nitt=101000, thin=10, burnin=1000))
predSMODsqmmE <- predict(SMODsqmmE)
resSMODsqmmE <- duckdat$sqmmE - predSMODsqmmE
predSMODareapperE <- predict(SMODareapperE)
resSMODareapperE <- duckdat$areapperE - predSMODareapperE
predSMODcolint <- predict(SMODcolint)
resSMODcolint <- duckdat$colint - predSMODcolint
predSMODmeanround <- predict(SMODmeanround)
resSMODmeanround <- duckdat$meanround - predSMODmeanround
predSMODinvOD <- predict(SMODinvOD)
resSMODinvOD <- duckdat$invOD - predSMODinvOD
pdf("residulsVSpred_space.pdf",height=3,width=12) #residual vs predicted values for space models
par(mfrow=c(1,5))
par(mar=c(4,2,2,0))
par(oma=c(0,2,0,0))
plot(resSMODsqmmE~predSMODsqmmE,xlab="",ylab="")
	mtext("Residuals",side=2,line=2)
	mtext("Area",side=3,line=0.5)
plot(resSMODinvOD~predSMODinvOD,xlab="",ylab="")
	mtext("Inv. Optical Density",side=3,line=0.5) #using inverse
plot(resSMODareapperE~predSMODareapperE,xlab="",ylab="")
	mtext("Aggregation",side=3,line=0.5)
	mtext("Predicted value",side=1,line=2.5)
plot(resSMODcolint~predSMODcolint,xlab="",ylab="")
	mtext("Color Intensity",side=3,line=0.5)
plot(resSMODmeanround~predSMODmeanround,xlab="",ylab="")
	mtext("Roundness",side=3,line=0.5)
dev.off()
pdf("predVSact_space.pdf",height=3,width=12)#predicted vs actual values for space models
par(mfrow=c(1,5))
par(mar=c(4,2,2,0))
par(oma=c(0,2,0,0))
plot(predSMODsqmmE~duckdat$sqmmE,xlab="",ylab="")
	mtext("Predicted value",side=2,line=2)
	mtext("Area",side=3,line=0.5)
plot(predSMODinvOD~duckdat$invOD,xlab="",ylab="")
	mtext("Inv. Optical Density",side=3,line=0.5)
plot(predSMODareapperE~duckdat$areapperE,xlab="",ylab="")
	mtext("Aggregation",side=3,line=0.5)
	mtext("Actual value",side=1,line=2.5)
plot(predSMODcolint~duckdat$colint,xlab="",ylab="")
	mtext("Color Intensity",side=3,line=0.5)
plot(predSMODmeanround~duckdat$meanround,xlab="",ylab="")
	mtext("Roundness",side=3,line=0.5)
dev.off()

############################################
#### Sources of variance (general); Figure 2
############################################

#Make a matrix of traits, fitness. 
#Then, make a duplicate matrix where optical density data is the measurement, rather than transformed scale
tfitmat <- data.frame(sqmmE=duckdat$sqmmE, invOD=duckdat$invOD, colint=duckdat$colint, agg = duckdat$areapperE , rnd=duckdat$meanround)
tfitmatregOD <- data.frame(sqmmE=duckdat$sqmmE, OD=duckdat$OD600, colint=duckdat$colint, agg = duckdat$areapperE , rnd=duckdat$meanround)

###Sums of Squares
ssbyvar <- function(response,category.vec){ #sums of squares function
	means <- tapply(response,category.vec,mean,na.rm=T) #take the means by category
	ssresid <- sum(sapply(sort(unique(category.vec)), function(z) sum( (response[category.vec==z] - means[names(means)==z])^2,na.rm=T ))) #square of difference of each datapoint from its associated treatment mean (residual variation)
	sstot <- sum((response-mean(response,na.rm=T))^2,na.rm=T) #square of difference of each datapoint from the grand mean (total variation)
	sst <- (sstot-ssresid) # total variation - residual variation = treatment variation
	return(sst/sstot) # treatment variance as a fraction of total variation
	}
pSSduck <- sapply(1:ncol(tfitmatregOD), function(z) ssbyvar(tfitmatregOD[,z],duckdat$duck))
pSSmicr <- sapply(1:ncol(tfitmatregOD), function(z) ssbyvar(tfitmatregOD[,z],duckdat$micr))
pSSzinc <- sapply(1:ncol(tfitmatregOD), function(z) ssbyvar(tfitmatregOD[,z],duckdat$zinc))
pSSduckmicr <- sapply(1:ncol(tfitmatregOD), function(z) ssbyvar(tfitmatregOD[,z],paste(duckdat$duck,duckdat$micr)))
pSSduckzinc <- sapply(1:ncol(tfitmatregOD), function(z) ssbyvar(tfitmatregOD[,z],paste(duckdat$duck,duckdat$zinc)))
pSSmicrzinc <- sapply(1:ncol(tfitmatregOD), function(z) ssbyvar(tfitmatregOD[,z],paste(duckdat$micr,duckdat$zinc)))
pSSduckmicrzinc <- sapply(1:ncol(tfitmatregOD), function(z) ssbyvar(tfitmatregOD[,z],paste(duckdat$duck,duckdat$micr,duckdat$zinc)))
pSS <- cbind(pSSduck, pSSmicr, pSSzinc, pSSduckmicr, pSSduckzinc, pSSmicrzinc, pSSduckmicrzinc)
pSSc <- cbind(pSSduck, pSSmicr, pSSzinc, 
			pSSduckmicr-pSSduck-pSSmicr, pSSduckzinc-pSSduck-pSSzinc, pSSmicrzinc-pSSmicr-pSSzinc, 
			pSSduckmicrzinc-pSSduck-pSSmicr-pSSzinc - (pSSduckmicr-pSSduck-pSSmicr) - (pSSduckzinc-pSSduck-pSSzinc) - (pSSmicrzinc-pSSmicr-pSSzinc) )

ylabv <- c( expression("G"[host]),expression("G"[microbe]),expression("E"),
		expression("G"[host]~X~"G"[microbe]),expression("G"[host]~X~"E"),expression("G"[microbe]~X~"E"),
expression("G"[host]~X~"G"[microbe]~X~"E"))
tnames<- c("Duckweed frond area","Microbial density","Color Intensity","Aggregation","Roundness")
pdf("bar-sources-of-var.pdf",height=4.75,width=3.75)
par(mar=c(9,4,1,1))
bar <-barplot(t(pSSc*100),beside=F,col=c(rgb(1,0,0),rgb(0,0,1),rgb(0.95,0.85,0),rgb(0.75,0,0.75),rgb(1,0.5,0),rgb(0,0.75,0),rgb(0.33,0.33,0.33)),ylim=c(0,50),xlim=c(0,6))
	axis(side=1,at=bar,labels=tnames,las=2)
	mtext(side=2,"%Variance Explained",line=2)
	legend(-0.4,51,ylabv[1:3], fill=c(rgb(1,0,0),rgb(0,0,1),rgb(0.95,0.85,0)),bty="n")
	legend(2,51,ylabv[4:7], fill=c(rgb(0.75,0,0.75),rgb(1,0.5,0),rgb(0,0.75,0),rgb(0.33,0.33,0.33)),bty="n")
dev.off()

###Heritability

Hmxe <- sapply(1:5, function(t) sapply(c(0,1), function (z) sapply(poptab$pop, function(m) 
						var(tapply(tfitmatregOD[duckdat$numzinc==z & duckdat$micr==m,t],duckdat$duck[duckdat$numzinc==z & duckdat$micr==m],mean,na.rm=T))/var(tfitmatregOD[duckdat$numzinc==z & duckdat$micr==m, t] ,na.rm=T) ) ) )
#Heritability (in the host/duckweed), calculated separately across zinc environment and microbiome treatments
##for each trait (1:5), each zinc level (0,1), and each microbiome: calculate the trait means for each duckweed genotype, the variance among these means, and divide by the total variance among datapoints within this section (trait, zinc level, microbiome treatment) of the data

Hgxe <- sapply(1:5, function(t) sapply(c(0,1), function (z) sapply(poptab$pop, function(d) 
						var(tapply(tfitmatregOD[duckdat$numzinc==z & duckdat$duck==d,t],duckdat$micr[duckdat$numzinc==z & duckdat$duck==d],mean,na.rm=T))/var(tfitmatregOD[duckdat$numzinc==z & duckdat$duck==d, t] ,na.rm=T) ) ) )
#'Heritability in the microbiome', calculated separately across zinc environment and microbiome treatments
##for each trait (1:5), each zinc level (0,1), and each duckweed genotype: calculate the trait means for each microbiome treatment, the variance among these means, and divide by the total variance among datapoints within this section (trait, zinc level, duckweed genotype) of the data

#calculate the same heritability and 'heritability in the microbiome' metrics, but on 10,000 randomized datasets
randHmxe <- array(NA,dim=c(20,5,10000))
randHgxe <- array(NA,dim=c(20,5,10000))
#takes some time.
set.seed(100)
for(i in 1:10000){
	RtfitmatregOD <- tfitmatregOD[sample(1:nrow(tfitmatregOD),nrow(tfitmatregOD)),]
	randHmxe[,,i] <- sapply(1:5, function(t) sapply(c(0,1), function (z) sapply(poptab$pop, function(m) 
						var(tapply(RtfitmatregOD[duckdat$numzinc==z & duckdat$micr==m,t],duckdat$duck[duckdat$numzinc==z & duckdat$micr==m],mean,na.rm=T))/var(RtfitmatregOD[duckdat$numzinc==z & duckdat$micr==m, t] ,na.rm=T) ) ) )
	randHgxe[,,i] <- sapply(1:5, function(t) sapply(c(0,1), function (z) sapply(poptab$pop, function(d) 
						var(tapply(RtfitmatregOD[duckdat$numzinc==z & duckdat$duck==d,t],duckdat$micr[duckdat$numzinc==z & duckdat$duck==d],mean,na.rm=T))/var(RtfitmatregOD[duckdat$numzinc==z & duckdat$duck==d, t] ,na.rm=T) ) ) )
} #columns are traits, rows are microbes, repeated 2x, first 10 rows low zinc, next 10 high zinc


cols <- c("red","orange","gold1","forestgreen","blue","blueviolet","magenta","brown","black","gray")

pdf("HeritabilityAcrossMZ.pdf",height=4,width=8)
layout(matrix(1:12,ncol=6,byrow=T),widths=c(1,1,1,1,1,1.8))
par(oma=c(4,5,2,0))
par(mar=c(0,0,1,0))
for(trt in 1:5){
plot(Hmxe[c(order(poptab$Nimperv),10+order(poptab$Nimperv)),trt]~rep(c(1,2),each=10),col=cols,pch=NA,xaxt="n",ylim=c(0,0.65),ylab="",xlab="",yaxt="n",xlim=c(0.5,2.5)) 
	abline(h=max(sapply(1:20, function(mz) HPDinterval(as.mcmc(randHmxe[mz,1,])))),lty=3 )
	points(Hmxe[c(order(poptab$Nimperv),10+order(poptab$Nimperv)),trt]~rep(c(1,2),each=10),col=cols,pch=16)
	arrows(1,Hmxe[order(poptab$Nimperv),trt],x1=2,y1=Hmxe[10+order(poptab$Nimperv),trt],length=0,col=cols)
	if(trt==1){mtext(side=3,"Duckweed",line=1.5)}
	if(trt==1){mtext(side=3,"frond area",line=0.5)}
	if(trt==2){mtext(side=3,"Microbial",line=1.5)}
	if(trt==2){mtext(side=3,"density",line=0.5)}
	if(trt==3){mtext(side=3,"Color",line=1.5)}
	if(trt==3){mtext(side=3,"intensity",line=0.5)}
	if(trt==4){mtext(side=3,"Aggregation",line=0.5)}
	if(trt==5){mtext(side=3,"Roundness",line=0.5)}
	if(trt==1){mtext(side=2,"in duckweed",line=2.25)}
	if(trt==1){mtext(side=2,"Heritability",line=3.5)}
	if(trt==1){axis(side=2)}
}
par(mar=c(0,0,0,0))
plot(1:10~c(1:10),pch=NA,xlab="",ylab="",bty="n",xaxt="n",yaxt="n")
	legend(1,9,
		c("Arsandco Pond","Broken Bridge","Cambelleville","Cedarvale Pond","Humber Bay Park East","Koffler Scientific Reserve","Moccasin Trail","Rattray","University of Toronto Missisauga","West Duffins")[order(poptab$Nimperv)],
		fill=cols,bty="n")
	mtext("Microbial community",side=3,at=1,line=-2,adj=0,cex=0.85)
par(mar=c(0,0,1,0))
for(trt in 1:5){
plot(Hgxe[c(order(poptab$Nimperv),10+order(poptab$Nimperv)),trt]~rep(c(1,2),each=10),col=cols,pch=NA,xaxt="n",ylim=c(0,0.65),ylab="",xlab="",yaxt="n",xlim=c(0.5,2.5)) 
	abline(h=max(sapply(1:20, function(mz) HPDinterval(as.mcmc(randHgxe[mz,1,])))),lty=3 )
	points(Hgxe[c(order(poptab$Nimperv),10+order(poptab$Nimperv)),trt]~rep(c(1,2),each=10),col=cols,pch=16)
	arrows(1,Hgxe[order(poptab$Nimperv),trt],x1=2,y1=Hgxe[10+order(poptab$Nimperv),trt],length=0,col=cols)
	axis(side=1,at=c(1,2),labels=c("Low","High"))
	if(trt==3){mtext(side=1,"Zinc level",line=2.5)}
	if(trt==1){mtext(side=2,"in the microbiome",line=2.25)}
	if(trt==1){mtext(side=2,"'Heritability'",line=3.5)}
	if(trt==1){axis(side=2)}
}
par(mar=c(0,0,0,0))
plot(1:10~c(1:10),pch=NA,xlab="",ylab="",bty="n",xaxt="n",yaxt="n")
	legend(1,9,
		c("Arsandco Pond","Broken Bridge","Cambelleville","Cedarvale Pond","Humber Bay Park East","Koffler Scientific Reserve","Moccasin Trail","Rattray","University of Toronto Missisauga","West Duffins")[order(poptab$Nimperv)],
		fill=cols,bty="n")
	mtext("Duckweed genotype",side=3,at=1,line=-2,adj=0,cex=0.85)
dev.off()



############################################
###test whether zince or source site explain Ghost, Gmicr, E, and interactive effects
###includes models evaluating effects for permeable surface area at site source (supplemental figures)
### and models evaluating potential local adaptation between duckweeds and microbes (Figure 3)
################################################

#for each trait/fitness dependent variable there are two models, one is the full model tested, the second is the model with n.s. terms removed (from stepwise procedure)
#see space models for example code of full stepwise procedure, results are repeatable.
summary(MCMCglmm(sqmmE ~ x + I(x^2) + y + I(y^2) + col + row + I(row^2) +
 duckNImpv + micrNImpv + numzinc + duckNImpv:numzinc + micrNImpv:numzinc + duckNImpv:micrNImpv + duckNImpv:micrNImpv:numzinc, 
 family="gaussian", data=duckdat,verbose=FALSE,pr=TRUE,nitt=11000, thin=10, burnin=1000)) 
summary(MCMCglmm(sqmmE ~ x + I(x^2) + y + I(y^2) + col + row + I(row^2) +
 duckNImpv + numzinc , 
 family="gaussian", data=duckdat,verbose=FALSE,pr=TRUE,nitt=11000, thin=10, burnin=1000)) 
summary(MCMCglmm(invOD ~ x + I(x^2) + y + I(col^2) + row + I(row^2) +
 duckNImpv + micrNImpv + numzinc + duckNImpv:numzinc + micrNImpv:numzinc + duckNImpv:micrNImpv + duckNImpv:micrNImpv:numzinc, 
 family="gaussian", data=duckdat,verbose=FALSE,pr=TRUE,nitt=11000, thin=10, burnin=1000))
 summary(MCMCglmm(invOD ~ x + I(x^2) + y + I(col^2) + row + I(row^2) +
 duckNImpv + micrNImpv  + duckNImpv:micrNImpv , 
 family="gaussian", data=duckdat,verbose=FALSE,pr=TRUE,nitt=11000, thin=10, burnin=1000))
summary(MCMCglmm(colint ~ x + I(x^2) + y + I(y^2) + col + I(col^2)+ row + I(row^2) +
 duckNImpv + micrNImpv + numzinc + duckNImpv:numzinc + micrNImpv:numzinc + duckNImpv:micrNImpv + duckNImpv:micrNImpv:numzinc, 
 family="gaussian", data=duckdat,verbose=FALSE,pr=TRUE,nitt=11000, thin=10, burnin=1000))
summary(MCMCglmm(colint ~ x + I(x^2) + y + I(y^2) + col + I(col^2)+ row + I(row^2) +
 duckNImpv + duckNImpv:micrNImpv , 
 family="gaussian", data=duckdat,verbose=FALSE,pr=TRUE,nitt=11000, thin=10, burnin=1000))
summary(MCMCglmm(areapperE ~ x + I(x^2) + y + I(y^2) + col + I(col^2)+ row + I(row^2)+
 duckNImpv + micrNImpv + numzinc + duckNImpv:numzinc + micrNImpv:numzinc + duckNImpv:micrNImpv + duckNImpv:micrNImpv:numzinc, 
 family="gaussian", data=duckdat,verbose=FALSE,pr=TRUE,nitt=11000, thin=10, burnin=1000))
summary(MCMCglmm(areapperE ~ x + I(x^2) + y + I(y^2) + col + I(col^2)+ row + I(row^2)+
 duckNImpv +  duckNImpv:micrNImpv + duckNImpv:micrNImpv:numzinc, 
 family="gaussian", data=duckdat,verbose=FALSE,pr=TRUE,nitt=11000, thin=10, burnin=1000))
summary(MCMCglmm(meanround ~ x + I(x^2) + y + I(y^2) + col + I(col^2) +
 duckNImpv + micrNImpv + numzinc + duckNImpv:numzinc + micrNImpv:numzinc + duckNImpv:micrNImpv + duckNImpv:micrNImpv:numzinc, 
 family="gaussian", data=duckdat,verbose=FALSE,pr=TRUE,nitt=11000, thin=10, burnin=1000))
 #none sig, leaving in duckNImpv improves DIC (negative slope), but n.s.
#####
###re-fit models with increased iterations for reporting results
set.seed(100) #run the following 4 models at once after setting seed for consistency
urbanduckmod <- MCMCglmm(sqmmE ~ x + I(x^2) + y + I(y^2) + col + row + I(row^2) + duckNImpv + numzinc, family="gaussian", data=duckdat,verbose=FALSE,pr=TRUE,nitt=101000, thin=10, burnin=1000)
urbanmicrmod <- MCMCglmm(invOD ~ x + I(x^2) + y + I(col^2)+ row + I(row^2) +  duckNImpv + micrNImpv + duckNImpv:micrNImpv, family="gaussian", data=duckdat,verbose=FALSE,pr=TRUE,nitt=101000, thin=10, burnin=1000)
urbanappEmod <- (MCMCglmm(areapperE ~ x + I(x^2) + y + I(y^2) + col + I(col^2)+ row + I(row^2)+ duckNImpv +  duckNImpv:micrNImpv + duckNImpv:micrNImpv:numzinc,  family="gaussian", data=duckdat,verbose=FALSE,pr=TRUE,nitt=101000, thin=10, burnin=1000))
urbancolintmod <- MCMCglmm(colint ~ x + I(x^2) + y + I(y^2) + col + I(col^2)+ row + I(row^2) + duckNImpv + duckNImpv:micrNImpv , family="gaussian", data=duckdat,verbose=FALSE,pr=TRUE,nitt=101000, thin=10, burnin=1000)
summary(urbanduckmod);summary(urbanmicrmod); summary(urbanappEmod); summary(urbancolintmod)
#results in supplementary tables
#roundness excluded because no treatment effects significant in this set of models

# create vectors of space effects for use in fitted solutions in figures, extracting predictions and HPDI for predictions.
#use average values for space effects, one model each for fitness and all traits but roundness
NI.s <- seq(from=min(poptab$Nimperv),to = max(poptab$Nimperv),length.out=1000)
NI.r <- range(poptab$Nimperv)
#in predictions, substitute average values for space effects
mx <- mean(duckdat$x)
mx2 <- mean(duckdat$x^2)
my <- mean(duckdat$y)
my2 <- mean(duckdat$y^2)
mcol <- mean(duckdat$col)
mcol2 <- mean(duckdat$col^2)
mrow <- mean(duckdat$row)
mrow2 <- mean(duckdat$row^2) 

duckUSol <- urbanduckmod$Sol
az.dd.ci <- sapply(1:1000, function(z) HPDinterval(as.mcmc(
		duckUSol[,1] + duckUSol[,2]*mx + duckUSol[,3]*mx2 + duckUSol[,4]*my + duckUSol[,5]*my2 + duckUSol[,6]*mcol + duckUSol[,7]*mrow + duckUSol[,8]*mrow2 +
				 duckUSol[,9]*NI.s[z] + duckUSol[,10]*0.5  ),prob=0.95) )
azpe.dd.prd <- data.frame(fit = sapply(1:1000, function(z) 	mean(duckUSol[,1] + duckUSol[,2]*mx + duckUSol[,3]*mx2 + duckUSol[,4]*my + duckUSol[,5]*my2 + duckUSol[,6]*mcol + duckUSol[,7]*mrow + duckUSol[,8]*mrow2 +
				 duckUSol[,9]*NI.s[z] + duckUSol[,10]*0.5 ) ), 
						lwr =  az.dd.ci[1,], upr = az.dd.ci[2,] )
						zincmn <- tapply(duckdat$sqmmE,duckdat$numzinc,mean,na.rm=T)
meanspacesqmm <- mean(	duckUSol[,1] + duckUSol[,2]*mx + duckUSol[,3]*mx2 + duckUSol[,4]*my + duckUSol[,5]*my2 + duckUSol[,6]*mcol + duckUSol[,7]*mrow + duckUSol[,8]*mrow2 + duckUSol[,10]*0.5 )
meanspaceSQ <- duckdat$sqmmE - 
			sapply(1:nrow(duckdat), function(z) mean(duckUSol[,1] + duckUSol[,2]*duckdat$x[z] + duckUSol[,3]*(duckdat$x[z]^2) + duckUSol[,4]*duckdat$y[z] + duckUSol[,5]*(duckdat$y[z]^2) + duckUSol[,6]*duckdat$col[z] + duckUSol[,7]*duckdat$row[z] + duckUSol[,8]*(duckdat$row[z]^2)  + duckUSol[,9]*duckdat$numzinc[z] )) +
			meanspacesqmm
azmns <- tapply(meanspaceSQ,duckdat$duckNImpv,mean,na.rm=T)
azses <- tapply(meanspaceSQ,duckdat$duckNImpv,std.error)
#
micrUSol <- urbanmicrmod$Sol
micr.dI.ci <- lapply(NI.r, function(micrP) sapply(NI.s, function(duckP)    1/HPDinterval(as.mcmc(
	micrUSol[,1] + micrUSol[,2]*mx + micrUSol[,3]*mx2 + micrUSol[,4]*my + micrUSol[,5]*mcol2+ micrUSol[,6]*mrow + micrUSol[,7]*mrow2 + 
 		micrUSol[,8]*duckP + micrUSol[,9]*micrP + micrUSol[,10]*duckP*micrP ),prob=0.95) )) 
micr.dI.mn <- lapply(NI.r, function(micrP) sapply(NI.s, function(duckP) 1/mean(
	micrUSol[,1] + micrUSol[,2]*mx + micrUSol[,3]*mx2 + micrUSol[,4]*my + micrUSol[,5]*mcol2+ micrUSol[,6]*mrow + micrUSol[,7]*mrow2 + 
 		micrUSol[,8]*duckP + micrUSol[,9]*micrP + micrUSol[,10]*duckP*micrP)  ) ) 
meanspaceODinv <- mean(	micrUSol[,1] + micrUSol[,2]*mx + micrUSol[,3]*mx2 + micrUSol[,4]*my + micrUSol[,5]*mcol2+ micrUSol[,6]*mrow + micrUSol[,7]*mrow2 )
meanspaceOD <- 1/(duckdat$invOD - 
			sapply(1:nrow(duckdat), function(z) mean(micrUSol[,1] + micrUSol[,2]*duckdat$x[z] + micrUSol[,3]*(duckdat$x[z]^2) + micrUSol[,4]*duckdat$y[z] + micrUSol[,5]*(duckdat$col[z]^2)+ micrUSol[,6]*duckdat$row[z] + micrUSol[,7]*(duckdat$row[z]^2) )) +
			meanspaceODinv)
odmn <- tapply(meanspaceOD,paste(duckdat$micrNImpv,duckdat$duckNImpv),mean,na.rm=T)
niMmn <- tapply(duckdat$micrNImpv,paste(duckdat$micrNImpv,duckdat$duckNImpv),mean,na.rm=T)
niDmn <- tapply(duckdat$duckNImpv,paste(duckdat$micrNImpv,duckdat$duckNImpv),mean,na.rm=T)
#
aggUSol <- urbanappEmod$Sol #note that for aggregation this extraction only predictions at high zinc. differences are n.s. at low zinc.
agg.dI.ci <- lapply(NI.r, function(micrP) sapply(NI.s, function(duckP)    HPDinterval(as.mcmc(
	aggUSol[,1] + aggUSol[,2]*mx + aggUSol[,3]*mx2 + aggUSol[,4]*my + aggUSol[,5]*my2 + aggUSol[,6]*mcol + aggUSol[,7]*mcol2+ aggUSol[,8]*mrow + aggUSol[,9]*mrow2 + 
 		aggUSol[,10]*duckP + aggUSol[,11]*duckP*micrP + aggUSol[,12]*duckP*micrP*1 ),prob=0.95) )) 
agg.dI.mn <- lapply(NI.r, function(micrP) sapply(NI.s, function(duckP) mean(
	aggUSol[,1] + aggUSol[,2]*mx + aggUSol[,3]*mx2 + aggUSol[,4]*my + aggUSol[,5]*my2 + aggUSol[,6]*mcol + aggUSol[,7]*mcol2+ aggUSol[,8]*mrow + aggUSol[,9]*mrow2 + 
 		aggUSol[,10]*duckP + aggUSol[,11]*duckP*micrP + aggUSol[,12]*duckP*micrP*1 ) )) 
meanspaceagg <- mean(	aggUSol[,1] + aggUSol[,2]*mx + aggUSol[,3]*mx2 + aggUSol[,4]*my + aggUSol[,5]*my2 + aggUSol[,6]*mcol + aggUSol[,7]*mcol2+ aggUSol[,8]*mrow + aggUSol[,9]*mrow2  )
meanspaceaggR <- duckdat$areapperE - 
			sapply(1:nrow(duckdat), function(z) mean(aggUSol[,1] + aggUSol[,2]*duckdat$x[z] + aggUSol[,3]*(duckdat$x[z]^2) + aggUSol[,4]*duckdat$y[z] + aggUSol[,5]*(duckdat$y[z]^2)+ aggUSol[,6]*duckdat$col[z] + aggUSol[,7]*(duckdat$col[z]^2)+  aggUSol[,8]*duckdat$row[z] + aggUSol[,9]*(duckdat$row[z]^2) )) +
			meanspaceagg
aggmn <- tapply(meanspaceaggR[duckdat$numzinc==1],paste(duckdat$micrNImpv[duckdat$numzinc==1],duckdat$duckNImpv[duckdat$numzinc==1]),mean,na.rm=T)
#
colUSol<- urbancolintmod$Sol
col.dI.ci <- lapply(NI.r, function(micrP) sapply(NI.s, function(duckP)    HPDinterval(as.mcmc(
	colUSol[,1] + colUSol[,2]*mx + colUSol[,3]*mx2 + colUSol[,4]*my + colUSol[,5]*my2 + colUSol[,6]*mcol + colUSol[,7]*mcol2+ colUSol[,8]*mrow + colUSol[,9]*mrow2 + 
 		colUSol[,10]*duckP + colUSol[,11]*duckP*micrP ),prob=0.95) )) 
col.dI.mn <- lapply(NI.r, function(micrP) sapply(NI.s, function(duckP) mean(
	colUSol[,1] + colUSol[,2]*mx + colUSol[,3]*mx2 + colUSol[,4]*my + colUSol[,5]*my2 + colUSol[,6]*mcol + colUSol[,7]*mcol2+ colUSol[,8]*mrow + colUSol[,9]*mrow2 + 
 		colUSol[,10]*duckP + colUSol[,11]*duckP*micrP ) )) 
meanspacecol <- mean(	colUSol[,1] + colUSol[,2]*mx + colUSol[,3]*mx2 + colUSol[,4]*my + colUSol[,5]*my2 + colUSol[,6]*mcol + colUSol[,7]*mcol2+ colUSol[,8]*mrow + colUSol[,9]*mrow2  )
meanspacecolR <- duckdat$colint - 
			sapply(1:nrow(duckdat), function(z) mean(colUSol[,1] + colUSol[,2]*duckdat$x[z] + colUSol[,3]*(duckdat$x[z]^2) + colUSol[,4]*duckdat$y[z] + colUSol[,5]*(duckdat$y[z]^2)+ colUSol[,6]*duckdat$col[z] + colUSol[,7]*(duckdat$col[z]^2)+  colUSol[,8]*duckdat$row[z] + colUSol[,9]*(duckdat$row[z]^2) )) +
			meanspacecol
colmn <- tapply(meanspacecolR,paste(duckdat$micrNImpv,duckdat$duckNImpv),mean,na.rm=T)

##Supplementary figure for significant effects in models
pdf("GxGxEurbanness.pdf",height=5.5,width=5) #points are treatment averaged space residuals (and zinc residuals for aggregation)
par(mar=c(0,3.5,0,1))
par(oma=c(7,1,2,0))
layout(matrix(1:4,nrow=2),widths=c(1,1))
suppressWarnings(plot(1~1,pch=NA,xlim=bufferX(NI.s,.05),ylim=range(c(azmns+azses,azmns-azses)),ylab="",xlab="",main="",xaxt="n" )) #R dislikes the formula 1~1
	polygon (x= c(NI.s, rev(NI.s) ), y= c( azpe.dd.prd[,2], rev(azpe.dd.prd[,3]) ) ,col=rgb(0,0,0,alpha=0.2),border=NA)
	points(azmns ~ sort(poptab$Nimperv),pch=16)
	lines(azpe.dd.prd[,1]~NI.s)
	arrows(x0 =sort(poptab$Nimperv), y0 = azmns-azses ,y1 = azmns+azses,length=0)
	mtext(expression("Duckweed frond area, mm"^2),side=2,line=2)
suppressWarnings(plot(1~1,pch=NA,xlim=bufferX(NI.s,0.05),ylim=bufferX(colmn,0.0),ylab="",xlab="",main="" )) #R dislikes the formula 1~1
	polygon (x= c(NI.s, rev(NI.s) ), y= c( col.dI.ci[[1]][1,], rev(col.dI.ci[[1]][2,]) ) ,col=rgb(0.75,0.55,0.35,alpha=0.2),border=NA)
  	lines(col.dI.mn[[1]]~NI.s,col=rgb(0.75,0.55,0.35))
 	polygon (x= c(NI.s, rev(NI.s) ), y= c( col.dI.ci[[2]][1,], rev(col.dI.ci[[2]][2,]) ) ,col=rgb(0,0.5,0,alpha=0.2),border=NA)
  	lines(col.dI.mn[[2]]~NI.s,col=rgb(0,0.5,0))
	mtext("Frond color intensity",side=2,line=2)
  	points(colmn[c(1:10,91:100)]~niDmn[c(1:10,91:100)],pch=16, col=rep(c(rgb(0.75,0.55,0.35),rgb(0,0.5,0)),each=10)  )
	axis(at=c(bufferX(NI.s,.075)), label=c("most urban","most rural"),line=1,tick=FALSE,las=2,side=1,cex.axis=0.75)
suppressWarnings(plot(1~1,pch=NA,xlim=bufferX(NI.s,0.05),ylim=bufferX(c(odmn),0.05),ylab="",xlab="",main="" ,xaxt="n")) #R dislikes the formula 1~1
	polygon (x= c(NI.s, rev(NI.s) ), y= c( micr.dI.ci[[1]][1,], rev(micr.dI.ci[[1]][2,]) ) ,col=rgb(0.75,0.55,0.35,alpha=0.2),border=NA)
  	lines(micr.dI.mn[[1]]~NI.s,col=rgb(0.75,0.55,0.35))
 	polygon (x= c(NI.s, rev(NI.s) ), y= c( micr.dI.ci[[2]][1,], rev(micr.dI.ci[[2]][2,]) ) ,col=rgb(0,0.5,0,alpha=0.2),border=NA)
  	lines(micr.dI.mn[[2]]~NI.s,col=rgb(0,0.5,0))
  	points(odmn[c(1:10,91:100)]~niDmn[c(1:10,91:100)],pch=16,col= rep(c(rgb(0.75,0.55,0.35), rgb(0,0.5,0)),each=10)  )
	mtext("Microbial density, OD 600 nm",side=2,line=2)
  	legend(62,0.133, c("Urban microbes", "Rural microbes"), fill= c(rgb(0.75,0.55,0.35),rgb(0,0.5,0)), bty="n")
suppressWarnings(plot(1~1,pch=NA,xlim=bufferX(NI.s,0.05),ylim=bufferX(aggmn,0.05),ylab="",xlab="",main="" )) #R dislikes the formula 1~1
	polygon (x= c(NI.s, rev(NI.s) ), y= c( agg.dI.ci[[1]][1,], rev(agg.dI.ci[[1]][2,]) ) ,col=rgb(0.75,0.55,0.35,alpha=0.2),border=NA)
  	lines(agg.dI.mn[[1]]~NI.s,col=rgb(0.75,0.55,0.35))
 	polygon (x= c(NI.s, rev(NI.s) ), y= c( agg.dI.ci[[2]][1,], rev(agg.dI.ci[[2]][2,]) ) ,col=rgb(0,0.5,0,alpha=0.2),border=NA)
  	lines(agg.dI.mn[[2]]~NI.s,col=rgb(0,0.5,0))###HIGH ZINC ONLY
  	points(aggmn[c(1:10,91:100)]~niDmn[c(1:10,91:100)],pch=16,col= rep(c(rgb(0.75,0.55,0.35), rgb(0,0.5,0)),each=10)  )
	mtext("Frond aggregation",side=2,line=2)
	axis(at=c(bufferX(NI.s,.075)), label=c("most urban","most rural"),line=1,tick=FALSE,las=2,side=1,cex.axis=0.75)
	mtext("%Permeable surface, duckweed site",side=1,line=5.5, at = 55)
dev.off()

##for supplementary figure of zinc effects by population.
zincpopmn <- tapply(duckdat$sqmmE,paste(duckdat$numzinc,duckdat$duck ),mean,na.rm=T)
zincpopse <- tapply(duckdat$sqmmE,paste(duckdat$numzinc,duckdat$duck),std.error)
zincMICRmn <- tapply(duckdat$invOD,paste(duckdat$numzinc,duckdat$micr ),mean,na.rm=T)
zincMICRse <- tapply(duckdat$invOD,paste(duckdat$numzinc,duckdat$micr),std.error)
pdf("zinc_eff_var.pdf",height=3.5,width=5)
par(mar=c(5,3.5,1,1))
par(oma=c(2,1,1,0))
layout(matrix(1:2,nrow=1),widths=c(1,1))
plot(zincpopmn~I(rep(poptab$Nimperv,times=2)+rep(c(0,0.25),each=10)),pch=16,col=rep(c(rgb(0,0,0), rgb(1,0,0)),each=10),ylab="",xlab="",ylim=c(2,9),xlim=bufferX(NI.s,.05))
	arrows( I(rep(poptab$Nimperv,times=2)+rep(c(0,0.25),each=10)),y0=zincpopmn-zincpopse,y1=zincpopmn+zincpopse,length=0,col= rep(c(rgb(0,0,0), rgb(1,0,0)),each=10) )
	mtext(expression("Duckweed frond area, mm"^2),side=2,line=2)
	mtext("most urban",side=2,line=0.25, at=-0.5,cex=0.75)
	mtext("most rural",side=4,line=-0.25, at=-0.5,cex=0.75)
	mtext("%Permeable surface,",side=1,line=3.5)
	mtext("duckweed site",side=1,line=4.5)
	mtext("a.",side=3,line=1,at=42)
plot(1/zincMICRmn~I(rep(poptab$Nimperv,times=2)+rep(c(0,0.25),each=10)),pch=16,col=rep(c(rgb(0,0,0), rgb(1,0,0)),each=10),ylab="",xlab="",ylim=c(0.06,0.085),xlim=bufferX(NI.s,.05))
	arrows( I(rep(poptab$Nimperv,times=2)+rep(c(0,0.25),each=10)),y0=1/(zincMICRmn-zincMICRse),y1=1/(zincMICRmn+zincMICRse),length=0,col= rep(c(rgb(0,0,0), rgb(1,0,0)),each=10) )
	mtext("Microbial density, OD 600 nm",side=2,line=2)
	mtext("most urban",side=2,line=0.25, at=0.0505,cex=0.75)
	mtext("most rural",side=4,line=-0.25, at=0.0505,cex=0.75)
	mtext("%Permeable surface,",side=1,line=3.5)
	mtext("microbe site",side=1,line=4.5)
	mtext("b.",side=3,line=1,at=42)
dev.off()


######################################################
### do sympatric combinations result in greater fitness and does contamination disrupt this?

#test by fitting models for duckweed area (fitness) and microbial density.
# fit sympatric effects for fitness. here must still retain minimum duck + micr random effects bc blanquart et al 2013 (see manuscript).
#turn up iterations if p-values borderline <.1
#reported models are full and best after iterative removal, see space models for example code of full stepwise procedure, results are repeatable.
summary(MCMCglmm(sqmmE ~ x + I(x^2) + y + I(y^2) + col + row + I(row^2) +
 SA + numzinc + SA:numzinc, random = ~ duck + micr,
 family="gaussian", data=duckdat,verbose=FALSE,pr=TRUE,nitt=11000, thin=10, burnin=1000)) 
summary(MCMCglmm(sqmmE ~ x + I(x^2) + y + I(y^2) + col + row + I(row^2) +
  numzinc + SA:numzinc, random = ~ duck + micr,
 family="gaussian", data=duckdat,verbose=FALSE,pr=TRUE,nitt=11000, thin=10, burnin=1000)) 
summary(MCMCglmm(sqmmE ~ x + I(x^2) + y + I(y^2) + col + row + I(row^2) +
 distItoF + numzinc + distItoF:numzinc, 
 family="gaussian", data=duckdat,verbose=FALSE,pr=TRUE,nitt=11000, thin=10, burnin=1000)) 
summary(MCMCglmm(sqmmE ~ x + I(x^2) + y + I(y^2) + col + row + I(row^2) +
 distItoF + numzinc , 
 family="gaussian", data=duckdat,verbose=FALSE,pr=TRUE,nitt=11000, thin=10, burnin=1000)) 
#distItoF  sig, and one other term must be removed. zinc & zincXdistItoF models apparently DIC equivalent. w/o zinc:distItoF and with main effect of zinc seems to have better behaviour
# sympatric microbes marginally costly at high zinc,  - SA effect is borderline; though model always fits better with it
summary(MCMCglmm(invOD ~ x + I(x^2) + y + I(col^2)+ row + I(row^2) + 
 SA + numzinc + SA:numzinc, random = ~ duck + micr,
 family="gaussian", data=duckdat,verbose=FALSE,pr=TRUE,nitt=11000, thin=10, burnin=1000)) 
summary(MCMCglmm(invOD ~ x + I(x^2) + y + I(col^2)+ row + I(row^2) + 
 SA:numzinc, random = ~ duck + micr, #only that  SA by NZ
 family="gaussian", data=duckdat,verbose=FALSE,pr=TRUE,nitt=11000, thin=10, burnin=1000)) 
summary(MCMCglmm(invOD ~ x + I(x^2) + y + I(col^2)+ row + I(row^2) + 
 distItoF + numzinc + distItoF:numzinc, 
 family="gaussian", data=duckdat,verbose=FALSE,pr=TRUE,nitt=11000, thin=10, burnin=1000)) 
summary(MCMCglmm(invOD ~ x + I(x^2) + y + I(col^2)+ row + I(row^2) + 
 distItoF, 
 family="gaussian", data=duckdat,verbose=FALSE,pr=TRUE,nitt=11000, thin=10, burnin=1000)) 
#only distItoF, still pos.

#refit best models 
set.seed(100) #run the following 4 models at once after setting seed for consistency
#models reported in supplement
AE.H2 <- MCMCglmm(sqmmE ~ x + I(x^2) + y + I(y^2) + col + row + I(row^2) + numzinc + SA:numzinc, random = ~ duck + micr, family="gaussian", data=duckdat,verbose=FALSE,pr=TRUE,nitt=101000, thin=10, burnin=1000)
invOD.H2 <- MCMCglmm(invOD ~ x + I(x^2) + y + I(col^2)+ row + I(row^2) +  SA:numzinc, random = ~ duck + micr, family="gaussian", data=duckdat,verbose=FALSE,pr=TRUE,nitt=101000, thin=10, burnin=1000)
AE.if <-(MCMCglmm(sqmmE ~ x + I(x^2) + y + I(y^2) + col + row + I(row^2) + distItoF + numzinc ,  family="gaussian", data=duckdat,verbose=FALSE,pr=TRUE,nitt=101000, thin=10, burnin=1000)) 
invOD.if <- MCMCglmm(invOD ~ x + I(x^2) + y + I(col^2)+ row + I(row^2) +  distItoF , family="gaussian", data=duckdat,verbose=FALSE,pr=TRUE,nitt=101000, thin=10, burnin=1000) 
summary(AE.H2); summary(invOD.H2); summary(AE.if); summary(invOD.if)

###Extract predictions from best models for reporting and figures
#area, SA
aeh2s <- AE.H2$Sol
dFSAinZ<- sapply(c(0,1), function(z) mean(
		aeh2s[,1] + aeh2s[,2]*mx + aeh2s[,3]*mx2 + aeh2s[,4]*my + aeh2s[,5]*my2 + aeh2s[,6]*mcol + aeh2s[,7]*mrow + aeh2s[,8]*mrow2 +
				 aeh2s[,9]*1 + + aeh2s[,10]*z*1 + rowMeans(aeh2s[,11:20])+ rowMeans(aeh2s[,21:30])  ) )
(dFSAinZ[1] - dFSAinZ[2]) / dFSAinZ[1]		 		 
dFSAinZ[1] ; dFSAinZ[2]
#OD, SA
oh2s <- invOD.H2$Sol
dOSAinZ<- sapply(c(0,1), function(z) mean(oh2s[,1] + oh2s[,2]*mx + oh2s[,3]*mx2 + oh2s[,4]*my + oh2s[,5]*mcol2+ oh2s[,6]*mrow + oh2s[,7]*mrow2 + 
 		oh2s[,8]*z*1 + rowMeans(oh2s[,9:18]) + rowMeans(oh2s[,19:28]) ))
1/(dOSAinZ[1]) ; 1/(dOSAinZ[2]) #/ dOSAinZ[1]		 		 
#area, distItoF
Aifsol <- AE.if$Sol
if.s <- seq(from=min(duckdat$distItoF),to=max(duckdat$distItoF),length.out=1000)
aif.ci <- sapply(if.s, function(z) HPDinterval(as.mcmc(
		Aifsol[,1] + Aifsol[,2]*mx + Aifsol[,3]*mx2 + Aifsol[,4]*my + Aifsol[,5]*my2 + Aifsol[,6]*mcol + Aifsol[,7]*mrow + Aifsol[,8]*mrow2 +
				 Aifsol[,9]*z + Aifsol[,10]*0.5  ),prob=0.95) )
aif.mn <-  sapply(if.s, function(z) 	mean(Aifsol[,1] + Aifsol[,2]*mx + Aifsol[,3]*mx2 + Aifsol[,4]*my + Aifsol[,5]*my2 + Aifsol[,6]*mcol + Aifsol[,7]*mrow + Aifsol[,8]*mrow2 +
				 Aifsol[,9]*z + Aifsol[,10]*0.5 ) )
(aif.mn[1]- aif.mn[1000])/aif.mn[1]
#OD, distItoF
ODifsol <- invOD.if$Sol
micr.if.ci <- sapply(if.s, function(z)    HPDinterval(as.mcmc(
	ODifsol[,1] + ODifsol[,2]*mx + ODifsol[,3]*mx2 + ODifsol[,4]*my + ODifsol[,5]*mcol2+ ODifsol[,6]*mrow + ODifsol[,7]*mrow2 + 
 		ODifsol[,8]*z ),prob=0.95) ) 
micr.if.mn <- sapply(if.s, function(z)    mean(
	ODifsol[,1] + ODifsol[,2]*mx + ODifsol[,3]*mx2 + ODifsol[,4]*my + ODifsol[,5]*mcol2+ ODifsol[,6]*mrow + ODifsol[,7]*mrow2 + 
 		ODifsol[,8]*z ) ) 
1/(micr.if.mn[1]); 1/(micr.if.mn[1000])

#calculate means for figure
dH2res.m  <- tapply(duckdat$sqmmE, paste(duckdat$numzinc,duckdat$SA),mean)
dH2res.se <- tapply(duckdat$sqmmE, paste(duckdat$numzinc,duckdat$SA),std.error)
mH2res.m  <- tapply(duckdat$invOD, paste(duckdat$numzinc,duckdat$SA),mean)
mH2res.se <- tapply(duckdat$invOD, paste(duckdat$numzinc,duckdat$SA),std.error)
mH2A.m  <- tapply(duckdat$sqmmE, paste(duckdat$duck,duckdat$micr,duckdat$numzinc),mean)
mH2A.se <- tapply(duckdat$sqmmE, paste(duckdat$duck,duckdat$micr,duckdat$numzinc),std.error)
mH2.m  <- tapply(duckdat$invOD, paste(duckdat$duck,duckdat$micr,duckdat$numzinc),mean)
mH2.se <- tapply(duckdat$invOD, paste(duckdat$duck,duckdat$micr, duckdat$numzinc),std.error)
mH2if <- tapply(duckdat$distItoF, paste(duckdat$duck,duckdat$micr, duckdat$numzinc),mean)
mH2sa <- tapply(duckdat$SA, paste(duckdat$duck,duckdat$micr, duckdat$numzinc),mean)
Ozincmn <- tapply(duckdat$invOD,duckdat$numzinc,mean,na.rm=T)
Ozincse <- tapply(duckdat$invOD,duckdat$numzinc,std.error)
zincmn <- tapply(duckdat$sqmmE,duckdat$numzinc,mean,na.rm=T)
zincse <- tapply(duckdat$sqmmE,duckdat$numzinc,std.error)
(zincmn[1]-zincmn[2])/ zincmn[1]

pdf("zincANDsymp_den.pdf",height=4.8,width=5.4)
par(mar=c(2.5,2.5,1.5,0))
par(oma=c(3.25,2.25,0.5,1.5))
layout(matrix(1:6,ncol=3),widths=c(0.75,1,1.5))
plot(zincmn~c(0,1),xlim=c(-0.5,1.5),pch=16,ylim=bufferX(c(zincmn+zincse,zincmn-zincse),0.1),xlab="",ylab="",xaxt="n")
	arrows(x0 =c(0,1), y0 = zincmn - zincse ,y1 = zincmn + zincse,length=0)
	axis(at=c(0,1),side=1,labels=c("Low","High"))
	mtext(expression("Duckweed frond area, mm"^2),side=2,line=3.1,cex=0.95)	
	mtext("a.",side=3,line=0.5,at=-1.3)
plot(1/Ozincmn~c(0,1),xlim=c(-0.5,1.5),pch=16,ylim=bufferX(1/(c(Ozincmn+Ozincse,Ozincmn-Ozincse)),0.1),xlab="",ylab="",xaxt="n")
	arrows(x0 =c(0,1), y0 = 1/(Ozincmn - Ozincse) ,y1 = 1/(Ozincmn + Ozincse),length=0)
	axis(at=c(0,1),side=1,labels=c("Low","High"))
	mtext("Microbial density, OD 600 nm",side=2,line=3.1,cex=0.95)	
	mtext("Zinc",side=1,line=2.5,cex=0.95)	
	mtext("d.",side=3,line=0.5,at=-1.3)
plot(dH2res.m~c(1,2,1,2),ylim=c(4.75,6.5),xlim=c(0.5,2.5),ylab="",xaxt="n",xlab="",pch=16,col=rep(c(rgb(0,0,0),rgb(1,0,0)),each=2),cex=1.25 )
	axis(side=1,at=c(1,2),labels=c("A","S"))
	mtext("b.",side=3,line=0.5,at=0)
	arrows(x0=c(1,2,1,2), y0 = dH2res.m-dH2res.se, y1= dH2res.m + dH2res.se,length=0,col=rep(c(rgb(0,0,0),rgb(1,0,0)),each=2) ,lwd=2)
plot(1/mH2res.m~c(1,2,1.1,2.1),ylim=c(0.0625,0.075),xlim=c(0.5,2.5),ylab="",xaxt="n",xlab="",pch=16,cex=1.25,col=rep(c(rgb(0,0,0),rgb(1,0,0)),each=2) )
	axis(side=1,at=c(1,2),labels=c("A","S"))
	arrows(x0=c(1,2,1.1,2.1), y0 = 1/(mH2res.m-mH2res.se), y1= 1/(mH2res.m + mH2res.se),length=0 ,lwd=2, col=rep(c(rgb(0,0,0),rgb(1,0,0)),each=2))
	mtext("Host-Community",side=1,line=2.5,cex=0.95)	
	mtext("Combination",side=1,line=3.75,cex=0.95)	
	mtext("e.",side=3,line=0.5,at=0)
	legend(0.35,0.0655,c("Low Zinc","High Zinc"),fill=c(rgb(0,0,0),rgb(1,0,0)),bty="n",cex=1.1)
suppressWarnings(plot(1~1,pch=NA,xlim=bufferX(if.s,0.05),ylim=bufferX(c(mH2A.m),0.1),ylab="",xlab="",main="" )) #R dislikes the formula 1~1
	polygon (x= c(if.s, rev(if.s) ), y= c( aif.ci[1,], rev(aif.ci[2,]) ) ,col=rgb(0,0,0,alpha=0.2),border=NA)
	points(mH2A.m ~ mH2if,pch=1,col=rgb(0,0,0),cex=0.5)
 	mtext("c.",side=3,line=0.5,at=0.15)
	lines(aif.mn~if.s)
suppressWarnings(plot(1~1,pch=NA,xlim=bufferX(if.s,0.05),ylim=bufferX(1/mH2.m,0.1),ylab="",xlab="",main="" )) #R dislikes the formula 1~1
	polygon (x= c(if.s, rev(if.s) ), y=1/ c( micr.if.ci[1,], rev(micr.if.ci[2,]) ) ,col=rgb(0,0,0,alpha=0.2),border=NA)
	points(1/mH2.m ~ mH2if,pch=1,col=rgb(0,0,0),cex=0.5)
 	mtext("f.",side=3,line=0.5,at=0.15)
	lines(1/micr.if.mn~if.s)
	mtext("Field-inocula Unifrac distance",side=1,line=2.5,cex=0.95)
dev.off()


#########################################################
########## DOES FITNESS CORRELATE WITH TRAIT EXPRESSSION
######### Figure 4
#######################################################

# mean scale all traits so that slopes are comparable estimates of the strength of "selection" (Hansen & Houle 2008)
# mean scale fitness at the full dataset ("global") level, and traits at the treatment ("local") level  (De Lisle & Svensson 2017)
# use raw OD, since the inverse transformation is not appropriate for this analysis (Hansen & Houle 2008)
tfitmat2 <- data.frame(sqmmE=duckdat$sqmmE, 
						untOD=duckdat$OD600, 
						colint=duckdat$colint, 
						agg = duckdat$areapperE,
						rnd=duckdat$meanround) 
localtrtsc <- as.data.frame(sapply(1:ncol(tfitmat2), function(trait) 
	sapply(unique(duckdat$micrNImpv), function(micr) 
		sapply(c(0,1), function(z) 
			tapply( tfitmat2[duckdat$micrNImpv==micr & duckdat$numzinc==z, trait ], duckdat$duck[duckdat$micrNImpv==micr & duckdat$numzinc==z ],mean, na.rm=T )
			 )  )   ))
colnames(localtrtsc) <- colnames(tfitmat2)  
localtrtsc$numzinc <- rep(rep(c(0,1),each=10),times=10)
localtrtsc$micrNImpv <- rep(unique(duckdat$micrNImpv),each=20)
#globally scale fitness, and locally scale traits
localtrtsc$sqmmE <- localtrtsc$sqmmE/mean(localtrtsc$sqmmE,na.rm=T)
localtrtsc$untOD <- sapply(1:nrow(localtrtsc), function(z)  localtrtsc$untOD[z] / mean(localtrtsc$untOD[ localtrtsc$numzinc==localtrtsc$numzinc[z] & localtrtsc$micrNImpv==localtrtsc$micrNImpv[z] ], na.rm=T   ) )
localtrtsc$colint <- sapply(1:nrow(localtrtsc), function(z)  localtrtsc$colint[z] / mean(localtrtsc$colint[ localtrtsc$numzinc==localtrtsc$numzinc[z] & localtrtsc$micrNImpv==localtrtsc$micrNImpv[z] ], na.rm=T   ) )
localtrtsc$agg <- sapply(1:nrow(localtrtsc), function(z)  localtrtsc$agg[z] / mean(localtrtsc$agg[ localtrtsc$numzinc==localtrtsc$numzinc[z] & localtrtsc$micrNImpv==localtrtsc$micrNImpv[z] ], na.rm=T   ) )
localtrtsc$rnd <- sapply(1:nrow(localtrtsc), function(z)  localtrtsc$rnd[z] / mean(localtrtsc$rnd[ localtrtsc$numzinc==localtrtsc$numzinc[z] & localtrtsc$micrNImpv==localtrtsc$micrNImpv[z] ], na.rm=T   ) )

##fit every model - so for EACH microbe, zinc, trait combo
utrts <- data.frame( trait = rep(2:5, each = 20) , zinc = rep(rep(0:1,each=10 ), times = 4),  micrNImpv = rep(poptab$Nimperv, times = 8) )
#vector of all unique microbe-zinc-trait combinations
linmod <- list()
intmod <- list()
set.seed(100) #set seed to produce identical objects as reported in MS if desired
for(i in 1:nrow(utrts)){ #for every  microbe-zinc treatment and trait
	rows <- which(localtrtsc$micrNImpv==utrts$micrNImpv[i] & localtrtsc$numzinc==utrts$zinc[i]) #identify rows of full dataset in  microbe-zinc treatment
	ymns <- localtrtsc[rows,1] #extract fitness, then traits (xmns)
	xmns <- localtrtsc[rows, utrts$trait[i] ]
	dfmns <- data.frame(y=ymns,x=xmns)
	intmod[[i]] <- MCMCglmm(y ~ 1,  data=dfmns,verbose=FALSE,nitt=51000, thin=10, burnin=1000)
	linmod[[i]] <- MCMCglmm(y ~ x,  data=dfmns,verbose=FALSE,nitt=51000, thin=10, burnin=1000) 
	#regress fitness on trait, or just intercept
}
modres <- data.frame ( #make dataframes including intercept model and linear model results from previous step
			linmodpx = unlist(lapply(linmod, function(z) summary(z)$solutions[2,5])),
			linmodx = unlist(lapply(linmod, function(z) summary(z)$solutions[2,1])),
			linDIC = unlist(lapply(linmod, function(z) z$DIC)),
			intDIC = unlist(lapply(intmod, function(z) z$DIC)) ) 
FitxTraitbyTreat<- as.data.frame(cbind(utrts,modres))
FitxTraitbyTreatDIC <- data.frame (			
			linDIC = unlist(lapply(linmod, function(z) z$DIC)),
			intDIC = unlist(lapply(intmod, function(z) z$DIC)) ) 
whichmod <- sapply(1:nrow(FitxTraitbyTreatDIC), function(z) which(FitxTraitbyTreatDIC[z,]==min(FitxTraitbyTreatDIC[z,])) )
table(whichmod) #evaluate whether trait or intercept model fits better for each unique  microbe-zinc treatment and trait combination
FitxTraitbyTreat[whichmod==1,]

tranges <- sapply(2:5,function(z) range(tapply(tfitmat2[,z],paste(duckdat$duck,duckdat$micr,duckdat$zinc),mean,na.rm=T)))

slopeHPDI <- t(sapply(1:length(linmod), function(z) HPDinterval(linmod[[z]]$Sol[,2])))
#extract HPDI for each slope of each linear model of fitness~trait (for each trait-microbe-zinc combination)

##functionalized version of the above
modelresults_FxTxT <- function(utrts, scaledtfitmat) {  #takes arguments of utrts, localstrtsc objects
	linmod <- list()
	intmod <- list()
	for(i in 1:nrow(utrts)){ 
		rows <- which(scaledtfitmat$micrNImpv==scaledtfitmat$micrNImpv[i] & scaledtfitmat$numzinc==utrts$zinc[i]) #identify rows of full dataset in  microbe-zinc treatment
		ymns <- scaledtfitmat[rows,1] #extract fitness, then traits (xmns)
		xmns <- scaledtfitmat[rows, utrts$trait[i] ]
		dfmns <- data.frame(y=ymns,x=xmns)
		intmod[[i]] <- MCMCglmm(y ~ 1,  data=dfmns,verbose=FALSE,nitt=51000, thin=10, burnin=1000)
		linmod[[i]] <- MCMCglmm(y ~ x,  data=dfmns,verbose=FALSE,nitt=51000, thin=10, burnin=1000) 
	}
	modres <- data.frame (
			linmodpx = unlist(lapply(linmod, function(z) summary(z)$solutions[2,5])),
			linmodx = unlist(lapply(linmod, function(z) summary(z)$solutions[2,1])),
			linDIC = unlist(lapply(linmod, function(z) z$DIC)),
			intDIC = unlist(lapply(intmod, function(z) z$DIC)) ) 
	return(modres)
}

###permutations, using permuted datasets and function above for calculating selection gradients, DIC.
######WARNING: Permutations may take a LONG time (especially on a personal laptop), suggest run overnight, reduce the iterations from 1000 to 100
#the below is commented out, with option to load the result from a previous run
# perm_FxTxT_winT <- array(dim=c(80,4,1000))
# umz <- data.frame(micr = rep(poptab$Nimperv,times=2), zinc=rep(c(0,1),each=10))
# set.seed(100) 
# for(i in 1:1000){	
# 	resamp_sctfit <- matrix(NA, ncol=4,nrow=nrow(localtrtsc))
# 	for(j in 1:nrow(umz)){
# 		trtrows <- which(localtrtsc$micrNImpv==umz$micr[j] & localtrtsc$numzinc == umz$zinc[j]) 
# 		resamp_sctfit[trtrows,] <- as.matrix(localtrtsc[ sample(trtrows,length(trtrows),replace=F) ,2:5 ])
# 	} 
# 	randtrt <- as.data.frame( cbind(localtrtsc[,1],resamp_sctfit,localtrtsc$numzinc,localtrtsc$micrNImpv) )
# 	#this randomization preserves y-values (and x-values) present in a given zinc/microbe treatment to allow DIC comparisons, but scrambles the trait-fit relationship
# 	colnames(randtrt) <- colnames(localtrtsc)
# 	perm_FxTxT_winT[,,i] <- as.matrix(modelresults_FxTxT(utrts,randtrt))
# print(i)
# }
# save(perm_FxTxT_winT,file="R inputs/trt_eff_perms.Rdata")
load(file="R inputs/trt_eff_perms.Rdata")

#compare DIC in real and permuted data.
DICsimP <- sapply(1:nrow(modres), function(z) findInterval(modres[z,3], sort(perm_FxTxT_winT[z,3,]) )/length(perm_FxTxT_winT[1,1,])) 
isDICsig <- DICsimP  < 0.05  & modres[,3]<modres[,4] #DIC outside permuted 95% bounds, and require also that the linear model fit better than an intercept model
which(isDICsig)

###plot main figure - Figure 4
betas <- data.frame(BOD= modres[1:20,2], Bagg = modres[41:60,2], Bcol = modres[21:40,2], Brnd = modres[61:80,2], 
					micr = utrts$micr[1:20], zinc = utrts$zinc[1:20], NImp = rep(poptab$Nimperv,times=2) )
pdf("treat_eff_sel.pdf",width=3.1,height=6.5)
par(mfrow=c(4,1))
par(mar=c(0,0,0,0))
par(oma=c(6,4,0.5,2))
plot(betas$BOD~I(as.numeric(as.factor(betas$NImp))+rep(c(-0.1,0.1),each=10)),col=rgb(betas$zinc,0,0),ylim=c(-12,11),xlab="",ylab="",xaxt="n")
	arrows(as.numeric(as.factor(betas$NImp))+rep(c(-0.1,0.1),each=10), y0 = slopeHPDI[1:20,1], y1=slopeHPDI[1:20,2],col=rgb(betas$zinc,0,0),length=0 )
	abline(h=0,lty=3)
	text(as.numeric(as.factor(betas$NImp))+rep(c(-0.1,0.1),each=10), -11.8, ifelse(isDICsig[1:20], "*",""), col=rgb(betas$zinc,0,0), cex=2 )
	legend(6,-1, c("Low zinc","High zinc"), fill=c(rgb(0,0,0),rgb(1,0,0)),bty="n",cex=1.25)
	text(5.5,10,"Microbial density",cex=1.4)
plot(betas$Bcol~I(as.numeric(as.factor(betas$NImp))+rep(c(-0.1,0.1),each=10)),col=rgb(betas$zinc,0,0),ylim=c(-12,11),xlab="",ylab="",xaxt="n")
	arrows(as.numeric(as.factor(betas$NImp))+rep(c(-0.1,0.1),each=10), y0 = slopeHPDI[21:40,1], y1=slopeHPDI[21:40,2],col=rgb(betas$zinc,0,0),length=0 )
	abline(h=0,lty=3)
	text(as.numeric(as.factor(betas$NImp))+rep(c(-0.1,0.1),each=10), -11.8, ifelse(isDICsig[21:40], "*",""), col=rgb(betas$zinc,0,0), cex=2 )
	text(5.5,10,"Color intensity",cex=1.4)
	mtext("Estimated selection gradient",side=2,line=2.5, at = -10)
plot(betas$Bagg~I(as.numeric(as.factor(betas$NImp))+rep(c(-0.1,0.1),each=10)),col=rgb(betas$zinc,0,0),ylim=c(-12,11),xlab="",ylab="",xaxt="n")
	arrows(as.numeric(as.factor(betas$NImp))+rep(c(-0.1,0.1),each=10), y0 = slopeHPDI[41:60,1], y1=slopeHPDI[41:60,2],col=rgb(betas$zinc,0,0),length=0 )
	abline(h=0,lty=3)
	text(as.numeric(as.factor(betas$NImp))+rep(c(-0.1,0.1),each=10), -11.8, ifelse(isDICsig[41:60], "*",""), col=rgb(betas$zinc,0,0), cex=2 )
	text(5.5,10,"Aggregation",cex=1.4)
plot(betas$Brnd~I(as.numeric(as.factor(betas$NImp))+rep(c(-0.1,0.1),each=10)),col=rgb(betas$zinc,0,0),ylim=c(-12,11),xlab="",ylab="",xaxt="n")
	arrows(as.numeric(as.factor(betas$NImp))+rep(c(-0.1,0.1),each=10), y0 = slopeHPDI[61:80,1], y1=slopeHPDI[61:80,2],col=rgb(betas$zinc,0,0),length=0 )
	abline(h=0,lty=3)
	text(as.numeric(as.factor(betas$NImp))+rep(c(-0.1,0.1),each=10), -11.8, ifelse(isDICsig[61:80], "*",""), col=rgb(betas$zinc,0,0), cex=2 )
	axis(at=1:10, labels = round(sort(poptab$Nimperv)),side=1)
	mtext("most urban",side=2,line=0, at=-21,cex=0.75)
	mtext("most rural",side=4,line=0, at=-20.5,cex=0.75)
	text(5.5,10,"Roundness",cex=1.4)
	mtext("%permeable surface area",side=1,line=3)
	mtext("of microbe site",side=1,line=5)
dev.off()

##make individual level predictions and HPDI from previously fitted models for selection gradients for figure
linpreds <- list()
x.seqs <- list()
linhpdi <- list()
for(i in 1:length(linmod)){
	rows <- which(duckdat$micrNImpv==utrts$micrNImpv[i] & duckdat$numzinc==utrts$zinc[i])
	xmns <- tapply(tfitmat2[ rows, utrts$trait[i] ]/mean(tfitmat2[ rows, utrts$trait[i] ],na.rm=T),duckdat$duck[rows],mean,na.rm=T)
	x.seqs[[i]] <- seq(from = range(xmns)[1], to = range(xmns)[2],length.out=1000)
	linpreds[[i]] <- sapply(x.seqs[[i]], function(z) mean(linmod[[i]]$Sol[,1] + linmod[[i]]$Sol[,2]*z) )
	linhpdi[[i]] <- sapply(x.seqs[[i]], function(z) HPDinterval(linmod[[i]]$Sol[,1] + linmod[[i]]$Sol[,2]*z,0.95) )
}
#supplemental figure
pdf("treat_eff_sel_all.pdf",width=6.75,height=8)
layout(matrix(1:40,ncol=4,byrow=F))#10 microbes X 5 traits (x 2 zinc)
par(mar=c(0,0,0,0))
par(oma=c(4,4,3,6))
for(z in 2:5){ #zth trait
	for(m in 1:10){ #mth microbe
		rowsH <- which(duckdat$micrNImpv==sort(poptab$Nimperv)[m] & duckdat$zinc=="H")
		rowsL <- which(duckdat$micrNImpv==sort(poptab$Nimperv)[m] & duckdat$zinc=="L")
		uindex <- which(utrts$micrNImpv==sort(poptab$Nimperv)[m] & utrts$trait==z)
		ymns <- localtrtsc[ localtrtsc$micrNImpv==sort(poptab$Nimperv)[m] , 1 ]
		xmns <- localtrtsc[ localtrtsc$micrNImpv==sort(poptab$Nimperv)[m] , z ]
		plot(ymns~xmns,xaxt="n",yaxt="n",ylim=c(-0.25,1.9),xlim=bufferX(range(x.seqs[which(utrts$trait==z)]),0.025),col=rgb(rep(c(0,1),each=10),0,0))
			lines(linpreds[[uindex[1]]]~x.seqs[[uindex[1]]], col=rgb(0,0,0) )
			lines(linpreds[[uindex[2]]]~x.seqs[[uindex[2]]], col=rgb(1,0,0) )
			polygon( c(x.seqs[[uindex[1]]],rev(x.seqs[[uindex[1]]])), 
				c(linhpdi[[uindex[1]]][1,],rev(linhpdi[[uindex[1]]][2,])),col=rgb(0,0,0,alpha=0.25),border=NA )
			polygon( c(x.seqs[[uindex[2]]],rev(x.seqs[[uindex[2]]])), 
				c(linhpdi[[uindex[2]]][1,],rev(linhpdi[[uindex[2]]][2,])),col=rgb(1,0,0,alpha=0.25),border=NA )
			text(range(x.seqs[which(utrts$trait==z)])[1],-0.2, ifelse(isDICsig[uindex[1]],"*",""),col=rgb(0,0,0),cex=1.5)
			text(range(x.seqs[which(utrts$trait==z)])[2],-0.2, ifelse(isDICsig[uindex[2]],"*",""),col=rgb(1,0,0),cex=1.5)
			if(z==2){axis(side=2)}
			if(m==10){axis(side=1)}
			if(m==1){mtext(tnames[z],side=3,line=0.5) }
			if(m==5 & z==2){mtext("Duckweed relative fitness",side=2,line=2,at=-1)}
			if(m==10 & z==3){mtext("Mean scaled trait value",side=1,line=2,at=range(x.seqs[which(utrts$trait==z)])[2])}
			if(z==5){mtext(round(sort(poptab$Nimperv)[m]),side=4,line=0.5)}
			if(z==5 & m==5){mtext("% permeable surface area, microbe site",side=4,line=2.5,at=-1) }
			if(z==5 & m==1){axis(at=1.95, label="most urban",line=NA,tick=FALSE,las=2,side=4)}
			if(z==5 & m==10){axis(at=-0.25, label="most rural",line=NA,tick=FALSE,las=2,side=4)}
	}
}
dev.off()

#check for patterns in betas (selection gradients) and zinc treatments or microbiome origin along urban to rural gradient.
#reported models are full and best after iterative removal, see space models for example code of full stepwise procedure, results are repeatable.
summary(MCMCglmm(BOD ~ zinc + NImp + zinc:NImp,data=betas, verbose=F,nitt=101000, burnin=1000,thin=10)) #best fit
summary(MCMCglmm(Bcol ~ zinc + NImp + zinc:NImp,data=betas, verbose=F,nitt=101000, burnin=1000,thin=10))
	summary(MCMCglmm(Bcol ~ 1,data=betas, verbose=F,nitt=101000, burnin=1000,thin=10)) #best fit
summary(MCMCglmm(Bagg ~ zinc + NImp + zinc:NImp,data=betas, verbose=F,nitt=101000, burnin=1000,thin=10))
	summary(MCMCglmm(Bagg ~ 1,data=betas, verbose=F,nitt=101000, burnin=1000,thin=10)) #best fit
summary(MCMCglmm(Brnd ~ zinc + NImp + zinc:NImp,data=betas, verbose=F,nitt=101000, burnin=1000,thin=10))
	summary(MCMCglmm(Brnd ~ zinc + zinc:NImp,data=betas, verbose=F,nitt=101000, burnin=1000,thin=10)) #best DIC but zinc main n.s., and interactive marginal
set.seed(100)
betaODNIZ <- (MCMCglmm(BOD ~ zinc + NImp + zinc:NImp,data=betas, verbose=F,nitt=1010000, burnin=10000,thin=100))
BodnizS <- betaODNIZ$Sol
bodnizH <- lapply(c(0,1), function(z) sapply(NI.s, function(n) HPDinterval(BodnizS[,1] + BodnizS[,2]*z + BodnizS[,3]*n + BodnizS[,4]*z*n ,0.95) ))
bodnizM <- lapply(c(0,1), function(z) sapply(NI.s, function(n) mean(BodnizS[,1] + BodnizS[,2]*z + BodnizS[,3]*n + BodnizS[,4]*z*n) ))

####supplemental figure
pdf("treat_eff_sel_ODNI.pdf",width=3,height=3.5)
par(mar=c(6,4,1,1))
plot(betas$BOD~betas$NImp,pch=16,ylim=c(-1,1),xlab="",ylab="",col=rgb(betas$zinc,0,0))
	abline(h=0,lty=3)
	text(5.5,9,"Microbial density",cex=1.5)
	polygon(c(NI.s,rev(NI.s)),c(bodnizH[[1]][1,], rev(bodnizH[[1]][2,])),col=rgb(0,0,0,alpha=0.25),border=NA )
	lines(bodnizM[[1]]~NI.s)
	polygon(c(NI.s,rev(NI.s)),c(bodnizH[[2]][1,], rev(bodnizH[[2]][2,])),col=rgb(1,0,0,alpha=0.25),border=NA )
	lines(bodnizM[[2]]~NI.s,col=rgb(1,0,0))
	mtext("%permeable surface area",side=1,line=3)
	mtext("of microbe site",side=1,line=4)
	mtext("Estimated selection gradient",side=2,line=2.5)
	mtext("most urban",side=2,line=0, at=-1.93,cex=0.75)
	mtext("most rural",side=4,line=0, at=-1.9,cex=0.75)
	legend(65,1, c("Low zinc","High zinc"), fill=c(rgb(0,0,0),rgb(1,0,0)),bty="n",cex=1.25)
dev.off()

##accounting for uncertainty in beta estimates
summary(MCMCglmm(sqmmE ~ untOD + untOD:micrNImpv + untOD:numzinc + untOD:numzinc:micrNImpv, data= localtrtsc, verbose=F ) ) #should be no main effs of zinc or microbe source, since they are scaled!
summary(MCMCglmm(sqmmE ~    untOD:numzinc:micrNImpv, data= localtrtsc, verbose=F ) )  #best
summary(MCMCglmm(sqmmE ~ colint + colint:micrNImpv + colint:numzinc + colint:numzinc:micrNImpv, data= localtrtsc, verbose=F ) )  #only trait main eff
summary(MCMCglmm(sqmmE ~ agg + agg:micrNImpv + agg:numzinc + agg:numzinc:micrNImpv, data= localtrtsc, verbose=F ) ) #trait main eff only; model with agg:numzinc fits better, but term is ns
summary(MCMCglmm(sqmmE ~ rnd + rnd:micrNImpv + rnd:numzinc + rnd:numzinc:micrNImpv, data= localtrtsc, verbose=F ) ) #only trait main eff
set.seed(100)
fitnessODbytrt <- MCMCglmm(sqmmE ~  untOD:numzinc:micrNImpv, data= localtrtsc, verbose=F, nitt=1010000, burnin=10000,thin=100) 


#########################################################
########## CONSISTENCY OF FITNESS, TRAIT EXPRESSSION
#######################################################

#calculate genotype means in each microbe (treatment means for each genotype in each microbe, but pooling across zinc levels) 
DmnsinM <- lapply(1:ncol(tfitmat), function(z) sapply(sort(unique(duckdat$micr)),function(m) 
		tapply(tfitmat[duckdat$micr== m,z],duckdat$duck[duckdat$micr==m],mean,na.rm=T)  ) )# duckweed are rows, microbes are columns each trait is own matrix
#calculate the correlation for every pairwise comparison of genotypic means across microbial community treatments
Bgb.Bu2 <- lapply(DmnsinM, function(z) cor(z)) #tables of duckweed genetic correlations across pairwise microbiome treatment, 1 matrix for each trait. rows and columns are microbe origins
#generally high, but some evidence of dissimilarity of expressed traits/fitness across duckweed genotypes when measured in different microbes -- esp for roundness.

##use the matrix to ask whether microbial community distance explains the variation in heritability (note for clones, heritability is covariance in genotype means)
sapply(1:5, function(z) cor(Bgb.Bu2[[z]][lower.tri(Bgb.Bu2[[z]])==T], wu.Mladpt[lower.tri(wu.Mladpt)==T]))
#rearrange unifrac distance and pairwise genotype mean correlations into a data frame
UnifracGCdat <- lapply(1:5, function(z) data.frame(GC = Bgb.Bu2[[z]][lower.tri(Bgb.Bu2[[z]])==T], Unifrac = wu.Mladpt[lower.tri(wu.Mladpt)==T]) )
#fit a model for each trait, correlation is dependent variable, independent is unifrac distance, extract predictions and HPDI
set.seed(100)
UnifracGCmods <-	lapply(1:5, function(z) MCMCglmm(GC ~ Unifrac, data=UnifracGCdat[[z]], burnin=1000 , nitt =101000 , thin=10, verbose=F ) )
#make predictions
UD.s <- seq(from = 0, to =max(wu.Mladpt), length.out=1000)
UnifracGCpreds <-	lapply(1:5, function(z) sapply(UD.s, function(D) 		mean(UnifracGCmods[[z]]$Sol[,1] + UnifracGCmods[[z]]$Sol[,2]*D) ) )
UnifracGCHPDIs <-	lapply(1:5, function(z) sapply(UD.s, function(D) HPDinterval(UnifracGCmods[[z]]$Sol[,1] + UnifracGCmods[[z]]$Sol[,2]*D,prob=0.95 ) ) )
lapply(UnifracGCmods,summary)
lapply(UnifracGCHPDIs , function(z) cbind(z[,1],z[,1000] ) )
lapply(UnifracGCpreds , function(z) cbind(z[1],z[1000] ) )

##Figure 5
##plot results 
tnames<- c("Duckweed area","Microbial density","Color Intensity","Aggregation","Roundness")
pdf("duckweed_pairwise_Vg_vs_microbeUnifrac_den.pdf",height=2.25,width=7.75)
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

##supplemental figure including plots of all pairwise genotype mean correlations
ijdf <- data.frame(i = unlist(sapply(10:2, function(z) rep(z,times=z-1)) ), j = unlist(sapply(10:2, function(z)  (z-1):1 )) )
ijdf$unid <- sapply(1:nrow(ijdf), function(z) wunfrc.M[ ijdf$i[z] , ijdf$j[z] ])
ijdf <- ijdf[order(ijdf$unid) ,]
pdf("scat_frondA.pdf",height=3,width=5)
par(mfrow=c(5,9))
par(mar=c(0,0,0,0))
par(oma=c(0,0,0,0))
for(k in 1:nrow(ijdf)){
	dat <- rbind( DmnsinM[[1]][ , ijdf$i[k] ], DmnsinM[[1]][ ,ijdf$j[k] ] )
	plot(dat[1,]~dat[2,],pch=1,xaxt="n",yaxt="n", ylim=bufferX(range(DmnsinM[[1]]),.1),xlim=bufferX(range(DmnsinM[[1]]),.1),lty=1,col=rgb(range01(ijdf$unid)[k],0,0))
	rho <- cor(dat[1,],dat[2,])
	if(rho > sort(UnifracGCdat[[1]]$GC,decreasing=T)[23]){ # >=.75
		coefs <-coef(lm(dat[1,]~dat[2,]))
		abline(a=coefs[1],b=coefs[2], col=rgb(range01(ijdf$unid)[k],0,0),lwd=3)#lty=ifelse(rho>0.85,1,2),lwd=3) 
		}
}
dev.off()
pdf("scat_micrD.pdf",height=3,width=5)
par(mfrow=c(5,9))
par(mar=c(0,0,0,0))
par(oma=c(0,0,0,0))
for(k in 1:nrow(ijdf)){
	dat <- rbind( DmnsinM[[2]][ , ijdf$i[k] ], DmnsinM[[2]][ ,ijdf$j[k] ] )
	plot(dat[1,]~dat[2,],pch=1,xaxt="n",yaxt="n", ylim=bufferX(range(DmnsinM[[2]]),.1),xlim=bufferX(range(DmnsinM[[2]]),.1),lty=1,col=rgb(range01(ijdf$unid)[k],0,0))
	rho <- cor(dat[1,],dat[2,])
	if(rho > sort(UnifracGCdat[[2]]$GC,decreasing=T)[23]){ 
		coefs <-coef(lm(dat[1,]~dat[2,]))
		abline(a=coefs[1],b=coefs[2], col=rgb(range01(ijdf$unid)[k],0,0),lwd=3)
		}
}
dev.off()
pdf("scat_col.pdf",height=3,width=5)
par(mfrow=c(5,9))
par(mar=c(0,0,0,0))
par(oma=c(0,0,0,0))
for(k in 1:nrow(ijdf)){
	dat <- rbind( DmnsinM[[3]][ , ijdf$i[k] ], DmnsinM[[3]][ ,ijdf$j[k] ] )
	plot(dat[1,]~dat[2,],pch=1,xaxt="n",yaxt="n", ylim=bufferX(range(DmnsinM[[3]]),.1),xlim=bufferX(range(DmnsinM[[3]]),.1),lty=1,col=rgb(range01(ijdf$unid)[k],0,0))
	rho <- cor(dat[1,],dat[2,])
	if(rho > sort(UnifracGCdat[[3]]$GC,decreasing=T)[23]){ 
		coefs <-coef(lm(dat[1,]~dat[2,]))
		abline(a=coefs[1],b=coefs[2], col=rgb(range01(ijdf$unid)[k],0,0),lwd=3)
		}
}
dev.off()
pdf("scat_agg.pdf",height=3,width=5)
par(mfrow=c(5,9))
par(mar=c(0,0,0,0))
par(oma=c(0,0,0,0))
for(k in 1:nrow(ijdf)){
	dat <- rbind( DmnsinM[[4]][ , ijdf$i[k] ], DmnsinM[[4]][ ,ijdf$j[k] ] )
	plot(dat[1,]~dat[2,],pch=1,xaxt="n",yaxt="n", ylim=bufferX(range(DmnsinM[[4]]),.1),xlim=bufferX(range(DmnsinM[[4]]),.1),lty=1,col=rgb(range01(ijdf$unid)[k],0,0))
	rho <- cor(dat[1,],dat[2,])
	if(rho > sort(UnifracGCdat[[4]]$GC,decreasing=T)[23]){ 
		coefs <-coef(lm(dat[1,]~dat[2,]))
		abline(a=coefs[1],b=coefs[2], col=rgb(range01(ijdf$unid)[k],0,0),lwd=3)
		}
}
dev.off()
pdf("scat_rnd.pdf",height=3,width=5)
par(mfrow=c(5,9))
par(mar=c(0,0,0,0))
par(oma=c(0,0,0,0))
for(k in 1:nrow(ijdf)){
	dat <- rbind( DmnsinM[[5]][ , ijdf$i[k] ], DmnsinM[[5]][ ,ijdf$j[k] ] )
	plot(dat[1,]~dat[2,],pch=1,xaxt="n",yaxt="n", ylim=bufferX(range(DmnsinM[[5]]),.1),xlim=bufferX(range(DmnsinM[[5]]),.1),lty=1,col=rgb(range01(ijdf$unid)[k],0,0))
	rho <- cor(dat[1,],dat[2,])
	if(rho > sort(UnifracGCdat[[5]]$GC,decreasing=T)[23]){ 
		coefs <-coef(lm(dat[1,]~dat[2,]))
		abline(a=coefs[1],b=coefs[2], col=rgb(range01(ijdf$unid)[k],0,0),lwd=3)
		}
}
dev.off()


###We do not have genotypic information for the duckweeds, but we might imagine they could follow an isolation by distance pattern.
#this small snippet checks if the effects of microbiomes on each genotype are more correlated for genotypes that are collected from closer locations
library(geosphere)
pairwiseD <- matrix(sapply(1:10, function(x) sapply(1:10, function(y) distHaversine(c(poptab$lon[x],poptab$lat[x]),c(poptab$lon[y],poptab$lat[y])))),nrow=10)
MmnsinD <- lapply(1:ncol(tfitmat), function(z) sapply(sort(unique(duckdat$duck)),function(d) 
		tapply(tfitmat[duckdat$duck== d,z],duckdat$micr[duckdat$duck==d],mean,na.rm=T)  ) )# microbes are rows, duckweeds are columns each trait is own matrix
MEffCor <- lapply(MmnsinD, function(z) cor(z)) #tables of microbial community effect correlations across pairwise duckweed genotype treatment, 1 matrix for each trait. rows and columns are duckweed origins
HaversineMEffCordat <- lapply(1:5, function(z) data.frame(ME = MEffCor[[z]][lower.tri(MEffCor[[z]])==T], distance = pairwiseD[lower.tri(pairwiseD)==T]) )
set.seed(100)
DistvMEffmods <-	lapply(1:5, function(z) MCMCglmm(ME ~ distance, data=HaversineMEffCordat[[z]], burnin=1000 , nitt =101000 , thin=10, verbose=F ) )
lapply(DistvMEffmods,summary) #non sig
sapply(1:5, function(z) cor(MEffCor[[z]][lower.tri(MEffCor[[z]])==T], pairwiseD[lower.tri(pairwiseD)==T])) #all correlations low




##############
###INVESTIGATE MICROBIOME COMPOSITION
############
##read in full taxonomic data
library(iNEXT)

feat.tab.qz 	<- read_qza("R inputs/Lemna_GxGxE_sepp_keptfeatures_table_denM.qza") #is filtered to rm streptophyta
feat.tabFM.qz 	<- read_qza("R inputs/Lemna_GxGxE_sepp_keptfeatures_table_den.qza") #is filtered to rm streptophyta
feat.tax.qz 	<- read_qza("R inputs/Lemna_GxGxE_taxonomy_den.qza")##NOT FILTERED TO kept features, or even to rm streptophyta (chloroplast, which are removed from the dataset, along with mitochondrial sequences)

Ltab.s	<- feat.tab.qz$data[,order(colnames(feat.tab.qz$data))] #inocula only
LtabFM.s	<- feat.tabFM.qz$data[,order(colnames(feat.tabFM.qz$data))] #both field and inocula
##sequences remaining after ALL filtering are just sums of LtabFM
 sum(LtabFM.s)#1,026,828
 range(colSums(LtabFM.s)) #19,239 84,116
 mean(colSums(LtabFM.s))# 51,341.4
 mean(colSums(LtabFM.s)[c(1,3,5,7,9)])# [1] 32,024 field
 mean(colSums(LtabFM.s)[c(2,4,6,8,10)]) #[1] 67,994.6 inocula
 
#split tax info into columns
 feat.tax 	<- feat.tax.qz$data ##ALL taxa - this is not filtered to taxa that occur in samples, or even to remove streptophyta
 taxlevs <- sapply(1:nrow(feat.tax), function(z) strsplit(as.character(feat.tax$Taxon[z]),split="; "))
 tax <- matrix(NA,ncol=7,nrow=nrow(feat.tax))
 for(i in 1:nrow(tax)){ tax[i,1:length(taxlevs[[i]])]<- taxlevs[[i]]}
 rownames(tax) <- as.character(feat.tax$Feature.ID)
##subset taxonomy ID to taxa that occur in samples (and that are not chloroplast, etc), for both inocula only and full set of communities
 Ltax.s <- 	 	tax[sapply(1:nrow(Ltab.s), function(z) which(rownames(tax)==rownames(Ltab.s)[z])),]
 LtaxFM.s <- 	 	tax[sapply(1:nrow(LtabFM.s), function(z) which(rownames(tax)==rownames(LtabFM.s)[z])),]#
 Ltax <- Ltax.s
 Ltab <- Ltab.s[rownames(Ltax.s),]
 rm(feat.tab.qz,feat.tax.qz, feat.tabFM.qz, feat.tax, taxlevs,tax)
#Set column sums to be 1, since 16s data is inherently proportional/relative
 Lprp <- Ltab
 for(i in 1:ncol(Ltab)){Lprp[,i] <- Ltab[,i]/sum(Ltab[,i])}
 LprpFM <- LtabFM.s
 for(i in 1:ncol(LtabFM.s)){LprpFM[,i] <- LtabFM.s[,i]/sum(LtabFM.s[,i])}
#Check that proportions and detailed taxonomy ID dataframes are organized the same
 sum(rownames(Ltax)==rownames(Lprp) )
 dim(Lprp)
#quick stats for ms
 colSums(sign(Lprp))
#   Ars.M.101   BrBr.M.105    Cam.M.107    Cdv.M.109 HumBPE.M.119    KSR.M.121   MocT.M.129    Rat.M.133    UTM.M.145  WestD.M.148 
#           33           25           24           38           21           26           38           48           41           24 

 ##coverage estimation
sapply(1:ncol(LtabFM.s), function(z) DataInfo(as.vector(LtabFM.s[ which(LtabFM.s[,z]>0) ,z]),datatype="abundance")$SC) 
# [1] 0.9994 1.0000 0.9999 1.0000 0.9994 1.0000 1.0000 0.9999 0.9999 1.0000 0.9998 1.0000 0.9998 1.0000 0.9997 1.0000 0.9998 1.0000 0.9997 1.0000


#field-inocula comparison
mtab <- LtabFM.s[,seq(from = 2, to = 20, by =2)]
ftab <- LtabFM.s[,seq(from = 1, to = 19, by =2)]
sum(sign(rowSums(mtab)))#232
sum(sign(rowSums(ftab)))#3949
sum( rowSums(ftab) >0 & rowSums(mtab)>0 ) #68 overlap from field and master
FcapanyM <- colSums(ftab[ rowSums(ftab) >0 & rowSums(mtab)>0,])/colSums(ftab[ ,])#percent of field community reads belonging to taxa in any inocula community, ranges from 1.1%-16%
FcapbyM <- sapply(1:10, function(z) sum(ftab[ (mtab[,z]) >0,z])/sum(ftab[ ,z])) #percent of field community reads belonging to taxa in the matched inocula community, 0-8.3%
MpresinanyF <- colSums(mtab[ rowSums(ftab) >0 ,])/colSums(mtab[ ,])#percent of reads in inocula communities present in any field community, ranges from 1.3-96%
colSums(sign(mtab[ rowSums(ftab) >0 ,]))/colSums(sign(mtab[ ,]))# percent of ASV found in any field sample
#   Ars.M.101   BrBr.M.105    Cam.M.107    Cdv.M.109 HumBPE.M.119    KSR.M.121   MocT.M.129    Rat.M.133    UTM.M.145  WestD.M.148 
#    0.3030303    0.3600000    0.4583333    0.1842105    0.5238095    0.3846154    0.5263158    0.5000000    0.4146341    0.5416667 
MpresinF <- sapply(1:10, function(z) sum(mtab[ (ftab[,z]) >0,z])/sum(mtab[ ,z])) # same, but for matched field site, ranges from 0-94%; 
sapply(1:10, function(z) sum(sign(mtab[ (ftab[,z]) >0,z]))/sum(sign(mtab[ ,z])))#now number of ASV also in matched site / total ASV again
# 0.06060606 0.00000000 0.12500000 0.00000000 0.23809524 0.07692308 0.31578947 0.33333333 0.07317073 0.33333333 

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
sapply(1:10, function(z) sum(unlist(fieldug[z])%in%unlist(fieldug[-z])) )# how many genera in one field site are found in any of the other 9 (note this still includes "g__" and "NA")
#  91  51  93  13  82  81 118 112  68  74
 sapply(1:10, function(z) sum(unlist(fieldug[z])%in%unlist(fieldug[-z])) )/colSums(sign(ftab))
 #   Ars.F.13   BrBr.F.20    Cam.F.22    Cdv.F.26 HumBPE.F.40    KSR.F.43   MocT.F.57    Rat.F.68    UTM.F.90  WestD.F.93 
#  0.08081705  0.18411552  0.14553991  0.30952381  0.20654912  0.12272727  0.15817694  0.13542926  0.15813953  0.14453125 

mean(MpresinF); mean(MpresinanyF); mean(MorGpresinF); mean(MorGpresinanyF)
# [1] 0.158025
# [1] 0.5232412
# [1] 0.5237108
# [1] 0.7959689

#continue with more in depth consideration of the taxa across communites
sum(rowSums(mtab)>0 ) #232 total in inocula communities
sapply(0:10, function(z) length(which(rowSums(sign(mtab)) == z))) # how many entries in mtab occur in 0, 1, ...all 10 of the inocula?
# [1] 3881  187   19   18    3    3    2    0    0    0    0
sum(length(which(rowSums(sign(mtab)) >= 2)) ) #45 occur in 2 or more inocula; 
sum(rowSums(sign(mtab)) >= 2 & rowSums(ftab)>0 )  # 31 of which occur in ftab
sum(length(which(rowSums(sign(mtab)) >= 3)) )  #26 occur in 3 or more inocula
sum(rowSums(sign(mtab)) >= 3 & rowSums(ftab)>0 )  # 19 of which occur in ftab
sum(length(which(rowSums(sign(mtab)) >= 4)) ) #8 occur in 4 or more inocul
sum(rowSums(sign(mtab)) >= 4 & rowSums(ftab)>0 )  # all but one of which occur in ftab
sapply(1:10, function(z) sum(ftab[,z]>0 & mtab[,z] > 0) ) # found in mtab column and corresponding field site
# [1]  2  0  3  0  5  2 12 16  3  8
sapply(1:10, function(z) sum(mtab[,z] > 0 & rowSums(ftab)>0) ) # found in mtab column and any field sample
# [1] 10  9 11  7 11 10 20 24 17 13
sapply(1:10, function(z) sum(mtab[,z] > 0) ) # found in mtab column
# [1] 33 25 24 38 21 26 38 48 41 24
unique(LtaxFM.s[rowSums(mtab)>0,6] [which(!LtaxFM.s[rowSums(mtab)>0,6]%in%LtaxFM.s[rowSums(ftab)>0,6])]) #genus names in inocula but not field
# [1] "g__Bacillus"          "g__Tolumonas"         "g__Rahnella"          "g__Achromobacter"     "g__Streptomyces"     
#  [6] "g__Gluconacetobacter" "g__Lysinibacillus"    "g__Ancylobacter"      "g__Averyella"         "g__Kaistia"          
# [11] "g__Serratia"         
unique(LtaxFM.s[rowSums(mtab)>0,6] [which(LtaxFM.s[rowSums(mtab)>0,6]%in%LtaxFM.s[rowSums(ftab)>0,6])]) # genus names in inocula and field
#  [1] "g__Shewanella"        "g__"                  "g__Rheinheimera"      NA                     "g__Azospirillum"     
#  [6] "g__Flavobacterium"    "g__Pseudomonas"       "g__Janthinobacterium" "g__Acidovorax"        "g__Stenotrophomonas" 
# [11] "g__Agrobacterium"     "g__Hydrogenophaga"    "g__Rhodoferax"        "g__Aeromonas"         "g__Paucibacter"      
# [16] "g__Variovorax"        "g__Flectobacillus"    "g__Rubrivivax"        "g__Kineosporia"       "g__Phenylobacterium" 
# [21] "g__Pedobacter"        "g__Caulobacter"       "g__Reyranella"        "g__Novosphingobium"   "g__Opitutus"         
# [26] "g__Chryseobacterium"  "g__Rhodococcus"       "g__Rhizobium"         "g__Sphingomonas"      "g__Dyadobacter"      
# [31] "g__Mycobacterium"     "g__Elstera"           "g__Mycoplana"         "g__Acinetobacter"     "g__Luteimonas"       
# [36] "g__Runella"           "g__Methylibium"       "g__Devosia"           "g__Dechloromonas"     "g__Pelomonas"        
# [41] "g__Pseudoxanthomonas" "g__Rhodobacter"       "g__Leptothrix"        "g__Exiguobacterium"  
length(unique(LtaxFM.s[rowSums(mtab)>0,6] [which(LtaxFM.s[rowSums(mtab)>0,6]%in% LtaxFM.s[rowSums(ftab)>0,6])]))-2 # same as above but count, here minus 2 for NA and "g__"
#42
length(unique(LtaxFM.s[rowSums(mtab)>0,6] [which(!LtaxFM.s[rowSums(mtab)>0,6]%in% LtaxFM.s[rowSums(ftab)>0,6])])) 
#so of 42 named genera in mtab, 31 are in in ftab
unique(LtaxFM.s[rowSums(mtab)>0,5] [which(!LtaxFM.s[rowSums(mtab)>0,5]%in% LtaxFM.s[rowSums(ftab)>0,5])])#same but for families
unique(LtaxFM.s[rowSums(mtab)>0,5] [which(LtaxFM.s[rowSums(mtab)>0,5]%in% LtaxFM.s[rowSums(ftab)>0,5])]) #38 families in inocula, 2 of those not in field (excluding "f__" and NA).
GnotinF <- unique(LtaxFM.s[rowSums(mtab)>0,6] [which(!LtaxFM.s[rowSums(mtab)>0,6]%in% LtaxFM.s[rowSums(ftab)>0,6])])
colSums(mtab[LtaxFM.s[,6]%in%GnotinF,]) / colSums(mtab) # what is the abundance, by sample, of genera in inocula that are not in any field community
#   Ars.M.101   BrBr.M.105    Cam.M.107    Cdv.M.109 HumBPE.M.119    KSR.M.121   MocT.M.129    Rat.M.133    UTM.M.145  WestD.M.148 
# 0.8133292120 0.0197958643 0.0000000000 0.1193128974 0.0001089963 0.0000000000 0.0000000000 0.0000473216 0.0610013331 0.0103336707 
GnotinFtab <- (LtabFM.s[LtaxFM.s[,6]%in%GnotinF,]) # we can see it is mainly 1 genus (and 1 ASV), in Arsandco and Cedarvale inocula
LtaxFM.s[rownames( GnotinFtab)[sapply(1:nrow(GnotinFtab), function(z) any(GnotinFtab[z,]/(colSums(mtab))>0.01) )], ] # There are 3 ASV over 1% where the genus is not present in the field, two are a bacillus


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
#pirellulaceae in field and inocula
 simplerfamsum[simplerfamsumnames=="Pirellulaceae",c(1,3,5,7,9,11,13,15,17,19)]
# [1] 0.041066710 0.019588384 0.078202730 0.008940174 0.025229042 0.175149775 0.040415779 0.083446712 0.059381864 0.245758696
 simplerfamsum[simplerfamsumnames=="Pirellulaceae",c(2,4,6,8,10,12,14,16,18,20)]
# [1] 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0002208341 0.0000000000 0.0000000000
# Rhodobacteriaceae
 simplerfamsum[simplerfamsumnames=="Rhodobacteraceae",c(1,3,5,7,9,11,13,15,17,19)]
# [1] 0.06846249 0.11918433 0.15529112 0.01237071 0.29589787 0.18148735 0.23534444 0.14110224 0.25556950 0.14482040
 simplerfamsum[simplerfamsumnames=="Rhodobacteraceae",c(2,4,6,8,10,12,14,16,18,20)]
# [1] 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0001892864 0.0000000000 0.0000000000
#Aeromonadaceae
 simplerfamsum[simplerfamsumnames=="Aeromonadaceae",c(1,3,5,7,9,11,13,15,17,19)]
#[1] 0.0001617861 0.0000000000 0.0003581804 0.0001559333 0.0003576636 0.0001237808 0.0050358319 0.0023858819 0.0004706621 0.0016428143
 simplerfamsum[simplerfamsumnames=="Aeromonadaceae",c(2,4,6,8,10,12,14,16,18,20)]
#[1] 7.133007e-05 8.380298e-01 1.694534e-01 0.000000e+00 0.000000e+00 9.482183e-01 1.667918e-04 9.312102e-01 1.106799e-01 9.474337e-01
#Pseudomonadaceae
 simplerfamsum[simplerfamsumnames=="Pseudomonadaceae",c(1,3,5,7,9,11,13,15,17,19)]
# [1] 0.0002966079 0.0016303963 0.0020694870 0.0022870212 0.0031639475 0.0012873199 0.0051972367 0.0073745440 0.0009021023 0.0013584810
 simplerfamsum[simplerfamsumnames=="Pseudomonadaceae",c(2,4,6,8,10,12,14,16,18,20)]
#  [1] 0.0032692948 0.1407035915 0.8280215600 0.2830407137 0.0001089963 0.0430150117 0.5947168710 0.0639472522 0.5823222693 0.0409640958

 #set color
library(Polychrome)
# newpal <- createPalette(nrow(simplerfamsum)+5, c(rgb(0.7,0,0),rgb(1,1,0),rgb(0,0,0.7)), M=1000)
siteorder <- c(order(poptab$Nimperv)*2-1,order(poptab$Nimperv)*2)
# save(newpal, file="colorpal.RData")
load(file="R inputs/colorpal.RData")

#interpreting unifrac
#pairwise distance for rel. abund of each family -- using only abundant families
famdist <- lapply(1:nrow(simplerfamsum), function(fam) sapply(1:ncol(simplerfamsum), function(samp) abs(simplerfamsum[fam,]-simplerfamsum[fam,samp] ) ) )
#compare that to the unifrac distance matrix
mant.fam.uni <- lapply(famdist, function(fam) mantel(wunfrc.B,fam))
mant.famuni.p <- unlist(lapply(mant.fam.uni, function(fam) fam$signif))
mant.famuni.r <- unlist(lapply(mant.fam.uni, function(fam) fam$statistic))
pdf("desc_unifrc_fam_den.pdf",height=4.5,width=6)
par(mar=c(10,4,2,1))
plot(mant.famuni.r[-c(3,20)]~c(1:18),ylab="",xlab="",xaxt="n",ylim=c(-0.2,1),pch=16)
	axis(at=c(1:18),labels=simplerfamsumnames[-c(3,20)],las=2,side=1)
	text((1:18)[ mant.famuni.p[-c(3,20)]<0.01 ],-0.15,"*")
	mtext("Family relative abundance difference vs. weighted unifrac distance", side =3, line=0.5)
	mtext("Correlation", side =2, line=2.5)
	abline(h=0,col=rgb(0.5,0.5,0.5,alpha=0.5),lty=3)
dev.off()

uniDdat <- data.frame(uniD=wuncdat[wuncdat$matched==1,1], nimp=poptab$Nimperv)
set.seed(10)
uniDm <- (MCMCglmm(uniD ~ nimp , family="gaussian", data=uniDdat,verbose=FALSE,pr=TRUE,nitt=101000, thin=10, burnin=1000)) #all kept
uniDs<-uniDm$Sol
n.s <- seq(from = min(uniDdat$nimp),to=max(uniDdat$nimp),length.out=1000)
uniD.p <- sapply(n.s, function(z) mean(uniDs[,1]+uniDs[,2]*z) )
uniD.ci <- sapply(n.s, function(z) HPDinterval(as.mcmc(uniDs[,1]+uniDs[,2]*z),prob=.95) )

##main microbiome composition figure
pdf("ladpt_microbiome_sim_den.pdf",height=6,width=9)
layout(matrix(c(1,1,2,3,4,2),nrow=2, byrow = TRUE),heights=c(3,3),widths=c(1.1,.75,.4))
par(oma=c(0,0,2,0))
par(mar=c(6,4,0,0))
bg <- barplot(simplerfamsum[,siteorder], col= newpal,ylab="",xlab="")#,names.arg=rep(round(poptab$Nimperv),times=2) )
axis(  at=c( (bg[5]+bg[6])/2, (bg[15]+bg[16])/2), side=3, labels = c("Field","Inocula") ,lty=0,line=-0.5,cex.axis=1.25)
axis( at = bg, side=1, labels = rep(sort(round(poptab$Nimperv)),times=2),lty=0, line=-1)
abline(v = (bg[10]+bg[11])/2, lty=1, lwd=3)
mtext("% permeable surface area",side=1,line=4.5)
mtext("a.",side=3, adj = -0.075, line = 0.5) #,adj=-0.12,line=-0.5)
mtext("Proportion",side=2,line=2)
par(mar=c(0,0,0,0))
plot(rep(1,times=3)~c(1:3), pch=NA,bty="n",xaxt="n",yaxt="n",ylab="",xlab="")
legend(x=0.9,y=1.4,legend=simplerfamsumnames,fill=newpal,bty="n",ncol=1,cex=1.1)
par(mar=c(6,4,1,2))
bg<- barplot(100*t(cbind(MpresinF,MorGpresinF,MpresinanyF,MorGpresinanyF))[,order(poptab$Nimperv)],beside=T,xlab="",
 		col=c("blue","gray"),density=c(100,100,50,50), ,names.arg=rep("",times=10),ylim=c(0,103))
axis( at = (bg[2,]+bg[3,])/2, side=1, labels = sort(round(poptab$Nimperv)),lty=0, line=-1)
mtext("% represented in field",side=2,line=2)
mtext("% permeable surface area",side=1,line=4.5)
mtext("b.",side=3,adj=-0.125,line = 0.25)
legend(x=-2,y=108,legend=c("ASV, matched site","Genus, matched","ASV, any site","Genus, any site"),fill=c("blue","gray","blue","gray"),density=c(100,100,50,50),bty="n",ncol=1,cex=1)
plot(wuncdat[wuncdat$matched==1,1]~poptab$Nimperv,pch=NA,ylab="",xlab="", ,cex.axis=1.25,ylim=c(0.25,0.9))#,
polygon(c(n.s,rev(n.s)),c(uniD.ci[1,],rev(uniD.ci[2,])),col=rgb(0,0,0,alpha=0.25),border=NA)
lines(uniD.p~n.s)
points(uniDdat$uniD~uniDdat$nimp,pch=16,cex=2)
 mtext("Field-inocula Unifrac distance",side = 2, line = 2.5,at=0.525)
 mtext("% permeable surface area", side =1, line = 4.5)
mtext("c.",side=3,adj=-0.2,line = 0.25)
 dev.off()


##############
##does microbiome composition explain any effects on fitness/traits?
##############
Mfamsum <- famsum[,c(2,4,6,8,10,12,14,16,18,20)]
Mfamsum <- Mfamsum[sign(rowSums(Mfamsum))>0,]
Mfams <- as.data.frame(t(Mfamsum))
colnames(Mfams) <- prettyfamnames[ sign(rowSums(famsum[,c(2,4,6,8,10,12,14,16,18,20)]))>0  ]

duckdat$aeros <- sapply(1:nrow(duckdat), function(z) Mfams$Aeromonadaceae[which(poptab$pop==duckdat$micr[z])])
duckdat$pseudos <- sapply(1:nrow(duckdat), function(z) Mfams$Pseudomonadaceae[which(poptab$pop==duckdat$micr[z])])

##placing into models of G, G, E, GG, GE, GE, GGE above. , seeing if best fit is improved in DIC from using micrNImpv.
###iterative term removal was performed, but not every model is pasted below, since there are many, but no terms of interest remain after removing n.s. terms iteratively.
##instead only the most complex model is printed, and the code user can perform iterative term removal themselves to see which model is selected
summary(MCMCglmm(sqmmE ~ x + I(x^2) + y + I(y^2) + col + row + I(row^2) +
 duckNImpv + aeros + numzinc + duckNImpv:numzinc + aeros:numzinc + duckNImpv:aeros + duckNImpv:aeros:numzinc, 
 family="gaussian", data=duckdat,verbose=FALSE,pr=TRUE,nitt=11000, thin=10, burnin=1000)) 
summary(MCMCglmm(sqmmE ~ x + I(x^2) + y + I(y^2) + col + row + I(row^2) +
 duckNImpv + pseudos + numzinc + duckNImpv:numzinc + pseudos:numzinc + duckNImpv:pseudos + duckNImpv:pseudos:numzinc, 
 family="gaussian", data=duckdat,verbose=FALSE,pr=TRUE,nitt=11000, thin=10, burnin=1000)) 
#all microbe terms drop out

summary(MCMCglmm(invOD ~ x + I(x^2) + y + I(col^2) + row + I(row^2) +
 duckNImpv + aeros + numzinc + duckNImpv:numzinc + aeros:numzinc + duckNImpv:aeros + duckNImpv:aeros:numzinc, 
 family="gaussian", data=duckdat,verbose=FALSE,pr=TRUE,nitt=11000, thin=10, burnin=1000)) 
#removing n.s. terms, some terms with aeros sig, but drop out eventually, and model fit worse than model with micrNImpv above
summary(MCMCglmm(invOD ~ x + I(x^2) + y + I(col^2) + row + I(row^2) +
 duckNImpv + pseudos + numzinc + duckNImpv:numzinc + pseudos:numzinc + duckNImpv:pseudos + duckNImpv:pseudos:numzinc, 
 family="gaussian", data=duckdat,verbose=FALSE,pr=TRUE,nitt=11000, thin=10, burnin=1000)) 
#pseudos terms all drop out

summary(MCMCglmm(colint ~ x + I(x^2) + y + I(y^2) + col + I(col^2) + row + I(row^2) +
 duckNImpv + aeros + numzinc + duckNImpv:numzinc + aeros:numzinc + duckNImpv:aeros + duckNImpv:aeros:numzinc, 
 family="gaussian", data=duckdat,verbose=FALSE,pr=TRUE,nitt=11000, thin=10, burnin=1000)) 
summary(MCMCglmm(colint ~ x + I(x^2) + y + I(y^2) + col + I(col^2) + row + I(row^2) +
 duckNImpv + pseudos + numzinc + duckNImpv:numzinc + pseudos:numzinc + duckNImpv:pseudos + duckNImpv:pseudos:numzinc, 
 family="gaussian", data=duckdat,verbose=FALSE,pr=TRUE,nitt=11000, thin=10, burnin=1000)) 
 #removing n.s. terms keeps some aeros and pseudos effects, but neither better DIC than micrNImpv model

summary(MCMCglmm(areapperE ~ x + I(x^2) + y + I(y^2) + col + I(col^2) + row + I(row^2) +
 duckNImpv + aeros + numzinc + duckNImpv:numzinc + aeros:numzinc + duckNImpv:aeros + duckNImpv:aeros:numzinc, 
 family="gaussian", data=duckdat,verbose=FALSE,pr=TRUE,nitt=11000, thin=10, burnin=1000)) 
summary(MCMCglmm(areapperE ~ x + I(x^2) + y + I(y^2) + col + I(col^2) + row + I(row^2) +
 duckNImpv + pseudos + numzinc + duckNImpv:numzinc + pseudos:numzinc + duckNImpv:pseudos + duckNImpv:pseudos:numzinc, 
 family="gaussian", data=duckdat,verbose=FALSE,pr=TRUE,nitt=11000, thin=10, burnin=1000))
 #removing n.s. terms keeps some aeros and pseudos effects, but neither better than micrNImpv model

summary(MCMCglmm(meanround ~ x + I(x^2) + y + I(y^2) + col + I(col^2) +
 duckNImpv + aeros + numzinc + duckNImpv:numzinc + aeros:numzinc + duckNImpv:aeros + duckNImpv:aeros:numzinc, 
 family="gaussian", data=duckdat,verbose=FALSE,pr=TRUE,nitt=11000, thin=10, burnin=1000)) 
#aeros terms all drop out
summary(MCMCglmm(meanround ~ x + I(x^2) + y + I(y^2) + col + I(col^2) +
 duckNImpv + pseudos + numzinc + duckNImpv:numzinc + pseudos:numzinc + duckNImpv:pseudos + duckNImpv:pseudos:numzinc, 
 family="gaussian", data=duckdat,verbose=FALSE,pr=TRUE,nitt=11000, thin=10, burnin=1000)) 
summary(MCMCglmm(meanround ~ x + I(x^2) + y + I(y^2) + col + I(col^2) +
  pseudos +   duckNImpv:pseudos , 
 family="gaussian", data=duckdat,verbose=FALSE,pr=TRUE,nitt=11000, thin=10, burnin=1000))
#this model with pseudos terms fits better than the model selected when using micrNImpv (a n.s. duckNImpv effect alone)

# refit single model that finds composition predicts effects
set.seed(100)
round.psd <- MCMCglmm(meanround ~ x + I(x^2) + y + I(y^2) + col + I(col^2) +
  pseudos +   duckNImpv:pseudos , 
 family="gaussian", data=duckdat,verbose=FALSE,pr=TRUE,nitt=101000, thin=10, burnin=1000)
summary(round.psd)

#calculate model predictions
Psd.s <- seq(from=min(duckdat$pseudos),to=max(duckdat$pseudos),length.out=1000)
Psd.r <- range(Mfams$Pseudomonadaceae)
rndPsol <- round.psd$Sol
rnd.dI.ci <- lapply(Psd.r, function(micrP) sapply(NI.s, function(duckP)    HPDinterval(as.mcmc(
	rndPsol[,1] + rndPsol[,2]*mx + rndPsol[,3]*mx2 + rndPsol[,4]*my + rndPsol[,5]*my2 + rndPsol[,6]*mcol + rndPsol[,7]*mcol2+ 
 		rndPsol[,8]*micrP + rndPsol[,9]*duckP*micrP ),prob=0.95) )) 
rnd.dI.mn <- lapply(Psd.r, function(micrP) sapply(NI.s, function(duckP) mean(
	rndPsol[,1] + rndPsol[,2]*mx + rndPsol[,3]*mx2 + rndPsol[,4]*my + rndPsol[,5]*my2 + rndPsol[,6]*mcol + rndPsol[,7]*mcol2+ 
 		rndPsol[,8]*micrP + rndPsol[,9]*duckP*micrP ) )) 
meanspacernd <- mean(	rndPsol[,1] + rndPsol[,2]*mx + rndPsol[,3]*mx2 + rndPsol[,4]*my + rndPsol[,5]*my2 + rndPsol[,6]*mcol + rndPsol[,7]*mcol2  )
meanspacerndR <- duckdat$meanround - 
			sapply(1:nrow(duckdat), function(z) mean(rndPsol[,1] + rndPsol[,2]*duckdat$x[z] + rndPsol[,3]*(duckdat$x[z]^2) + rndPsol[,4]*duckdat$y[z] + rndPsol[,5]*(duckdat$y[z]^2)+ rndPsol[,6]*duckdat$col[z] + rndPsol[,7]*(duckdat$col[z]^2) )) +
			meanspacernd
rndmn <- tapply(meanspacerndR,paste(duckdat$pseudos,duckdat$duckNImpv),mean,na.rm=T)

#plot figure
pdf("roundness_pseudos_den.pdf",height=4,width=3.5)
par(mar=c(5.5,3.5,0.5,0.5))
suppressWarnings(plot(1~1,pch=NA,xlim=bufferX(NI.s,0.05),ylim=bufferX(rndmn,0.0),ylab="",xlab="",main="" ))#R dislikes formula 1~1
	polygon (x= c(NI.s, rev(NI.s) ), y= c( rnd.dI.ci[[1]][1,], rev(rnd.dI.ci[[1]][2,]) ) ,col=rgb(0.5,0.5,0.5,alpha=0.2),border=NA)
  	lines(rnd.dI.mn[[1]]~NI.s,col=rgb(0.5,0.5,0.5))
 	polygon (x= c(NI.s, rev(NI.s) ), y= c( rnd.dI.ci[[2]][1,], rev(rnd.dI.ci[[2]][2,]) ) ,col=rgb(0,0,1,alpha=0.2),border=NA)
  	lines(rnd.dI.mn[[2]]~NI.s,col=rgb(0,0,1))
	mtext("Frond roundness",side=2,line=2)
	mtext("%Permeable surface, duckweed site",side=1,line=5.5, at = 55)
  	points(rndmn[c(1:10,91:100)]~niDmn[c(1:10,91:100)],pch=16, col=rep(c(rgb(0.5,0.5,0.5),rgb(0,0,1)),each=10)  ) #now color is pseudomonadaceae amount
	axis(at=c(bufferX(NI.s,.075)), label=c("most urban","most rural"),line=1,tick=FALSE,las=2,side=1,cex.axis=0.75)
	legend(62,0.50, c("low %Pseudomonadaceae","high %Pseudomonadaceae"), fill = c(rgb(0.5,0.5,0.5),rgb(0,0,1)) , bty="n")
dev.off()