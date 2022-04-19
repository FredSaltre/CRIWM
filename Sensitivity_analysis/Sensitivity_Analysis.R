##########################################################
## Simulate time series & apply extinction time algorithms
##########################################################
rm(list=ls(all=TRUE))

Tseries <- function(nb, Text, maxage, prob.dist) # GENERATING TIME SERIES DEPENDING ON A GIVEN DISTRIBUTION  
{
  if (prob.dist == "uniform") {
    TS_uni <- sort(unique(round(runif(nb,Text,maxage)))) # random uniform
    nTS_uni <- length(TS_uni)
    #plot(TS_uni, rep(0,nTS_uni),pch=4,xlim=c(Text,maxage),ylab="",xlab="",cex=1.5,yaxt="n")
    #abline(v=Text,lty=1,lwd=4)
    out<-TS_uni
    }
  
  if (prob.dist == "linear") {
    ## calculate the probability of distribution for linear +=====
    nbTS<-maxage-Text;timevec<-seq(Text,maxage,1);prlin.vec <- seq(0,1,1/nbTS)#timevec = all age of the time series incremented by 1
    problinear <- sort(prlin.vec/sum(prlin.vec)) # probability distribution linear
    ##====================================
    TS_lin <- sort(sample(timevec,nb,replace=F,prob=problinear)) 
    nTS_lin <- length(TS_lin)
    #plot(TS_lin, rep(0,nTS_lin),pch=4,xlim=c(Text,maxage),ylab="",xlab="",cex=1.5,yaxt="n")
    #abline(v=Text,lty=1,lwd=4)
    out<-TS_lin
    }
  
  if (prob.dist == "sigm") {
    ## calculate the probability of distribution for sigmoidal +=====
    nbTS<-maxage-Text;timevec<-seq(Text,maxage,1);
    sig.st <- (timevec - mean(timevec))/(1.2*nbTS/10)
    prsig.vec <- 1/(1+exp(-sig.st))
    probs.sran <- prsig.vec - min(prsig.vec); probs.sran <- probs.sran/sum(probs.sran)
    ##====================================
    TS_sig <- sort(sample(timevec,nb,replace=F,prob=probs.sran)) # probability distribution sigmoidal (half linear)
    nTS_sig <- length(TS_sig)
    plot(TS_sig, rep(0,nTS_sig),pch=4,xlim=c(Text,maxage),ylab="",xlab="",cex=1.5, yaxt="n")
    abline(v=Text,lty=1,lwd=4)
    out<-TS_sig
    }
  
  if (prob.dist == "expon") {
    ## calculate the probability of distribution for exponential +=====
    nbTS<-maxage-Text;timevec<-seq(Text,maxage,1);
    prlin.vec <- seq(0,1,1/nbTS);fit.exp <- lm(prlin.vec~log(timevec,base=20))
    prexp.vec <- (20^(predict.lm(fit.exp)))/sum(20^(predict.lm(fit.exp)))
    probs.eran <- prexp.vec - min(prexp.vec);probs.eran <- probs.eran/sum(probs.eran)
    ##====================================
    TS_exp <- sort(sample(timevec,nb,replace=F,prob=probs.eran))
    nTS_exp <- length(TS_exp)
    #plot(TS_exp, rep(0,nTS_exp),pch=4,xlim=c(Text,maxage),ylab="",xlab="",cex=1.5, yaxt="n")
    #abline(v=Text,lty=1,lwd=4)
    out<-TS_exp
    }
  
  if (prob.dist == "log") {
    ## calculate the probability of distribution for logarithmic +=====
    nbTS<-maxage-Text;timevec<-seq(Text,maxage,1);prlin.vec <- seq(0,1,1/nbTS)
    fit.lg <- lm(prlin.vec~exp(timevec/sum(timevec)))
    prlg.vec <- (log(sum(timevec)*predict.lm(fit.lg)))/sum(log(sum(timevec)*predict.lm(fit.lg)))
    probs.lgran <- prlg.vec - min(prlg.vec);probs.lgran <- probs.lgran/sum(probs.lgran)
    ##====================================
    TS_log <- sort(sample(timevec,nb,replace=F,prob=probs.lgran))
    nTS_log <- length(TS_log)
    #plot(TS_log, rep(0,nTS_log),pch=4,xlim=c(Text,maxage),xlab="YBP",ylab="",cex=1.5,yaxt="n")
    #abline(v=Text,lty=1,lwd=4)
    out<-TS_log
    }
  return(out);
  rm(list=ls(all=TRUE)); 
}

TSerror<-function(dat,opt1,opt2)# ADDING ERRORS TO THE AGES
{
  if (opt1 == "proportional") {
    #the size of the error will increase from 10% on the first age to 60% on the last age
    ndat<-length(dat)
    dat2<-rbind(dat[1],dat[ndat]);edat2<-rbind(dat[1]*0.1,dat[ndat]*0.6)
    dat3<-cbind(dat2,edat2);colnames(dat3)<-c("time","error")
    dat3<-as.data.frame(dat3)
    fit<-lm(formula = error ~ time, data = dat3)
    edat2<-dat*fit$coefficients[2]+fit$coefficients[1]
    if (opt2 == "symetrical") {
      #same error both side of the age
      sd<-runif(ndat)*edat2
      sd<-cbind(sd,sd)
    }
    if (opt2 == "nonsymetrical") {
      #different error each side of the age
      sd1<-runif(ndat)*edat2;sd2<-runif(ndat)*edat2;
      sd<-cbind(sd1,sd2)
    }
  }
  if (opt1 == "random") {
    #the size of the error will random between 10% and 60% of the ages
    ndat<-length(dat)
    edat<-runif(ndat,min=10,max=60);
    edat2<-dat*(edat/100);
    if (opt2 == "symetrical") {
      #same error both side of the age
      sd<-runif(ndat)*edat2
      sd<-cbind(sd,sd)
    }
    if (opt2 == "nonsymetrical") {
      #different error each side of the age
      sd1<-runif(ndat)*edat2;sd2<-runif(ndat)*edat2;
      sd<-cbind(sd1,sd2)
    }
  }
  out<-cbind(dat,sd);colnames(out)<-c("age","rightSD","leftSD")
  return(out)
  rm(list=ls(all=TRUE)); 
}

TStat<-function(dat)#calculate the statistical chgaracterisitc of the time series
  # we focus on Parameters î, σ2i, ε, and σ2ε (Table 2 in Saltré et al. 2015)
{
  inter<-diff(dat[,1]);Allsd<-c(dat[,2],dat[,3]);
  MeanInt<-mean(inter);VarInt<-var(inter)
  MeanAllsd<-mean(Allsd);VarAllsd<-var(Allsd);
  out<-rbind(MeanInt,VarInt,MeanAllsd,VarAllsd)
  rownames(out)<-c("IntervalAverage","IntervalVariance","ErrorAverage","ErrorVAriance")
  return(out)
  rm(list=ls(all=TRUE)); 
}

##=====================================================================
## MAIN PROGRAM
##=====================================================================
##Simulate fossil records of 5,000 years each following Bradshaw et al. 2012 (original description of GRIWM). 
##I used this approach to assess modality in Appendix C.
##Parameters î, σ2i, ε, and σ2ε (Table 2 in Saltré et al. 2015) will be absolute whether you express them as percentages or not.
##the start year for any simulated 5,000-year-long fossil record should be placed randomly at 100 to 9,707 years BP

source("griwm_beta_RUN.R")
source("criwm_beta_RUN.R")
library(rlist)
library(lhs)

capage <- 14706  # oldest record considered = time series capped at 14,706 BP
maxnb<-100 #max number of sample
Windnb<-c(2,maxnb);WindText<-c(200,9707);
nm<-1000;npar<-4;simid<-randomLHS(nm, npar)
prob.dist <- c("uniform","linear","sigm","expon","log")
error.opt1 <- c("proportional","random")
error.opt2 <- c("symetrical","nonsymetrical")
KnownExt <- numeric(nm) 

##START OF THE LOOP OVER SIMULATION
for (i in 1:nm) { 
  #print((i/nm)*100)
  # Parameters setting
  nb <- round(simid[i,1]*(Windnb[2]-Windnb[1]),0)+Windnb[1] # number of records to simulate +> will be between 2 and 100
  Text <- round(simid[i,2]*(WindText[2]-WindText[1]),0)+WindText[1] # 'true' extinction date => between 100BP and 9,707
  KnownExt[i] <- Text; #save all the know extinction for later analysis
  maxage <- round(runif(1,min = Text, max = capage-maxnb),0)+maxnb  # oldest record considered
  
  pid<-round(simid[i,3]*(length(prob.dist)-1))+1
  per<-round(simid[i,4]*(length(error.opt1)-1))+1

  input<-cbind(nb,Text,maxage,prob.dist[pid],error.opt1[per])
  colnames(input)<-c("nb","Extinction","MaxTS","distribution","errortype")

  dat<-Tseries(nb, Text, maxage, prob.dist[pid])#generating the time series
  dat2<-TSerror(dat, error.opt1[per],error.opt2[1])#generating error associated with age
  StatTimeseries<-TStat(dat2)
  dat3<-cbind(dat2[,1],dat2[,2])
  colnames(dat3)<-c("age","sd")

  write.table(dat3, 'transitdata.txt', append = FALSE, sep = " ", dec = ".",row.names = FALSE, col.names = FALSE)

  #run GRIWM Biased version
  fredGRIWM1 <- griwm(fossil_data = "transitdata.txt",radiocarbon = "all", calibra = TRUE, cal_curve = "shcal20", upper14C = 55000,
              signor_lipps = "ext", biased = TRUE, alpha = 0.05, resample = 10000, cal_save = FALSE, resample_save = FALSE, griwm_save = FALSE)
  #run GRIWM unBiased version
  fredGRIWM2 <- griwm(fossil_data = "transitdata.txt",radiocarbon = "all", calibra = TRUE, cal_curve = "shcal20", upper14C = 55000,
                    signor_lipps = "ext", biased = FALSE, alpha = 0.05, resample = 10000, cal_save = FALSE, resample_save = FALSE, griwm_save = FALSE)
  #run CRIWM Biased version
  fredCRIWM1 <- criwm(fossil_data = "transitdata.txt", radiocarbon = "all", calibra = TRUE, cal_curve = "intcal20", upper14C = 55000,
      signor_lipps = "ext", biased = TRUE, alpha = 0.05, resample = 10000, cal_save = FALSE, resample_save = FALSE, criwm_save = FALSE)

  #run CRIWM unBiased version
  fredCRIWM2 <- criwm(fossil_data = "transitdata.txt", radiocarbon = "all", calibra = TRUE, cal_curve = "intcal20", upper14C = 55000,
                    signor_lipps = "ext", biased = FALSE, alpha = 0.05, resample = 10000,cal_save = FALSE, resample_save = FALSE, criwm_save = FALSE)

  out<-list(input,StatTimeseries,fredGRIWM1$griwm,fredGRIWM2$griwm,fredCRIWM1$criwm,fredCRIWM2$criwm)
  names(out) <- c("Inputs", "StatsTS","OutputsGRIWMbiased","OutputsGRIWMunbiased","OutputsCRIWMbiased","OutputsCRIWMunbiased")

  list.save(out, paste0("Simul_outputs_", i, ".rds"))

  rm(out,fredGRIWM1,fredGRIWM2,fredCRIWM1,fredCRIWM2,dat,dat2,dat3,StatTimeseries,nb,Text,per,pid,maxage,input)
}
write.csv(KnownExt,"KnownExtinction_uncalibrated.csv")

##=====================================================================
## IMPORTANT TRANSITION STEP !!!!!!!!
##=====================================================================
## The dates in 'KnownExtinction_uncalibrated.csv' need to be calibrated with a calibration curve
# once calibrated the file takes the name 'KnownExtinction.csv'

##=====================================================================
## ANALYSIS OF OUTPUTS BASED on Hilbert-Schmidt Independence Criterion (HSIC)
##=====================================================================
rm(list=ls(all=TRUE))
library(rlist)
library(boot)
library(sensitivity)
library(ggplot2)

data <- read.csv("KnownExtinction.csv",header=T,sep=",",dec=".")
nm<-dim(data)
mat<-matrix(NA, nrow = nm[1], ncol = 15);
colnames(mat)<-c("ID","Text","nb","IntervalAverage","IntervalVariance","ErrorAverage","ErrorVariance","GRIWMbiased_Accuracy",
                 "GRIWMbiased_Precision","GRIWMunbiased_Accuracy","GRIWMunbiased_Precision","CRIWMbiased_Accuracy",
                 "CRIWMbiased_Precision","CRIWMunbiased_Accuracy","CRIWMunbiased_Precision")
for (i in 1:nm[1]) { 
  print((i/nm[1])*100)
  new4<-list.load(paste0("Simul_outputs_", i, ".rds")) #load the results
  mat[i,1:7]<-c(i,data$age.cal.[i],as.numeric(new4$Inputs[1]),new4$StatsTS[1],new4$StatsTS[2],new4$StatsTS[3],new4$StatsTS[4]) #adding the inputs
  mat[i,8:15]<-c(abs(data$age.cal.[i]-new4$OutputsGRIWMbiased[2,2]),abs(new4$OutputsGRIWMbiased[2,3]-new4$OutputsGRIWMbiased[2,1]), 
                 abs(data$age.cal.[i]-new4$OutputsGRIWMunbiased[2,2]),abs(new4$OutputsGRIWMunbiased[2,3]-new4$OutputsGRIWMunbiased[2,1]),                                         
                 abs(data$age.cal.[i]-new4$OutputsCRIWMbiased[2,2]),abs(new4$OutputsCRIWMbiased[2,3]-new4$OutputsCRIWMbiased[2,1]), 
                 abs(data$age.cal.[i]-new4$OutputsCRIWMunbiased[2,2]),abs(new4$OutputsCRIWMunbiased[2,3]-new4$OutputsCRIWMunbiased[2,1])) #adding the inputs
  rm(new4)
}
mat2<-as.data.frame(mat);id<-which(is.na(mat2$IntervalVariance));mat3<-mat2[-c(id),]
write.csv(mat3,file = "GRIWM-CRIWM_sensitivity.csv",row.names = TRUE,col.names = TRUE) 


n <- 100;
X <- as.data.frame(mat3[,3:7]);#predictors

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sensitivity analysis on GRIWM
# biased estinates => on accuracy
HSIC_Gr1 <- sensiHSIC(model = NULL, X,  kernelX = "rbf", kernelY = "rbf",nboot = 1000, conf = 0.95, estimator.type = "U-stat", test.method = "Permutation",B = 10, y = mat3$GRIWMbiased_Accuracy)
yGr1 <- mat3$GRIWMbiased_Accuracy
HSIC_Gr1_out <- tell(HSIC_Gr1,yGr1)
ggplot(HSIC_Gr1_out);

# trends
HSIC_Gr1_LowA<-lowess(X$ErrorAverage,yGr1)
HSIC_Gr1_LowV<-lowess(X$ErrorVariance,yGr1)

plot(X$ErrorAverage,yGr1);lines(lowess(X$ErrorAverage,yGr1),col=2)
plot(X$ErrorVariance,yGr1);lines(lowess(X$ErrorVariance,yGr1),col=2)

outGr1<-data.frame(X$ErrorAverage,X$ErrorVariance,yGr1,HSIC_Gr1_LowA$x,HSIC_Gr1_LowA$y,HSIC_Gr1_LowV$x,HSIC_Gr1_LowV$y)
colnames(outGr1)<-c('ErrorAverage','ErrorVariance','HSICindex','XLowessAvg','YLowessAvg','XLowessVar','YLowessVar')

write.csv(outGr1,file = "HSIC_GRIWMbiased_accuracy.csv",row.names = F,col.names = TRUE) 

#======================================
# biased estinates => on precision
HSIC_Gr2 <- sensiHSIC(model = NULL, X,  kernelX = "rbf", kernelY = "rbf",nboot = 1000, conf = 0.95, estimator.type = "U-stat", test.method = "Permutation",B = 10, y = mat3$GRIWMbiased_Precision)
yGr2 <- mat3$GRIWMbiased_Precision
HSIC_Gr2_out <- tell(HSIC_Gr2,yGr2)
ggplot(HSIC_Gr2_out)

# trends
HSIC_Gr2_LowA<-lowess(X$ErrorAverage,yGr2)
HSIC_Gr2_LowV<-lowess(X$ErrorVariance,yGr2)

plot(X$ErrorAverage,yGr2);lines(lowess(X$ErrorAverage,yGr2),col=2)
plot(X$ErrorVariance,yGr2);lines(lowess(X$ErrorVariance,yGr2),col=2)

outGr2<-data.frame(X$ErrorAverage,X$ErrorVariance,yGr2,HSIC_Gr2_LowA$x,HSIC_Gr2_LowA$y,HSIC_Gr2_LowV$x,HSIC_Gr2_LowV$y)
colnames(outGr2)<-c('ErrorAverage','ErrorVariance','HSICindex','XLowessAvg','YLowessAvg','XLowessVar','YLowessVar')

write.csv(outGr2,file = "HSIC_GRIWMbiased_precision.csv",row.names = F,col.names = TRUE) 

#======================================
# unbiased estinates => on accuracy
HSIC_Gr3 <- sensiHSIC(model = NULL, X,  kernelX = "rbf", kernelY = "rbf",nboot = 1000, conf = 0.95, estimator.type = "U-stat", test.method = "Permutation",B = 10, y = mat3$GRIWMunbiased_Accuracy)
yGr3 <- mat3$GRIWMunbiased_Accuracy
HSIC_Gr3_out <- tell(HSIC_Gr3,yGr3)
ggplot(HSIC_Gr3_out)

# trends
HSIC_Gr3_LowA<-lowess(X$ErrorAverage,yGr3)
HSIC_Gr3_LowV<-lowess(X$ErrorVariance,yGr3)

plot(X$ErrorAverage,yGr3);lines(lowess(X$ErrorAverage,yGr3),col=2)
plot(X$ErrorVariance,yGr3);lines(lowess(X$ErrorVariance,yGr3),col=2)

outGr3<-data.frame(X$ErrorAverage,X$ErrorVariance,yGr3,HSIC_Gr3_LowA$x,HSIC_Gr3_LowA$y,HSIC_Gr3_LowV$x,HSIC_Gr3_LowV$y)
colnames(outGr3)<-c('ErrorAverage','ErrorVariance','HSICindex','XLowessAvg','YLowessAvg','XLowessVar','YLowessVar')

write.csv(outGr3,file = "HSIC_GRIWMunbiased_accuracy.csv",row.names = F,col.names = TRUE) 

#======================================
# unbiased estinates => on precision
HSIC_Gr4 <- sensiHSIC(model = NULL, X,  kernelX = "rbf", kernelY = "rbf",nboot = 1000, conf = 0.95, estimator.type = "U-stat", test.method = "Permutation",B = 10, y = mat3$GRIWMunbiased_Precision)
yGr4 <- mat3$GRIWMunbiased_Precision
HSIC_Gr4_out <- tell(HSIC_Gr4,yGr4)
ggplot(HSIC_Gr4_out)

# trends
HSIC_Gr4_LowA<-lowess(X$ErrorAverage,yGr4)
HSIC_Gr4_LowV<-lowess(X$ErrorVariance,yGr4)

plot(X$ErrorAverage,yGr4);lines(lowess(X$ErrorAverage,yGr4),col=2)
plot(X$ErrorVariance,yGr4);lines(lowess(X$ErrorVariance,yGr4),col=2)

outGr4<-data.frame(X$ErrorAverage,X$ErrorVariance,yGr4,HSIC_Gr4_LowA$x,HSIC_Gr4_LowA$y,HSIC_Gr4_LowV$x,HSIC_Gr4_LowV$y)
colnames(outGr4)<-c('ErrorAverage','ErrorVariance','HSICindex','XLowessAvg','YLowessAvg','XLowessVar','YLowessVar')

write.csv(outGr4,file = "HSIC_GRIWMunbiased_precision.csv",row.names = F,col.names = TRUE) 


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sensitivity analysis on CRIWM
# biased estinates
HSIC_Cr1 <- sensiHSIC(model = NULL, X,  kernelX = "rbf", kernelY = "rbf",nboot = 1000, conf = 0.95, estimator.type = "U-stat", test.method = "Permutation",B = 10, y = mat3$CRIWMbiased_Accuracy)
yCr1 <- mat3$CRIWMbiased_Accuracy
HSIC_Cr1_out <- tell(HSIC_Cr1,yCr1)
ggplot(HSIC_Cr1_out)

# trends
HSIC_Cr1_LowA<-lowess(X$ErrorAverage,yCr1)
HSIC_Cr1_LowV<-lowess(X$ErrorVariance,yCr1)

plot(X$ErrorAverage,yCr1);lines(lowess(X$ErrorAverage,yCr1),col=2)
plot(X$ErrorVariance,yCr1);lines(lowess(X$ErrorVariance,yCr1),col=2)

outCr1<-data.frame(X$ErrorAverage,X$ErrorVariance,yCr1,HSIC_Cr1_LowA$x,HSIC_Cr1_LowA$y,HSIC_Cr1_LowV$x,HSIC_Cr1_LowV$y)
colnames(outCr1)<-c('ErrorAverage','ErrorVariance','HSICindex','XLowessAvg','YLowessAvg','XLowessVar','YLowessVar')

write.csv(outCr1,file = "HSIC_CRIWMbiased_accuracy.csv",row.names = F,col.names = TRUE) 

#======================================
# biased estinates => on precision
HSIC_Cr2 <- sensiHSIC(model = NULL, X,  kernelX = "rbf", kernelY = "rbf",nboot = 1000, conf = 0.95, estimator.type = "U-stat", test.method = "Permutation",B = 10, y = mat3$CRIWMbiased_Precision)
yCr2 <- mat3$CRIWMbiased_Precision
HSIC_Cr2_out <- tell(HSIC_Cr2,yCr2)
ggplot(HSIC_Cr2_out)

# trends
HSIC_Cr2_LowA<-lowess(X$ErrorAverage,yCr2)
HSIC_Cr2_LowV<-lowess(X$ErrorVariance,yCr2)

plot(X$ErrorAverage,yCr2);lines(lowess(X$ErrorAverage,yCr2),col=2)
plot(X$ErrorVariance,yCr2);lines(lowess(X$ErrorVariance,yCr2),col=2)

outCr2<-data.frame(X$ErrorAverage,X$ErrorVariance,yCr2,HSIC_Cr2_LowA$x,HSIC_Cr2_LowA$y,HSIC_Cr2_LowV$x,HSIC_Cr2_LowV$y)
colnames(outCr2)<-c('ErrorAverage','ErrorVariance','HSICindex','XLowessAvg','YLowessAvg','XLowessVar','YLowessVar')

write.csv(outCr2,file = "HSIC_CRIWMbiased_precision.csv",row.names = F,col.names = TRUE) 

#======================================
# unbiased estinates => on accuracy
HSIC_Cr3 <- sensiHSIC(model = NULL, X,  kernelX = "rbf", kernelY = "rbf",nboot = 1000, conf = 0.95, estimator.type = "U-stat", test.method = "Permutation",B = 10, y = mat3$CRIWMunbiased_Accuracy)
yCr3 <- mat3$CRIWMunbiased_Accuracy
HSIC_Cr3_out <- tell(HSIC_Cr3,yCr3)
ggplot(HSIC_Cr3_out)

# trends
HSIC_Cr3_LowA<-lowess(X$ErrorAverage,yCr3)
HSIC_Cr3_LowV<-lowess(X$ErrorVariance,yCr3)

plot(X$ErrorAverage,yCr3);lines(lowess(X$ErrorAverage,yCr3),col=2)
plot(X$ErrorVariance,yCr3);lines(lowess(X$ErrorVariance,yCr3),col=2)

outCr3<-data.frame(X$ErrorAverage,X$ErrorVariance,yCr3,HSIC_Cr3_LowA$x,HSIC_Cr3_LowA$y,HSIC_Cr3_LowV$x,HSIC_Cr3_LowV$y)
colnames(outCr3)<-c('ErrorAverage','ErrorVariance','HSICindex','XLowessAvg','YLowessAvg','XLowessVar','YLowessVar')

write.csv(outCr3,file = "HSIC_CRIWMunbiased_accuracy.csv",row.names = F,col.names = TRUE) 

#======================================
# unbiased estinates => on precision
HSIC_Cr4 <- sensiHSIC(model = NULL, X,  kernelX = "rbf", kernelY = "rbf",nboot = 1000, conf = 0.95, estimator.type = "U-stat", test.method = "Permutation",B = 10, y = mat3$CRIWMunbiased_Precision)
yCr4 <- mat3$CRIWMunbiased_Precision
HSIC_Cr4_out <- tell(HSIC_Cr4,yCr4)
ggplot(HSIC_Cr4_out)

# trends
HSIC_Cr4_LowA<-lowess(X$ErrorAverage,yCr4)
HSIC_Cr4_LowV<-lowess(X$ErrorVariance,yCr4)

plot(X$ErrorAverage,yCr4);lines(lowess(X$ErrorAverage,yCr4),col=2)
plot(X$ErrorVariance,yCr4);lines(lowess(X$ErrorVariance,yCr4),col=2)

outCr4<-data.frame(X$ErrorAverage,X$ErrorVariance,yCr4,HSIC_Cr4_LowA$x,HSIC_Cr4_LowA$y,HSIC_Cr4_LowV$x,HSIC_Cr4_LowV$y)
colnames(outCr4)<-c('ErrorAverage','ErrorVariance','HSICindex','XLowessAvg','YLowessAvg','XLowessVar','YLowessVar')

write.csv(outCr4,file = "HSIC_CRIWMunbiased_precision.csv",row.names = F,col.names = TRUE) 


##================================ saving outputs ============================
# Griwm & Crim biased ACCURACY
GC1_Accuracy_biased=data.frame(HSIC_Gr1$S,HSIC_Cr1$S) #save ouputs
colnames(GC1_Accuracy_biased)<-c("GRIWM_original","GRIWM_bias","GRIWM_stderror","GRIWM_minCI","GRIWM_maxCI",
                                 "CRIWM_original","CRIWM_bias","CRIWM_stderror","CRIWM_minCI","CRIWM_maxCI")
rownames(GC1_Accuracy_biased)<-c("nb","IntervalAverage","IntervalVariance","ErrorAverage","ErrorVariance")
write.csv(GC1_Accuracy_biased,file = "HSIC_GRIWM-CRIWM_biased_accuracy.csv",row.names = TRUE,col.names = TRUE) 

# Griwm & Crim biased PRECISION
GC2_Precision_biased=data.frame(HSIC_Gr2$S,HSIC_Cr2$S) #save ouputs
colnames(GC2_Precision_biased)<-c("GRIWM_original","GRIWM_bias","GRIWM_stderror","GRIWM_minCI","GRIWM_maxCI",
                                  "CRIWM_original","CRIWM_bias","CRIWM_stderror","CRIWM_minCI","CRIWM_maxCI")
rownames(GC2_Precision_biased)<-c("nb","IntervalAverage","IntervalVariance","ErrorAverage","ErrorVariance")
write.csv(GC2_Precision_biased,file = "HSIC_GRIWM-CRIWM_biased_precision.csv",row.names = TRUE,col.names = TRUE) 

# Griwm & Crim UNbiased ACCURACY
GC3_Accuracy_unbiased=data.frame(HSIC_Gr3$S,HSIC_Cr3$S) #save ouputs
colnames(GC3_Accuracy_unbiased)<-c("GRIWM_original","GRIWM_bias","GRIWM_stderror","GRIWM_minCI","GRIWM_maxCI",
                                   "CRIWM_original","CRIWM_bias","CRIWM_stderror","CRIWM_minCI","CRIWM_maxCI")
rownames(GC3_Accuracy_unbiased)<-c("nb","IntervalAverage","IntervalVariance","ErrorAverage","ErrorVariance")
write.csv(GC3_Accuracy_unbiased,file = "HSIC_GRIWM-CRIWM_unbiased_accuracy.csv",row.names = TRUE,col.names = TRUE) 

# Griwm & Crim UNbiased PRECISION
GC4_Precision_unbiased=data.frame(HSIC_Gr4$S,HSIC_Cr4$S) #save ouputs
colnames(GC4_Precision_unbiased)<-c("GRIWM_original","GRIWM_bias","GRIWM_stderror","GRIWM_minCI","GRIWM_maxCI",
                                    "CRIWM_original","CRIWM_bias","CRIWM_stderror","CRIWM_minCI","CRIWM_maxCI")
rownames(GC4_Precision_unbiased)<-c("nb","IntervalAverage","IntervalVariance","ErrorAverage","ErrorVariance")
write.csv(GC4_Precision_unbiased,file = "HSIC_GRIWM-CRIWM_unbiased_precision.csv",row.names = TRUE,col.names = TRUE) 




