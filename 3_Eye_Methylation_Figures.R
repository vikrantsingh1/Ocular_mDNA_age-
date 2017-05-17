setwd("/Users/vikrants/Desktop/DNAmage/hannum_result/")
library(WGCNA)
library(sqldf)
library(impute)
library(lattice)
dataout <- read.csv("betavalues.output.csv")
#######################################HORVATH_Method_PLOTS##########################################
signif(dataout$DNAmAge,2)
datSample<-read.csv("sampleidannotat.csv")
DNAmAge1=dataout$DNAmAge
medianAbsDev=function(x,y) median(abs(x-y),na.rm=TRUE)
medianAbsDev1=signif(medianAbsDev(DNAmAge1, datSample$Age),2)
par(mfrow=c(1,3))
#verboseScatterplot(datSample$Age,DNAmAge,xlab="Chronological Age",ylab="DNAm Age",ylim=c(10,100), pch =c(1,2,3,4)[as.numeric(datSample$SampleLabel)], cex = 1.3, cex.axis = 1, cex.lab = 1, cex.main = 1, col= rainbow(10)[round(datSample$Age/max(datSample$Age)*10)], main=paste("All, err=", medianAbsDev1) ) ;abline(0,1) 

#verboseScatterplot(datSample$Age,DNAmAge,xlab="Chronological Age",ylab="DNAm Age",xlim=c(10,80),ylim=c(10,80), pch =c(1,2,3,4)[as.numeric(datSample$SampleLabel)], cex = 1.8, cex.axis = 1, cex.lab = 1,cex.main=1,col.main="white",col= c("black"));abline(0,1) 

verboseScatterplot(datSample$Age,DNAmAge1,xlab="Chronological Age",ylab="DNAm Age Horvath",xlim=c(10,80), ylim=c(10,80), pch =c(1,2,3,4)[as.numeric(datSample$SampleLabel)], cex = 2.3, cex.axis = 1.3, cex.lab = 1.3, cex.main = 1.3, col= c("black"), main=paste("Horvath") ) ;abline(0,1) 
legend("topleft", legend= c("Blood","Neurosensory Retina","Optic Nerve","RPE/Choriod"),pch=c(1,2,3,4),bty = "n",text.width = NULL,y.intersp =1,x.intersp = 1 ,cex = 1.4)

AgeaccelerationDiff1 = DNAmAge1 - datSample$Age
samlabel1 <- factor(c(datSample$SampleLabel) )
samlabel1 <- ordered(datSample$SampleLabel, levels=c("Blood","Optic Nerve","RPE/Choroid","Neurosensory Retina") )

verboseBarplot(AgeaccelerationDiff1,samlabel1, xlab= "Tissues", ylab = "Difference in epigenetic and chronological age(years)",ylim=c(-60,20), cex = 1, cex.axis = 1, cex.lab = 1, cex.main =1,col.main="white", numberStandardErrors = 1)
#verboseBarplot(AgeaccelerationDiff, datSample$CauseofDeath,main = "Cause of Death", xlab ="Cause of Death", ylab = "Age Acceleration Diff",ylim=c(-70,70), cex =0.7, cex.axis =1,cex.lab =1, cex.main =1, numberStandardErrors = 1)
verboseScatterplot(datSample$PreservationTimeInterval,AgeaccelerationDiff1,xlab="Post Mortem Interval",xlim=c(2,12),ylim=c(-60,20), pch =c(1,2,3,4),cex = 1.8, cex.axis =1,cex.lab =1,cex.main =1,main ="C",col.main="white")

#########################################################################################
##########################HANNUM_Method_PLOTS###########################################
signif(dataout$BioAge4HA,2)
datSample<-read.csv("sampleidannotat.csv")
DNAmAge=dataout$BioAge4HA
medianAbsDev=function(x,y) median(abs(x-y),na.rm=TRUE)
medianAbsDev1=signif(medianAbsDev(DNAmAge, datSample$Age),2)
par(mfrow=c(1,3))
#verboseScatterplot(datSample$Age,DNAmAge,xlab="Chronological Age",ylab="DNAm Age",ylim=c(10,100), pch =c(1,2,3,4)[as.numeric(datSample$SampleLabel)], cex = 1.3, cex.axis = 1, cex.lab = 1, cex.main = 1, col= rainbow(10)[round(datSample$Age/max(datSample$Age)*10)], main=paste("All, err=", medianAbsDev1) ) ;abline(0,1) 

#verboseScatterplot(datSample$Age,DNAmAge,xlab="Chronological Age",ylab="DNAm Age",xlim=c(10,80),ylim=c(10,80), pch =c(1,2,3,4)[as.numeric(datSample$SampleLabel)], cex = 1.8, cex.axis = 1, cex.lab = 1,cex.main=1,col.main="white",col= c("black"));abline(0,1) 

verboseScatterplot(datSample$Age,DNAmAge,xlab="Chronological Age",ylab="DNAm Age Hannum",xlim=c(10,80), ylim=c(10,140), pch =c(1,2,3,4)[as.numeric(datSample$SampleLabel)], cex = 2.3, cex.axis = 1.3, cex.lab = 1.3, cex.main = 1.3, col= c("black"), main=paste("Hannum") ) ;abline(0,1) 
legend("topleft", legend= c("Blood","Neurosensory Retina","Optic Nerve","RPE/Choriod"),pch=c(1,2,3,4),bty = "n",text.width = NULL,y.intersp =1,x.intersp = 1 ,cex = 1.4)

AgeaccelerationDiff = DNAmAge - datSample$Age
samlabel <- factor(c(datSample$SampleLabel) )
samlabel <- ordered(datSample$SampleLabel, levels=c("Blood","Optic Nerve","RPE/Choroid","Neurosensory Retina") )

verboseBarplot(AgeaccelerationDiff,samlabel, xlab= "Tissues", ylab = "Difference in epigenetic and chronological age(years)",ylim=c(-70,70), cex = 1, cex.axis = 1, cex.lab = 1, cex.main =1,col.main="white", numberStandardErrors = 1)
#verboseBarplot(AgeaccelerationDiff, datSample$CauseofDeath,main = "Cause of Death", xlab ="Cause of Death", ylab = "Age Acceleration Diff",ylim=c(-70,70), cex =0.7, cex.axis =1,cex.lab =1, cex.main =1, numberStandardErrors = 1)
verboseScatterplot(datSample$PreservationTimeInterval,AgeaccelerationDiff,xlab="Post Mortem Interval",xlim=c(2,12),ylim=c(-70,70), pch =c(1,2,3,4),cex = 1.8, cex.axis =1,cex.lab =1,cex.main =1,main ="C",col.main="white")
########################################################################################################