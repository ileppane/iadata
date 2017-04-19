library(lme4)
library(lmerTest)
library(lattice)
library(effects)
library(reshape2)
library(plyr)
library(scales)
library(gridExtra)
library(ggplot2)
library(arm)
rm(list = ls())

# READ DATA FROM CSV-FILES INTO DATA FRAMES
mydata <- read.table("iadata.csv", header=TRUE, sep=",")
bdata <- read.table("baselinedata.csv", header=TRUE, sep=",")
bindata <- read.table("bindata.csv", header=TRUE, sep=",")
# ADD LEVELS
mydata$person <- as.factor(mydata$person)
mydata$treatm <- relevel(mydata$treatm, "pas") # MAKE PASSIVE TREATMENT THE REFERENCE
mydata$iscr <- log(1+mydata$iscr) # LOG-TRANSFORM ISCR

# HERE WE FORM THE BINDATA FOR THE EMG
# FURR BROW: IF CORR IS THE ONLY ACTIVE IN THE BIN
mydata$furrbrow <- ifelse(mydata$bcorr < bindata$bincs1 & mydata$bzygo > bindata$binzm1 & mydata$borbi > bindata$binoo1,1,0) +
 			 ifelse(mydata$bcorr < bindata$bincs2 & mydata$bzygo > bindata$binzm2 & mydata$borbi > bindata$binoo2,1,0) +
			 ifelse(mydata$bcorr < bindata$bincs3 & mydata$bzygo > bindata$binzm3 & mydata$borbi > bindata$binoo3,1,0) +
			 ifelse(mydata$bcorr < bindata$bincs4 & mydata$bzygo > bindata$binzm4 & mydata$borbi > bindata$binoo4,1,0) +
			 ifelse(mydata$bcorr < bindata$bincs5 & mydata$bzygo > bindata$binzm5 & mydata$borbi > bindata$binoo5,1,0) +
			 ifelse(mydata$bcorr < bindata$bincs6 & mydata$bzygo > bindata$binzm6 & mydata$borbi > bindata$binoo6,1,0)
# DUCHENNE: IF ZYGO & ORBI ARE ACTIVE AND CORR IS NOT
mydata$duchenne <- ifelse(mydata$bzygo < bindata$binzm1 & mydata$borbi < bindata$binoo1 & mydata$bcorr > bindata$bincs1,1,0) +
		 	 ifelse(mydata$bzygo < bindata$binzm2 & mydata$borbi < bindata$binoo2 & mydata$bcorr > bindata$bincs2,1,0) +
		 	 ifelse(mydata$bzygo < bindata$binzm3 & mydata$borbi < bindata$binoo3 & mydata$bcorr > bindata$bincs3,1,0) +
			 ifelse(mydata$bzygo < bindata$binzm4 & mydata$borbi < bindata$binoo4 & mydata$bcorr > bindata$bincs4,1,0) +
			 ifelse(mydata$bzygo < bindata$binzm5 & mydata$borbi < bindata$binoo5 & mydata$bcorr > bindata$bincs5,1,0) +
			 ifelse(mydata$bzygo < bindata$binzm6 & mydata$borbi < bindata$binoo6 & mydata$bcorr > bindata$bincs6,1,0)
# NONDUCHENNE: IF ONLY ZYGO IS ACTIVE
mydata$nonduchenne <- ifelse(mydata$bzygo < bindata$binzm1 & mydata$borbi > bindata$binoo1 & mydata$bcorr > bindata$bincs1,1,0) +
			    ifelse(mydata$bzygo < bindata$binzm2 & mydata$borbi > bindata$binoo2 & mydata$bcorr > bindata$bincs2,1,0) +
			    ifelse(mydata$bzygo < bindata$binzm3 & mydata$borbi > bindata$binoo3 & mydata$bcorr > bindata$bincs3,1,0) +
			    ifelse(mydata$bzygo < bindata$binzm4 & mydata$borbi > bindata$binoo4 & mydata$bcorr > bindata$bincs4,1,0) +
			    ifelse(mydata$bzygo < bindata$binzm5 & mydata$borbi > bindata$binoo5 & mydata$bcorr > bindata$bincs5,1,0) +
			    ifelse(mydata$bzygo < bindata$binzm6 & mydata$borbi > bindata$binoo6 & mydata$bcorr > bindata$bincs6,1,0)
mydata$time <- rep(1:78,33) # NOT USED

# CENTERING
mydata$cenduchenne <- mydata$duchenne - rep(sapply(split(mydata$duchenne,mydata$person),mean),each=78)
mydata$cennonduchenne <- mydata$nonduchenne - rep(sapply(split(mydata$nonduchenne,mydata$person),mean),each=78)
mydata$cenfurrbrow <- mydata$furrbrow - rep(sapply(split(mydata$furrbrow,mydata$person),mean),each=78)
mydata$cenempathy <- mydata$empathy-mean(mydata$empathy)
mydata$cenbrcorr <- mydata$brcorr - rep(sapply(split(mydata$brcorr,mydata$person),mean),each=78)
mydata$cenbrzygo <- mydata$brzygo - rep(sapply(split(mydata$brzygo,mydata$person),mean),each=78)
mydata$cenbrorbi <- mydata$brorbi - rep(sapply(split(mydata$brorbi,mydata$person),mean),each=78)


##
# MAIN RESULTS HERE!
#
# EMGs ARE COUNT DATA, ISCR IS LOGARITHMIZED AND CONTINUOUS
# DUMMY CODING: INTERCEPT REFERS TO THE CELL MEAN OF THE REF GROUP ("pas")
summary(m2du<-lmer(duchenne~treatm+(1|person), mydata, REML=F))
summary(m2nd<-lmer(nonduchenne~treatm+(1|person), mydata, REML=F))
summary(m2fb<-lmer(furrbrow~treatm+(1|person), mydata, REML=F))
summary(m2ar<-lmer(iscr~treatm+(1|person), mydata, REML=F))
summary(m2ar2 <- lmer(iscr~duchenne+nonduchenne+furrbrow+(1|person), mydata, REML=F))

##
# TABLE 1
# in interaction effects the variables need to be subject mean centered
mydata$cenduchenne<-mydata$cenduchenne
mydata$cennonduchenne<-mydata$cennonduchenne
mydata$cenfurrbrow<-mydata$cenfurrbrow
summary(m2ar22 <- lmer(iscr~cenduchenne+cennonduchenne+cenfurrbrow+(cenduchenne+cennonduchenne+cenfurrbrow):treatm+(1|person), mydata, REML=F))
plot(allEffects(m2ar22))

##
# FOCUSED CONTRASTS
#
# NON-DUCHENNE: BETWEEN INQ AND ADV
summary(lmer(nonduchenne~treatm+ (1|person), subset(mydata,treatm!="pas"), REML=F))
# AROUSAL: BETWEEN INQ AND ADV
summary(lmer(iscr~treatm+ (1|person), subset(mydata,treatm!="pas"), REML=F))
# FURRBROWS AS FUNCTION OF DUCHENNES AND NON-DUCHENNES
summary(lmer(furrbrow~duchenne+nonduchenne+(1|person), mydata, REML=F))
# INQUIRY TREATMENT: IS DUCHENNE GREATER THAN NON-DUCHENNE?
t.test(mydata$duchenne[mydata$treatm=="inq"], mydata$nonduchenne[mydata$treatm=="inq"], paired=T, alternative="greater")



#############################
# FIGURE 1: main treatment effects
# SCALE REF LEVEL (PASSIVE) TO ZERO!
plf <- data.frame(fit=c(fixef(m2du),fixef(m2nd),fixef(m2fb),fixef(m2ar)), 
                  se= c(se.fixef(m2du),se.fixef(m2nd),se.fixef(m2fb),se.fixef(m2ar)),
                  tr=c("du","du","du","nd","nd","nd","fb","fb","fb","ar","ar","ar"),
			row=c(1,"Advocacy","Inquiry",1,"Advocacy","Inquiry",1,"Advocacy","Inquiry",1,"Advocacy","Inquiry"))
plf <- subset(plf,row!=1)
plfplot1 <- ggplot(subset(plf,tr=="du"), aes(y=fit,x=row)) + 
            geom_bar(stat="identity", fill="#ababab") +
            geom_errorbar(aes(ymin=fit-se, ymax=fit+se),width=0,size=1) +
            annotate("text", x = 1, y = subset(plf,tr=="du")$fit[1]+subset(plf,tr=="du")$se[1], label = "", size=10) +
            annotate("text", x = 2, y = subset(plf,tr=="du")$fit[2]+subset(plf,tr=="du")$se[2]+.01, label = "***", size=10) +
            theme_classic()+ylim(c(-.15,.65)) +
            ylab("Mean EMG bin count") + ggtitle("Duchenne") + xlab("")
plfplot2 <- ggplot(subset(plf,tr=="nd"), aes(y=fit,x=row)) + 
            geom_bar(stat="identity", fill="#ababab") +
            geom_errorbar(aes(ymin=fit-se, ymax=fit+se),width=0,size=1) +
            annotate("text", x = 1, y = subset(plf,tr=="nd")$fit[1]+subset(plf,tr=="nd")$se[1]+.005, label = "*", size=10) +
            annotate("text", x = 2, y = subset(plf,tr=="nd")$fit[2]+subset(plf,tr=="nd")$se[2]+.005, label = "*", size=10) +
            theme_classic()+ylim(c(0,.15)) +
            ylab("Mean EMG bin count") + ggtitle("non-Duchenne") + xlab("")
plfplot3 <- ggplot(subset(plf,tr=="fb"), aes(y=fit,x=row)) + 
            geom_bar(stat="identity", fill="#ababab", position=position_dodge()) +
            geom_errorbar(aes(ymin=fit-se, ymax=fit+se),width=0,size=1) +
            annotate("text", x = 1, y = subset(plf,tr=="fb")$fit[1]+subset(plf,tr=="fb")$se[1]+.015, label = "*", size=10) +
            annotate("text", x = 2, y = -1, label = "***", size=10) +
            theme_classic()+ylim(c(-1,.4)) +
            ylab("Mean EMG bin count") + ggtitle("Furrowed brows") + xlab("")
plfplot4 <- ggplot(subset(plf,tr=="ar"), aes(y=fit,x=row)) + 
            geom_bar(stat="identity", fill="#ababab") +
            geom_errorbar(aes(ymin=fit-se, ymax=fit+se),width=0,size=1) +
            annotate("text", x = 1, y = subset(plf,tr=="ar")$fit[1]+subset(plf,tr=="ar")$se[1]+.01, label = "**", size=10) +
            annotate("text", x = 2, y = subset(plf,tr=="ar")$fit[2]+subset(plf,tr=="ar")$se[2]+.01, label = "***", size=10) +
            theme_classic()+ylim(c(0,0.55)) + 
            ylab("Mean ISCR (log uSs)") + ggtitle("Arousal") + xlab("")
grid.arrange(plfplot1,plfplot2,plfplot3,plfplot4, ncol=2, nrow=2)

##
# EMPATHY LMMs (TABLE 2)
summary(lmer(duchenne~cenempathy*treatm+(1|person), mydata, REML=F))
summary(lmer(nonduchenne~cenempathy*treatm + (1|person), mydata, REML=F))
summary(lmer(furrbrow~cenempathy*treatm + (1|person), mydata, REML=F))
summary(lmer(iscr~cenempathy*treatm+(1|person), mydata, REML=F))
