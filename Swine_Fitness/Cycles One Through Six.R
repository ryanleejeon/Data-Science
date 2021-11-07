################################################################################
#####################     WUR Single Trait Analysis     ########################
################################################################################



##===========    Library     =============##

library(lme4qtl)
library(lme4)
library(emmeans)
library(tidyverse)
library(car)
library(carData)
library(dplyr)
library(effects)
library(qpcR)
library(readr)


##========    Loading in Data     ==========##
setwd("/Volumes/RYAN JEON/past/(1) WUR Traits/R Codes (Research)")
load("/Volumes/RYAN JEON/past/(1) WUR Traits/R Codes (Research)/Data_WUR.Rdata")
prod <- read_csv("/Volumes/RYAN JEON/past/(1) WUR Traits/R Codes (Research)/prod_KyuSang.csv")
G <- read.table("~/Desktop/Backup/past/(1) WUR Traits/R Codes (Research)/asreml/G.txt", header=FALSE)
IDs<-read.table("~/Desktop/Backup/past/(1) WUR Traits/R Codes (Research)/asreml/IDs_orig.txt")
ALL <- read_delim("ALL CYCLE_KyuSang.txt", 
                  "\t", escape_double = FALSE, trim_ws = TRUE)



#============   G Matrix    =============#

# Filtering, then changing Matrix Row and Column Names 

a<-filter(data.wur,PassedClean=="TRUE") # this code turns dimensions from 2163 by 3 to 2140 by 3
newID<-a$ID # extracting the IDs from the filtered data.wur file
colnames(G)<-newID # changing column names of my G matrix to that of the new ID
row.names(G)<-newID # changing row names of my G matrix to that of the new ID
GG<-data.matrix(G, rownames.force = NA) # important to change to a matrix





#===========    Class-Changing Variables and adding WUR Data to the main production data set    ===========#

data.production <- left_join(ALL, data.wur) 

data.production$SowID<-as.factor(prod$SowID) # sow ID is a factor 
data.production$Nur2ADG <- as.numeric(data.production$Nur2ADG) # Nur2ADG is numeric
data.production$Pen<-as.factor(prod$NurPenBatch) # Pen is a factor
data.production$FinisherPen<-as.factor(prod$FinPenBatch)
data.production$WURFac <- as.factor(data.production$WUR)
data.production$BatchNumber<- as.factor(data.production$BatchNum)
data.production$Wt2AgeNum<-as.numeric(data.production$Wt2Age)
data.production$Trt <- prod$nTrtsPer180

# Data.Production to include Binary Response (Death)
data.production$family<-prod$Died #where 1 means they died, 0 means alive

# set 5's to missing values
data.production$WUR[data.production$WUR > 4.9 & data.production$WUR < 5.1] <- NA

#Response Variables:
data.production$Trt <- prod$nTrtsPer180
data.production$AFI <-prod$ADFI
data.production$Fin <-prod$FinADG
data.production$DBF <-prod$DBF
data.production$RMSE <-prod$FI_RMSE
data.production$QRP <-prod$FI_QRP
data.production$DLD <-prod$DLD
data.production$RFI <-prod$RFI

# data production has 2273 pigs... so... 

DATA<-subset(data.production,ID %in% a$ID) 
# (above) DATA is a subset of data production where ID's that do not match ID of the filtered data set are removed
is.na(DATA) <- DATA=="."



####################  Models of Each Response Variable ########

## Model for Nur2ADG
model<-relmatLmer(Nur2ADG~BatchNumber+Wt2AgeNum+WURFac + SowID + (1|Pen)+ (1|ID), DATA, relmat=list(ID=GG))
summary(model)
anova(model)
emmeans(model, "WURFac", contr="pairwise")
plot(allEffects(model))
plot(model)


## Model for Treatment
modelTrtPEN<-relmatLmer(Trt~BatchNumber+Wt2AgeNum+WURFac+ SowID + (1|Pen)+ (1|ID), DATA, relmat=list(ID=GG))
summary(modelTrtPEN)
anova(modelTrtPEN)
emmeans(modelTrtPEN, "WURFac", contr="pairwise")
plot(allEffects(modelTrtPEN))

modelTrtFin<-relmatLmer(Trt~BatchNumber+Wt2AgeNum+WURFac+ SowID + (1|FinisherPen)+ (1|ID), DATA, relmat=list(ID=GG))
summary(modelTrtFin)
anova(modelTrtFin)
emmeans(modelTrtFin, "WURFac", contr="pairwise")
plot(allEffects(modelTrtFin))

######  Average Feed Intake

modelAFI<-relmatLmer(AFI~BatchNumber+Wt2AgeNum+WURFac+ SowID + (1|FinisherPen)+ (1|ID), DATA, relmat=list(ID=GG))
summary(modelAFI)
anova(modelAFI)
emmeans(modelAFI, "WURFac", contr="pairwise")
plot(allEffects(modelAFI))

###### Finishing ADG
modelFin<-relmatLmer(Fin~BatchNumber+Wt2AgeNum+WURFac+ SowID + (1|FinisherPen)+ (1|ID), DATA, relmat=list(ID=GG))
summary(modelFin)
anova(modelFin)
emmeans(modelFin, "WURFac", contr="pairwise")
plot(allEffects(modelFin))


###### Back Fat
modelBF<-relmatLmer(DBF~BatchNumber+Wt2AgeNum+WURFac+ SowID + (1|FinisherPen)+ (1|ID), DATA, relmat=list(ID=GG))
summary(modelBF)
anova(modelBF)
emmeans(modelBF, "WURFac", contr="pairwise")
plot(allEffects(modelBF))


##### RMSE
modelRMSE<-relmatLmer(RMSE~BatchNumber+Wt2AgeNum+WURFac+ SowID +(1|FinisherPen)+ (1|ID), DATA, relmat=list(ID=GG))
summary(modelRMSE)
anova(modelRMSE)
emmeans(modelRMSE, "WURFac", contr="pairwise")
plot(allEffects(modelRMSE))


##### QRP
modelQRP<-relmatLmer(QRP~BatchNumber+Wt2AgeNum+WURFac+ SowID + (1|FinisherPen)+ (1|ID), DATA, relmat=list(ID=GG))
summary(modelQRP)
anova(modelQRP)
emmeans(modelQRP, "WURFac", contr="pairwise")
plot(allEffects(modelQRP))

#### RFI
modelRFI<-relmatLmer(RFI~BatchNumber+Wt2AgeNum+WURFac+ SowID + (1|FinisherPen)+ (1|ID), DATA, relmat=list(ID=GG))
summary(modelRFI)
anova(modelRFI)
emmeans(modelRFI, "WURFac", contr="pairwise")
plot(allEffects(modelRFI))

#### DLD
modelDLD<-relmatLmer(DLD~BatchNumber+Wt2AgeNum+WURFac+ SowID + (1|FinisherPen)+ (1|ID), DATA, relmat=list(ID=GG))
summary(modelDLD)
anova(modelDLD)
emmeans(modelDLD, "WURFac", contr="pairwise")
plot(allEffects(modelDLD))





############################ Binary Response
## Model for Nur2ADG
model_binary<-relmatGlmer(family~BatchNumber+Wt2AgeNum+WURFac +(1|SowID)+ (1|Pen)+ (1|ID), DATA, relmat=list(ID=GG), family=binomial)
summary(model_binary)
anova(model_binary)
emmeans(model_binary, "WURFac", contr="pairwise")
plot(allEffects(model_binary))



