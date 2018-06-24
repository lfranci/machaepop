#Luciana Franci 24th Fev 2015
#LM, LMM, GLM and GLM for growth and mortality of Machaerium cuspidatum

# in this script I used the forest uderstory in 2009 in the models
# the understory classes were combined: d+l = d; h+o+s = s

library(lme4)
alldata <- read.csv("merged.data_150225.csv", header = T, sep = ";")
names(alldata)

#Calculating relative growth rate (RGR):  RGR = (ln W2 €“ ln W1)/(t- â€“ t1)
alldata$RGR <- ((log(alldata$dia.5a)-log(alldata$dia.2a))/(2009-1998))
alldata$transect.x <- as.factor(alldata$transect.x)
alldata$region <- as.factor(alldata$region)

seedling <- subset(alldata, sz.grp.2a == "1" | sz.grp.2a == "2")
nonseedling <- subset(alldata, sz.grp.2a == "3" | sz.grp.2a == "4" | sz.grp.2a == "5" |sz.grp.2a == "6" | 
                        sz.grp.2a == "7" | sz.grp.2a == "8" | sz.grp.2a == "9")
notclimb <- subset(nonseedling, gr.fm.2a == "u" | gr.fm.2a == "g" | gr.fm.2a == "a")
climb <- subset(nonseedling, gr.fm.2a == "h" | gr.fm.2a == "c")

#Removing Nas
seedling <- seedling[!is.na(seedling$RGR),]
seedling <- seedling[!is.na(seedling$mean.can.max.2a.5a.),]
seedling <- seedling[!is.na(seedling$region),]
seedling <- seedling[!is.na(seedling$undergr.5a.comb)]

notclimb <- notclimb[!is.na(notclimb$RGR),]
notclimb <- notclimb[!is.na(notclimb$mean.can.max.2a.5a.),]
notclimb <- notclimb[!is.na(notclimb$region),]
notclimb <- notclimb[!is.na(notclimb$undergr.5a.comb)]

climb <- climb[!is.na(climb$RGR),]
climb <- climb[!is.na(climb$mean.can.max.2a.5a.),]
climb <- climb[!is.na(climb$region),]
climb <- climb[!is.na(climb$undergr.5a.comb)]

seedling <- subset(seedling, !inclin.5a.y == "-9")
notclimb <- subset(notclimb, !inclin.5a.y == "-9")
climb <- subset(climb, !inclin.5a.y == "-9")
seedling <- subset(seedling, !undergr.5a.comb == "x")
notclimb <- subset(notclimb, !undergr.5a.comb == "x")
climb <- subset(climb, !undergr.5a.comb == "x")

seedling$undergr.5a.comb <- as.factor(as.character(seedling$undergr.5a.comb))
notclimb$undergr.5a.comb <- as.factor(as.character(notclimb$undergr.5a.comb))
climb$undergr.5a.comb <- as.factor(as.character(climb$undergr.5a.comb))

##################################################################################################################
#Relative growth rate - Linear Mixed Effects Models

#################################################################################################################
##Seedling

#Logit transformation to seedling$RGR
library(psych)
seedling$RGR <- logit(seedling$RGR)

#Selecting the model - LMM for Relative growth rate
#First, include all the relevant fixed effects. Use the AIC to select the best model

#Without random effect
s.lm.1 <- lm(RGR ~ mean.can.max.2a.5a. + inclin.5a.y + palms_PC2 + soil_PC1 + undergr.5a.comb, data = seedling)
summary(s.lm.1)

#LMM with random effect

#Random effect: transect
s.lmm.1 <- lmer(RGR ~ mean.can.max.2a.5a. + inclin.5a.y + palms_PC2 + soil_PC1 + undergr.5a.comb + (1|transect.x), 
                control = lmerControl(optimizer = "bobyqa"), REML = T, data = seedling)
summary(s.lmm.1)
#Random effect: region
s.lmm.2 <- lmer(RGR ~ mean.can.max.2a.5a. + inclin.5a.y + palms_PC2 + soil_PC1 + undergr.5a.comb + (1|region), 
                control = lmerControl(optimizer = "bobyqa"), REML = T, data = seedling)
summary(s.lmm.2)
#Random effect: transect nested within region
s.lmm.3 <- lmer(RGR ~ mean.can.max.2a.5a. + inclin.5a.y + palms_PC2 + soil_PC1 + undergr.5a.comb + (1|region/transect.x), 
                control = lmerControl(optimizer = "bobyqa"), REML = T, data = seedling)
summary(s.lmm.3)

#Checking the normality of the residuals
par(mfrow=c(2,2), oma = c(0.3, 0.2, 2, 0), mgp = c(2.5,1,0))
hist(seedling$RGR - predict(s.lm.1))
hist(seedling$RGR - predict(s.lmm.1))
hist(seedling$RGR - predict(s.lmm.2))
hist(seedling$RGR - predict(s.lmm.3))

AIC(s.lm.1, s.lmm.1, s.lmm.2, s.lmm.3)
#According to AIC the best model is the s.lm.1
s.lm.1 <- lm(RGR ~ mean.can.max.2a.5a. + inclin.5a.y + undergr.5a.comb +  palms_PC2 + soil_PC1, data = seedling)
summary(s.lm.1)
#Using the p-value as exclusing criteria.
#Without the canopy
s.lm.1.1 <- lm(RGR ~ inclin.5a.y + undergr.5a.comb +  palms_PC2 + soil_PC1, data = seedling)
summary(s.lm.1.1)
#Without the understory
s.lm.1.2 <- lm(RGR ~ inclin.5a.y + palms_PC2 + soil_PC1, data = seedling)
summary(s.lm.1.2)
#Without the steepness
s.lm.1.3 <- lm(RGR ~ palms_PC2 + soil_PC1, data = seedling)
summary(s.lm.1.3)
#Without the soil
s.lm.1.4 <- lm(RGR ~ soil_PC1, data = seedling)
summary(s.lm.1.4)

AIC(s.lm.1, s.lm.1.1, s.lm.1.2, s.lm.1.3, s.lm.1.4)

#Validating the model
par(mfrow=c(2,2), oma = c(0.3, 0.2, 2, 0), mgp = c(2.5,1,0))
plot(s.lm.1.3)

#################################################################################################################
#Non-climber
#Without random effect
nc.lm.1 <- lm(RGR ~ mean.can.max.2a.5a. + inclin.5a.y + undergr.5a.comb + palms_PC2 + soil_PC1, data = notclimb)
summary(nc.lm.1)
#Random effect: transect
nc.lmm.1 <- lmer(RGR ~ mean.can.max.2a.5a. + inclin.5a.y + undergr.5a.comb + palms_PC2 + soil_PC1 + (1|transect.x), 
                 control = lmerControl(optimizer = "bobyqa"), REML = T, data = notclimb)
summary(nc.lmm.1)
#Random effect: region
nc.lmm.2 <- lmer(RGR ~ mean.can.max.2a.5a. + inclin.5a.y + undergr.5a.comb + palms_PC2 + soil_PC1 + (1|region), 
                 control = lmerControl(optimizer = "bobyqa"), REML = T, data = notclimb)
summary(nc.lmm.2)
#Random effect: transect nested within region
nc.lmm.3 <- lmer(RGR ~ mean.can.max.2a.5a. + inclin.5a.y + undergr.5a.comb + palms_PC2 + soil_PC1 + (1|region/transect.x), 
                 control = lmerControl(optimizer = "bobyqa"), REML = T, data = notclimb)
summary(nc.lmm.3)

par(mfrow=c(2,2), oma = c(0.3, 0.2, 2, 0), mgp = c(2.5,1,0))
hist(notclimb$RGR - predict(nc.lm.1))
hist(notclimb$RGR - predict(nc.lmm.1))
hist(notclimb$RGR - predict(nc.lmm.2))
hist(notclimb$RGR - predict(nc.lmm.3))

#Comparing the models
AIC(nc.lm.1, nc.lmm.1, nc.lmm.2, nc.lmm.3)

#According to AIC the best model is the nc.lm.1
nc.lm.1 <- lm(RGR ~ mean.can.max.2a.5a. + inclin.5a.y + undergr.5a.comb + palms_PC2 + soil_PC1, data = notclimb)
summary(nc.lm.1)

#Using p-value as exclusion criteria
#Without canopy
nc.lm.1.1 <- lm(RGR ~ inclin.5a.y + undergr.5a.comb + palms_PC2 + soil_PC1, data = notclimb)
summary(nc.lm.1.1)

#Without soil
nc.lm.1.2 <- lm(RGR ~ inclin.5a.y + undergr.5a.comb + palms_PC2, data = notclimb)
summary(nc.lm.1.2)

#Without palms
nc.lm.1.3 <- lm(RGR ~ inclin.5a.y + undergr.5a.comb, data = notclimb)
summary(nc.lm.1.3)

#Without understory
nc.lm.1.4 <- lm(RGR ~ inclin.5a.y, data = notclimb)
summary(nc.lm.1.4)

#Selecting the model
AIC(nc.lm.1, nc.lm.1.1, nc.lm.1a, nc.lm.1.2, nc.lm.1.3)

#Validating the model
par(mfrow=c(2,2), oma = c(0.3, 0.2, 2, 0), mgp = c(2.5,1,0))
plot(nc.lm.1.2)

#################################################################################################################
##Climber
#LM without random effect
c.lm.1 <- lm(RGR ~ mean.can.max.2a.5a. + inclin.5a.y + undergr.5a.comb + palms_PC2 + soil_PC1, data = climb)
summary(c.lm.1)
#Random effect: transect
c.lmm.1 <- lmer(RGR ~ mean.can.max.2a.5a. + inclin.5a.y + undergr.5a.comb + palms_PC2 + soil_PC1 + (1|transect.x), 
                control = lmerControl(optimizer = "bobyqa"), data = climb)
summary(c.lmm.1)
#Random effect: region
c.lmm.2 <- lmer(RGR ~ mean.can.max.2a.5a. + inclin.5a.y + undergr.5a.comb + palms_PC2 + soil_PC1 + (1|region), 
                control = lmerControl(optimizer = "bobyqa"), data = climb)
summary(c.lmm.2)
#Random effect: transect nested within region
c.lmm.3 <- lmer(RGR ~ mean.can.max.2a.5a. + inclin.5a.y + undergr.5a.comb + palms_PC2 + soil_PC1 + (1|region/transect.x), 
                control = lmerControl(optimizer = "bobyqa"), data = climb)
summary(c.lmm.3)

#Comparing the models
AIC(c.lm.1, c.lmm.1, c.lmm.2, c.lmm.3)
par(mfrow=c(2,2), oma = c(0.3, 0.2, 2, 0), mgp = c(2.5,1,0))
hist(climb$RGR - predict(c.lm.1))
hist(climb$RGR - predict(c.lmm.1))
hist(climb$RGR - predict(c.lmm.2))
hist(climb$RGR - predict(c.lmm.3))

#According to AIC the best model is the c.lm.1
c.lm.1 <- lm(RGR ~ mean.can.max.2a.5a. + inclin.5a.y + undergr.5a.comb + palms_PC2 + soil_PC1, data = climb)
summary(c.lm.1)

#Without canopy
c.lm.1.1 <- lm(RGR ~ inclin.5a.y + undergr.5a.comb + palms_PC2 + soil_PC1, data = climb)
summary(c.lm.1.1)

#Without palms
c.lm.1.2 <- lm(RGR ~ inclin.5a.y + undergr.5a.comb + soil_PC1, data = climb)
summary(c.lm.1.2)

#Without steepness
c.lm.1.3 <- lm(RGR ~ undergr.5a.comb + soil_PC1, data = climb)
summary(c.lm.1.3)

#Without soil
c.lm.1.4 <- lm(RGR ~ undergr.5a.comb, data = climb)
summary(c.lm.1.4)

#Selecting the model
AIC(c.lm.1, c.lm.1.1, c.lm.1.2, c.lm.1.3, c.lm.1.4)

#Validating the model
par(mfrow=c(2,2), oma = c(0.3, 0.2, 2, 0), mgp = c(2.5,1,0))
plot(c.lm.1.4)

##################################################################################################################
######Mortality - Generalized Linear Mixed Effects Models

# 0 = dead; 1 = alive

seedling1 <- subset(alldata, sz.grp.2a == "1" | sz.grp.2a == "2")
nonseedling1 <- subset(alldata, sz.grp.2a == "3" | sz.grp.2a == "4" | sz.grp.2a == "5" |sz.grp.2a == "6" | 
                         sz.grp.2a == "7" | sz.grp.2a == "8" | sz.grp.2a == "9")
notclimb1 <- subset(nonseedling1, gr.fm.2a == "u" | gr.fm.2a == "g" | gr.fm.2a == "a")
climb1 <- subset(nonseedling1, gr.fm.2a == "h" | gr.fm.2a == "c")

#Removing Nas
seedling1 <- seedling1[!is.na(seedling1$mean.can.max.2a.5a.),]
notclimb1 <- notclimb1[!is.na(notclimb1$mean.can.max.2a.5a.),]
climb1 <- climb1[!is.na(climb1$mean.can.max.2a.5a.),]

seedling1 <- seedling1[!is.na(seedling1$undergr.5a.comb),]
notclimb1 <- notclimb1[!is.na(notclimb1$undergr.5a.comb),]
climb1 <- climb1[!is.na(climb1$undergr.5a.comb),]

seedling1 <- subset(seedling1, !inclin.5a.y == "-9")
notclimb1 <- subset(notclimb1, !inclin.5a.y == "-9")
climb1 <- subset(climb1, !inclin.5a.y == "-9")

seedling1 <- subset(seedling1, !undergr.5a.comb == "x")
notclimb1 <- subset(notclimb1, !undergr.5a.comb == "x")
climb1 <- subset(climb1, !undergr.5a.comb == "x")

seedling1 <- seedling1[!is.na(seedling1$region),]
notclimb1 <- notclimb1[!is.na(notclimb1$region),]
climb1 <- climb1[!is.na(climb1$region),]

seedling1$deads.2a.5a <- as.factor(seedling1$deads.2a.5a)
table(seedling1$deads.2a.5a)
notclimb1$deads.2a.5a <- as.factor(notclimb1$deads.2a.5a)
table(notclimb1$deads.2a.5a)
climb1$deads.2a.5a <- as.factor(climb1$deads.2a.5a)
table(climb1$deads.2a.5a)

seedling1$undergr.5a.comb <- as.factor(as.character(seedling1$undergr.5a.comb))
table(seedling1$undergr.5a.comb)
notclimb1$undergr.5a.comb <- as.factor(as.character(notclimb1$undergr.5a.comb))
table(notclimb1$undergr.5a.comb)
climb1$undergr.5a.comb <- as.factor(as.character(climb1$undergr.5a.comb))
table(climb1$undergr.5a.comb)

##########################################################################################################
##Seedling

#GLM without random effect
s.glm.1 <- glm(deads.2a.5a ~ mean.can.max.2a.5a. + inclin.5a.y + undergr.5a.comb + palms_PC2 + soil_PC1, 
               family = binomial, data = seedling1)
summary(s.glm.1)

#GLMM with random effect
#Random effect: transect
s.glmm.1 <- glmer(deads.2a.5a ~ mean.can.max.2a.5a. + inclin.5a.y + undergr.5a.comb + palms_PC2 + soil_PC1 + 
                    (1|transect.x), family = binomial, control = glmerControl(optimizer = "bobyqa"), data = seedling1)
summary(s.glmm.1)
#Random effect: region
s.glmm.2 <- glmer(deads.2a.5a ~ mean.can.max.2a.5a. + inclin.5a.y + undergr.5a.comb + palms_PC2 + soil_PC1 + (1|region), 
                  family = binomial, control = glmerControl(optimizer = "bobyqa"), data = seedling1)
summary(s.glmm.2)
#Random effect: transect nested within region
s.glmm.3 <- glmer(deads.2a.5a ~ mean.can.max.2a.5a. + inclin.5a.y + undergr.5a.comb + palms_PC2 + soil_PC1 + 
                    (1|region/transect.x), family = binomial, control = glmerControl(optimizer = "bobyqa"), 
                  data = seedling1)
summary(s.glmm.3)

#How to compare the models?! How do I know which one is the best one?
AIC(s.glm.1, s.glmm.1, s.glmm.2, s.glmm.3)

#According to AIC the best model is the one without the random effect
s.glm.1 <- glm(deads.2a.5a ~ mean.can.max.2a.5a. + inclin.5a.y + undergr.5a.comb + palms_PC2 + soil_PC1, 
               family = binomial, data = seedling1)
summary(s.glm.1)

#Without understory
s.glm.1.1 <- glm(deads.2a.5a ~ mean.can.max.2a.5a. + inclin.5a.y + palms_PC2 + soil_PC1, 
                 family = binomial, data = seedling1)
summary(s.glm.1.1)

#Without soil
s.glm.1.2 <- glm(deads.2a.5a ~ mean.can.max.2a.5a. + inclin.5a.y + palms_PC2, family = binomial, 
                 data = seedling1)
summary(s.glm.1.2)

#Without palms
s.glm.1.3 <- glm(deads.2a.5a ~ mean.can.max.2a.5a. + inclin.5a.y, family = binomial, 
                 data = seedling1)
summary(s.glm.1.3)

#Without canopy
s.glm.1.4 <- glm(deads.2a.5a ~ inclin.5a.y, family = binomial, 
                 data = seedling1)
summary(s.glm.1.4)


AIC(s.glm.1, s.glm.1.1, s.glm.1.2, s.glm.1.3, s.glm.1.4)

##########################################################################################################
#Non-climber

#GLM without random effect
nc.glm.1 <- glm(deads.2a.5a ~ mean.can.max.2a.5a. + inclin.5a.y + undergr.5a.comb + palms_PC2 + soil_PC1, 
                family = binomial, data = notclimb1)
summary(nc.glm.1)

#GLMM with random effect
#Random effect: transect
nc.glmm.1 <- glmer(deads.2a.5a ~ mean.can.max.2a.5a. + inclin.5a.y + undergr.5a.comb + palms_PC2 + soil_PC1 + 
                     (1|transect.x), family = binomial, control = glmerControl(optimizer = "bobyqa"), data = notclimb1)
summary(nc.glmm.1)
#Random effect: region
nc.glmm.2 <- glmer(deads.2a.5a ~ mean.can.max.2a.5a. + inclin.5a.y + undergr.5a.comb + palms_PC2  + soil_PC1+ (1|region), 
                   family = binomial, control = glmerControl(optimizer = "bobyqa"), data = notclimb1)
summary(nc.glmm.2)
#Random effect: transect nested within region
nc.glmm.3 <- glmer(deads.2a.5a ~ mean.can.max.2a.5a. + inclin.5a.y + undergr.5a.comb + palms_PC2 + soil_PC1 + 
                     (1|region/transect.x), family = binomial, control = glmerControl(optimizer = "bobyqa"), 
                   data = notclimb1)
summary(nc.glmm.3)

AIC(nc.glm.1, nc.glmm.1, nc.glmm.2, nc.glmm.3)

#According to AIC the best model is the one with transects as random effects
nc.glmm.1 <- glmer(deads.2a.5a ~ mean.can.max.2a.5a. + inclin.5a.y + undergr.5a.comb + palms_PC2 + soil_PC1 + 
                     (1|transect.x), family = binomial, control = glmerControl(optimizer = "bobyqa"), data = notclimb1)
summary(nc.glmm.1)

#Without understory
nc.glmm.1.1 <- glmer(deads.2a.5a ~ mean.can.max.2a.5a. + inclin.5a.y + palms_PC2 + soil_PC1 + (1|transect.x), 
                     family = binomial, control = glmerControl(optimizer = "bobyqa"), data = notclimb1)
summary(nc.glmm.1.1)

#Without soil
nc.glmm.1.2 <- glmer(deads.2a.5a ~ mean.can.max.2a.5a. + inclin.5a.y + palms_PC2 + (1|transect.x), 
                     family = binomial, control = glmerControl(optimizer = "bobyqa"), data = notclimb1)
summary(nc.glmm.1.2)

#Without palms
nc.glmm.1.3 <- glmer(deads.2a.5a ~ mean.can.max.2a.5a. + inclin.5a.y + (1|transect.x), 
                     family = binomial, control = glmerControl(optimizer = "bobyqa"), data = notclimb1)
summary(nc.glmm.1.3)

#Without canopy
nc.glmm.1.4 <- glmer(deads.2a.5a ~ inclin.5a.y + (1|transect.x), 
                     family = binomial, control = glmerControl(optimizer = "bobyqa"), data = notclimb1)
summary(nc.glmm.1.4)


AIC(nc.glmm.1, nc.glmm.1.1, nc.glmm.1.2, nc.glmm.1.3, nc.glmm.1.4)

##########################################################################################################
##Climber

#GLM without random effect
c.glm.1 <- glm(deads.2a.5a ~ mean.can.max.2a.5a. + inclin.5a.y + undergr.5a.comb + palms_PC2 + soil_PC1, 
               family = binomial, data = climb1)
summary(c.glm.1)

#GLMM with random effect
#Random effect: transect
c.glmm.1 <- glmer(deads.2a.5a ~ mean.can.max.2a.5a. + inclin.5a.y + undergr.5a.comb + palms_PC2 + soil_PC1 + 
                    (1|transect.x), family = binomial, control = glmerControl(optimizer = "bobyqa"), data = climb1)
summary(c.glmm.1)
#Random effect: region
c.glmm.2 <- glmer(deads.2a.5a ~ mean.can.max.2a.5a. + inclin.5a.y + undergr.5a.comb + palms_PC2 + soil_PC1 + (1|region), 
                  family = binomial, control = glmerControl(optimizer = "bobyqa"), data = climb1)
summary(c.glmm.2)
#Random effect: transect nested within region
c.glmm.3 <- glmer(deads.2a.5a ~ mean.can.max.2a.5a. + inclin.5a.y + undergr.5a.comb + palms_PC2 + soil_PC1 + 
                    (1|region/transect.x), family = binomial, control = glmerControl(optimizer = "bobyqa"), 
                  data = climb1)
summary(c.glmm.3)

AIC(c.glm.1, c.glmm.1, c.glmm.2, c.glmm.3)

#According to AIC the best model is the one with transects as random effects
c.glmm.1 <- glmer(deads.2a.5a ~ mean.can.max.2a.5a. + inclin.5a.y + undergr.5a.comb + palms_PC2 + soil_PC1 + 
                    (1|transect.x), family = binomial, control = glmerControl(optimizer = "bobyqa"), data = climb1)
summary(c.glmm.1)

#without palms
c.glmm.1.1 <- glmer(deads.2a.5a ~ mean.can.max.2a.5a. + inclin.5a.y + undergr.5a.comb + soil_PC1 + 
                      (1|transect.x), family = binomial, control = glmerControl(optimizer = "bobyqa"), data = climb1)
summary(c.glmm.1.1)
#Without soil
c.glmm.1.2 <- glmer(deads.2a.5a ~ mean.can.max.2a.5a. + inclin.5a.y + undergr.5a.comb + (1|transect.x), 
                    family = binomial, control = glmerControl(optimizer = "bobyqa"), data = climb1)
summary(c.glmm.1.2)

#Without understory
c.glmm.1.3 <- glmer(deads.2a.5a ~ mean.can.max.2a.5a. + inclin.5a.y + (1|transect.x), 
                    family = binomial, control = glmerControl(optimizer = "bobyqa"), data = climb1)
summary(c.glmm.1.3)

#Without steepness
c.glmm.1.4 <- glmer(deads.2a.5a ~ mean.can.max.2a.5a. + (1|transect.x), 
                    family = binomial, control = glmerControl(optimizer = "bobyqa"), data = climb1)
summary(c.glmm.1.4)

AIC(c.glmm.1, c.glmm.1.1,c.glmm.1.2, c.glmm.1.3, c.glmm.1.4)

