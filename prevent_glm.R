library(fitdistrplus)
library(logspline)
library(glmmTMB)
library(lme4)
library(ggplot2)
library(AICcmodavg)
library(gamlss)

setwd("/Users/shine/Documents/MSc/Neuro\ Research/PREVENT_Study/dump/")

dataset = read.csv("wb_ctrb_aging.csv")
dataset = read.csv("processing_speed_ctrb.csv")

dataset = read.csv("verbal_memory_ctrb.csv",header=T, stringsAsFactors=T)

dataset = read.csv("vis_memory_ctrb.csv")

#======== Determine outcome variable distribution 
descdist(dataset$wbModalCtrb, discrete = FALSE, boot=1000)
descdist(dataset$ModalCtrb, discrete = FALSE, boot=1000)
descdist(dataset$AvgCtrb, discrete = FALSE, boot=1000)

descdist(dataset$CogScore, discrete = TRUE, boot=1000)

#======== fit outcome variable & test 
fit.normal <- fitdist(dataset$wbModalCtrb, "normal")
fit.lognormal <- fitdist(dataset$ModalCtrb, "lnorm")
fit.beta <- fitdist(dataset$ModalCtrb, "beta")
fit.gamma <- fitdist(dataset$wbAvgCtrb, "gamma")

plot(fit.normal)

shapiro.test(dataset$wbModalCtrb)
shapiro.test(dataset$ModalCtrb)

plot(density(resid(lme.wbModalCtrb)))
qqnorm(resid(lme.wbModalCtrb))
qqline(resid(lme.wbModalCtrb))

#======== Construct LME model 
lme.wbModalCtrb <- lmer(wbModalCtrb ~ age + group + gender + education + (1|time), data=dataset)

lme.ModalCtrbREML <- lmer(log(ModalCtrb) ~ Age + Group + BioSex + EducationYrs + (1|Time), data=dataset)
lme_model <- lmer(log(ModalCtrb) ~ Age + Group + BioSex + EducationYrs + (1|Time), data=dataset, REML = FALSE)
null_wo_group = lmer(log(ModalCtrb) ~ Age + BioSex + EducationYrs + (1|Time), data=dataset, REML = FALSE)
null_wo_age = lmer(log(ModalCtrb) ~ Group + BioSex + EducationYrs + (1|Time), data=dataset, REML = FALSE)

anova(null_wo_group, lme_model)
anova(null_wo_age, lme_model)


summary(lme.ModalCtrb)

#======== Construct GLMM model w/ Binomial distribution
glmm.1 <- glmer(CogScore/15 ~ FPN + AUD + FPN*AUD + Group + BioSex + EducationYrs + Iq + (1|Time), data=dataset, family = "binomial")


#======== Construct GLMM model w/ Lognormal distribution
glmm.1 <- glmmTMB(CogScore ~ FPN + CON + FPN*CON + Group + BioSex + EducationYrs + Iq + (1|Time), data=dataset, family = lognormal())

#======== Construct GLMM model w/ Weibull distribution
glmm.weibull <- glmer(ModalCtrb ~ group + age + gender + education + (1|time), data=dataset, family = weibull())

#======== Construct GLMM model w/ Beta distribution
glmm.beta <- glmmTMB(ModalCtrb ~ Group + Age + BioSex + EducationYrs + (1|Time), data=dataset, family = beta_family())

#======== Construct GLMM model w/ Gamma distribution
gamma.4 <- glmmTMB(CogScore ~ FPN + Group + BioSex + EducationYrs + (1|Time), data=dataset, family = Gamma(link="log"))

#======== Skill mack test 
library('Skillings.Mack')

sink("ski-mack-result.csv")
Ski.Mack(as.matrix(dataset$wbAvgCtrb),groups = as.matrix(dataset$Group), blocks = as.matrix(dataset$Time))
sink()


