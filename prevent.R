library(lme4)
library(ggplot2)
library(MuMIn)
setwd("/Users/shine/Documents/MSc/Neuro\ Research/PREVENT_Study/dump/")
library(AICcmodavg)

dataset = read.csv("ctrb_aging.csv")

# ---- Whole brain average controllability with age
# (1) Selecting a model that best fits the data 
lm_1 <- lmer(data = dataset, wbAvgCtrb ~ age*group + gender + education + (1|time), REML = FALSE)
lm_2 <- lmer(data = dataset, wbAvgCtrb ~ age + group + gender + education + (1|time), REML = FALSE)
#lm_3 <- lmer(data = dataset, wbAvgCtrb ~ age + gender + education + (1|time), REML = FALSE)
lm_4 <- lmer(data = dataset, wbAvgCtrb ~ age + gender + education + (1|time) + (1|group), REML = FALSE)

# using AIC to select best model for the data. The lower the AIC score, the better the model. 
models <- list(lm_1, lm_2, lm_4)
model.names <- c('model1', 'model2', 'model4')
aictab(cand.set = models, modnames = model.names)

# Result: lm_3 is the best (group parameter left out)
# However, if we require group to be a co-variate/parameter, then the best model is lm_2. 

# (2) Q: Is whole-brain average ctrb dependent on age and C/T group? 
avg_ctrb.model = lmer(wbAvgCtrb ~ age + group + gender + education + (1|time), data=dataset)
avg_ctrb.model_2 = lmer(wbAvgCtrb ~ age + group + gender + education + (1|time), data=dataset, REML = FALSE)
avg_ctrb.null = lmer(wbAvgCtrb ~ group + gender + education + (1|time), data=dataset, REML = FALSE)

# test model against null hypothesis - do we use ANOVA? Or something else? I think something else... 
anova(avg_ctrb.null, avg_ctrb.model_2)

# Results: 
# 1. No sig. difference between control and TIA group 
# 2. Age significantly correlates with whole-brain average ctrb. We find that as an individual age, their 
# whole-brain average controllability increases with age (p = 0.003605)

# (3) plots
lme_mod = lmer(data = dataset, wbAvgCtrb ~ age + group + gender + education + (1|time))
plot(wbAvgCtrb ~ age, data = dataset)
points(dataset$age[order(dataset$age)], fitted(lme_mod)[order(dataset$age)], pch = 19)


## Spaghetti Plot of the fit longitudinal data
lme_mod = lmer(data = dataset, wbAvgCtrb ~ age + group + gender + education + (1|time))

pfit <-
  ggplot(data = dataset, aes(
    x = age,
    y = fitted(lme_mod),
    group = time)
  )
pfit + geom_point() + geom_line() 

# ====================================================================================
# ---- Whole brain MODAL controllability with age
dataset = read.csv("ctrb_aging.csv")

# (1) Selecting a model that best fits the data 
lm_1 <- lmer(data = dataset, wbModalCtrb ~ age*group + gender + education + (1|time), REML = FALSE)
lm_2 <- lmer(data = dataset, wbModalCtrb ~ age + group + gender + education + (1|time), REML = FALSE)
lm_3 <- lmer(data = dataset, wbModalCtrb ~ age + gender + education + (1|time), REML = FALSE)
lm_4 <- lmer(data = dataset, wbModalCtrb ~ age + gender + education + (1|time) + (1|group), REML = FALSE)

# using AIC to select best model for the data. The lower the AIC score, the better the model. 
models <- list(lm_1, lm_2, lm_3, lm_4)
model.names <- c('model1', 'model2', 'model3', 'model4')
aictab(cand.set = models, modnames = model.names)

# Result: lm_2 is the best!! This is interesting since it means having the group included is a better model than it removed. 
# This is very different from the avg. ctrb. model where group is better to be removed as a parameter!!!

# (2) Q: Is whole-brain modal ctrb dependent on age and C/T group? 
mod_ctrb.model = lmer(wbModalCtrb ~ age + group + gender + education + (1|time), data=dataset)
mod_ctrb.model_2 = lmer(wbModalCtrb ~ age + group + gender + education + (1|time), data=dataset, REML = FALSE)
mod_ctrb.null_wo_group = lmer(wbModalCtrb ~ age + gender + education + (1|time), data=dataset, REML = FALSE)
mod_ctrb.null_wo_age = lmer(wbModalCtrb ~ group + gender + education + (1|time), data=dataset, REML = FALSE)

# test model against null hypothesis - do we use ANOVA? Or something else? I think something else... 
anova(mod_ctrb.null_wo_group, mod_ctrb.model_2)
anova(mod_ctrb.null_wo_age, mod_ctrb.model_2)

# Results: 
# 1. Group (C/TIA) makes a significant "impact" on the model (p = 0.008861)
# 2. Age does NOT make a significant "impact" on the model !
# 3. Don't understand how to read model output, but it seems like it's saying "TIA" group has lower w-b modal ctrb. 

# (3) plots -- doesn't work, need to separate dataset by group
lme_mod = lmer(data = dataset, wbModalCtrb ~ age + group + gender + education + (1|time))
plot(wbModalCtrb ~ age, data = dataset)
points(dataset$age[order(dataset$age)], fitted(lme_mod)[order(dataset$age)], pch = 19)


## Spaghetti Plot of the fit longitudinal data -- doesn't work, need to separate dataset by group
lme_mod = lmer(data = dataset, wbModalCtrb ~ age + group + gender + education + (1|time))

pfit <-
  ggplot(data = dataset, aes(
    x = age,
    y = fitted(lme_mod),
    group = time)
  )
pfit + geom_point() + geom_line() 


# ====================================================================================
# ---- score vs age
dataset = read.csv("ctrb_aging_cogscore.csv")

# (1) Selecting a model that best fits the data 
lm_1 <- lmer(data = dataset, score ~ age*group + gender + education + iq + (1|time), REML = FALSE)
lm_2 <- lmer(data = dataset, score ~ age + group + gender + education + iq + (1|time), REML = FALSE)
lm_3 <- lmer(data = dataset, score ~ age + gender + education + iq + (1|time), REML = FALSE)
lm_4 <- lmer(data = dataset, score ~ age + gender + education + iq + (1|time) + (1|group), REML = FALSE)

models <- list(lm_1, lm_2, lm_3, lm_4)
model.names <- c('model1', 'model2', 'model3', 'model4')
aictab(cand.set = models, modnames = model.names)

# Model 2 is the best one, followed by model 3

# (2) Q: Is cog. score dependent on age and C/T group? 
cog.model = lmer(score ~ age + group + gender + education + iq + (1|time), data=dataset)
cog.model_2 = lmer(score ~ age + group + gender + education + iq + (1|time), data=dataset, REML = FALSE)
cog.null_wo_group = lmer(score ~ age + gender + education + iq + (1|time), data=dataset, REML = FALSE)
cog.null_wo_age = lmer(score ~ group + gender + education + iq + (1|time), data=dataset, REML = FALSE)

# test model against null hypothesis - do we use ANOVA? Or something else? I think something else... 
anova(cog.null_wo_group, cog.model_2)
anova(cog.null_wo_age, cog.model_2)


# ====================================================================================
# ---- executive function score vs modal ctrb
dataset = read.csv("ctrb_aging_cogscore.csv")

# (1) Selecting a model that best fits the data 
# Since modal ctrb depends on group, including both modal ctrb and group could cause overfitting and dependency between dependent variables. 
# Therefore the model for this should not include group but age. 
lm_1 <- lmer(data = dataset, score ~ wbModalCtrb*age + gender + education + iq + (1|time), REML = FALSE)
lm_2 <- lmer(data = dataset, score ~ wbModalCtrb + age + gender + education + iq + (1|time), REML = FALSE)
lm_3 <- lmer(data = dataset, score ~ wbModalCtrb + gender + education + iq + (1|time), REML = FALSE)
lm_4 <- lmer(data = dataset, score ~ wbModalCtrb + gender + education + iq + (1|time) + (1|age), REML = FALSE)

models <- list(lm_1, lm_2, lm_3, lm_4)
model.names <- c('model1', 'model2', 'model3', 'model4')
aictab(cand.set = models, modnames = model.names)

# Model 1 is the best one, followed by model 2

# (2) Q: Is cog. score dependent on modal controllability
cog.model = lmer(score ~ wbModalCtrb*age + gender + education + iq + (1|time), data=dataset)
cog.model_2 = lmer(score ~ wbModalCtrb*age + gender + education + iq + (1|time), data=dataset, REML = FALSE)
cog.null_wo_age = lmer(score ~ wbModalCtrb + gender + education + iq + (1|time), data=dataset, REML = FALSE)
cog.null_wo_ctrb = lmer(score ~ age + gender + education + iq + (1|time), data=dataset, REML = FALSE)

# test model against null hypothesis - do we use ANOVA? Or something else? I think something else... 
anova(cog.null_wo_ctrb, cog.model_2)
anova(cog.null_wo_age, cog.model_2)


# ====================================================================================
# ---- Functional subnetwork controllability with age
dataset = read.csv("Visual_ctrb_aging.csv")

lm_1 <- lmer(data = dataset, AvgCtrb ~ Age*Group + Group + Age + BioSex + EducationYrs + (1|Time), REML = FALSE)
lm_2 <- lmer(data = dataset, AvgCtrb ~ Group + Age + BioSex + EducationYrs + (1|Time), REML = FALSE)
lm_3 <- lmer(data = dataset, AvgCtrb ~ Age + BioSex + EducationYrs + (1|Time), REML = FALSE)
lm_4 <- lmer(data = dataset, AvgCtrb ~ Group + BioSex + EducationYrs + (1|Time), REML = FALSE)

models <- list(lm_1, lm_2, lm_3, lm_4)
model.names <- c('model1', 'model2', 'model3', 'model4')
aictab(cand.set = models, modnames = model.names)

best_model = lmer(AvgCtrb ~ Group + Age + BioSex + EducationYrs + (1|Time), data=dataset)
best_model_2 = lmer(AvgCtrb ~ Group + Age + BioSex + EducationYrs + (1|Time), data=dataset, REML = FALSE)
null_model_wo_group = lmer(AvgCtrb ~ Age + BioSex + EducationYrs + (1|Time), data=dataset, REML = FALSE)
null_model_wo_age = lmer(AvgCtrb ~ Group + BioSex + EducationYrs + (1|Time), data=dataset, REML = FALSE)

anova(null_model_wo_group, best_model_2)
anova(null_model_wo_age, best_model_2)

summary(best_model)

# ---------------------------------------------------------------------------------
# ---- modal ctrb
lm_1 <- lmer(data = dataset, ModalCtrb ~ Age*Group + Group + Age + BioSex + EducationYrs + (1|Time), REML = FALSE)
lm_2 <- lmer(data = dataset, ModalCtrb ~ Group + Age + BioSex + EducationYrs + (1|Time), REML = FALSE)
lm_3 <- lmer(data = dataset, ModalCtrb ~ Age + BioSex + EducationYrs + (1|Time), REML = FALSE)
lm_4 <- lmer(data = dataset, ModalCtrb ~ Group + BioSex + EducationYrs + (1|Time), REML = FALSE)

models <- list(lm_1, lm_2, lm_3, lm_4)
model.names <- c('model1', 'model2', 'model3', 'model4')
aictab(cand.set = models, modnames = model.names)

best_model = lmer(ModalCtrb ~ Group + Age + BioSex + EducationYrs + (1|Time), data=dataset)
best_model_2 = lmer(ModalCtrb ~ Group + Age + BioSex + EducationYrs + (1|Time), data=dataset, REML = FALSE)
null_model_wo_group = lmer(ModalCtrb ~ Age + BioSex + EducationYrs + (1|Time), data=dataset, REML = FALSE)
null_model_wo_age = lmer(ModalCtrb ~ Group + BioSex + EducationYrs + (1|Time), data=dataset, REML = FALSE)

anova(null_model_wo_group, best_model_2)
anova(null_model_wo_age, best_model_2)

# AvgCrb - age and grou pboth significantly impact the mean FPN avg ctrb 
# modalCtrb - only age significantly impact mean FPN modal ctrb. 


# dataset = read.csv("Ventral attention_ctrb_aging.csv")

# ====================================================================================
# ---- Hub nodes vs age
dataset = read.csv("avg_ctrb_hub_102_aging.csv")
dataset = read.csv("modal_ctrb_hub_101_aging.csv")
dataset = read.csv("modal_ctrb_hub_7_aging.csv")


lm_1 <- lmer(data = dataset, AvgCtrb ~ Age*Group + Group + Age + BioSex + EducationYrs + (1|Time), REML = FALSE)
lm_2 <- lmer(data = dataset, AvgCtrb ~ Group + Age + BioSex + EducationYrs + (1|Time), REML = FALSE)
lm_3 <- lmer(data = dataset, AvgCtrb ~ Age + BioSex + EducationYrs + (1|Time), REML = FALSE)
lm_4 <- lmer(data = dataset, AvgCtrb ~ Group + BioSex + EducationYrs + (1|Time), REML = FALSE)

models <- list(lm_1, lm_2, lm_3, lm_4)
model.names <- c('model1', 'model2', 'model3', 'model4')
aictab(cand.set = models, modnames = model.names)



lm_1 <- lmer(data = dataset, ModalCtrb ~ Age*Group + Group + Age + BioSex + EducationYrs + (1|Time), REML = FALSE)
lm_2 <- lmer(data = dataset, ModalCtrb ~ Group + Age + BioSex + EducationYrs + (1|Time), REML = FALSE)
lm_3 <- lmer(data = dataset, ModalCtrb ~ Age + BioSex + EducationYrs + (1|Time), REML = FALSE)
lm_4 <- lmer(data = dataset, ModalCtrb ~ Group + BioSex + EducationYrs + (1|Time), REML = FALSE)

models <- list(lm_1, lm_2, lm_3, lm_4)
model.names <- c('model1', 'model2', 'model3', 'model4')
aictab(cand.set = models, modnames = model.names)

mod_ctrb.model = lmer(ModalCtrb ~ Age + BioSex + EducationYrs + (1|Time), data=dataset)
mod_ctrb.model_2 = lmer(ModalCtrb ~ Age + BioSex + EducationYrs + (1|Time), data=dataset, REML = FALSE)
mod_ctrb.null = lmer(ModalCtrb ~ BioSex + EducationYrs + (1|Time), data=dataset, REML = FALSE)

anova(mod_ctrb.null, mod_ctrb.model_2)

# ====================================================================================
# ---- Cog score and Ctrb
dataset = read.csv("processing_speed_avg_ctrb.csv")
dataset = read.csv("executive_function_wais_ctrb.csv")

lm_1 <- lmer(data = dataset, CogScore ~ VIS*FPN + VIS + FPN + Group + BioSex + EducationYrs + Iq + (1|Time), REML = FALSE)
lm_2 <- lmer(data = dataset, CogScore ~ FPN + Group + BioSex + EducationYrs + Iq + (1|Time), REML = FALSE)
lm_3 <- lmer(data = dataset, CogScore ~ VIS + Group + BioSex + EducationYrs + Iq + (1|Time), REML = FALSE)
lm_4 <- lmer(data = dataset, CogScore ~ Group + BioSex + EducationYrs + Iq + (1|Time), REML = FALSE)

models <- list(lm_1, lm_2, lm_3, lm_4)
model.names <- c('model1', 'model2', 'model3', 'model4')
aictab(cand.set = models, modnames = model.names)

best_model = lmer(CogScore ~ FPN + Group + BioSex + EducationYrs + Iq + (1|Time), data = dataset)
best_model_2 = lmer(CogScore ~ FPN + Group + BioSex + EducationYrs + Iq + (1|Time), data = dataset, REML=FALSE)
null_model_wo_ctrb = lmer(CogScore ~ Group + BioSex + EducationYrs + Iq + (1|Time), data = dataset, REML=FALSE)
null_model_wo_group = lmer(CogScore ~ FPN + BioSex + EducationYrs + Iq + (1|Time), data = dataset, REML=FALSE)

anova(null_model_wo_ctrb, best_model_2)
anova(null_model_wo_group, best_model_2)

summary(best_model)

# -------
lm_1 <- lmer(data = dataset, CogScore ~ CON*FPN + CON + FPN + Group + BioSex + EducationYrs + Iq + (1|Time), REML = FALSE)
lm_2 <- lmer(data = dataset, CogScore ~ FPN + Group + BioSex + EducationYrs + Iq + (1|Time), REML = FALSE)
lm_3 <- lmer(data = dataset, CogScore ~ CON + Group + BioSex + EducationYrs + Iq + (1|Time), REML = FALSE)
lm_4 <- lmer(data = dataset, CogScore ~ Group + BioSex + EducationYrs + Iq + (1|Time), REML = FALSE)

models <- list(lm_1, lm_2, lm_3, lm_4)
model.names <- c('model1', 'model2', 'model3', 'model4')
aictab(cand.set = models, modnames = model.names)

best_model = lmer(CogScore ~ CON*FPN + CON + FPN + Group + BioSex + EducationYrs + Iq + (1|Time), data = dataset)
best_model_2 = lmer(CogScore ~ CON*FPN + CON + FPN + Group + BioSex + EducationYrs + Iq + (1|Time), data = dataset, REML=FALSE)
null_model_wo_con_ctrb = lmer(CogScore ~ FPN + Group + BioSex + EducationYrs + Iq + (1|Time), data = dataset, REML=FALSE)
null_model_wo_fpn_ctrb = lmer(CogScore ~ CON + Group + BioSex + EducationYrs + Iq + (1|Time), data = dataset, REML=FALSE)
null_model_wo_group = lmer(CogScore ~ CON*FPN + CON + FPN + BioSex + EducationYrs + Iq + (1|Time), data = dataset, REML=FALSE)

anova(null_model_wo_con_ctrb, best_model_2)
anova(null_model_wo_fpn_ctrb, best_model_2)
anova(null_model_wo_group, best_model_2)



