# Analyst Confidence

library(tidyverse)
library(patchwork)
library(car)  #for the Anova function
library(emmeans)  #for the emmeans function = estimated marginal means
library(DHARMa)  # for testing LM assumptions
library(glmmTMB)
library(performance) #for getting r-squared value from glmmTMB
  #installer: install.packages(c("insight", "performance"), repos = 'https://easystats.r-universe.dev')
rm(list=ls())

options(contrasts=c("contr.sum","contr.poly"))
options(scipen = 999)     #remove scientific notation from plots

#### Read-in joined data file and make it easier to work with ####

d = read_csv("combined data deidentified.csv")

# Change names of operation tasks to be easily identifiable (all qx.7_0_n)
d = rename(d, 
           backgroundSubtract = qx.7_0_5_rank,
           flatField = qx.7_0_14_rank,
           lowPassFilter = qx.7_0_12_rank,
           deconvolution = qx.7_0_6_rank,
           zIntensityCorrect = qx.7_0_7_rank,
           thresholdManual = qx.7_0_1_rank,
           thresholdGlobal = qx.7_0_2_rank,
           thresholdLocal = qx.7_0_3_rank,
           segmentAlgorithm = qx.7_0_8_rank,
           machineLearning = qx.7_0_9_rank,
           binaryFilter = qx.7_0_4_rank,
           binaryOperation = qx.7_0_11_rank,
           outlierExclusion = qx.7_0_10_rank,
           other1 = qx.7_0_13_rank,
           other2 = qx.7_0_16_rank)

# make column for grouping by set type (nuclei or fish)
d = mutate(d, setType = str_sub(setName, end=-2) )



#### analyst confidence ####

# correlation between analyst confidence and LSA

nuclsa_confidence = ggplot( filter(d, setType=="nuclei"), aes(x=qx.5_1, y=lsa_mse_transformed))+
  geom_jitter(alpha=0.25, width=0.1, height=0)+
  scale_y_log10()+
  scale_x_continuous(breaks=c(0,2,4,6,8,10), limits=c(-0.1,10.1))+
  xlab("Confidence")+
  ylab("LSA MSE")

mod.nuclsa_confidence = glmmTMB( round(log10(lsa_mse_transformed*1000)*100) ~ qx.5_1 +(1|id), data=filter(d, setType=="nuclei"), family="nbinom1")
plot( simulateResiduals( mod.nuclsa_confidence)) 
Anova(mod.nuclsa_confidence, type=3)
r2(mod.nuclsa_confidence)


fishlsa_confidence = ggplot( filter(d, setType=="fish"), aes(x=qx.5_1, y=lsa_mse_transformed))+
  geom_jitter(alpha=0.25, width=0.1, height=0)+
  scale_y_log10()+
  scale_x_continuous(breaks=c(0,2,4,6,8,10), limits=c(-0.1,10.1))+
  xlab("Confidence")+
  ylab("LSA MSE")

mod.fishlsa_confidence = glmmTMB( round(log10(lsa_mse_transformed*1000)*100) ~ qx.5_1 +(1|id), data=filter(d, setType=="fish"), family="nbinom1")
plot( simulateResiduals( mod.fishlsa_confidence)) 
Anova(mod.fishlsa_confidence, type=3)
r2(mod.fishlsa_confidence)


fishjac_confidence = ggplot(d, aes(x=qx.5_1, y=jac_transformed))+
  geom_jitter(alpha=0.25, width=0.1, height=0)+
  scale_x_continuous(breaks=c(0,2,4,6,8,10), limits=c(-0.1,10.1))+
  xlab("Confidence")+
  ylab("Jaccard Index")

mod.fishjac_confidence = glmmTMB( log10(1-jac_transformed) ~ sqrt(10-qx.5_1) + (1|id), data=filter(d, setType=="fish"), family="gaussian" )
plot( simulateResiduals( mod.fishjac_confidence)) 
Anova(mod.fishjac_confidence, type=3)
r2(mod.fishjac_confidence)



#### analyst self-expressed expertise ####

nuclsa_expertise = ggplot( filter(d, setType=="nuclei"), aes(x=q4.2_1, y=lsa_mse_transformed))+
  geom_jitter(alpha=0.25, width=0.1, height=0)+
  scale_y_log10()+
  scale_x_continuous(breaks=c(0,2,4,6,8,10), limits=c(-0.1,10.1))+
  xlab("Expertise")+
  ylab("LSA MSE")

mod.nuclsa_expertise = glmmTMB( round(log10(lsa_mse_transformed*1000)*100) ~ q4.2_1 +(1|id), data=filter(d, setType=="nuclei"), family="nbinom1")
plot( simulateResiduals( mod.nuclsa_expertise)) 
Anova(mod.nuclsa_expertise, type=3)
r2(mod.nuclsa_expertise)


fishlsa_expertise = ggplot( filter(d, setType=="fish"), aes(x=q4.2_1, y=lsa_mse_transformed))+
  geom_jitter(alpha=0.25, width=0.1, height=0)+
  scale_y_log10()+
  scale_x_continuous(breaks=c(0,2,4,6,8,10), limits=c(-0.1,10.1))+
  xlab("Expertise")+
  ylab("LSA MSE")

mod.fishlsa_expertise = glmmTMB( round(log10(lsa_mse_transformed*1000)*100) ~ q4.2_1 +(1|id), data=filter(d, setType=="fish"), family="nbinom1")
plot( simulateResiduals( mod.fishlsa_expertise)) 
Anova(mod.fishlsa_expertise, type=3)
r2(mod.fishlsa_expertise)


fishjac_expertise = ggplot( filter(d, setType=="fish"), aes(x=q4.2_1, y=jac_transformed))+
  geom_jitter(alpha=0.25, width=0.1, height=0)+
  scale_x_continuous(breaks=c(0,2,4,6,8,10), limits=c(-0.1,10.1))+
  xlab("Expertise")+
  ylab("Jaccard Index")

mod.fishjac_expertise = glmmTMB( log10(1-jac_transformed) ~ sqrt(10-q4.2_1) + (1|id), data=filter(d, setType=="fish"), family="gaussian" )
plot( simulateResiduals( mod.fishjac_expertise)) 
Anova(mod.fishjac_expertise, type=3)
r2(mod.fishjac_expertise)


#### Correlation betweeen expertise and confidence ####

mod.expert = glmmTMB( sqrt(10-qx.5_1) ~ sqrt(10-q4.2_1) + (1|id), data=d, family="gaussian")
plot( simulateResiduals( mod.expert))
Anova(mod.expert, type=3)
r2(mod.expert)
  #conditional R2 is for both fixed and random effects
  #marginal R2 is fixed effects only

correlation = ggplot(d, aes(q4.2_1, qx.5_1) )+
  geom_jitter(height=0.2, width=0.2, alpha=0.25)+
  xlab("Expertise")+
  ylab("Confidence")+
  #coord_fixed(ratio=1) +
  scale_x_continuous(breaks=c(0,2,4,6,8,10))+
  scale_y_continuous(breaks=c(0,2,4,6,8,10))
  


#### Plot ####
# default plot theme
theme_set( theme_classic()+
             theme(
               axis.title.x = element_text(size=12, face="bold"),
               axis.text.x = element_text(size=10, color="black"),
               axis.title.y = element_text(size=12, face="bold"),
               axis.text.y = element_text(size=10, color="black")
             ))

(nuclsa_confidence | nuclsa_expertise) / 
  (fishlsa_confidence | fishlsa_expertise) /
  (fishjac_confidence | fishjac_expertise) /
  (correlation + plot_spacer()) +
  plot_annotation(tag_levels='A')+
  plot_layout(widths=1)&
  theme(plot.tag=element_text(size=14, face="bold"))




