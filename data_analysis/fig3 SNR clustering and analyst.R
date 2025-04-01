# Analyze combined effects of SNR and clustering

library(tidyverse)
library(car)  #for the Anova function
library(emmeans)  #for the emmeans function = estimated marginal means
library(DHARMa)  # for testing LM assumptions
library(glmmTMB)
library(patchwork)
rm(list=ls())

# Use sum-to-zero contrasts (to match Type III sums of squares)
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


#### label variants by clustering and SNR ####
d = mutate(d, clustering = if_else(str_detect(setName,"2|4"), TRUE, FALSE) )
d = mutate(d, snr = if_else(str_detect(setName,"3|4"), TRUE, FALSE) )


#### By-analyst normalization ####
#  Was not able to run stats including analyst as fixed effect - had to normalize instead

# normalize to sample1, LSA MSE
lsaNorm = d %>%
  group_by(responseid, setType) %>%
  mutate(norm = lsa_mse_transformed/mean(lsa_mse_transformed, na.rm=TRUE) )
  
# check the distributions
ggplot(filter(lsaNorm, setType=="nuclei"), aes(sqrt(norm)))+
  geom_histogram(bins=20) + facet_wrap(~setName)

ggplot(filter(lsaNorm, setType=="fish"), aes(sqrt(norm)))+
  geom_histogram(bins=20) + facet_wrap(~setName)



#### Nuclei Stats ####
# gaussian model with clustering and snr interaction; sqrt transformation
mod.nuclei.norm = glmmTMB( sqrt(norm)~clustering+snr, data=filter(lsaNorm, setType=="nuclei"))
plot( simulateResiduals( mod.nuclei.norm))
Anova(mod.nuclei.norm, type=3)

lsaNorm %>% group_by(setType, clustering, snr) %>% summarize(n=n())


nuc_clust = ggplot(filter(d, setType=="nuclei"), aes(x=clustering, y=lsa_mse_transformed))+
  geom_boxplot(outlier.shape=NA)+
  geom_jitter(width=0.2, alpha=0.5)+
  scale_y_log10()+
  scale_x_discrete(limits=c("FALSE","TRUE"), labels=c("low","high"))+
  ylab("LSA MSE")+
  xlab("Object Clustering")

nuc_snr = ggplot(filter(d, setType=="nuclei"), aes(x=snr, y=lsa_mse_transformed))+
  geom_boxplot(outlier.shape=NA)+
  geom_jitter(width=0.2, alpha=0.5)+
  scale_y_log10()+
  scale_x_discrete(limits=c("FALSE","TRUE"), labels=c("high","low"))+
  ylab("LSA MSE")+
  xlab("Signal-to-Noise Ratio")


#### FISH LSA Stats ####

mod.fish.norm = glmmTMB( sqrt(norm)~clustering+snr, data=filter(lsaNorm, setType=="fish"))
plot( simulateResiduals( mod.fish.norm))
Anova(mod.fish.norm, type=3)

d %>% filter(setType=="fish") %>% group_by(clustering, snr) %>% summarize(n=n())

fishlsa_clust = ggplot(filter(d, setType=="fish"), aes(x=clustering, y=lsa_mse_transformed))+
  geom_boxplot(outlier.shape=NA)+
  geom_jitter(width=0.2, alpha=0.5)+
  scale_y_log10()+
  scale_x_discrete(limits=c("FALSE","TRUE"), labels=c("low","high"))+
  ylab("LSA MSE")+
  xlab("Object Clustering")

fishlsa_snr = ggplot(filter(d, setType=="fish"), aes(x=snr, y=lsa_mse_transformed))+
  geom_boxplot(outlier.shape=NA)+
  geom_jitter(width=0.2, alpha=0.5)+
  scale_y_log10()+
  scale_x_discrete(limits=c("FALSE","TRUE"), labels=c("high","low"))+
  ylab("LSA MSE")+
  xlab("Signal-to-Noise Ratio")


#### FISH Jaccard Stats ####
mod.fish.jac = glmmTMB( ( log10(1-jac_transformed) ) ~clustering+snr+(1|id), data=filter(d, setType=="fish") )
plot( simulateResiduals( mod.fish.jac))
Anova(mod.fish.jac, type=3)

fishjac_clust = ggplot(filter(d, setType=="fish"), aes(x=clustering, y=jac_transformed))+
  geom_boxplot(outlier.shape=NA)+
  geom_jitter(width=0.2, alpha=0.5)+
  scale_x_discrete(limits=c("FALSE","TRUE"), labels=c("low","high"))+
  ylab("Jaccard")+
  xlab("Object Clustering")

fishjac_snr = ggplot(filter(d, setType=="fish"), aes(x=snr, y=jac_transformed))+
  geom_boxplot(outlier.shape=NA)+
  geom_jitter(width=0.2, alpha=0.5)+
  scale_x_discrete(limits=c("FALSE","TRUE"), labels=c("high","low"))+
  ylab("Jaccard")+
  xlab("Signal-to-Noise Ratio")



#### Plot ####
# default plot theme
theme_set( theme_classic()+
             theme(
               axis.title.x = element_text(size=12, face="bold"),
               axis.text.x = element_text(size=10, color="black"),
               axis.title.y = element_text(size=12, face="bold"),
               axis.text.y = element_text(size=10, color="black")
             ))

(nuc_clust | nuc_snr)/(fishlsa_clust | fishlsa_snr)/(fishjac_clust | fishjac_snr) +
  plot_annotation(tag_levels='A') &
  theme(plot.tag=element_text(size=14, face="bold"))




