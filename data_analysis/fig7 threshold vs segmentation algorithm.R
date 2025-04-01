# Thresholding vs. segmentation algorithm

library(tidyverse)
library(patchwork)
library(car)  #for the Anova function
library(emmeans)  #for the emmeans function = estimated marginal means
library(DHARMa)  # for testing LM assumptions
library(glmmTMB)
rm(list=ls())

options(contrasts=c("contr.sum","contr.poly"))
options(scipen = 999)


#### Read-in joined data file and make it easier to work with ####

# Joined Data including coded "other" steps
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

# Convert to long-form by operation
operationOrder = pivot_longer(d, cols=backgroundSubtract:other2, names_to="operation", values_to="rank", values_drop_na=TRUE)

# replace "other1" and "other2" with actual operation names (from the coded columns)
operationOrder = mutate(operationOrder, operation = case_when(
  operation=="other1" ~ qx.7_13_text_coded,
  operation=="other2" ~ qx.7_16_text_coded,
  .default = operation
))

# Group operations into higher-order (more general) categories
#view(unique(operationOrder$operation)) # list all operations
operationOrder = mutate( operationOrder, operationCategory = case_when(
  operation %in% c("binaryFilter","connectedComponents","splitTouchingObjects") ~"binaryOperation",
  operation %in% c("backgroundSubtract","edgeDetection","flatField","lowPassFilter","sharpen") ~ "filter",
  operation %in% c("3DobjectDetection","distanceAssignmentClustering","spotDetection","surfaceCreationImaris") ~"segmentAlgorithm",
  str_detect(operation, "threshold") ~ "threshold",
  .default = operation
))



#### Basic Strategy: thresholding vs. segmentation ####

strategy = operationOrder %>%
  group_by(id, setName, setType, lsa_mse_transformed, jac_transformed) %>%
  summarize(thresh = any(operationCategory=="threshold"),
            seg = any(operationCategory %in% c("segmentAlgorithm","machineLearning")),
            threshNoSeg = thresh & !seg,  #used thresholding only (no seg)
            other = !thresh & !seg) #people who didn't include either ("same as previous" or likely error)


#mod.fish = glmmTMB( log10(lsa_mse_transformed)~threshNoSeg+(1|id), data=filter(strategy, !other, setType=="fish", jac_transformed>0.5))
mod.fish = glmmTMB( round(log10(lsa_mse_transformed*1000)*1000) ~threshNoSeg+(1|id), data=filter(strategy, !other, setType=="fish", jac_transformed>0.5), family="nbinom1")
plot( simulateResiduals( mod.fish))
Anova(mod.fish, type=3)
filter(strategy, !other, setType=="fish", jac_transformed>0.5) %>% 
  group_by(threshNoSeg) %>%
  summarize(n=n())

fishlsa = ggplot( filter(strategy, !other, setType=="fish", jac_transformed>0.5) , aes(x=threshNoSeg, y=lsa_mse_transformed))+
  geom_boxplot(outlier.shape=NA)+
  geom_jitter(alpha=0.5)+
  scale_y_log10()+
  scale_x_discrete(breaks=c("FALSE","TRUE"), labels=c("segmentation\nalgorithm","threshold"))+
  xlab("")+
  ylab("LSA MSE")


#mod.nuclei = glmmTMB( log10(lsa_mse_transformed)~threshNoSeg+(1|id), data=filter(strategy, !other, setType=="nuclei"))
mod.nuclei = glmmTMB( round(log10(lsa_mse_transformed*1000)*1000) ~threshNoSeg+(1|id), data=filter(strategy, !other, setType=="nuclei"), family="nbinom1")
plot( simulateResiduals( mod.nuclei)) 
Anova(mod.nuclei, type=3)
filter(strategy, !other, setType=="nuclei") %>%
  group_by(threshNoSeg) %>%
  summarize(n=n())

nuclsa = ggplot( filter(strategy, !other, setType=="nuclei") , aes(x=threshNoSeg, y=lsa_mse_transformed))+
  geom_boxplot(outlier.shape=NA)+
  geom_jitter(alpha=0.5)+
  scale_y_log10()+
  scale_x_discrete(breaks=c("FALSE","TRUE"), labels=c("segmentation\nalgorithm","threshold"))+
  xlab("")+
  ylab("LSA MSE")


mod.fish.jac = glmmTMB( log10(1-jac_transformed)~threshNoSeg+(1|id), data=filter(strategy, !other, setType=="fish"))
plot( simulateResiduals( mod.fish.jac))
Anova(mod.fish.jac, type=3)
filter(strategy, !other, setType=="fish") %>%
  group_by(threshNoSeg) %>%
  summarize(n=n())

fishjac = ggplot( filter(strategy, !other, setType=="fish") , aes(x=threshNoSeg, y=jac_transformed))+
  geom_boxplot(outlier.shape=NA)+
  geom_jitter(alpha=0.5)+
  scale_x_discrete(breaks=c("FALSE","TRUE"), labels=c("segmentation\nalgorithm","threshold"))+
  xlab("")+
  ylab("Jaccard")


# default plot theme
theme_set( theme_classic()+
             theme(
               axis.title.x = element_text(size=12, face="bold"),
               axis.text.x = element_text(size=10, color="black", angle=45, vjust=1, hjust=1),
               axis.title.y = element_text(size=12, face="bold"),
               axis.text.y = element_text(size=10, color="black")
             ))

# Plot
(nuclsa | fishlsa | fishjac) +
  plot_annotation(tag_levels='A') &
  theme(plot.tag=element_text(size=14, face="bold"))
