# software usage

library(tidyverse)
library(car)  #for the Anova function
library(emmeans)  #for the emmeans function = estimated marginal means
library(DHARMa)  # for testing LM assumptions
library(glmmTMB)
library(patchwork)
rm(list=ls())

options(contrasts=c("contr.sum","contr.poly"))
options(scipen = 999)     #remove scientific notation from plots


#### Read-in joined data file and make it easier to work with ####
d = read_csv("combined data deidentified.csv")

# make column for grouping by set type (nuclei or fish)
d = mutate(d, setType = str_sub(setName, end=-2) )

# create summary of qx.4 (did you write your own code etc.) to be binary (rather than 3 options)
d = mutate(d, wroteCode = !(qx.4=="No"))

# clean-up software (qx.2)
d = mutate(d, software = case_when(
  str_detect(qx.2, ",") ~ "multiple",  #if multiple software were used
  qx.2=="NIS Elements (Nikon)" ~ "Elements",
  qx.2=="Imaris (Bitplane)" ~ "Imaris",
  qx.2=="Image Pro Plus (Media Cybernetics)" ~ "Image-Pro",
  qx.2=="Cell Profiler" ~ "CellProfiler",
  .default = qx.2
))


#### Software usage ####
# each software only counted once per analyst (per set)
# so if an analyst used imagej for all images in nuclei, that is counted as 1
# but if an analyst used imagej for one and python for the other three, then that counts as 1 for imagej and 1 for python

useByImage = d %>% group_by(id, setType, software) %>% summarize(n=n())
useByAnalyst = useByImage %>% group_by(setType, software) %>% summarize(n=n())

usage_plot = ggplot(useByAnalyst, aes(x=n, y=software))+
  facet_grid(cols=vars(setType), labeller=labeller(setType=c("fish"="FISH","nuclei"="Nuclei")))+
  geom_col()+
  scale_y_discrete(limits=rev)+
  ylab("")+
  xlab("Number of Analysts")


#### Segmentation trends by software ####

# group low-usage sofware into "other"
d = mutate(d, software2 = case_when(
  software %in% c("Icy","Image-Pro","Elements") ~ "other",
  .default = software
))

nuclsa_full = ggplot(filter(d, setType=="nuclei", software2!="other"), aes(x=software2, y=lsa_mse_transformed))+
  geom_boxplot(outlier.shape=NA)+
  geom_jitter(alpha=0.5, width=0.2)+
  scale_y_log10()+
  ylab("Nuclei LSA MSE")+
  xlab("")+
  theme(
    axis.text.x = element_text(angle=90)
  )

fishlsa_full = ggplot(filter(d, setType=="fish", jac_transformed>0.5, software2!="other"), aes(x=software2, y=lsa_mse_transformed))+
  geom_boxplot(outlier.shape=NA)+
  geom_jitter(alpha=0.5, width=0.2)+
  scale_y_log10()+
  ylab("FISH LSA MSE")+
  xlab("")+
  theme(
    axis.text.x = element_text(angle=90)
  )

fishjac_full = ggplot(filter(d, setType=="fish", software2!="other"), aes(x=software2, y=jac_transformed))+
  geom_boxplot(outlier.shape=NA)+
  geom_jitter(alpha=0.5, width=0.2)+
  ylab("FISH Jaccard")+
  xlab("")+
  theme(
    axis.text.x = element_text(angle=90)
  )

# what are the multiple?
mult = d %>% distinct(id, setType, .keep_all=TRUE) %>% 
  filter(software=="multiple") %>% 
  select(qx.2, setType)



#### wrote own code ####

# create summary of qx.4 (did you write your own code etc.) to be binary (rather than 3 options)
d = mutate(d, wroteCode = !(qx.4=="No"))

# each wrote? only counted once per analyst per software (per set)
wroteByImage = d %>% group_by(id, setType, software2, wroteCode) %>% summarize(n=n())
wroteByAnalyst = wroteByImage %>% group_by(setType, software2, wroteCode) %>% summarize(n=n())

ggplot(filter(wroteByAnalyst, software2!="other"), aes(x=n, y=software2, fill=wroteCode))+
  facet_grid(cols=vars(setType))+
  geom_col(position=position_dodge2(preserve="single"))



#plot defaults
theme_set( theme_classic()+
             theme(
               axis.title.x = element_text(size=12, face="bold"),
               axis.text.x = element_text(size=10, color="black"),
               axis.title.y = element_text(size=12, face="bold"),
               axis.text.y = element_text(size=10, color="black"),
               strip.text = element_text(size=10, color="black")
             ))

# plot
(usage_plot | nuclsa_full) /
(fishlsa_full | fishjac_full) +
  plot_annotation(tag_levels='A') &
  theme(plot.tag=element_text(size=14, face="bold"))
