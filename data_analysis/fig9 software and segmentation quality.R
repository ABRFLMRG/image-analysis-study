# software and segmentation quality

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

# Group software by "type" 
#  Type based on whether there were multiple software used
#   and whether they are more programming/free-form or pre-made macro-based
d = mutate(d, softType = case_when(
  str_detect(qx.2, ",") ~ "multiple",  #if multiple software were used
  qx.2 %in% c("Python", "FIJI/ImageJ", "Cell Profiler") ~ "python/imagej/cellProfiler",
  qx.2 %in% c("Imaris (Bitplane)","Arivis","Aivia") ~ "imaris/arivis/aivia",
  .default = "other" #Elements, Icy, ImagePro...
))



#### Software and segmentation quality ####
mod.nuclei = glmmTMB( round(log10(lsa_mse_transformed*1000)*1000) ~softType+(1|id), data=filter(d, setType=="nuclei", softType != "other", softType != "multiple"), family="nbinom1")
plot( simulateResiduals( mod.nuclei)) 
Anova(mod.nuclei, type=3)
filter(d, setType=="nuclei", softType != "other", softType != "multiple") %>%
  group_by(softType) %>%
  summarize(n=n())

nuclsa_plot = ggplot(filter(d, setType=="nuclei", softType %in% c("python/imagej/cellProfiler","imaris/arivis/aivia")), 
                    aes(x=softType, y=lsa_mse_transformed))+
  geom_boxplot(outlier.shape=NA)+
  geom_jitter(alpha=0.5, width=0.2)+
  scale_y_log10()+
  scale_x_discrete(limits=c("imaris/arivis/aivia","python/imagej/cellProfiler"), 
                   labels=c("Imaris\nArivis\nAivia","Python\nImagej\nCellProfiler"))+
  ylab("LSA MSE")+
  xlab("Software type")

mod.fish.lsa = glmmTMB( round(log10(lsa_mse_transformed*1000)*100) ~softType+(1|setName)+(1|id), 
                        data=filter(d, setType=="fish", softType != "other", softType != "multiple", jac_transformed>0.5), 
                        family="nbinom1")
plot( simulateResiduals( mod.fish.lsa)) 
Anova(mod.fish.lsa, type=3)
filter(d, setType=="fish", softType != "other", softType != "multiple", jac_transformed>0.5) %>%
  group_by(softType) %>%
  summarize(n=n())


fishlsa_plot = ggplot(filter(d, setType=="fish", jac_transformed>0.5, softType %in% c("python/imagej/cellProfiler","imaris/arivis/aivia")), 
                     aes(x=softType, y=lsa_mse_transformed))+
  geom_boxplot(outlier.shape=NA)+
  geom_jitter(alpha=0.5, width=0.2)+
  scale_y_log10()+
  scale_x_discrete(limits=c("imaris/arivis/aivia","python/imagej/cellProfiler"), 
                   labels=c("Imaris\nArivis\nAivia","Python\nImagej\nCellProfiler"))+
  ylab("LSA MSE")+
  xlab("Software type")

mod.fish.jac = glmmTMB( (log10(1-jac_transformed)) ~softType+(1|setName)+(1|id), data=filter(d, setType=="fish", softType != "other", softType != "multiple"), family="gaussian")
plot( simulateResiduals( mod.fish.jac)) 
Anova(mod.fish.jac, type=3)
filter(d, setType=="fish", softType != "other", softType != "multiple") %>%
  group_by(softType) %>%
  summarize(n=n())


fishjac_plot = ggplot(filter(d, setType=="fish", softType %in% c("python/imagej/cellProfiler","imaris/arivis/aivia")), 
                     aes(x=softType, y=jac_transformed))+
  geom_boxplot(outlier.shape=NA)+
  geom_jitter(width=0.2, alpha=0.5)+
  scale_x_discrete(limits=c("imaris/arivis/aivia","python/imagej/cellProfiler"), 
                   labels=c("Imaris\nArivis\nAivia","Python\nImagej\nCellProfiler"))+
  ylab("Jaccard Index")+
  xlab("Software type")



#### Outliers and software ####
# Nuclei outliers defined by me as LSA MSE > 10
# FISH outliers defined by me as Jaccard < 0.5
#  (FISH LSA MSE didn't have well-defined outliers)

dnuc = d %>% filter(setType=="nuclei", !is.na(lsa_mean_transformed)) %>% mutate(outlier = lsa_mse_transformed >10)
dfish = d %>% filter(setType=="fish", !is.na(jac_transformed)) %>% mutate(outlier = jac_transformed <0.5)

nuclsa_out = ggplot(dnuc, aes(x=softType, group=outlier, fill=outlier))+
  geom_bar(position = position_stack(reverse = TRUE))+
  xlab("Software type")+
  scale_x_discrete(limits=c("imaris/arivis/aivia", "multiple", "other", "python/imagej/cellProfiler"), 
                   labels=c("Imaris\nArivis\nAivia","multiple","other","Python\nImagej\nCellProfiler"))+
  scale_fill_grey(breaks=c("FALSE","TRUE"), labels=c("normal", "outlier"), start=0.3, end=0.6)+
  ylab("Count")

dnuc %>% filter(outlier) %>% select(id, setName, softType, qx.2)
#The "multiple" used cellprofiler and python, #11 was only one to use imaris


fishjac_out = ggplot(dfish, aes(x=softType, group=outlier, fill=outlier))+
  geom_bar(position = position_stack(reverse = TRUE))+
  xlab("Software type")+
  scale_x_discrete(limits=c("imaris/arivis/aivia", "multiple", "other", "python/imagej/cellProfiler"), 
                   labels=c("Imaris\nArivis\nAivia","multiple","other","Python\nImagej\nCellProfiler"))+
  scale_fill_grey(breaks=c("FALSE","TRUE"), labels=c("normal", "outlier"), start=0.3, end=0.6)+
  ylab("Count")

dfish %>% filter(outlier) %>% select(id, setName, softType, qx.2)
# other = Icy, multiple = python + cell profiler



#### Figure ####

#plot defaults
theme_set( theme_classic()+
             theme(
               axis.title.x = element_text(size=12, face="bold"),
               axis.text.x = element_text(size=10, color="black"),
               axis.title.y = element_text(size=12, face="bold"),
               axis.text.y = element_text(size=10, color="black"),
               strip.text = element_text(size=10, color="black"),
               legend.title = element_blank()
             ))

# plot
(nuclsa_plot + nuclsa_out) / (fishlsa_plot + fishjac_out) / (fishjac_plot + plot_spacer()) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag=element_text(size=14, face="bold"))

((nuclsa_plot / fishlsa_plot / fishjac_plot) | (nuclsa_out / fishjac_out / plot_spacer())) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag=element_text(size=14, face="bold")) 
  