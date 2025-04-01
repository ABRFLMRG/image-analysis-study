# Basic participant demographic data

library(tidyverse)
library(patchwork)
rm(list=ls())


#### Read-in joined data file and make it easier to work with ####

d = read_csv("combined data deidentified.csv")

# make column for grouping by set type (nuclei or fish)
d = mutate(d, setType = str_sub(setName, end=-2) )

# make reduced dataset with just one entry per participant and demographic data only
demo = d %>% select(recordeddate:q4.4, q3.10coded)
demo = distinct(demo, id, .keep_all=TRUE)



#### n for image sets ####
d %>% group_by(setName) %>% summarize(n=n())



#### Participant Expertise ####

# number of analysts of each degree
degree = demo %>% group_by(q3.9) %>% summarize(n=n())
plotDegree = ggplot(degree, aes(x=n, y=q3.9))+
  geom_point(size=3)+
  geom_segment( aes(x=0, xend=n, y=q3.9, yend=q3.9), linewidth=1)+
  scale_y_discrete(limits=c("Doctorate","Masters","Bachelors"))+
  xlab("Count")+
  ylab("Highest Degree")


# number of analysts of each expertise
demo = demo %>% mutate( q3.4mod = case_when(  
  str_detect(q3.4, "Core") ~ "Core Facility",             #group all core facility together
  q3.4=="Other" ~ "Industry",                             #all the "others" were industry-related
  q3.4=="Faculty/Principal Investigator" ~ "Faculty/PI",  #rename
  .default = q3.4
))
position = demo %>% group_by(q3.4mod) %>% summarize(n=n())
plotPosition = ggplot(position, aes(x=n, y=q3.4mod))+
  geom_point(size=3)+
  geom_segment( aes(x=0, xend=n, y=q3.4mod, yend=q3.4mod), size=1)+
  xlab("Count")+
  ylab("Position")


# number of analysts of each degree area
field = demo %>% group_by(q3.10coded) %>% summarize(n=n())
plotField = ggplot(field, aes(x=n, y=q3.10coded))+
  geom_point(size=3)+
  geom_segment( aes(x=0, xend=n, y=q3.10coded, yend=q3.10coded), size=1)+
  scale_y_discrete(limits=c("biology","bioengineering","biophysics","chemistry","computerEngineering","electricalEngineering","math","medicalImageProcessing"),
                   labels=c("Biology","Bioengineering","Biophysics","Chemistry","Computer Eng.","Electrical Eng.","Math","Image Processing"))+
  xlab("Count")+
  ylab("Field of Highest Degree")



#### Location ####

# Institution
institution = demo %>% group_by(q3.7) %>% summarize(n=n())
plotInstitution = ggplot(institution, aes(x=n, y=q3.7))+
  geom_point(size=3)+
  geom_segment( aes(x=0, xend=n, y=q3.7, yend=q3.7), size=1)+
  scale_y_discrete(limits=c("Academic - University","Academic - not University","Business","Government"),
                   labels=c("University","Other Academic","Business","Government"))+
  xlab("Count")+
  ylab("Type of Institution")


# Country
country = demo %>% group_by(q3.8) %>% summarize( n = n() )
plotCountry = ggplot( filter(country, !is.na(q3.8)), aes(x=n, y=q3.8))+
  geom_point(size=3)+
  geom_segment( aes(x=0, xend=n, y=q3.8, yend=q3.8), size=1)+
  scale_y_discrete(limits=c("USA","United Kingdom","Germany","Chile","Canada","Belgium"))+
  xlab("Count")+
  ylab("Country")



#### Self-expressed image analysis expertise ####
expertise = demo %>% group_by(q4.2_1) %>% summarize(n=n())
plotExpertise = ggplot( expertise, aes(x=n, y=q4.2_1))+
  geom_point(size=3)+
  geom_segment( aes(x=0, xend=n, y=q4.2_1, yend=q4.2_1), size=1)+
  scale_y_continuous(breaks=c(2,4,6,8,10))+
  xlab("Count")+
  ylab("Self-rated Expertise")



#### Plots ####
# default plot theme
theme_set( theme_classic()+
             theme(
               axis.title.x = element_text(size=12, face="bold"),
               axis.text.x = element_text(size=10, color="black"),
               axis.title.y = element_text(size=12, face="bold"),
               axis.text.y = element_text(size=10, color="black")
             ))

#layout
design = "AB
          CD
          EF"

plotField + plotDegree + plotPosition + plotExpertise + plotInstitution + plotCountry +
  plot_annotation(tag_levels = 'A')+
  plot_layout(design=design)&
  theme(plot.tag=element_text(size=14, face="bold"))
                                      