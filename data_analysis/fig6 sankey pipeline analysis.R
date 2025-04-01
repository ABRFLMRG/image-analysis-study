# Sankey Pipeline Operation Order Analysis

#If ggsankey isn't already installed
#install.packages("remotes") #ggsankey isn't on CRAN so need to download from github
#remotes::install_github("davidsjoberg/ggsankey")
library(tidyverse)
library(ggsankey)
library(patchwork)
rm(list=ls())



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

#anonymous IDs actually start at 10, but it looks weird when we plot, so simply shifting to start at 1
d = mutate(d, id = id-10)  

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


#### Sankey Diagrams ####

# There isn't much space on the Sankey Diagrams for labels, so abbreviate all operation names.
#  Yes, I know it would have been easier/better to do this earlier in the code, but hindsight is 20-20.
operationOrder = mutate( operationOrder, operationAbb = case_when(
  operationCategory=="deconvolution" ~ "D",
  operationCategory=="machineLearning" ~ "ML",
  operationCategory=="binaryOperation" ~ "BO",
  operationCategory=="segmentAlgorithm" ~ "SA",
  operationCategory=="manualCorrection" ~ "ML",
  operationCategory=="masking" ~ "M",
  operationCategory=="threshold" ~ "T",
  operationCategory=="maximaDetection" ~ "MD",
  operationCategory=="normalization" ~ "N",
  operationCategory=="outlierExclusion" ~ "OE",
  operationCategory=="resizeOrRescale" ~ "RR",
  operationCategory=="filter" ~ "F"
))

#converts back to wide, but with columns as ranks
operationOrder_wide = pivot_wider(operationOrder, id_cols=c("responseid","setName", "setType"),
                                  names_from="rank", values_from="operationAbb", names_prefix="step", names_sort=TRUE)

# convert to format needed for ggsankey; there are only 7 ranks (no one used >7 steps)
sankeyNuclei = operationOrder_wide %>% 
  filter(setName=="nuclei1") %>%
  make_long(step1, step2, step3, step4, step5, step6, step7)

sankeyNuclei = filter(sankeyNuclei, !is.na(node)) #removes NA

nuc_sankey_plot = ggplot(sankeyNuclei, aes(x=x, next_x=next_x, node=node, next_node=next_node, fill=factor(node), label=node))+
  geom_sankey(flow.alpha = 0.5, node.color = "black", show.legend = FALSE)+
  geom_sankey_label(size = 3, color = "black", fill= "white", alpha=0.5)+
  scale_x_discrete(limits=c("step1","step2","step3","step4","step5","step6","step7"),
                   labels=c("1","2","3","4","5","6","7"))+
  #scale_fill_viridis_d(option="turbo") +
  scale_fill_viridis_d(option="turbo", limits=c("D","ML","BO","SA","ML","M","T","MD","N","OE","RR","F"))+
  #theme_sankey(base_size = 12)+
  labs(x="Step", y="Nuclei Image #1 Pipelines")

# convert to format needed for ggsankey; no one used > 5 steps for the FISH
sankeyFish = operationOrder_wide %>% 
  filter(setName=="fish1") %>%
  make_long(step1, step2, step3, step4, step5)

sankeyFish = filter(sankeyFish, !is.na(node)) #removes NA

fish_sankey_plot = ggplot(sankeyFish, aes(x=x, next_x=next_x, node=node, next_node=next_node, fill=factor(node), label=node))+
  geom_sankey(flow.alpha = 0.5, node.color = "black", show.legend = FALSE)+
  geom_sankey_label(size = 3, color = "black", fill= "white", alpha=0.5)+
  scale_x_discrete(limits=c("step1","step2","step3","step4","step5"),
                   labels=c("1","2","3","4","5"))+
  scale_fill_viridis_d(option="turbo", limits=c("D","ML","BO","SA","ML","M","T","MD","N","OE","RR","F"))+
  #theme_sankey(base_size = 12)+
  labs(x="Step", y="FISH Image #1 Pipelines")


# default plot theme
theme_set( theme_sankey(base_size = 12)+
             theme(
               axis.title.x = element_text(size=12, face="bold"),
               axis.text.x = element_text(size=10, color="black"),
               axis.title.y = element_text(size=12, face="bold")
             ))
#plot
nuc_sankey_plot / fish_sankey_plot + 
  plot_annotation(tag_levels='A') &
  theme(plot.tag=element_text(size=14, face="bold"))


