# Pipeline Operation Order Analysis


library(tidyverse)
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




#### Pipeline Order Comparison ####

pipe_nuc = ggplot( filter(operationOrder, setType=="nuclei"), aes(x=rank, y=setName))+
  facet_grid(rows=vars(id)) +
  geom_tile(aes(fill=operationCategory), color="black")+
  scale_fill_viridis_d(option="turbo", 
                       limits=c("deconvolution","machineLearning","binaryOperation","segmentAlgorithm","manualCorrection","masking","threshold","maximaDetection","normalization","outlierExclusion","resizeOrRescale","filter"),
                       labels=c("deconvolution","machine learning","binary operation","segmentation algorithm","manual correction","masking","threshold","maxima detection","normalization","outlier exclusion","resize/rescale","filter"))+
  scale_x_discrete(limits=factor(c(1,2,3,4,5,6,7)), expand=c(0,0))+
  xlab("Order in Pipeline")+
  ylab("Nuclei Sets by Analyst")

pipe_fish = ggplot( filter(operationOrder, setType=="fish"), aes(x=rank, y=setName))+
  facet_grid(rows=vars(id)) +
  geom_tile(aes(fill=operationCategory), color="black")+
  scale_fill_viridis_d(option="turbo", 
                       limits=c("deconvolution","machineLearning","binaryOperation","segmentAlgorithm","manualCorrection","masking","threshold","maximaDetection","normalization","outlierExclusion","resizeOrRescale","filter"),
                       labels=c("deconvolution","machine learning","binary operation","segmentation algorithm","manual correction","masking","threshold","maxima detection","normalization","outlier exclusion","resize/rescale","filter"))+
  scale_x_discrete(limits=factor(c(1,2,3,4,5,6,7)), expand=c(0,0))+
  xlab("Order in Pipeline")+
  ylab("FISH Sets by Analyst")


# default plot theme
theme_set( theme_gray()+
             theme(
               panel.grid.minor=element_blank(),
               panel.grid.major=element_blank(),
               axis.text.y=element_blank(),
               axis.ticks.y=element_blank(),
               legend.title=element_blank(),
               axis.title.x = element_text(size=12, face="bold"),
               axis.text.x = element_text(size=10, color="black"),
               axis.title.y = element_text(size=12, face="bold"),
               strip.text = element_text(size=10)
             ))
#plot
pipe_nuc + 
  pipe_fish + 
  plot_annotation(tag_levels='A')+
  plot_layout(guides="collect") & 
  theme(legend.position="bottom", plot.tag=element_text(size=14, face="bold")) 



