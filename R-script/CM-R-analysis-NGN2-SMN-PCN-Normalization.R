
#load required libraries
library(tidyverse)
library(cowplot)
library(ggrepel)
library(factoextra)
library(DescTools)
library(Rtsne)
theme_set(theme_cowplot())

#clear workspace and set options 
rm(list=ls(all=T))
pdf.options(useDingbats=F, width=8, height=8)
setwd("/Data-Experiments/Analysis-R/R-SMN-NGN2-PCN/Normalization")

### To load in all csv files you want to normalize and combine into one (data1)
  #Example with STING Area per cell - can be aplied to any other parameter
temp = list.files(pattern = ".*.csv")
files = lapply(temp, read.csv)
data1 = do.call(bind_rows, files)

pdf("SMN-NGN2-PCN-Normalization-AreaPerCell.pdf")

data2 <- data1 %>% group_by(genotype,line,well,replicate) %>% summarize(
  AreaNucl_well=mean(area_STING_total_site),
  TUJ_well=mean(tuj_site),
  AreaPerCell=AreaNucl_well/TUJ_well)

# calculate mean of control wells within each replicate
data3 <- data2 %>% filter(genotype == "wt") %>% group_by(replicate) %>%  summarize(
  mean_area_rep = mean(AreaPerCell)
)

## normalizing to each respective replicate, see data3 for values (example values are added below)
data5 <-data2 %>% mutate(norm = case_when(
  replicate == "1" ~ 7,
  replicate == "2" ~ 56), norm = as.numeric(norm))

data5 <- data5 %>% mutate(norm_effect = AreaPerCell/norm)

data5$replicate = as.factor(data5$replicate)

# plotting all replicates 
ggplot(data5, aes(x=line,
             y=norm_effect)) +
  geom_jitter(aes(color=line, shape = replicate),height=0,width=0.1,size=3,alpha=0.7) + theme_bw() +
  theme(strip.text.x = element_text(size = 8, colour = "black", angle = 90))+
  ggtitle("Normalized Total Area of STING \n(per cell per site per well)") +
  ylim(0,5) + ylab("Normalized Total Area of STING \n(per cell per site per well)") + 
  xlab("") + theme(axis.title = element_text(size = 10)) + #+
  stat_summary(fun.data=mean_se,geom="crossbar",color="darkgrey",fill="grey10",alpha=0.2) +
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = 'bold')) +
  theme(axis.text.x=element_text(angle=30,hjust=1)) +
  scale_color_manual(values=c("grey30","deepskyblue")) 

dev.off()

### run a t-test - value example below
t.test(norm_effect ~ line, data = data5)
# P = 0.015



