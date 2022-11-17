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
setwd("/Data-Experiments/Analysis-R/SMN-NGN2-PCN/example-STING/Data")

#create platemap dataframe with sample information -example
platemap<-data.frame(well=c("E03","E07",
                            "E04","E08",
                            "E05","E09",
                            "E06","E10"),
                     genotype=c("wt","mut",
                                "wt","mut",
                                "wt","mut",
                                "wt","mut"),
                     age=c("D35","D35",
                           "D35","D35",
                           "D35","D35",
                           "D35","D35"),
                     iso=c("isoTDP43","isoTDP43",
                            "isoTDP43","isoTDP43",
                            "isoTDP43","isoTDP43",
                            "isoTDP43","isoTDP43"),
                     line=c("TDP43+/+","TDP43+/G298S",
                           "TDP43+/+","TDP43+/G298S",
                           "TDP43+/+","TDP43+/G298S",
                           "TDP43+/+","TDP43+/G298S"),
                     stringsAsFactors=F)

#identify files and pull metadata from file names
files<-list.files(pattern=".csv",recursive=T,full.names=F)
files_live<-files[which((files %>% str_detect("hoescht_count")))]

files_Areabody<-files[which((files %>% str_detect("Area_Tuj_body_count")))]
files_AreaNeurites<-files[which((files %>% str_detect("Area_Tuj_neurite_count")))]

files_STINGtotalcell<-files[which((files %>% str_detect("STING_Tuj_Total_intensity_count")))]
files_STINGtotalPeriNucl<-files[which((files %>% str_detect("STING_Tuj_Peri_Nuclear_intensity_count")))]
files_STINGtotalneurite<-files[which((files %>% str_detect("STING_Tuj_Neurite_intensity_count")))]
files_STINGTujcellparticle<-files[which((files %>% str_detect("STING_Tuj_Total_particle_count")))]


files<-files[which((files %>% str_detect("_count")))]
meta<-files %>% str_split(pattern="/",simplify=T) %>% as.data.frame(stringsAsFactors=F)
meta<-cbind(meta,meta$V4 %>% str_sub(1,10),stringsAsFactors=F) 
colnames(meta)<-c('sample','data','folder','subfold','key')
meta<-meta %>% select(sample,data,key)
meta$file<-files %>% str_sub(end=-5L)
meta<-cbind(meta,meta$key %>% str_split(pattern="_",simplify=T) %>% as.data.frame(stringsAsFactors=F))
colnames(meta)<-c('sample','data','key','file','well','channel','site')
meta<-meta %>% select(file,sample,key,well,site)
meta<-left_join(meta,platemap,by="well")
meta$key2<-str_c(meta$well,meta$site,sep="_")

# cells 
#load data from cell bodies count
data_live<-lapply(files_live,function(x) {read.csv(x,header=T,stringsAsFactors=F)})
names(data_live)<-files_live %>% str_split(patter="/",simplify=T) %>% .[,4] %>% str_split(pattern="\\.",simplify=T) %>% .[,1]
for(i in seq_along(data_live)) {
  data_live[[i]]$key<-rep(meta %>% filter(str_detect(file,"_hoescht_count")) %>% pull(key) %>% .[i],nrow(data_live[[i]])) }  # line of code to copy past
DT<-do.call(rbind,data_live)
colnames(DT)<-c("name","object_tuj","area_tuj","Average_tuj_Size","Percentage_tuj_Area","key")
DT<-cbind(DT,DT$key %>% str_split(pattern="_",simplify=T),stringsAsFactors=F)
colnames(DT)<-c("name","object_tuj","area_tuj","Average_tuj_Size","Percentage_tuj_Area","key","well","channel","site")
DT<-DT %>% select(key,object_tuj,area_tuj,Average_tuj_Size,Percentage_tuj_Area,well,channel,site)

##summarize cell count data
dt_summ<-DT %>% group_by(well)%>%
  summarize(tuj_well=mean(object_tuj),sites_well=length(unique(key))) %>%
  ungroup()
dt_summ2<-DT %>% group_by(well,site) %>%
  summarize(tuj_site=object_tuj,avg_area_tuj=Average_tuj_Size) %>%
  ungroup()


#### STING TOTAL

#total signal
#load data from intensity of STING in all cell count
data_STING_full_cell<-lapply(files_STINGtotalcell,function(x) {read.csv(x,header=T,stringsAsFactors=F)})
names(data_STING_full_cell)<-files_STINGtotalcell %>% str_split(patter="/",simplify=T) %>% .[,4] %>% str_split(pattern="\\.",simplify=T) %>% .[,1]
for(i in seq_along(data_STING_full_cell)) {
  data_STING_full_cell[[i]]$key<-rep(meta %>% filter(str_detect(file,"STING_Tuj_Total_intensity")) %>% pull(key) %>% .[i],nrow(data_STING_full_cell[[i]])) }  # line of code to copy past
DLs<-do.call(rbind,data_STING_full_cell)
colnames(DLs)<-c("object_STING_total","area_STING_total","MeanGrey_STING_total","StdDev_STING_total","MinGrey_STING_total",
                 "MaxGrey_STING_total","Feret_STING_total","IntDen_STING_total","RawIntDen_STING_total","FeretX_STING_total","FeretY_STING_total",
                 "Feretangle_STING_total","MinFeret_STING_total","MinThr_STING_total","MaxThr_STING_total","key")
DLs<-cbind(DLs,DLs$key %>% str_split(pattern="_",simplify=T),stringsAsFactors=F)
colnames(DLs)<-c("object_STING_total","area_STING_total","MeanGrey_STING_total","StdDev_STING_total","MinGrey_STING_total",
                 "MaxGrey_STING_total","Feret_STING_total","IntDen_STING_total","RawIntDen_STING_total","FeretX_STING_total","FeretY_STING_total",
                 "Feretangle_STING_total","MinFeret_STING_total","MinThr_STING_total","MaxThr_STING_total","key","well","channel","site")
DLs<-DLs %>% select(key,area_STING_total,MeanGrey_STING_total,IntDen_STING_total,RawIntDen_STING_total,well,channel,site)
DLs<-DLs %>% filter_all(all_vars(!is.na(.)))

##summarize dead cell data
DLs_summ<-DLs %>% group_by(well)%>%
  summarize(area_STING_total_well=mean(area_STING_total),IntDen_STING_total_well=mean(IntDen_STING_total),MeanGrey_STING_total_well=mean(MeanGrey_STING_total)) %>%
  ungroup()
DLs_summ2<-DLs %>% group_by(well,site) %>%
  summarize(area_STING_total_site=sum(area_STING_total),IntDen_STING_total_site=sum(IntDen_STING_total),MeanGrey_STING_total_site=sum(MeanGrey_STING_total)) %>%
  ungroup()


#periNuclear signal
#load data from intensity of STING in all cell count
data_STING_Peri_Nucl<-lapply(files_STINGtotalPeriNucl,function(x) {read.csv(x,header=T,stringsAsFactors=F)})
names(data_STING_Peri_Nucl)<-files_STINGtotalPeriNucl %>% str_split(patter="/",simplify=T) %>% .[,4] %>% str_split(pattern="\\.",simplify=T) %>% .[,1]
for(i in seq_along(data_STING_Peri_Nucl)) {
  data_STING_Peri_Nucl[[i]]$key<-rep(meta %>% filter(str_detect(file,"STING_Tuj_Peri_Nuclear_intensity")) %>% pull(key) %>% .[i],nrow(data_STING_Peri_Nucl[[i]])) }  # line of code to copy past
DLspn<-do.call(rbind,data_STING_Peri_Nucl)
colnames(DLspn)<-c("object_STING_Peri_Nucl","area_STING_Peri_Nucl","MeanGrey_STING_Peri_Nucl","StdDev_STING_Peri_Nucl","MinGrey_STING_Peri_Nucl",
                   "MaxGrey_STING_Peri_Nucl","Feret_STING_Peri_Nucl","IntDen_STING_Peri_Nucl","RawIntDen_STING_Peri_Nucl","FeretX_STING_Peri_Nucl","FeretY_STING_Peri_Nucl",
                   "Feretangle_STING_Peri_Nucl","MinFeret_STING_Peri_Nucl","MinThr_STING_Peri_Nucl","MaxThr_STING_Peri_Nucl","key")
DLspn<-cbind(DLspn,DLspn$key %>% str_split(pattern="_",simplify=T),stringsAsFactors=F)
colnames(DLspn)<-c("object_STING_Peri_Nucl","area_STING_Peri_Nucl","MeanGrey_STING_Peri_Nucl","StdDev_STING_Peri_Nucl","MinGrey_STING_Peri_Nucl",
                   "MaxGrey_STING_Peri_Nucl","Feret_STING_Peri_Nucl","IntDen_STING_Peri_Nucl","RawIntDen_STING_Peri_Nucl","FeretX_STING_Peri_Nucl","FeretY_STING_Peri_Nucl",
                   "Feretangle_STING_Peri_Nucl","MinFeret_STING_Peri_Nucl","MinThr_STING_Peri_Nucl","MaxThr_STING_Peri_Nucl","key","well","channel","site")
DLspn<-DLspn %>% select(key,object_STING_Peri_Nucl,area_STING_Peri_Nucl,MeanGrey_STING_Peri_Nucl,IntDen_STING_Peri_Nucl,RawIntDen_STING_Peri_Nucl,well,channel,site)
DLspn<-DLspn %>% filter_all(all_vars(!is.na(.)))

##summarize dead cell data
DLspn_summ<-DLspn %>% group_by(well)%>%
  summarize(PeriNucl_well=length((key)),area_STING_Peri_Nucl_well=mean(area_STING_Peri_Nucl),
            IntDen_STING_Peri_Nucl_well=mean(IntDen_STING_Peri_Nucl),
            MeanGrey_STING_Peri_Nucl_well=mean(MeanGrey_STING_Peri_Nucl)) %>%
  ungroup()

DLspn_summ2<-DLspn %>% group_by(well,site) %>%
  summarize(PeriNucl_site=length((key)),
            area_STING_Peri_Nucl_site=sum(area_STING_Peri_Nucl),IntDen_STING_Peri_Nucl_site=sum(IntDen_STING_Peri_Nucl),
            MeanGrey_STING_Peri_Nucl_site=sum(MeanGrey_STING_Peri_Nucl)) %>%
  ungroup()

#####particle STING
##############load data from STING in all FULL TUJ cells
data_STING_particle_total_TUJ<-lapply(files_STINGTujcellparticle,function(x) {read.csv(x,header=T,stringsAsFactors=F)})
names(data_STING_particle_total_TUJ)<-files_STINGTujcellparticle %>% str_split(patter="/",simplify=T) %>% .[,4] %>% str_split(pattern="\\.",simplify=T) %>% .[,1]
for(i in seq_along(data_STING_particle_total_TUJ)) {
  data_STING_particle_total_TUJ[[i]]$key<-rep(meta %>% filter(str_detect(file,"STING_Tuj_Total_particle")) %>% pull(key) %>% .[i],nrow(data_STING_particle_total_TUJ[[i]])) }  # line of code to copy past
DTTpartS<-do.call(rbind,data_STING_particle_total_TUJ)
colnames(DTTpartS)<-c("slice_STING_particle_total_TUJ","object_STING_particle_total_TUJ","area_total_STING_particle_total_TUJ","avg_Size_STING_particle_total_TUJ","Percent_area_STING_particle_total_TUJ",
                      "Mean_STING_particle_total_TUJ","Feret_STING_particle_total_TUJ","FeretX_particle_total_TUJ","FeretY_STING_particle_total_TUJ","FeretAngle_STING_particle_total_TUJ","MinFeret_STING_particle_total_TUJ",
                      "IntDen_STING_particle_total_TUJ","key")
DTTpartS<-cbind(DTTpartS,DTTpartS$key %>% str_split(pattern="_",simplify=T),stringsAsFactors=F)
colnames(DTTpartS)<-c("slice_STING_particle_total_TUJ","object_STING_particle_total_TUJ","area_total_STING_particle_total_TUJ","avg_Size_STING_particle_total_TUJ","Percent_area_STING_particle_total_TUJ",
                      "Mean_STING_particle_total_TUJ","Feret_STING_particle_total_TUJ","FeretX_particle_total_TUJ","FeretY_STING_particle_total_TUJ","FeretAngle_STING_particle_total_TUJ","MinFeret_STING_particle_total_TUJ",
                      "IntDen_STING_particle_total_TUJ","key","well","channel","site")
DTTpartS<-DTTpartS %>% select(key,object_STING_particle_total_TUJ,area_total_STING_particle_total_TUJ,avg_Size_STING_particle_total_TUJ,Mean_STING_particle_total_TUJ,IntDen_STING_particle_total_TUJ,well,channel,site)
DTTpartS<-DTTpartS %>% filter_all(all_vars(!is.na(.)))

##summarize STING particles cell data
DTTpartS_summ<-DTTpartS %>% group_by(well)%>%
  summarize(object_STING_particle_total_TUJ_well=mean(object_STING_particle_total_TUJ),area_total_STING_particle_total_TUJ_well=mean(area_total_STING_particle_total_TUJ),
            avg_Size_STING_particle_total_TUJ_well=mean(avg_Size_STING_particle_total_TUJ),IntDen_STING_particle_total_TUJ_well=mean(IntDen_STING_particle_total_TUJ)) %>%
  ungroup()
DTTpartS_summ2<-DTTpartS %>% group_by(well,site) %>%
  summarize(object_STING_particle_total_TUJ_site=(object_STING_particle_total_TUJ),area_total_STING_particle_total_TUJ_site=(area_total_STING_particle_total_TUJ),
            avg_Size_STING_particle_total_TUJ_site=(avg_Size_STING_particle_total_TUJ),IntDen_STING_particle_total_TUJ_site=(IntDen_STING_particle_total_TUJ)) %>%
  ungroup()

## build STING dataframe
DT_test<-DT
DT_test<-left_join(platemap, DT_test %>% unique(),by=c('well'))
DT_test<-DT_test %>% ungroup() %>% select(
  key,well,site,genotype,age,line,iso)

meta1<-DT_test

Data_summary_STING <- meta1

Data_summary_STING<-left_join(Data_summary_STING,dt_summ %>% unique(),by=c('well'))
Data_summary_STING<-left_join(Data_summary_STING,dt_summ2 %>% unique(),by=c('well','site'))

Data_summary_STING<-left_join(Data_summary_STING,DLs_summ %>% unique(),by=c('well'))
Data_summary_STING<-left_join(Data_summary_STING,DLs_summ2 %>% unique(),by=c('well','site')) 

Data_summary_STING<-left_join(Data_summary_STING,DLspn_summ %>% unique(),by=c('well'))
Data_summary_STING<-left_join(Data_summary_STING,DLspn_summ2 %>% unique(),by=c('well','site')) 

Data_summary_STING<-left_join(Data_summary_STING,DTTpartS_summ %>% unique(),by=c('well'))
Data_summary_STING<-left_join(Data_summary_STING,DTTpartS_summ2 %>% unique(),by=c('well','site'))


#STING data 
Data_summary_STING<-Data_summary_STING %>% ungroup() %>% select(
  key,well,site,genotype,age,mitotic,line,iso,
  tuj_well,sites_well,tuj_site,avg_area_tuj,
  
  area_STING_total_well,IntDen_STING_total_well,MeanGrey_STING_total_well,
  area_STING_total_site,IntDen_STING_total_site,MeanGrey_STING_total_site,
  
  PeriNucl_well,area_STING_Peri_Nucl_well,IntDen_STING_Peri_Nucl_well,MeanGrey_STING_Peri_Nucl_well,
  PeriNucl_site,area_STING_Peri_Nucl_site,IntDen_STING_Peri_Nucl_site,MeanGrey_STING_Peri_Nucl_site,
  
  object_STING_particle_total_TUJ_well,area_total_STING_particle_total_TUJ_well,avg_Size_STING_particle_total_TUJ_well,IntDen_STING_particle_total_TUJ_well,
  object_STING_particle_total_TUJ_site,area_total_STING_particle_total_TUJ_site,avg_Size_STING_particle_total_TUJ_site,IntDen_STING_particle_total_TUJ_site)
  
 
colnames(Data_summary_STING)<-c("key","well","site","genotype","age","mitotic","line","iso",
                                "tuj_well","sites_well","tuj_site","avg_area_tuj",
                                "area_STING_total_well","IntDen_STING_total_well","MeanGrey_STING_total_well",
                                "area_STING_total_site","IntDen_STING_total_site","MeanGrey_STING_total_site",
                                "PeriNucl_well","area_STING_Peri_Nucl_well","IntDen_STING_Peri_Nucl_well","MeanGrey_STING_Peri_Nucl_well",
                                "PeriNucl_site","area_STING_Peri_Nucl_site","IntDen_STING_Peri_Nucl_site","MeanGrey_STING_Peri_Nucl_site",
                                "object_STING_particle_total_TUJ_well","area_total_STING_particle_total_TUJ_well","avg_Size_STING_particle_total_TUJ_well","IntDen_STING_particle_total_TUJ_well",
                                "object_STING_particle_total_TUJ_site","area_total_STING_particle_total_TUJ_site","avg_Size_STING_particle_total_TUJ_site","IntDen_STING_particle_total_TUJ_site")
                               
Data_summary_STING[is.na(Data_summary_STING)]<-0

#STING
Data_summary_STING$line<-as.character(Data_summary_STING$line)
Data_summary_STING$line<-factor(Data_summary_STING$line,levels = unique(Data_summary_STING$line))

Data_summary_STING$age<-as.character(Data_summary_STING$age)
Data_summary_STING$age<-factor(Data_summary_STING$age,levels = unique(Data_summary_STING$age))

Data_summary_STING$mitotic<-as.character(Data_summary_STING$mitotic)
Data_summary_STING$mitotic<-factor(Data_summary_STING$mitotic,levels = unique(Data_summary_STING$mitotic))

Data_summary_STING$genotype<-as.character(Data_summary_STING$genotype)
Data_summary_STING$genotype<-factor(Data_summary_STING$genotype,levels = unique(Data_summary_STING$genotype))

#to order the treatment in the order we want
Data_summary_H2AX$line<-factor(Data_summary_H2AX$line,level=c("TDP43+/+","TDP43+/G298S"))
Data_summary_STING$line<-factor(Data_summary_STING$line,level=c("TDP43+/+","TDP43+/G298S"))

pdf("iPSC-SMN-NGN2-PCN-STING-plots.pdf")



####################################Cell bodies ######################################################

#total Hoescht

### QC
## raw data (all wells)
ggplot(Data_summary_STING %>% select(key, line,tuj_site,well,genotype) %>% unique(),
       aes(x = well, y=tuj_site, color = line)) +
  scale_color_manual(values=c("grey50","royalblue4")) + theme_bw() + 
  theme(strip.text.x = element_text(size = 8, colour = "black", angle = 0))+
  ggtitle("DAY 35 Raw Data: Hoechst+ cells \n(per site)") + 
  ylim(0,NA) + ylab("Hoechst+ cells \n(per site)") + xlab("") + theme(axis.title = element_text(size = 10)) + 
  stat_summary(fun.data=mean_cl_boot,geom="crossbar",fill="grey50",alpha=0.2) +
  # geom_boxplot() +
  #geom_point() +
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = 'bold')) +
  geom_jitter(height=0,width=0.1,size=2,alpha=0.7) +
  theme(axis.text.x=element_text(angle=30,hjust=1))

## remove wells with little survival, seeD from plottiDg 
Data_summary_STING1 <- Data_summary_STING
Data_summary_STING2 <- Data_summary_STING1 %>% filter(tuj_site > 2)

## raw data (removing wells with cells < 3)
ggplot(Data_summary_STING2 %>% select(key, line,tuj_site,well,genotype)  %>% unique(),
       aes(x = well, y=tuj_site, color = line)) +
  scale_color_manual(values=c("grey50","royalblue4")) + theme_bw() + 
  theme(strip.text.x = element_text(size = 8, colour = "black", angle = 0))+
  ggtitle("Raw Data (Removing wells with few cells): Hoechst+ cells \n(per site)") + 
  ylim(0,NA) + ylab("Hoechst+ cells \n(per site)") + xlab("") + theme(axis.title = element_text(size = 10)) + 
  stat_summary(fun.data=mean_cl_boot,geom="crossbar",fill="grey50",alpha=0.2) +
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = 'bold')) +
  geom_jitter(height=0,width=0.1,size=2,alpha=0.7) +
  theme(axis.text.x=element_text(angle=30,hjust=1))

hoeschtcount_datasum <- Data_summary_STING2 %>% group_by(well) %>% summarize(mean_tujcount = mean(tuj_site),
                                                                             n_site    = n(),
                                                                             sd_tujcount   = 1.5*sd(tuj_site))

hoeschtcount_datasum %>% summarize(mean_mean = mean(mean_tujcount),
                                   mean_sd = mean(sd_tujcount))
# 15.7+13.3 = 29 -> 30
# 15.7-13.3 = 2.4 -> 3

## plot range determined by mean number of cells per site per well +/- 1.5 standard deviation 
ggplot(Data_summary_STING2 %>% select(key, line,tuj_site,well,genotype) %>% filter(tuj_site > 2 & tuj_site < 31),
       aes(x = well, y=tuj_site, color = line)) +
  scale_color_manual(values=c("grey50","royalblue4")) + theme_bw() + 
  theme(strip.text.x = element_text(size = 8, colour = "black", angle = 0))+
  ggtitle("Cleaned Data: Hoechst+ cells \n(per site)") + 
  ylim(0,30) + ylab("Hoechst+ cells \n(per site)") + xlab("") + theme(axis.title = element_text(size = 10)) + 
  stat_summary(fun.data=mean_cl_boot,geom="crossbar",fill="grey50",alpha=0.2) +
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = 'bold')) +
  geom_jitter(height=0,width=0.1,size=2,alpha=0.7) +
  theme(axis.text.x=element_text(angle=30,hjust=1))

## assigned cleaned df with absolute range
Data_summary_STING_clean2 <- Data_summary_STING2 %>% filter(tuj_site > 2 & tuj_site < 31)

# check sample size or number of sites per well. 
Data_summary_STING_clean2 %>% group_by(well) %>% summarize(n = n())
#tuj_count_eq_3 <- Data_summary_STING_clean2 %>% filter(tuj_site == 3)

# plot cleaned data PER SITE
ggplot(Data_summary_STING_clean2 %>% select(key,line,tuj_site,well,genotype) %>% unique(),
       aes(x=genotype,y=(tuj_site),color=line)) + 
  scale_color_manual(values=c("grey50","royalblue4")) + theme_bw() + 
  theme(strip.text.x = element_text(size = 8, colour = "black", angle = 0))+
  ggtitle("Cleaned Data: Hoescht_live+ cells \n(per site)") + 
  ylim(0,NA) + ylab("Hoescht_live+ cells \n(per site)") + xlab("") + theme(axis.title = element_text(size = 10)) + 
  stat_summary(fun.data=mean_cl_boot,geom="crossbar",fill="grey50",alpha=0.2) +
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = 'bold')) +
  geom_jitter(height=0,width=0.1,size=2,alpha=0.7) +
  theme(axis.text.x=element_text(angle=30,hjust=1)) 

#plot cleaned data PER WELL
ggplot(Data_summary_STING_clean2 %>% group_by(genotype,line,well) %>% summarize(
  TUJ_well=mean(tuj_site)),
  aes(x=genotype,y=(TUJ_well),color=line)) + theme_bw() + 
  theme(strip.text.x = element_text(size = 8, colour = "black", angle = 0))+
  ggtitle("Cleaned Data: Hoechst+ cells \n(per site per well)") + 
  ylim(0,NA) + ylab("Hoechst+ cells \n(per site per well)") + xlab("") + theme(axis.title = element_text(size = 10)) + 
  stat_summary(fun.data=mean_cl_boot,geom="crossbar",fill="grey50",alpha=0.2) +
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = 'bold')) +
  geom_jitter(height=0,width=0.1,size=2,alpha=0.7) +
  theme(axis.text.x=element_text(angle=30,hjust=1)) +
  scale_color_manual(values=c("grey50","royalblue4"))


####### STING
# integrated intensity 
#total

### QC
## check for outliers
ggplot(Data_summary_STING_clean2 %>% select(IntDen_STING_total_site,MeanGrey_STING_total_site,area_STING_total_site,well,
                                            line,genotype,key) %>% unique(),
       aes(x=well,y=(IntDen_STING_total_site),color=line)) + 
  scale_color_manual(values=c("grey50","royalblue4")) + theme_bw() + 
  # facet_grid(~iso) + 
  theme(strip.text.x = element_text(size = 8, colour = "black", angle = 0))+
  ggtitle("Integrated Intensity of total STING\n(site ) \n(remove outliers >6.5e5)") + 
  ylim(0,NA) + ylab("Integrated Intensity of total STING\n(per site )") + xlab("") + theme(axis.title = element_text(size = 10)) + 
  stat_summary(fun.data=mean_sdl,fun.args=list(mult=2),geom="crossbar",fill="grey50",alpha=0.2) +
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = 'bold')) +
  geom_jitter(height=0,width=0.1,size=2,alpha=0.7) +
  theme(axis.text.x=element_text(angle=30,hjust=1)) 

# identify and remove outliers outside of 2SD
Data_summary_STING_clean_int2 <- Data_summary_STING_clean2 %>% group_by(well) %>% filter(
  IntDen_STING_total_site <= mean(IntDen_STING_total_site) + 2*sd(IntDen_STING_total_site))  %>% ungroup()

## check number of sites/well to ensure outliers removed
n_sting <- Data_summary_STING_clean2 %>% group_by(well) %>% summarize(n = n())
n_sting_clean <- Data_summary_STING_clean_int2 %>% group_by(well) %>% summarize(n = n())


#PER SITE
#total STING 
## strict absolute range 
ggplot(Data_summary_STING_clean_int2 %>% select(IntDen_STING_total_site,MeanGrey_STING_total_site,area_STING_total_site,
                                            line,genotype,key) %>% unique(),
       aes(x=genotype,y=(IntDen_STING_total_site),color=line)) + 
  scale_color_manual(values=c("grey50","royalblue4")) + theme_bw() + 
  # facet_grid(~iso) + 
  theme(strip.text.x = element_text(size = 8, colour = "black", angle = 0))+
  ggtitle("Clean, Removing Outliers: \n Integrated Intensity of total STING\n(site )") + 
  ylim(0,NA) + ylab("Integrated Intensity of total STING\n(per site )") + xlab("") + theme(axis.title = element_text(size = 10)) + 
  stat_summary(fun.data=mean_cl_boot,geom="crossbar",fill="grey50",alpha=0.2) +
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = 'bold')) +
  geom_jitter(height=0,width=0.1,size=2,alpha=0.7) +
  theme(axis.text.x=element_text(angle=30,hjust=1)) 

#PER WELL
#total STING 
# strict abs range (removing large outlier)
ggplot(Data_summary_STING_clean_int2 %>% group_by(genotype,line,well) %>% summarize(
  IntDenNucl_well=mean(IntDen_STING_total_site)),
  aes(x=genotype,y=(IntDenNucl_well),color=line)) + theme_bw() + 
  #facet_grid(~well) + 
  theme(strip.text.x = element_text(size = 8, colour = "black", angle = 0))+
  ggtitle("Clean, Removing Outliers: \n Integrated Intensity of total STING\n(well )") + 
  ylim(0,NA) + ylab("Integrated Intensity of total STING\n(well )") + xlab("") + theme(axis.title = element_text(size = 10)) + 
  stat_summary(fun.data=mean_cl_boot,geom="crossbar",fill="grey50",alpha=0.2) +
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = 'bold')) +
  geom_jitter(height=0,width=0.1,size=2,alpha=0.7) +
  theme(axis.text.x=element_text(angle=30,hjust=1)) +
  scale_color_manual(values=c("grey50","royalblue4")) 


#########PER CELL
#STING total
# strict abs range 
ggplot(Data_summary_STING_clean_int2 %>% group_by(genotype,line,well) %>% summarize(
  IntDenNucl_well=mean(IntDen_STING_total_site),TUJ_well=mean(tuj_site)),
  aes(x=genotype,y=((IntDenNucl_well)/(TUJ_well)),color=line)) + theme_bw() + 
  # facet_grid(~well) + 
  theme(strip.text.x = element_text(size = 8, colour = "black", angle = 0))+
  ggtitle("Clean, Removing Outliers: \n Integrated Intensity of total STING per cell\n(well )") + 
  ylim(0,20000) + ylab("Integrated Intensity of total STING per cell\n(well )") + xlab("") + theme(axis.title = element_text(size = 10)) + 
  stat_summary(fun.data = mean_se, geom = "crossbar", fill="grey50",alpha=0.2) +
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = 'bold')) +
  geom_jitter(height=0,width=0.1,size=2,alpha=0.7) +
  theme(axis.text.x=element_text(angle=30,hjust=1)) +
  scale_color_manual(values=c("grey50","royalblue4")) 


####### STING
# Area
#total

### QC
## check for outliers
ggplot(Data_summary_STING_clean2 %>% select(IntDen_STING_total_site,MeanGrey_STING_total_site,area_STING_total_site,well,
                                            line,genotype,key) %>% unique(),
       aes(x=well,y=(area_STING_total_site),color=line)) + 
  scale_color_manual(values=c("grey50","royalblue4")) + theme_bw() + 
  # facet_grid(~iso) + 
  theme(strip.text.x = element_text(size = 8, colour = "black", angle = 0))+
  ggtitle("Integrated Intensity of total STING\n(site ) \n(remove outliers >6.5e5)") + 
  ylim(0,NA) + ylab("Integrated Intensity of total STING\n(per site )") + xlab("") + theme(axis.title = element_text(size = 10)) + 
  stat_summary(fun.data=mean_sdl,fun.args=list(mult=2),geom="crossbar",fill="grey50",alpha=0.2) +
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = 'bold')) +
  geom_jitter(height=0,width=0.1,size=2,alpha=0.7) +
  theme(axis.text.x=element_text(angle=30,hjust=1)) 

# identify and remove outliers outside of 2SD
Data_summary_STING_clean_area2 <- Data_summary_STING_clean2 %>% group_by(well) %>% filter(
  area_STING_total_site <= mean(area_STING_total_site) + 2*sd(area_STING_total_site))  %>% ungroup()


#########PER CELL
#STING total
# strict abs range 
ggplot(Data_summary_STING_clean_area2 %>% group_by(genotype,line,well) %>% summarize(
  AreaNucl_well=mean(area_STING_total_site),TUJ_well=mean(tuj_site)),
  aes(x=genotype,y=((AreaNucl_well)/(TUJ_well)),color=line)) + theme_bw() + 
  # facet_grid(~well) + 
  theme(strip.text.x = element_text(size = 8, colour = "black", angle = 0))+
  ggtitle("Clean, Removing Outliers: \n Area of total STING per cell\n(well )") + 
  ylim(0,NA) + ylab("Area of total STING per cell\n(well )") + xlab("") + theme(axis.title = element_text(size = 10)) + 
  stat_summary(fun.data = mean_se, geom = "crossbar", fill="grey50",alpha=0.2) +
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = 'bold')) +
  geom_jitter(height=0,width=0.1,size=2,alpha=0.7) +
  theme(axis.text.x=element_text(angle=30,hjust=1)) +
  scale_color_manual(values=c("grey50","royalblue4"))


#### STING POSITIVE CELLS / PERINUCLEAR STING ANALYSIS 

##analysis of STING positive cells 
meta2<-meta1
meta2<-left_join(meta2, DLspn %>% unique(), by=c('key'))
meta22 <- left_join(meta2, DT %>% unique(), by="key") 

# filter out sites not within cell count range (to be consistent with previous data cleaning)
meta23 <- meta22 %>% filter(object_tuj > 2 & object_tuj < 31)

meta23<-meta23%>% ungroup() %>% select(
  key,well.x,site.x,genotype,age,mitotic,line,iso,object_STING_Peri_Nucl,area_STING_Peri_Nucl,
  MeanGrey_STING_Peri_Nucl,IntDen_STING_Peri_Nucl,RawIntDen_STING_Peri_Nucl)

colnames(meta23) <- c("key","well","site","genotype","age","mitotic","line","iso","object_Nuclei","area_STING_Peri_Nucl","
                     MeanGrey_STING_Peri_Nucl","IntDen_STING_Peri_Nucl","RawIntDen_STING_Peri_Nucl")

# plot integrated intensity of perinuclear sting vs. log(area of peri sting), set xlim(-4,5)
meta23 %>% ggplot() + 
  aes(x = log(area_STING_Peri_Nucl), y = IntDen_STING_Peri_Nucl, color = line) + 
  scale_color_manual(values=c("grey50","royalblue4")) + theme_bw() + 
  # facet_grid(~iso) + 
  theme(strip.text.x = element_text(size = 8, colour = "black", angle = 0))+
  # ggtitle("") + 
  ylim(0,NA) + xlim(-4,7) +
  ylab("Integrated Intensity of perinuclear STING") + xlab("log(Area of perinuclear STING)") + theme(axis.title = element_text(size = 10)) + 
  #stat_summary(fun.data=mean_cl_boot,geom="crossbar",fill="grey50",alpha=0.2) +
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = 'bold')) +
  geom_jitter(height=0,width=0.1,size=2,alpha=0.7) +
  theme(axis.text.x=element_text(angle=30,hjust=1)) 
  
## set threshold above which perinuclear sting intensity is analyzed
meta3 <- meta23 %>% filter(log(area_STING_Peri_Nucl)>=2.5) 

##summarize dead cell data
meta3_summ<-meta3 %>% group_by(well)%>%
  summarize(PeriNucl_well=length((key)),area_STING_Peri_Nucl_well=mean(area_STING_Peri_Nucl),
            IntDen_STING_Peri_Nucl_well=mean(IntDen_STING_Peri_Nucl)) %>%
  ungroup()

meta3_summ2<-meta3 %>% group_by(well,site,key) %>%
  summarize(PeriNucl_site=length((key)),
            area_STING_Peri_Nucl_site=sum(area_STING_Peri_Nucl),IntDen_STING_Peri_Nucl_site=sum(IntDen_STING_Peri_Nucl)) %>%
  ungroup()

meta5<-left_join(meta3_summ2,meta1 %>% unique(),by=c('key'))
#replace NA by 0
meta5[is.na(meta5)]<-0

meta5<-meta5 %>% ungroup() %>% select(
  key,well.x,site.x,genotype,age,mitotic,line,iso,PeriNucl_site,area_STING_Peri_Nucl_site,
  IntDen_STING_Peri_Nucl_site)
colnames(meta5) <- c("key","well","site","genotype","age","mitotic","line","iso","PeriNucl_site","area_STING_Peri_Nucl_site",
                     "IntDen_STING_Peri_Nucl_site")

meta5<-left_join(meta5,dt_summ %>% unique(),by=c('well'))
meta5<-left_join(meta5,dt_summ2 %>% unique(),by=c('well','site'))

meta5$line<-factor(meta5$line,level=c("TDP43+/+","TDP43+/G298S"))


## total number of STING+ cells
## check for outliers
ggplot(meta5 %>% select(PeriNucl_site,IntDen_STING_Peri_Nucl_site,well,
                        line,genotype,key) %>% unique(),
       aes(x=well,y=(PeriNucl_site),color=line)) + 
  scale_color_manual(values=c("grey50","royalblue4")) + theme_bw() + 
  # facet_grid(~iso) + 
  theme(strip.text.x = element_text(size = 8, colour = "black", angle = 0))+
  ggtitle("STING positive cells\n(site )") + 
  ylim(0,NA) + ylab("STING positive cells\n(per site ) \n check for outliers") + xlab("") + theme(axis.title = element_text(size = 10)) + 
  stat_summary(fun.data=mean_sdl,geom="crossbar",fill="grey50",alpha=0.2) +
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = 'bold')) +
  geom_jitter(height=0,width=0.1,size=2,alpha=0.7) +
  theme(axis.text.x=element_text(angle=30,hjust=1)) 

###PER SITE
ggplot(meta5 %>% select(PeriNucl_site,IntDen_STING_Peri_Nucl_site,line,genotype,key) %>% unique(),
       aes(x=genotype,y=(PeriNucl_site),color=line)) + 
  scale_color_manual(values=c("grey50","royalblue4")) + theme_bw() + 
  # facet_grid(~iso) + 
  theme(strip.text.x = element_text(size = 8, colour = "black", angle = 0))+
  ggtitle("STING positive cells\n(site )") + 
  ylim(0,NA) + ylab("STING positive cells\n(per site )") + xlab("") + theme(axis.title = element_text(size = 10)) + 
  stat_summary(fun.data=mean_cl_boot,geom="crossbar",fill="grey50",alpha=0.2) +
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = 'bold')) +
  geom_jitter(height=0,width=0.1,size=2,alpha=0.7) +
  theme(axis.text.x=element_text(angle=30,hjust=1)) 


#PER WELL
#total STING 
ggplot(meta5 %>% group_by(genotype,line,well) %>% summarize(
  PeriNucl_well1=sum(PeriNucl_site)),
  aes(x=line,y=(PeriNucl_well1),color=line)) + theme_bw() + 
  #facet_grid(~well) + 
  theme(strip.text.x = element_text(size = 8, colour = "black", angle = 0))+
  ggtitle("STING positive cells\n(well )") + 
  ylim(0,NA) + ylab("STING positive cells\n(well )") + xlab("") + theme(axis.title = element_text(size = 10)) + 
  stat_summary(fun.data=mean_cl_boot,geom="crossbar",fill="grey50",alpha=0.2) +
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = 'bold')) +
  geom_jitter(height=0,width=0.1,size=2,alpha=0.7) +
  theme(axis.text.x=element_text(angle=30,hjust=1)) +
  scale_color_manual(values=c("grey50","royalblue4")) 


#########PER CELL
#STING total 
ggplot(meta5 %>% group_by(genotype,line,well) %>% summarize(
  PeriNucl_well1=sum(PeriNucl_site),TUJ_well=sum(tuj_site)),
  aes(x=line,y=((PeriNucl_well1)/(TUJ_well))*100,color=line)) + theme_bw() + 
  #facet_grid(~well) + 
  theme(strip.text.x = element_text(size = 8, colour = "black", angle = 0))+
  ggtitle("Percentage of STING positive cells\n(well)") + 
  ylim(0,NA) + ylab("Percentage of STING positive cells\n(well)") + xlab("") + theme(axis.title = element_text(size = 10)) + 
  stat_summary(fun.data = mean_se, geom = "crossbar", fill="grey50",alpha=0.2) +
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = 'bold')) +
  geom_jitter(height=0,width=0.1,size=2,alpha=0.7) +
  theme(axis.text.x=element_text(angle=30,hjust=1)) +
  scale_color_manual(values=c("grey50","royalblue4")) 



### integrated intensity- use dataset without removing areas under a certain threshold

meta23_summ<-meta23 %>% group_by(well)%>%
  summarize(PeriNucl_well=length((key)),area_STING_Peri_Nucl_well=mean(area_STING_Peri_Nucl),
            IntDen_STING_Peri_Nucl_well=mean(IntDen_STING_Peri_Nucl)) %>%
  ungroup()

meta23_summ2<-meta23 %>% group_by(well,site,key) %>%
  summarize(PeriNucl_site=length((key)),
            area_STING_Peri_Nucl_site=sum(area_STING_Peri_Nucl),IntDen_STING_Peri_Nucl_site=sum(IntDen_STING_Peri_Nucl)) %>%
  ungroup()

meta24<-left_join(meta23_summ2,meta1 %>% unique(),by=c('key'))
#replace NA by 0
meta24[is.na(meta24)]<-0

meta24<-meta24 %>% ungroup() %>% select(
  key,well.x,site.x,genotype,age,mitotic,line,iso,PeriNucl_site,area_STING_Peri_Nucl_site,
  IntDen_STING_Peri_Nucl_site)
colnames(meta24) <- c("key","well","site","genotype","age","mitotic","line","iso","PeriNucl_site","area_STING_Peri_Nucl_site",
                     "IntDen_STING_Peri_Nucl_site")

meta24<-left_join(meta24,dt_summ %>% unique(),by=c('well'))
meta24<-left_join(meta24,dt_summ2 %>% unique(),by=c('well','site'))

meta24$line<-factor(meta24$line,level=c("TDP43+/+","TDP43+/G298S"))


## check for outliers
ggplot(meta24 %>% select(PeriNucl_site,IntDen_STING_Peri_Nucl_site,well,
                        line,genotype,key) %>% unique(),
       aes(x=well,y=(IntDen_STING_Peri_Nucl_site),color=line)) + 
  scale_color_manual(values=c("grey50","royalblue4")) + theme_bw() + 
  # facet_grid(~iso) + 
  theme(strip.text.x = element_text(size = 8, colour = "black", angle = 0))+
  ggtitle("Integrated Intensity of cytoplasmic STING\n(site )") + 
  ylim(0,NA) + ylab("Integrated Intensity of cytoplasmic STING\n(per site )") + xlab("") + theme(axis.title = element_text(size = 10)) + 
  stat_summary(fun.data=mean_sdl,geom="crossbar",fill="grey50",alpha=0.2) +
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = 'bold')) +
  geom_jitter(height=0,width=0.1,size=2,alpha=0.7) +
  theme(axis.text.x=element_text(angle=30,hjust=1)) 

# identify and remove outliers
meta25 <- meta24 %>% group_by(well) %>% 
  filter(IntDen_STING_Peri_Nucl_site <= mean(IntDen_STING_Peri_Nucl_site) + 2*sd(IntDen_STING_Peri_Nucl_site)) %>% ungroup()


#PER SITE
#total STING 
ggplot(meta25 %>% select(PeriNucl_site,IntDen_STING_Peri_Nucl_site,
                        line,genotype,key) %>% unique(),
       aes(x=line,y=(IntDen_STING_Peri_Nucl_site),color=line)) + 
  scale_color_manual(values=c("grey50","royalblue4")) + theme_bw() + 
  # facet_grid(~iso) + 
  theme(strip.text.x = element_text(size = 8, colour = "black", angle = 0))+
  ggtitle("Integrated Intensity of cytoplasmic STING\n(site )") + 
  ylim(0,NA) + ylab("Integrated Intensity of cytoplasmic STING\n(per site )") + xlab("") + theme(axis.title = element_text(size = 10)) + 
  stat_summary(fun.data=mean_cl_boot,geom="crossbar",fill="grey50",alpha=0.2) +
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = 'bold')) +
  geom_jitter(height=0,width=0.1,size=2,alpha=0.7) +
  theme(axis.text.x=element_text(angle=30,hjust=1)) 

#PER WELL
#total STING 
ggplot(meta25 %>% group_by(genotype,line,well) %>% summarize(
  IntDen_STING_Peri_Nucl_well=mean(IntDen_STING_Peri_Nucl_site)),
  aes(x=line,y=(IntDen_STING_Peri_Nucl_well),color=line)) + theme_bw() + 
  #facet_grid(~well) + 
  theme(strip.text.x = element_text(size = 8, colour = "black", angle = 0))+
  ggtitle("Integrated Intensity of perinuclear STING\n(well)") + 
  ylim(0,NA) + ylab("Integrated Intensity of perinuclear STING\n(well)") + xlab("") + theme(axis.title = element_text(size = 10)) + 
  stat_summary(fun.data=mean_cl_boot,geom="crossbar",fill="grey50",alpha=0.2) +
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = 'bold')) +
  geom_jitter(height=0,width=0.1,size=2,alpha=0.7) +
  theme(axis.text.x=element_text(angle=30,hjust=1)) +
  scale_color_manual(values=c("grey50","royalblue4")) 


#########PER CELL
#integrated intensity of perinuclear STING (set upper bound to 6000)
ggplot(meta25 %>% group_by(genotype,line,well) %>% summarize(
  IntDenNucl_well=mean(IntDen_STING_Peri_Nucl_site),TUJ_well=mean(tuj_site)),
  aes(x=line,y=((IntDenNucl_well)/(TUJ_well)),color=line)) + theme_bw() + 
  #facet_grid(~well) + 
  theme(strip.text.x = element_text(size = 8, colour = "black", angle = 0))+
  ggtitle("Integrated Intensity of perinuclear STING per cell\n(well)") + 
  ylim(0,NA) + ylab("Integrated Intensity of perinuclear STING per cell\n(well)") + xlab("") + theme(axis.title = element_text(size = 10)) + 
  stat_summary(fun.data = mean_se, geom = "crossbar", fill="grey50",alpha=0.2) +
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = 'bold')) +
  geom_jitter(height=0,width=0.1,size=2,alpha=0.7) +
  theme(axis.text.x=element_text(angle=30,hjust=1)) +
  scale_color_manual(values=c("grey50","royalblue4")) 

## check for outliers (Area)
ggplot(meta24 %>% select(PeriNucl_site,area_STING_Peri_Nucl_site,well,
                         line,genotype,key) %>% unique(),
       aes(x=well,y=(area_STING_Peri_Nucl_site),color=line)) + 
  scale_color_manual(values=c("grey50","royalblue4")) + theme_bw() + 
  # facet_grid(~iso) + 
  theme(strip.text.x = element_text(size = 8, colour = "black", angle = 0))+
  ggtitle("Integrated Intensity of cytoplasmic STING\n(site)") + 
  ylim(0,NA) + ylab("Integrated Intensity of cytoplasmic STING\n(per site)") + xlab("") + theme(axis.title = element_text(size = 10)) + 
  stat_summary(fun.data=mean_sdl,geom="crossbar",fill="grey50",alpha=0.2) +
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = 'bold')) +
  geom_jitter(height=0,width=0.1,size=2,alpha=0.7) +
  theme(axis.text.x=element_text(angle=30,hjust=1)) 

# identify and remove outliers
meta26 <- meta24 %>% group_by(well) %>% 
  filter(area_STING_Peri_Nucl_site <= mean(area_STING_Peri_Nucl_site) + 2*sd(area_STING_Peri_Nucl_site)) %>% ungroup()


yl <- expression("Area of perinuclear STING per cell (well, um"^"2"~")")

#########PER CELL
#Area perinuclear sting (set upper bound to 4)
ggplot(meta26 %>% group_by(genotype,line,well) %>% summarize(
  areaNucl_well=sum(area_STING_Peri_Nucl_site),TUJ_well=sum(tuj_site)),
  aes(x=line,y=((areaNucl_well)/(TUJ_well)),color=line)) + theme_bw() + 
  #facet_grid(~well) + 
  theme(strip.text.x = element_text(size = 8, colour = "black", angle = 0))+
  ggtitle(yl) + 
  ylim(0,NA) + ylab(yl) + xlab("") + theme(axis.title = element_text(size = 10)) + 
  stat_summary(fun.data = mean_se, geom = "crossbar", fill="grey50",alpha=0.2) +
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = 'bold')) +
  geom_jitter(height=0,width=0.1,size=2,alpha=0.7) +
  theme(axis.text.x=element_text(angle=30,hjust=1)) +
  scale_color_manual(values=c("grey50","royalblue4")) 


dev.off()

Data_summary_STING_clean_int2 <- Data_summary_STING_clean_int2 %>% mutate(replicate = "X")
Data_summary_STING_clean_area2 <- Data_summary_STING_clean_area2 %>% mutate(replicate = "X")
meta5 <- meta5 %>% mutate(replicate = "X")
meta25 <- meta25 %>% mutate(replicate = "X")
meta26<- meta26 %>% mutate(replicate = "X")

write.csv(Data_summary_STING_clean_int2, "iPSC-example-STING-intpercell.csv")
write.csv(Data_summary_STING_clean_area2, "iPSC-example-STING-areapercell.csv")
write.csv(meta5, "iPSC-example-STING-pctposcell.csv")
write.csv(meta25, "iPSC-example-STING-peri-intpercell.csv")
write.csv(meta26, "iPSC-example-STING-peri-areapercell.csv")



