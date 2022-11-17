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
setwd("/Data-Experiments/Analysis-R/pIRF3/Data") 


#create platemap dataframe with sample information
platemap<-data.frame(well=c('B08','B09','B19','B11','B12','B13','B14','B15','B16','B17',
                            'C08','C09','C19','C11','C12','C13','C14','C15','C16','C17'),
                     
                     genotype=c('wt','wt','wt','wt','wt','wt','wt','wt','wt','wt',
                                'wt','wt','wt','wt','wt','wt','wt','wt','wt','wt'),
                     
                     line=c('PCN','PCN','PCN','PCN','PCN','PCN','PCN','PCN','PCN','PCN',
                            'PCN','PCN','PCN','PCN','PCN','PCN','PCN','PCN','PCN','PCN'),
                     
                     treatment=c('DMSO','DMSO','DMSO','DMSO','DMSO','DMSO','DMSO','DMSO','DMSO','DMSO',
                                 'DMSO','DMSO','DMSO','DMSO','DMSO','DMSO','DMSO','DMSO','DMSO','DMSO',
                                 'etoposide','etoposide','etoposide','etoposide','etoposide','etoposide','etoposide','etoposide','etoposide','etoposide',
                                 'etoposide','etoposide','etoposide','etoposide','etoposide','etoposide','etoposide','etoposide','etoposide','etoposide'),
                     
                     duration=c('1H','1H','1H','1H','1H','1H','1H','1H','1H','1H',
                                '1H','1H','1H','1H','1H','1H','1H','1H','1H','1H'),stringsAsFactors=F)



#identify files and pull metadata from file names
files<-list.files(pattern=".csv",recursive=T,full.names=F)

#cells Live
files_live<-files[which((files %>% str_detect("hoescht_count")))]
files_ctip2_pos<-files[which((files %>% str_detect("ctip2_POS_nucl_count")))]
files_ctip2_neg<-files[which((files %>% str_detect("ctip2_NEG_nucl_count")))]

#intensity
files_pIRF3_ctip2_neg_cyto_intensity<-files[which((files %>% str_detect("pirf3_ctip2_NEGbody_intensity_count")))]
files_pIRF3_ctip2_neg_nucl_intensity<-files[which((files %>% str_detect("pirf3_ctip2_NEGnucl_intensity_count")))]

files_pIRF3_ctip2_pos_cyto_intensity<-files[which((files %>% str_detect("pirf3_ctip2_POSbody_intensity_count")))]
files_pIRF3_ctip2_pos_nucl_intensity<-files[which((files %>% str_detect("pirf3_ctip2_POSnucl_intensity_count")))]

# pull metadata 
files<-files[which((files %>% str_detect("_count")))]
meta<-files %>% str_split(pattern="/",simplify=T) %>% as.data.frame(stringsAsFactors=F)
meta<-cbind(meta,meta$V4 %>% str_sub(1,10),stringsAsFactors=F) #changed from V4 to V2 

colnames(meta)<-c('sample','data','folder','subfold','key')
meta<-meta %>% select(sample,data,key)
meta$file<-files %>% str_sub(end=-5L)
meta<-cbind(meta,meta$key %>% str_split(pattern="_",simplify=T) %>% as.data.frame(stringsAsFactors=F))
colnames(meta)<-c('sample','data','key','file','well','channel','site')
meta<-meta %>% select(file,sample,key,well,site)
meta<-left_join(meta,platemap %>% unique(),by="well")
meta$key2<-str_c(meta$well,meta$site,sep="_")


# load data from hoescht count (all nuclei)
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

#load data from cell bodies count - ctip2 positive
data_ctip2_pos<-lapply(files_ctip2_pos,function(x) {read.csv(x,header=T,stringsAsFactors=F)})
names(data_ctip2_pos)<-files_ctip2_pos %>% str_split(patter="/",simplify=T) %>% .[,4] %>% str_split(pattern="\\.",simplify=T) %>% .[,1]
for(i in seq_along(data_ctip2_pos)) {
  data_ctip2_pos[[i]]$key<-rep(meta$key[i],nrow(data_ctip2_pos[[i]])) }
DCP<-do.call(rbind,data_ctip2_pos)
colnames(DCP)<-c("name","object_ctip2_pos","area_ctip2_pos","Average_ctip2_pos_Size","Percentage_ctip2_pos_Area","key")
DCP<-cbind(DCP,DCP$key %>% str_split(pattern="_",simplify=T),stringsAsFactors=F)
colnames(DCP)<-c("name","object_ctip2_pos","area_ctip2_pos","Average_ctip2_pos_Size","Percentage_ctip2_pos_Area","key","well","channel","site")
DCP<-DCP %>% select(key,object_ctip2_pos,area_ctip2_pos,Average_ctip2_pos_Size,Percentage_ctip2_pos_Area,well,channel,site)

##summarize ctip2+ cell count data
dcp_summ<-DCP %>% group_by(well)%>%
  summarize(ctip2_pos_well=sum(object_ctip2_pos),sites_well=length(unique(key))) %>%
  ungroup()
dcp_summ2<-DCP %>% group_by(well,site) %>%
  summarize(ctip2_pos_site=object_ctip2_pos,avg_area_ctip2_pos=Average_ctip2_pos_Size) %>%
  ungroup()


#load data from cell bodies count - ctip2 negative
data_ctip2_neg<-lapply(files_ctip2_neg,function(x) {read.csv(x,header=T,stringsAsFactors=F)})
names(data_ctip2_neg)<-files_ctip2_neg %>% str_split(patter="/",simplify=T) %>% .[,4] %>% str_split(pattern="\\.",simplify=T) %>% .[,1]
for(i in seq_along(data_ctip2_neg)) {
  data_ctip2_neg[[i]]$key<-rep(meta$key[i],nrow(data_ctip2_neg[[i]])) }
DCN<-do.call(rbind,data_ctip2_neg)
colnames(DCN)<-c("name","object_ctip2_neg","area_ctip2_neg","Average_ctip2_neg_Size","Percentage_ctip2_neg_Area","key")
DCN<-cbind(DCN,DCN$key %>% str_split(pattern="_",simplify=T),stringsAsFactors=F)
colnames(DCN)<-c("name","object_ctip2_neg","area_ctip2_neg","Average_ctip2_neg_Size","Percentage_ctip2_neg_Area","key","well","channel","site")
DCN<-DCN %>% select(key,object_ctip2_neg,area_ctip2_neg,Average_ctip2_neg_Size,Percentage_ctip2_neg_Area,well,channel,site)

##summarize ctip2- cell count data
dcn_summ<-DCN %>% group_by(well)%>%
  summarize(ctip2_neg_well=sum(object_ctip2_neg),sites_well=length(unique(key))) %>%
  ungroup()
dcn_summ2<-DCN %>% group_by(well,site) %>%
  summarize(ctip2_neg_site=object_ctip2_neg,avg_area_ctip2_neg=Average_ctip2_neg_Size) %>%
  ungroup()


#######CTIP2 POSITIVE

#CTIP2 POSITIVE cytoplasm
#load intensity data from pIRF3 in ctip2+ cell count (cytobody)
data_pIRF3_ctip2_pos_cyto_intensity<-lapply(files_pIRF3_ctip2_pos_cyto_intensity,function(x) {read.csv(x,header=T,stringsAsFactors=F)})
names(data_pIRF3_ctip2_pos_cyto_intensity)<-files_pIRF3_ctip2_pos_cyto_intensity %>% str_split(patter="/",simplify=T) %>% .[,4] %>% str_split(pattern="\\.",simplify=T) %>% .[,1]
for(i in seq_along(data_pIRF3_ctip2_pos_cyto_intensity)) {
  data_pIRF3_ctip2_pos_cyto_intensity[[i]]$key<-rep(meta$key[i],nrow(data_pIRF3_ctip2_pos_cyto_intensity[[i]])) }
DCPIcyto<-do.call(rbind,data_pIRF3_ctip2_pos_cyto_intensity)
colnames(DCPIcyto)<-c("object_pIRF3_ctip2_pos_cyto_intensity","area_pIRF3_ctip2_pos_cyto_intensity","MeanGrey_pIRF3_ctip2_pos_cyto_intensity","StdDev_pIRF3_ctip2_pos_cyto_intensity","MinGrey_pIRF3_ctip2_pos_cyto_intensity",
                      "MaxGrey_pIRF3_ctip2_pos_cyto_intensity","Feret_pIRF3_ctip2_pos_cyto_intensity","IntDen_pIRF3_ctip2_pos_cyto_intensity","RawIntDen_pIRF3_ctip2_pos_cyto_intensity","FeretX_pIRF3_ctip2_pos_cyto_intensity","FeretY_pIRF3_ctip2_pos_cyto_intensity",
                      "Feretangle_pIRF3_ctip2_pos_cyto_intensity","MinFeret_pIRF3_ctip2_pos_cyto_intensity","MinThr_pIRF3_ctip2_pos_cyto_intensity","MaxThr_pIRF3_ctip2_pos_cyto_intensity","key")
DCPIcyto<-cbind(DCPIcyto,DCPIcyto$key %>% str_split(pattern="_",simplify=T),stringsAsFactors=F)
colnames(DCPIcyto)<-c("object_pIRF3_ctip2_pos_cyto_intensity","area_pIRF3_ctip2_pos_cyto_intensity","MeanGrey_pIRF3_ctip2_pos_cyto_intensity","StdDev_pIRF3_ctip2_pos_cyto_intensity","MinGrey_pIRF3_ctip2_pos_cyto_intensity",
                      "MaxGrey_pIRF3_ctip2_pos_cyto_intensity","Feret_pIRF3_ctip2_pos_cyto_intensity","IntDen_pIRF3_ctip2_pos_cyto_intensity","RawIntDen_pIRF3_ctip2_pos_cyto_intensity","FeretX_pIRF3_ctip2_pos_cyto_intensity","FeretY_pIRF3_ctip2_pos_cyto_intensity",
                      "Feretangle_pIRF3_ctip2_pos_cyto_intensity","MinFeret_pIRF3_ctip2_pos_cyto_intensity","MinThr_pIRF3_ctip2_pos_cyto_intensity","MaxThr_pIRF3_ctip2_pos_cyto_intensity","key","well","channel","site")
DCPIcyto<-DCPIcyto %>% select(key,object_pIRF3_ctip2_pos_cyto_intensity,area_pIRF3_ctip2_pos_cyto_intensity,MeanGrey_pIRF3_ctip2_pos_cyto_intensity,IntDen_pIRF3_ctip2_pos_cyto_intensity,RawIntDen_pIRF3_ctip2_pos_cyto_intensity,well,channel,site)
DCPIcyto<-DCPIcyto %>% filter_all(all_vars(!is.na(.)))

##summarize CTIP2 POSITIVE cytobody data
dcpicyto_summ2<-DCPIcyto %>% group_by(well,site)%>%
  summarize(object_pIRF3_ctip2_pos_cyto_intensity_site=length(unique(object_pIRF3_ctip2_pos_cyto_intensity)),area_pIRF3_ctip2_pos_cyto_intensity_site=sum(area_pIRF3_ctip2_pos_cyto_intensity),
            IntDen_pIRF3_ctip2_pos_cyto_intensity_site=sum(IntDen_pIRF3_ctip2_pos_cyto_intensity),MeanGrey_pIRF3_ctip2_pos_cyto_intensity_site=sum(MeanGrey_pIRF3_ctip2_pos_cyto_intensity)) %>%
  ungroup()

dcpicyto_summ<-dcpicyto_summ2 %>% group_by(well)%>%
  summarize(object_pIRF3_ctip2_pos_cyto_intensity_well=sum(object_pIRF3_ctip2_pos_cyto_intensity_site),area_pIRF3_ctip2_pos_cyto_intensity_well=mean(area_pIRF3_ctip2_pos_cyto_intensity_site),
            IntDen_pIRF3_ctip2_pos_cyto_intensity_well=mean(IntDen_pIRF3_ctip2_pos_cyto_intensity_site),MeanGrey_pIRF3_ctip2_pos_cyto_intensity_well=mean(MeanGrey_pIRF3_ctip2_pos_cyto_intensity_site)) %>%
  ungroup()


#CTIP2 POSITIVE nucleus
#load intensity data from pIRF3 in ctip2+ nuclei 
data_pIRF3_ctip2_pos_nucl_intensity<-lapply(files_pIRF3_ctip2_pos_nucl_intensity,function(x) {read.csv(x,header=T,stringsAsFactors=F)})
names(data_pIRF3_ctip2_pos_nucl_intensity)<-files_pIRF3_ctip2_pos_nucl_intensity %>% str_split(patter="/",simplify=T) %>% .[,4] %>% str_split(pattern="\\.",simplify=T) %>% .[,1]
for(i in seq_along(data_pIRF3_ctip2_pos_nucl_intensity)) {
  data_pIRF3_ctip2_pos_nucl_intensity[[i]]$key<-rep(meta$key[i],nrow(data_pIRF3_ctip2_pos_nucl_intensity[[i]])) }
DCPInucl<-do.call(rbind,data_pIRF3_ctip2_pos_nucl_intensity)
colnames(DCPInucl)<-c("object_pIRF3_ctip2_pos_nucl_intensity","area_pIRF3_ctip2_pos_nucl_intensity","MeanGrey_pIRF3_ctip2_pos_nucl_intensity","StdDev_pIRF3_ctip2_pos_nucl_intensity","MinGrey_pIRF3_ctip2_pos_nucl_intensity",
                      "MaxGrey_pIRF3_ctip2_pos_nucl_intensity","Feret_pIRF3_ctip2_pos_nucl_intensity","IntDen_pIRF3_ctip2_pos_nucl_intensity","RawIntDen_pIRF3_ctip2_pos_nucl_intensity","FeretX_pIRF3_ctip2_pos_nucl_intensity","FeretY_pIRF3_ctip2_pos_nucl_intensity",
                      "Feretangle_pIRF3_ctip2_pos_nucl_intensity","MinFeret_pIRF3_ctip2_pos_nucl_intensity","MinThr_pIRF3_ctip2_pos_nucl_intensity","MaxThr_pIRF3_ctip2_pos_nucl_intensity","key")
DCPInucl<-cbind(DCPInucl,DCPInucl$key %>% str_split(pattern="_",simplify=T),stringsAsFactors=F)
colnames(DCPInucl)<-c("object_pIRF3_ctip2_pos_nucl_intensity","area_pIRF3_ctip2_pos_nucl_intensity","MeanGrey_pIRF3_ctip2_pos_nucl_intensity","StdDev_pIRF3_ctip2_pos_nucl_intensity","MinGrey_pIRF3_ctip2_pos_nucl_intensity",
                      "MaxGrey_pIRF3_ctip2_pos_nucl_intensity","Feret_pIRF3_ctip2_pos_nucl_intensity","IntDen_pIRF3_ctip2_pos_nucl_intensity","RawIntDen_pIRF3_ctip2_pos_nucl_intensity","FeretX_pIRF3_ctip2_pos_nucl_intensity","FeretY_pIRF3_ctip2_pos_nucl_intensity",
                      "Feretangle_pIRF3_ctip2_pos_nucl_intensity","MinFeret_pIRF3_ctip2_pos_nucl_intensity","MinThr_pIRF3_ctip2_pos_nucl_intensity","MaxThr_pIRF3_ctip2_pos_nucl_intensity","key","well","channel","site")
DCPInucl<-DCPInucl %>% select(key,object_pIRF3_ctip2_pos_nucl_intensity,area_pIRF3_ctip2_pos_nucl_intensity,MeanGrey_pIRF3_ctip2_pos_nucl_intensity,IntDen_pIRF3_ctip2_pos_nucl_intensity,RawIntDen_pIRF3_ctip2_pos_nucl_intensity,well,channel,site)
DCPInucl<- DCPInucl %>% filter_all(all_vars(!is.na(.)))

##summarize pIRF3 in CTIP2+ positive nuclei cell data
dcpinucl_summ2<-DCPInucl %>% group_by(well,site)%>%
  summarize(object_pIRF3_ctip2_pos_nucl_intensity_site=length(unique(object_pIRF3_ctip2_pos_nucl_intensity)),area_pIRF3_ctip2_pos_nucl_intensity_site=sum(area_pIRF3_ctip2_pos_nucl_intensity),
            IntDen_pIRF3_ctip2_pos_nucl_intensity_site=sum(IntDen_pIRF3_ctip2_pos_nucl_intensity),MeanGrey_pIRF3_ctip2_pos_nucl_intensity_site=sum(MeanGrey_pIRF3_ctip2_pos_nucl_intensity)) %>%
  ungroup()

dcpinucl_summ<-dcpinucl_summ2 %>% group_by(well)%>%
  summarize(object_pIRF3_ctip2_pos_nucl_intensity_well=sum(object_pIRF3_ctip2_pos_nucl_intensity_site),area_pIRF3_ctip2_pos_nucl_intensity_well=mean(area_pIRF3_ctip2_pos_nucl_intensity_site),
            IntDen_pIRF3_ctip2_pos_nucl_intensity_well=mean(IntDen_pIRF3_ctip2_pos_nucl_intensity_site),MeanGrey_pIRF3_ctip2_pos_nucl_intensity_well=mean(MeanGrey_pIRF3_ctip2_pos_nucl_intensity_site)) %>%
  ungroup()


### CTIP2 negative cells 

#CTIP2 negative cytoplasm
#load intensity data from pIRF3 in ctip2- cytobody count
data_pIRF3_ctip2_neg_cyto_intensity<-lapply(files_pIRF3_ctip2_neg_cyto_intensity,function(x) {read.csv(x,header=T,stringsAsFactors=F)})
names(data_pIRF3_ctip2_neg_cyto_intensity)<-files_pIRF3_ctip2_neg_cyto_intensity %>% str_split(patter="/",simplify=T) %>% .[,4] %>% str_split(pattern="\\.",simplify=T) %>% .[,1]
for(i in seq_along(data_pIRF3_ctip2_neg_cyto_intensity)) {
  data_pIRF3_ctip2_neg_cyto_intensity[[i]]$key<-rep(meta$key[i],nrow(data_pIRF3_ctip2_neg_cyto_intensity[[i]])) }
DCNIcyto<-do.call(rbind,data_pIRF3_ctip2_neg_cyto_intensity)
colnames(DCNIcyto)<-c("object_pIRF3_ctip2_neg_cyto_intensity","area_pIRF3_ctip2_neg_cyto_intensity","MeanGrey_pIRF3_ctip2_neg_cyto_intensity","StdDev_pIRF3_ctip2_neg_cyto_intensity","MinGrey_pIRF3_ctip2_neg_cyto_intensity",
                      "MaxGrey_pIRF3_ctip2_neg_cyto_intensity","Feret_pIRF3_ctip2_neg_cyto_intensity","IntDen_pIRF3_ctip2_neg_cyto_intensity","RawIntDen_pIRF3_ctip2_neg_cyto_intensity","FeretX_pIRF3_ctip2_neg_cyto_intensity","FeretY_pIRF3_ctip2_neg_cyto_intensity",
                      "Feretangle_pIRF3_ctip2_neg_cyto_intensity","MinFeret_pIRF3_ctip2_neg_cyto_intensity","MinThr_pIRF3_ctip2_neg_cyto_intensity","MaxThr_pIRF3_ctip2_neg_cyto_intensity","key")
DCNIcyto<-cbind(DCNIcyto,DCNIcyto$key %>% str_split(pattern="_",simplify=T),stringsAsFactors=F)
colnames(DCNIcyto)<-c("object_pIRF3_ctip2_neg_cyto_intensity","area_pIRF3_ctip2_neg_cyto_intensity","MeanGrey_pIRF3_ctip2_neg_cyto_intensity","StdDev_pIRF3_ctip2_neg_cyto_intensity","MinGrey_pIRF3_ctip2_neg_cyto_intensity",
                      "MaxGrey_pIRF3_ctip2_neg_cyto_intensity","Feret_pIRF3_ctip2_neg_cyto_intensity","IntDen_pIRF3_ctip2_neg_cyto_intensity","RawIntDen_pIRF3_ctip2_neg_cyto_intensity","FeretX_pIRF3_ctip2_neg_cyto_intensity","FeretY_pIRF3_ctip2_neg_cyto_intensity",
                      "Feretangle_pIRF3_ctip2_neg_cyto_intensity","MinFeret_pIRF3_ctip2_neg_cyto_intensity","MinThr_pIRF3_ctip2_neg_cyto_intensity","MaxThr_pIRF3_ctip2_neg_cyto_intensity","key","well","channel","site")
DCNIcyto<-DCNIcyto %>% select(key,object_pIRF3_ctip2_neg_cyto_intensity,area_pIRF3_ctip2_neg_cyto_intensity,MeanGrey_pIRF3_ctip2_neg_cyto_intensity,IntDen_pIRF3_ctip2_neg_cyto_intensity,RawIntDen_pIRF3_ctip2_neg_cyto_intensity,well,channel,site)
DCNIcyto<-DCNIcyto %>% filter_all(all_vars(!is.na(.)))

##summarize pIRF3 in ctip2- cytobody data
dcnicyto_summ2<-DCNIcyto %>% group_by(well,site)%>%
  summarize(object_pIRF3_ctip2_neg_cyto_intensity_site=length(unique(object_pIRF3_ctip2_neg_cyto_intensity)),area_pIRF3_ctip2_neg_cyto_intensity_site=sum(area_pIRF3_ctip2_neg_cyto_intensity),
            IntDen_pIRF3_ctip2_neg_cyto_intensity_site=sum(IntDen_pIRF3_ctip2_neg_cyto_intensity),MeanGrey_pIRF3_ctip2_neg_cyto_intensity_site=sum(MeanGrey_pIRF3_ctip2_neg_cyto_intensity)) %>%
  ungroup()

dcnicyto_summ<-dcnicyto_summ2 %>% group_by(well)%>%
  summarize(object_pIRF3_ctip2_neg_cyto_intensity_well=sum(object_pIRF3_ctip2_neg_cyto_intensity_site),area_pIRF3_ctip2_neg_cyto_intensity_well=mean(area_pIRF3_ctip2_neg_cyto_intensity_site),
            IntDen_pIRF3_ctip2_neg_cyto_intensity_well=mean(IntDen_pIRF3_ctip2_neg_cyto_intensity_site),MeanGrey_pIRF3_ctip2_neg_cyto_intensity_well=mean(MeanGrey_pIRF3_ctip2_neg_cyto_intensity_site)) %>%
  ungroup()


#CTIP2 NEGATIVE nucleus
#load intensity data from pIRF3 in ctip2- nuclei
data_pIRF3_ctip2_neg_nucl_intensity<-lapply(files_pIRF3_ctip2_neg_nucl_intensity,function(x) {read.csv(x,header=T,stringsAsFactors=F)})
names(data_pIRF3_ctip2_neg_nucl_intensity)<-files_pIRF3_ctip2_neg_nucl_intensity %>% str_split(patter="/",simplify=T) %>% .[,4] %>% str_split(pattern="\\.",simplify=T) %>% .[,1]
for(i in seq_along(data_pIRF3_ctip2_neg_nucl_intensity)) {
  data_pIRF3_ctip2_neg_nucl_intensity[[i]]$key<-rep(meta$key[i],nrow(data_pIRF3_ctip2_neg_nucl_intensity[[i]])) }
DCNInucl<-do.call(rbind,data_pIRF3_ctip2_neg_nucl_intensity)
colnames(DCNInucl)<-c("object_pIRF3_ctip2_neg_nucl_intensity","area_pIRF3_ctip2_neg_nucl_intensity","MeanGrey_pIRF3_ctip2_neg_nucl_intensity","StdDev_pIRF3_ctip2_neg_nucl_intensity","MinGrey_pIRF3_ctip2_neg_nucl_intensity",
                      "MaxGrey_pIRF3_ctip2_neg_nucl_intensity","Feret_pIRF3_ctip2_neg_nucl_intensity","IntDen_pIRF3_ctip2_neg_nucl_intensity","RawIntDen_pIRF3_ctip2_neg_nucl_intensity","FeretX_pIRF3_ctip2_neg_nucl_intensity","FeretY_pIRF3_ctip2_neg_nucl_intensity",
                      "Feretangle_pIRF3_ctip2_neg_nucl_intensity","MinFeret_pIRF3_ctip2_neg_nucl_intensity","MinThr_pIRF3_ctip2_neg_nucl_intensity","MaxThr_pIRF3_ctip2_neg_nucl_intensity","key")
DCNInucl<-cbind(DCNInucl,DCNInucl$key %>% str_split(pattern="_",simplify=T),stringsAsFactors=F)
colnames(DCNInucl)<-c("object_pIRF3_ctip2_neg_nucl_intensity","area_pIRF3_ctip2_neg_nucl_intensity","MeanGrey_pIRF3_ctip2_neg_nucl_intensity","StdDev_pIRF3_ctip2_neg_nucl_intensity","MinGrey_pIRF3_ctip2_neg_nucl_intensity",
                      "MaxGrey_pIRF3_ctip2_neg_nucl_intensity","Feret_pIRF3_ctip2_neg_nucl_intensity","IntDen_pIRF3_ctip2_neg_nucl_intensity","RawIntDen_pIRF3_ctip2_neg_nucl_intensity","FeretX_pIRF3_ctip2_neg_nucl_intensity","FeretY_pIRF3_ctip2_neg_nucl_intensity",
                      "Feretangle_pIRF3_ctip2_neg_nucl_intensity","MinFeret_pIRF3_ctip2_neg_nucl_intensity","MinThr_pIRF3_ctip2_neg_nucl_intensity","MaxThr_pIRF3_ctip2_neg_nucl_intensity","key","well","channel","site")

DCNInucl<-DCNInucl %>% select(key,object_pIRF3_ctip2_neg_nucl_intensity,area_pIRF3_ctip2_neg_nucl_intensity,MeanGrey_pIRF3_ctip2_neg_nucl_intensity,IntDen_pIRF3_ctip2_neg_nucl_intensity,RawIntDen_pIRF3_ctip2_neg_nucl_intensity,well,channel,site)
DCNInucl<- DCNInucl %>% filter_all(all_vars(!is.na(.)))
#DFNt %>% na.omit()
#DFNt[is.na(DFNt)]<-0

##summarize pIRF3 in ctip2- nuclei cell data
DCNInucl_summ2<-DCNInucl %>% group_by(well,site)%>%
  summarize(object_pIRF3_ctip2_neg_nucl_intensity_site=length(unique(object_pIRF3_ctip2_neg_nucl_intensity)),area_pIRF3_ctip2_neg_nucl_intensity_site=sum(area_pIRF3_ctip2_neg_nucl_intensity),
            IntDen_pIRF3_ctip2_neg_nucl_intensity_site=sum(IntDen_pIRF3_ctip2_neg_nucl_intensity),MeanGrey_pIRF3_ctip2_neg_nucl_intensity_site=sum(MeanGrey_pIRF3_ctip2_neg_nucl_intensity)) %>%
  ungroup()

DCNInucl_summ<-DCNInucl_summ2 %>% group_by(well)%>%
  summarize(object_pIRF3_ctip2_neg_nucl_intensity_well=sum(object_pIRF3_ctip2_neg_nucl_intensity_site),area_pIRF3_ctip2_neg_nucl_intensity_well=mean(area_pIRF3_ctip2_neg_nucl_intensity_site),
            IntDen_pIRF3_ctip2_neg_nucl_intensity_well=mean(IntDen_pIRF3_ctip2_neg_nucl_intensity_site),MeanGrey_pIRF3_ctip2_neg_nucl_intensity_well=mean(MeanGrey_pIRF3_ctip2_neg_nucl_intensity_site)) %>%
  ungroup()


## make final DF
## build H2AX and STING DFs

DT_test<-DT
DT_test<-left_join(DT_test,platemap %>% unique(),by=c('well'))
DT_test<-DT_test %>% ungroup() %>% select(
  key,well,site,genotype,treatment,line,duration)

meta1<-DT_test

Data_summary <- meta1

##join with cell body information

#cells

Data_summary<-left_join(Data_summary %>% unique(),dcp_summ %>% unique(),by=c('well'))
Data_summary<-left_join(Data_summary,dcp_summ2,by=c('well','site'))

Data_summary<-left_join(Data_summary,dcn_summ,by=c('well'))
Data_summary<-left_join(Data_summary,dcn_summ2,by=c('well','site'))



#intensity 
Data_summary<-left_join(Data_summary,dcpinucl_summ,by=c('well'))
Data_summary<-left_join(Data_summary,dcpinucl_summ2,by=c('well','site'))

Data_summary<-left_join(Data_summary,dcpicyto_summ,by=c('well'))
Data_summary<-left_join(Data_summary,dcpicyto_summ2,by=c('well','site'))

Data_summary<-left_join(Data_summary,DCNInucl_summ,by=c('well'))
Data_summary<-left_join(Data_summary,DCNInucl_summ2,by=c('well','site'))

Data_summary<-left_join(Data_summary,dcnicyto_summ,by=c('well'))
Data_summary<-left_join(Data_summary,dcnicyto_summ2,by=c('well','site'))


# file, sample, key2
Data_summary<-Data_summary %>% ungroup() %>% select(
  key,well,site,genotype,line,treatment,duration,
  ctip2_pos_well,sites_well.x,ctip2_pos_site,avg_area_ctip2_pos,
  ctip2_neg_well,ctip2_neg_site,avg_area_ctip2_neg,
  
  #intensity
  object_pIRF3_ctip2_pos_nucl_intensity_well,area_pIRF3_ctip2_pos_nucl_intensity_well,IntDen_pIRF3_ctip2_pos_nucl_intensity_well,MeanGrey_pIRF3_ctip2_pos_nucl_intensity_well,
  object_pIRF3_ctip2_pos_nucl_intensity_site,area_pIRF3_ctip2_pos_nucl_intensity_site,IntDen_pIRF3_ctip2_pos_nucl_intensity_site,MeanGrey_pIRF3_ctip2_pos_nucl_intensity_site,
   
  object_pIRF3_ctip2_pos_cyto_intensity_well,area_pIRF3_ctip2_pos_cyto_intensity_well,IntDen_pIRF3_ctip2_pos_cyto_intensity_well,MeanGrey_pIRF3_ctip2_pos_cyto_intensity_well,
  object_pIRF3_ctip2_pos_cyto_intensity_site,area_pIRF3_ctip2_pos_cyto_intensity_site,IntDen_pIRF3_ctip2_pos_cyto_intensity_site,MeanGrey_pIRF3_ctip2_pos_cyto_intensity_site,
  
  object_pIRF3_ctip2_neg_nucl_intensity_well,area_pIRF3_ctip2_neg_nucl_intensity_well,IntDen_pIRF3_ctip2_neg_nucl_intensity_well,MeanGrey_pIRF3_ctip2_neg_nucl_intensity_well,
  object_pIRF3_ctip2_neg_nucl_intensity_site,area_pIRF3_ctip2_neg_nucl_intensity_site,IntDen_pIRF3_ctip2_neg_nucl_intensity_site,MeanGrey_pIRF3_ctip2_neg_nucl_intensity_site,
  
  object_pIRF3_ctip2_neg_cyto_intensity_well,area_pIRF3_ctip2_neg_cyto_intensity_well,IntDen_pIRF3_ctip2_neg_cyto_intensity_well,MeanGrey_pIRF3_ctip2_neg_cyto_intensity_well,
  object_pIRF3_ctip2_neg_cyto_intensity_site,area_pIRF3_ctip2_neg_cyto_intensity_site,IntDen_pIRF3_ctip2_neg_cyto_intensity_site,MeanGrey_pIRF3_ctip2_neg_cyto_intensity_site)


colnames(Data_summary)<-c("key","well","site","genotype","line","treatment","duration",
                          "ctip2_pos_well","sites_well","ctip2_pos_site","avg_area_ctip2_pos",
                          "ctip2_neg_well","ctip2_neg_site","avg_area_ctip2_neg",
                          
                          # intensity
                          "object_pIRF3_ctip2_pos_nucl_intensity_well","area_pIRF3_ctip2_pos_nucl_intensity_well","IntDen_pIRF3_ctip2_pos_nucl_intensity_well","MeanGrey_pIRF3_ctip2_pos_nucl_intensity_well",
                          "object_pIRF3_ctip2_pos_nucl_intensity_site","area_pIRF3_ctip2_pos_nucl_intensity_site","IntDen_pIRF3_ctip2_pos_nucl_intensity_site","MeanGrey_pIRF3_ctip2_pos_nucl_intensity_site",
                          
                          "object_pIRF3_ctip2_pos_cyto_intensity_well","area_pIRF3_ctip2_pos_cyto_intensity_well","IntDen_pIRF3_ctip2_pos_cyto_intensity_well","MeanGrey_pIRF3_ctip2_pos_cyto_intensity_well",
                          "object_pIRF3_ctip2_pos_cyto_intensity_site","area_pIRF3_ctip2_pos_cyto_intensity_site","IntDen_pIRF3_ctip2_pos_cyto_intensity_site","MeanGrey_pIRF3_ctip2_pos_cyto_intensity_site",
                          
                          "object_pIRF3_ctip2_neg_nucl_intensity_well","area_pIRF3_ctip2_neg_nucl_intensity_well","IntDen_pIRF3_ctip2_neg_nucl_intensity_well","MeanGrey_pIRF3_ctip2_neg_nucl_intensity_well",
                          "object_pIRF3_ctip2_neg_nucl_intensity_site","area_pIRF3_ctip2_neg_nucl_intensity_site","IntDen_pIRF3_ctip2_neg_nucl_intensity_site","MeanGrey_pIRF3_ctip2_neg_nucl_intensity_site",
                    
                          "object_pIRF3_ctip2_neg_cyto_intensity_well","area_pIRF3_ctip2_neg_cyto_intensity_well","IntDen_pIRF3_ctip2_neg_cyto_intensity_well","MeanGrey_pIRF3_ctip2_neg_cyto_intensity_well",
                          "object_pIRF3_ctip2_neg_cyto_intensity_site","area_pIRF3_ctip2_neg_cyto_intensity_site","IntDen_pIRF3_ctip2_neg_cyto_intensity_site","MeanGrey_pIRF3_ctip2_neg_cyto_intensity_site")

Data_summary2<-Data_summary

Data_summary2$ctip2_pos_site[which(Data_summary2$ctip2_pos_site%>%is.na())]<-0
Data_summary2$ctip2_neg_site[which(Data_summary2$ctip2_neg_site%>%is.na())]<-0

Data_summary2$object_pIRF3_ctip2_pos_nucl_intensity_site[which(Data_summary2$object_pIRF3_ctip2_pos_nucl_intensity_site%>%is.na())]<-0
Data_summary2$object_pIRF3_ctip2_pos_cyto_intensity_site[which(Data_summary2$object_pIRF3_ctip2_pos_cyto_intensity_site%>%is.na())]<-0
Data_summary2$object_pIRF3_ctip2_neg_nucl_intensity_site[which(Data_summary2$object_pIRF3_ctip2_neg_nucl_intensity_site%>%is.na())]<-0
Data_summary2$object_pIRF3_ctip2_neg_cyto_intensity_site[which(Data_summary2$object_pIRF3_ctip2_neg_cyto_intensity_site%>%is.na())]<-0

Data_summary2[is.na(Data_summary2)]<-0

# add a few columns corresponding to avg ctip2+/- cells/site and total cells/site (adding ctip2+/- cells), total cyto/nucl intensity
Data_summary2$avg_ctip2_pos_site<-Data_summary2$ctip2_pos_well/Data_summary2$sites_well
Data_summary2$avg_ctip2_neg_site<-Data_summary2$ctip2_neg_well/Data_summary2$sites_well
Data_summary2$Total_tuj_site<-Data_summary2$ctip2_pos_site+Data_summary2$ctip2_neg_site
Data_summary2$Total_nucl_intensity<-Data_summary2$IntDen_pIRF3_ctip2_pos_nucl_intensity_site+Data_summary2$IntDen_pIRF3_ctip2_neg_nucl_intensity_site
Data_summary2$Total_cyto_intensity<-Data_summary2$IntDen_pIRF3_ctip2_pos_cyto_intensity_site+Data_summary2$IntDen_pIRF3_ctip2_neg_cyto_intensity_site
Data_summary2$Total_nucl_area<-Data_summary2$area_pIRF3_ctip2_pos_nucl_intensity_site+Data_summary2$area_pIRF3_ctip2_neg_nucl_intensity_site



## reorder factor levels as desired 
Data_summary2$line<-as.character(Data_summary2$line)
Data_summary2$line<-factor(Data_summary2$line,levels = unique(Data_summary2$line))

Data_summary2$treatment<-as.character(Data_summary2$treatment)
Data_summary2$treatment<-factor(Data_summary2$treatment,levels = unique(Data_summary2$treatment))

Data_summary2$genotype<-as.character(Data_summary2$genotype)
Data_summary2$genotype<-factor(Data_summary2$genotype,levels = unique(Data_summary2$genotype))

Data_summary2$duration<-as.character(Data_summary2$duration)
Data_summary2$duration<-factor(Data_summary2$duration,levels = unique(Data_summary2$duration))

Data_summary2$duration<-factor(Data_summary2$duration,level=c("1H"))



pdf("pIRF3-plots.pdf")

########################################## pIRF3 ANALYSIS #############################################################################


############################# DMSO - ETOPOSIDE
##### cell number

#total cells 
ggplot(Data_summary2 %>% filter(treatment %in% c('DMSO','etoposide')) %>% group_by(treatment,genotype,line,well,duration) %>% summarize(
  IntDensCTIP2pos_well=mean(IntDen_pIRF3_ctip2_pos_cyto_intensity_site),tuj_well=mean(Total_tuj_site),ctip2_pos_well=mean(ctip2_pos_site),ctip2_neg_well=mean(ctip2_neg_site)),
  aes(x=duration,y=((tuj_well)),color=treatment)) +theme_bw() + 
  facet_grid(~treatment) + theme(strip.text.x = element_text(size = 8, colour = "black", angle = 90))+
  ggtitle("Mean number of TUJ+ neurons (per site)") + 
  ylim(0,80) + ylab("Mean number of TUJ+ neurons (per site)") + xlab("") + theme(axis.title = element_text(size = 10)) + 
  stat_summary(fun.data=mean_cl_boot,geom="crossbar",fill="grey50",alpha=0.2) +
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = 'bold')) +
  geom_jitter(height=0,width=0.1,size=2,alpha=0.7) +
  theme(axis.text.x=element_text(angle=30,hjust=1)) +
  scale_color_manual(values=c("grey20","blue2")) 

#total cells per site per well 
ggplot(Data_summary2 %>% filter(treatment %in% c('DMSO','etoposide')),
       aes(x=well,y=Total_tuj_site, color=treatment)) +theme_bw() + 
  # facet_grid(~treatment) + 
  theme(strip.text.x = element_text(size = 8, colour = "black", angle = 90))+
  ggtitle("Mean number of TUJ+ neurons (per site)") + 
  ylim(0,NA) + ylab("Mean number of TUJ+ neurons (per site)") + xlab("") + theme(axis.title = element_text(size = 10)) + 
  stat_summary(fun.data=mean_sdl,fun.args=list(mult=1.5),geom="crossbar",fill="grey50",alpha=0.2) +
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = 'bold')) +
  geom_jitter(height=0,width=0.1,size=2,alpha=0.7) +
  theme(axis.text.x=element_text(angle=30,hjust=1)) +
  scale_color_manual(values=c("grey20","blue2")) 

hoeschtcount_ds <- Data_summary2 %>% filter(treatment %in% c('DMSO','etoposide')) %>% group_by(well) %>% summarize(tuj_site = mean(Total_tuj_site),
                                                                                                                             n_site    = n(),
                                                                                                                             sd_tuj_site   = 1.5*sd(Total_tuj_site))

hoeschtcount_ds %>% summarize(mean_mean = mean(tuj_site),
                              mean_sd = mean(sd_tuj_site))

# calculations (all)
# 33.5+28.5 = 62 -> 62
# 33.5-28.5 = 5 -> 5


## absolute range (strict 8-57)
Data_summary_clean_etop <- Data_summary2 %>% filter(treatment %in% c('DMSO','etoposide')) %>% filter(Total_tuj_site > 4 & Total_tuj_site < 63)

# check sample size or number of sites per well. for D15, lowest n=5 (J13), highest n=18 (E12)
DS_clean_etop_sum <- Data_summary_clean_etop %>% group_by(well) %>% summarize(n = n())

#plot cleaned cell count per site per well
ggplot(Data_summary_clean_etop %>% filter(treatment %in% c('DMSO','etoposide')),
       aes(x=well,y=Total_tuj_site, color=treatment)) +theme_bw() + 
  # facet_grid(~treatment) + 
  theme(strip.text.x = element_text(size = 8, colour = "black", angle = 90))+
  ggtitle("Mean number of TUJ+ neurons (per site)") + 
  ylim(0,NA) + ylab("Mean number of TUJ+ neurons (per site)") + xlab("") + theme(axis.title = element_text(size = 10)) + 
  stat_summary(fun.data=mean_se,geom="crossbar",fill="grey50",alpha=0.2) +
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = 'bold')) +
  geom_jitter(height=0,width=0.1,size=2,alpha=0.7) +
  theme(axis.text.x=element_text(angle=30,hjust=1)) +
  scale_color_manual(values=c("grey20","blue2")) 


#ctip2 positive
ggplot(Data_summary_clean_etop %>% filter(treatment %in% c('DMSO','etoposide')) %>% group_by(treatment,genotype,line,well,duration) %>% summarize(
  IntDensCTIP2pos_well=mean(IntDen_pIRF3_ctip2_pos_cyto_intensity_site),ctip2_pos_well=mean(ctip2_pos_site)),
  aes(x=duration,y=((ctip2_pos_well)),color=treatment)) +theme_bw() + 
  facet_grid(~treatment) + theme(strip.text.x = element_text(size = 8, colour = "black", angle = 90))+
  ggtitle("Mean number of CTIP2+ cortical neurons (per site)") + 
  ylim(0,60) + ylab("Mean number of CTIP2+ cortical neurons (per site)") + xlab("") + theme(axis.title = element_text(size = 10)) + 
  stat_summary(fun.data=mean_cl_boot,geom="crossbar",fill="grey50",alpha=0.2) +
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = 'bold')) +
  geom_jitter(height=0,width=0.1,size=2,alpha=0.7) +
  theme(axis.text.x=element_text(angle=30,hjust=1)) +
  scale_color_manual(values=c("grey20","blue2")) 


#ctip2 positive - normalized to all cells
ggplot(Data_summary_clean_etop %>% filter(treatment %in% c('DMSO','etoposide')) %>% group_by(treatment,genotype,line,well,duration) %>% summarize(
  IntDensCTIP2pos_well=mean(IntDen_pIRF3_ctip2_pos_cyto_intensity_site),tuj_well=mean(Total_tuj_site),ctip2_pos_well=mean(ctip2_pos_site)),
  aes(x=duration,y=((ctip2_pos_well)/(tuj_well))*100,color=treatment)) +theme_bw() + 
  facet_grid(~treatment) + theme(strip.text.x = element_text(size = 8, colour = "black", angle = 90))+
  ggtitle("Percentage of CTIP2+ cortical neurons (per site)") + 
  ylim(0,NA) + ylab("Percentage of CTIP2+ cortical neurons (per site)") + xlab("") + theme(axis.title = element_text(size = 10)) + 
  stat_summary(fun.data=mean_cl_boot,geom="crossbar",fill="grey50",alpha=0.2) +
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = 'bold')) +
  geom_jitter(height=0,width=0.1,size=2,alpha=0.7) +
  theme(axis.text.x=element_text(angle=30,hjust=1)) +
  scale_color_manual(values=c("grey20","blue2")) 

#ctip2 negative 
ggplot(Data_summary_clean_etop %>% filter(treatment %in% c('DMSO','etoposide')) %>% group_by(treatment,genotype,line,well,duration) %>% summarize(
  IntDensCTIP2neg_well=mean(IntDen_pIRF3_ctip2_neg_cyto_intensity_site),tuj_well=mean(Total_tuj_site),ctip2_neg_well=mean(ctip2_neg_site)),
  aes(x=duration,y=((ctip2_neg_well)),color=treatment)) +theme_bw() + 
  facet_grid(~treatment) + theme(strip.text.x = element_text(size = 8, colour = "black", angle = 90))+
  ggtitle("Mean number of CTIP2- cortical neurons (per site)") + 
  ylim(0,60) + ylab("Mean number of CTIP2- cortical neurons (per site)") + xlab("") + theme(axis.title = element_text(size = 10)) + 
  stat_summary(fun.data=mean_cl_boot,geom="crossbar",fill="grey50",alpha=0.2) +
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = 'bold')) +
  geom_jitter(height=0,width=0.1,size=2,alpha=0.7) +
  theme(axis.text.x=element_text(angle=30,hjust=1)) +
  scale_color_manual(values=c("grey20","blue2")) 

#ctip2 negative PER CELL - normalized to all cells
ggplot(Data_summary_clean_etop %>% filter(treatment %in% c('DMSO','etoposide')) %>% group_by(treatment,genotype,line,well,duration) %>% summarize(
  IntDensCTIP2neg_well=mean(IntDen_pIRF3_ctip2_neg_cyto_intensity_site),tuj_well=mean(Total_tuj_site),ctip2_neg_well=mean(ctip2_neg_site)),
  aes(x=duration,y=((ctip2_neg_well)/(tuj_well)*100),color=treatment)) +theme_bw() + 
  facet_grid(~treatment) + theme(strip.text.x = element_text(size = 8, colour = "black", angle = 90))+
  ggtitle("Percentage of CTIP2- cortical neurons (per site)") + 
  ylim(0,NA) + ylab("Percentage of CTIP2- cortical neurons (per site)") + xlab("") + theme(axis.title = element_text(size = 10)) + 
  stat_summary(fun.data=mean_cl_boot,geom="crossbar",fill="grey50",alpha=0.2) +
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = 'bold')) +
  geom_jitter(height=0,width=0.1,size=2,alpha=0.7) +
  theme(axis.text.x=element_text(angle=30,hjust=1)) +
  scale_color_manual(values=c("grey20","blue2")) 



#PER CELL - pIRF3 nuclear
##integrated intensity 


### QC (CTIP2+ neurons )
ggplot(Data_summary_clean_etop %>% filter(treatment %in% c('DMSO','etoposide')),
       aes(x=well,y=IntDen_pIRF3_ctip2_pos_nucl_intensity_site,color=treatment)) +theme_bw() + 
  # facet_grid(~treatment) + 
  theme(strip.text.x = element_text(size = 8, colour = "black", angle = 90))+
  ggtitle("QC: Integrated Intensity of nuclear pIRF3 in \n CTIP2+ cortical neurons (per site per well)") + 
  ylim(0,NA) + ylab("Integrated Intensity of nuclear pIRF3 in \n CTIP2+ cortical neurons (per site per well)") + xlab("") + theme(axis.title = element_text(size = 10)) + 
  stat_summary(fun.data=mean_sdl,geom="crossbar",fill="grey50",alpha=0.2) +
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = 'bold')) +
  geom_jitter(height=0,width=0.1,size=2,alpha=0.7) +
  theme(axis.text.x=element_text(angle=30,hjust=1)) +
  scale_color_manual(values=c("grey20","blue2")) 

# identify and remove outliers outside of 2SD
Data_summary_clean_etop_POSint <- Data_summary_clean_etop %>% group_by(well) %>% filter(
  IntDen_pIRF3_ctip2_pos_nucl_intensity_site <= mean(IntDen_pIRF3_ctip2_pos_nucl_intensity_site) + 2*sd(IntDen_pIRF3_ctip2_pos_nucl_intensity_site))  %>% ungroup()

### QC (CTIP2- neurons)
ggplot(Data_summary_clean_etop %>% filter(treatment %in% c('DMSO','etoposide')),
       aes(x=well,y=IntDen_pIRF3_ctip2_neg_nucl_intensity_site,color=treatment)) +theme_bw() + 
  # facet_grid(~treatment) + 
  theme(strip.text.x = element_text(size = 8, colour = "black", angle = 90))+
  ggtitle("QC: Integrated Intensity of nuclear pIRF3 in \n CTIP2- cortical neurons (per site per well)") + 
  ylim(0,NA) + ylab("Integrated Intensity of nuclear pIRF3 in \n CTIP2- cortical neurons (per site per well)") + xlab("") + theme(axis.title = element_text(size = 10)) + 
  stat_summary(fun.data=mean_sdl,geom="crossbar",fill="grey50",alpha=0.2) +
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = 'bold')) +
  geom_jitter(height=0,width=0.1,size=2,alpha=0.7) +
  theme(axis.text.x=element_text(angle=30,hjust=1)) +
  scale_color_manual(values=c("grey20","blue2")) 

# identify and remove outliers outside of 2SD
Data_summary_clean_etop_NEGint <- Data_summary_clean_etop %>% group_by(well) %>% filter(
  IntDen_pIRF3_ctip2_neg_nucl_intensity_site <= mean(IntDen_pIRF3_ctip2_neg_nucl_intensity_site) + 2*sd(IntDen_pIRF3_ctip2_neg_nucl_intensity_site))  %>% ungroup()



########pIRF3 -ctip2 positive
#pIRF3 nucl - integrated intensity
ggplot(Data_summary_clean_etop_POSint %>% filter(treatment %in% c('DMSO','etoposide')) %>% group_by(treatment,genotype,line,well,duration) %>% summarize(
  IntDensCTIP2pos_well=mean(IntDen_pIRF3_ctip2_pos_nucl_intensity_site),ctip2_pos_well=mean(ctip2_pos_site),tuj_well=mean(Total_tuj_site),ctip2_neg_well=mean(ctip2_neg_site)),
  aes(x=duration,y=(IntDensCTIP2pos_well),color=treatment)) +theme_bw() + 
  facet_grid(~treatment) + theme(strip.text.x = element_text(size = 8, colour = "black", angle = 90))+
  ggtitle("Average Integrated Intensity of nuclear pIRF3 in \n CTIP2+ cortical neurons (per site per well)") + 
  ylim(0,NA) + ylab("Average Integrated Intensity of nuclear pIRF3 in \n CTIP2+ cortical neurons (per site per well)") + xlab("") + theme(axis.title = element_text(size = 10)) + 
  stat_summary(fun.data=mean_se,geom="crossbar",fill="grey50",alpha=0.2) +
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = 'bold')) +
  geom_jitter(height=0,width=0.1,size=2,alpha=0.7) +
  theme(axis.text.x=element_text(angle=30,hjust=1)) +
  scale_color_manual(values=c("grey20","blue2")) 

#pIRF3 nucl - integrated intensity - PER ctip2+ CELL
ggplot(Data_summary_clean_etop_POSint %>% filter(treatment %in% c('DMSO','etoposide')) %>% group_by(treatment,genotype,line,well,duration) %>% summarize(
  IntDensCTIP2pos_well=mean(IntDen_pIRF3_ctip2_pos_nucl_intensity_site),ctip2_pos_well=mean(ctip2_pos_site),tuj_well=mean(Total_tuj_site),ctip2_neg_well=mean(ctip2_neg_site)),
  aes(x=duration,y=((IntDensCTIP2pos_well)/(ctip2_pos_well)),color=treatment)) +theme_bw() + 
  facet_grid(~treatment) + theme(strip.text.x = element_text(size = 8, colour = "black", angle = 90))+
  ggtitle("Integrated Intensity of nuclear pIRF3 in \n CTIP2+ cortical neurons (per cell per site per well)") + 
  ylim(0,550) + ylab("Integrated Intensity of nuclear pIRF3 in \n CTIP2+ cortical neurons (per cell)") + xlab("") + theme(axis.title = element_text(size = 10)) + 
  stat_summary(fun.data=mean_se,geom="crossbar",fill="grey50",alpha=0.2) +
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = 'bold')) +
  geom_jitter(height=0,width=0.1,size=2,alpha=0.7) +
  theme(axis.text.x=element_text(angle=30,hjust=1)) +
  scale_color_manual(values=c("grey20","blue2")) 

########pIRF3 -ctip2 negative
#pIRF3 nucl - integrated intensity 
ggplot(Data_summary_clean_etop_NEGint %>% filter(treatment %in% c('DMSO','etoposide')) %>% group_by(treatment,genotype,line,well,duration) %>% summarize(
  IntDensCTIP2neg_well=mean(IntDen_pIRF3_ctip2_neg_nucl_intensity_site),ctip2_pos_well=mean(ctip2_pos_site),tuj_well=mean(Total_tuj_site),ctip2_neg_well=mean(ctip2_neg_site)),
  aes(x=duration,y=(IntDensCTIP2neg_well),color=treatment)) +theme_bw() + 
  facet_grid(~treatment) + theme(strip.text.x = element_text(size = 8, colour = "black", angle = 90))+
  ggtitle("Average Integrated Intensity of nuclear pIRF3 \n in CTIP2- cortical neurons (per site per well)") + 
  ylim(0,NA) + ylab("Average Integrated Intensity of nuclear pIRF3 \n in CTIP2- cortical neurons (per site per well)") + xlab("") + theme(axis.title = element_text(size = 10)) + 
  stat_summary(fun.data=mean_se,geom="crossbar",fill="grey50",alpha=0.2) +
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = 'bold')) +
  geom_jitter(height=0,width=0.1,size=2,alpha=0.7) +
  theme(axis.text.x=element_text(angle=30,hjust=1)) +
  scale_color_manual(values=c("grey20","blue2")) 

#pIRF3 nucl - integrated intensity - PER ctip2- CELL
ggplot(Data_summary_clean_etop_NEGint %>% filter(treatment %in% c('DMSO','etoposide')) %>% group_by(treatment,genotype,line,well,duration) %>% summarize(
  IntDensCTIP2neg_well=mean(IntDen_pIRF3_ctip2_neg_nucl_intensity_site),ctip2_pos_well=mean(ctip2_pos_site),tuj_well=mean(Total_tuj_site),ctip2_neg_well=mean(ctip2_neg_site)),
  aes(x=duration,y=((IntDensCTIP2neg_well)/(ctip2_neg_well)),color=treatment)) +theme_bw() + 
  facet_grid(~treatment) + theme(strip.text.x = element_text(size = 8, colour = "black", angle = 90))+
  ggtitle("Integrated Intensity of nuclear pIRF3 \n in CTIP2- cortical neurons\n(per cell per site per well)") + 
  ylim(0,550) + ylab("Integrated Intensity of nuclear pIRF3 in CTIP2- cortical neurons\n(per cell per site per well)") + xlab("") + theme(axis.title = element_text(size = 10)) + 
  stat_summary(fun.data=mean_se,geom="crossbar",fill="grey50",alpha=0.2) +
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = 'bold')) +
  geom_jitter(height=0,width=0.1,size=2,alpha=0.7) +
  theme(axis.text.x=element_text(angle=30,hjust=1)) +
  scale_color_manual(values=c("grey20","blue2")) 



####### pIRF3 nuclear AREA
##pIRF3 AREA


### QC (CTIP2+ neurons )
ggplot(Data_summary_clean_etop %>% filter(treatment %in% c('DMSO','etoposide')),
       aes(x=well,y=area_pIRF3_ctip2_pos_nucl_intensity_site,color=treatment)) +theme_bw() + 
  # facet_grid(~treatment) + 
  theme(strip.text.x = element_text(size = 8, colour = "black", angle = 90))+
  ggtitle("QC: Area of nuclear pIRF3 in \n CTIP2+ cortical neurons (per site per well)") + 
  ylim(0,NA) + ylab("Area of nuclear pIRF3 in \n CTIP2+ cortical neurons (per site per well)") + xlab("") + theme(axis.title = element_text(size = 10)) + 
  stat_summary(fun.data=mean_sdl,geom="crossbar",fill="grey50",alpha=0.2) +
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = 'bold')) +
  geom_jitter(height=0,width=0.1,size=2,alpha=0.7) +
  theme(axis.text.x=element_text(angle=30,hjust=1)) +
  scale_color_manual(values=c("grey20","blue2")) 

# identify and remove outliers outside of 2SD
Data_summary_clean_etop_POSarea <- Data_summary_clean_etop %>% group_by(well) %>% filter(
  area_pIRF3_ctip2_pos_nucl_intensity_site <= mean(area_pIRF3_ctip2_pos_nucl_intensity_site) + 2*sd(area_pIRF3_ctip2_pos_nucl_intensity_site))  %>% ungroup()

### QC (CTIP2- neurons)
ggplot(Data_summary_clean_etop %>% filter(treatment %in% c('DMSO','etoposide')),
       aes(x=well,y=area_pIRF3_ctip2_neg_nucl_intensity_site,color=treatment)) +theme_bw() + 
  # facet_grid(~treatment) + 
  theme(strip.text.x = element_text(size = 8, colour = "black", angle = 90))+
  ggtitle("QC: Area of nuclear pIRF3 in \n CTIP2- cortical neurons (per site per well)") + 
  ylim(0,NA) + ylab("Area of nuclear pIRF3 in \n CTIP2- cortical neurons (per site per well)") + xlab("") + theme(axis.title = element_text(size = 10)) + 
  stat_summary(fun.data=mean_sdl,geom="crossbar",fill="grey50",alpha=0.2) +
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = 'bold')) +
  geom_jitter(height=0,width=0.1,size=2,alpha=0.7) +
  theme(axis.text.x=element_text(angle=30,hjust=1)) +
  scale_color_manual(values=c("grey20","blue2")) 

# identify and remove outliers outside of 2SD
Data_summary_clean_etop_NEGarea <- Data_summary_clean_etop %>% group_by(well) %>% filter(
  area_pIRF3_ctip2_neg_nucl_intensity_site <= mean(area_pIRF3_ctip2_neg_nucl_intensity_site) + 2*sd(area_pIRF3_ctip2_neg_nucl_intensity_site))  %>% ungroup()


########pIRF3 AREA -ctip2 positive
ggplot(Data_summary_clean_etop_POSarea %>% filter(treatment %in% c('DMSO','etoposide')) %>% group_by(treatment,genotype,line,well,duration) %>% summarize(
  IntDensCTIP2pos_well=mean(area_pIRF3_ctip2_pos_nucl_intensity_site),ctip2_pos_well=mean(ctip2_pos_site),tuj_well=mean(Total_tuj_site),ctip2_neg_well=mean(ctip2_neg_site)),
  aes(x=duration,y=(IntDensCTIP2pos_well),color=treatment)) +theme_bw() +
  facet_grid(~treatment) + theme(strip.text.x = element_text(size = 8, colour = "black", angle = 90))+
  ggtitle("Average area of nuclear pIRF3 in \n CTIP2+ cortical neurons (per site per well)") +
  ylim(0,NA) + ylab("Average area of nuclear pIRF3 in \n CTIP2+ cortical neurons(per site per well)") + xlab("") + theme(axis.title = element_text(size = 10)) +
  stat_summary(fun.data=mean_se,geom="crossbar",fill="grey50",alpha=0.2) +
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = 'bold')) +
  geom_jitter(height=0,width=0.1,size=2,alpha=0.7) +
  theme(axis.text.x=element_text(angle=30,hjust=1)) +
  scale_color_manual(values=c("grey20","blue2"))

#pIRF3 nucl particule - integrated intensity - PER ctip2+ CELL
ggplot(Data_summary_clean_etop_POSarea %>% filter(treatment %in% c('DMSO','etoposide')) %>% group_by(treatment,genotype,line,well,duration) %>% summarize(
  IntDensCTIP2pos_well=mean(area_pIRF3_ctip2_pos_nucl_intensity_site),ctip2_pos_well=mean(ctip2_pos_site),tuj_well=mean(Total_tuj_site),ctip2_neg_well=mean(ctip2_neg_site)),
  aes(x=duration,y=((IntDensCTIP2pos_well)/(ctip2_pos_well)),color=treatment)) +theme_bw() +
  facet_grid(~treatment) + theme(strip.text.x = element_text(size = 8, colour = "black", angle = 90))+
  ggtitle("Area of nuclear pIRF3 in CTIP2+ cortical neurons\n(per cell per site per well)") +
  ylim(0,0.4) + ylab("Area of nuclear pIRF3 in CTIP2+ cortical neuronsl\n(per cell per site per well)") + xlab("") + theme(axis.title = element_text(size = 10)) +
  stat_summary(fun.data=mean_se,geom="crossbar",fill="grey50",alpha=0.2) +
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = 'bold')) +
  geom_jitter(height=0,width=0.1,size=2,alpha=0.7) +
  theme(axis.text.x=element_text(angle=30,hjust=1)) +
  scale_color_manual(values=c("grey20","blue2"))

########pIRF3 -ctip2 negative
#pIRF3 nucl - integrated intensity 
ggplot(Data_summary_clean_etop_NEGarea %>% filter(treatment %in% c('DMSO','etoposide')) %>% group_by(treatment,genotype,line,well,duration) %>% summarize(
  IntDensCTIP2neg_well=mean(area_pIRF3_ctip2_neg_nucl_intensity_site),ctip2_pos_well=mean(ctip2_pos_site),tuj_well=mean(Total_tuj_site),ctip2_neg_well=mean(ctip2_neg_site)),
  aes(x=duration,y=(IntDensCTIP2neg_well),color=treatment)) +theme_bw() +
  facet_grid(~treatment) + theme(strip.text.x = element_text(size = 8, colour = "black", angle = 90))+
  ggtitle("Average area of nuclear pIRF3 in \n CTIP2- cortical neurons (per site per well)") +
  ylim(0,NA) + ylab("Average area of nuclear pIRF3 in \n CTIP2- cortical neurons (per site per well)") + xlab("") + theme(axis.title = element_text(size = 10)) +
  stat_summary(fun.data=mean_se,geom="crossbar",fill="grey50",alpha=0.2) +
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = 'bold')) +
  geom_jitter(height=0,width=0.1,size=2,alpha=0.7) +
  theme(axis.text.x=element_text(angle=30,hjust=1)) +
  scale_color_manual(values=c("grey20","blue2"))

#pIRF3 nucl - integrated intensity - PER ctip2- CELL
ggplot(Data_summary_clean_etop_NEGarea %>% filter(treatment %in% c('DMSO','etoposide')) %>% group_by(treatment,genotype,line,well,duration) %>% summarize(
  IntDensCTIP2neg_well=mean(area_pIRF3_ctip2_neg_nucl_intensity_site),ctip2_pos_well=mean(ctip2_pos_site),tuj_well=mean(Total_tuj_site),ctip2_neg_well=mean(ctip2_neg_site)),
  aes(x=duration,y=((IntDensCTIP2neg_well)/(ctip2_neg_well)),color=treatment)) +theme_bw() +
  facet_grid(~treatment) + theme(strip.text.x = element_text(size = 8, colour = "black", angle = 90))+
  ggtitle("Area of nuclear pIRF3 in \n CTIP2- cortical neurons (per cell per site per well)") +
  ylim(0,0.4) + ylab("Area of nuclear pIRF3 in \n CTIP2- cortical neurons (per cell per site per well)") + xlab("") + theme(axis.title = element_text(size = 10)) +
  stat_summary(fun.data=mean_se,geom="crossbar",fill="grey50",alpha=0.2) +
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = 'bold')) +
  geom_jitter(height=0,width=0.1,size=2,alpha=0.7) +
  theme(axis.text.x=element_text(angle=30,hjust=1)) +
  scale_color_manual(values=c("grey20","blue2"))

# run t-tests ETOPOSIDE
etop_POSarea_ttest <- Data_summary_clean_etop_POSarea %>% filter(treatment %in% c('DMSO','etoposide')) %>% filter(duration == "1H") %>% group_by(treatment,genotype,line,well,duration) %>% summarize(
  IntDensCTIP2pos_well=mean(area_pIRF3_ctip2_pos_nucl_intensity_site),ctip2_pos_well=mean(ctip2_pos_site),y = IntDensCTIP2pos_well/ctip2_pos_well)

etop_NEGarea_ttest <- Data_summary_clean_etop_NEGarea %>% filter(treatment %in% c('DMSO','etoposide')) %>% filter(duration == "1H") %>% group_by(treatment,genotype,line,well,duration) %>% summarize(
  IntDensCTIP2neg_well=mean(area_pIRF3_ctip2_neg_nucl_intensity_site),ctip2_neg_well=mean(ctip2_neg_site),y = IntDensCTIP2neg_well/ctip2_neg_well)

t.test(y ~ treatment, data = etop_POSarea_ttest)
t.test(y ~ treatment, data = etop_NEGarea_ttest)



# ALL CELL ANALYSIS 

### etoposide
## QC (all neurons)# 
########pIRF3 -all cells
ggplot(Data_summary_clean_etop %>% filter(treatment %in% c('DMSO','etoposide')),
       aes(x=well,y=Total_nucl_intensity,color=treatment)) +theme_bw() + 
  # facet_grid(~treatment) + 
  theme(strip.text.x = element_text(size = 8, colour = "black", angle = 90))+
  ggtitle("QC: Integrated Intensity of nuclear pIRF3 in \n all cortical neurons (per site per well)") + 
  ylim(0,NA) + ylab("Integrated Intensity of nuclear pIRF3 in \n all cortical neurons (per site per well)") + xlab("") + theme(axis.title = element_text(size = 10)) + 
  stat_summary(fun.data=mean_sdl,geom="crossbar",fill="grey50",alpha=0.2) +
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = 'bold')) +
  geom_jitter(height=0,width=0.1,size=2,alpha=0.7) +
  theme(axis.text.x=element_text(angle=30,hjust=1)) +
  scale_color_manual(values=c("grey20","blue2")) 

# identify and remove outliers outside of 2SD
Data_summary_clean_etop_int <- Data_summary_clean_etop %>% group_by(well) %>% filter(
  Total_nucl_intensity <= mean(Total_nucl_intensity) + 2*sd(Total_nucl_intensity))  %>% ungroup()

########pIRF3 -ALL cells
#pIRF3 nucl - integrated intensity 
ggplot(Data_summary_clean_etop_int %>% filter(treatment %in% c('DMSO','etoposide')) %>% group_by(treatment,genotype,line,well,duration) %>% summarize(
  IntDens_well=mean(Total_nucl_intensity),ctip2_pos_well=mean(ctip2_pos_site),tuj_well=mean(Total_tuj_site),ctip2_neg_well=mean(ctip2_neg_site)),
  aes(x=duration,y=(IntDens_well),color=treatment)) +theme_bw() +
  facet_grid(~treatment) + theme(strip.text.x = element_text(size = 8, colour = "black", angle = 90))+
  ggtitle("Total Integrated Intensity of nuclear pIRF3 in \n Primary Cortical neurons (per site per well)") +
  ylim(0,NA) + ylab("Total Integrated Intensity of nuclear pIRF3 in \n Primary Cortical neurons (per site per well)") + xlab("") + theme(axis.title = element_text(size = 10)) +
  stat_summary(fun.data=mean_cl_boot,geom="crossbar",fill="grey50",alpha=0.2) +
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = 'bold')) +
  geom_jitter(height=0,width=0.1,size=2,alpha=0.7) +
  theme(axis.text.x=element_text(angle=30,hjust=1)) +
  scale_color_manual(values=c("grey20","blue2"))

#pIRF3 nucl - integrated intensity - per ALL cells
ggplot(Data_summary_clean_etop_int %>% filter(treatment %in% c('DMSO','etoposide')) %>% group_by(treatment,genotype,line,well,duration) %>% summarize(
  IntDens_well=mean(Total_nucl_intensity),tuj_well=mean(Total_tuj_site)),
  aes(x=duration,y=((IntDens_well)/(tuj_well)),color=treatment)) +theme_bw() + 
  facet_grid(~treatment) + theme(strip.text.x = element_text(size = 8, colour = "black", angle = 90))+
  ggtitle("Total Integrated Intensity of nuclear pIRF3 \n in Primary Cortical Neurons\n(per cell per site per well)") + 
  ylim(0,NA) + ylab("Total Intensity of nuclear pIRF3 in Primary Cortical Neurons\n(per cell per site per well)") + xlab("") + theme(axis.title = element_text(size = 10)) + 
  stat_summary(fun.data=mean_se,geom="crossbar",fill="grey50",alpha=0.2) +
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = 'bold')) +
  geom_jitter(height=0,width=0.1,size=2,alpha=0.7) +
  theme(axis.text.x=element_text(angle=30,hjust=1)) +
  scale_color_manual(values=c("grey20","blue2")) 


####### pIRF3 nuclear AREA
##pIRF3 AREA

### QC (ALL neurons )
ggplot(Data_summary_clean_etop %>% filter(treatment %in% c('DMSO','etoposide')),
       aes(x=well,y=Total_nucl_area,color=treatment)) +theme_bw() + 
  # facet_grid(~treatment) + 
  theme(strip.text.x = element_text(size = 8, colour = "black", angle = 90))+
  ggtitle("QC: Area of nuclear pIRF3 in \n all cortical neurons (per site per well)") + 
  ylim(0,NA) + ylab("Area of nuclear pIRF3 in \n all cortical neurons (per site per well)") + xlab("") + theme(axis.title = element_text(size = 10)) + 
  stat_summary(fun.data=mean_sdl,geom="crossbar",fill="grey50",alpha=0.2) +
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = 'bold')) +
  geom_jitter(height=0,width=0.1,size=2,alpha=0.7) +
  theme(axis.text.x=element_text(angle=30,hjust=1)) +
  scale_color_manual(values=c("grey20","blue2")) 

Data_summary_clean_etop_area <- Data_summary_clean_etop %>% group_by(well) %>% filter(
  Total_nucl_area <= mean(Total_nucl_area) + 2*sd(Total_nucl_area))  %>% ungroup()

#pIRF3 nucl - area - PER site
ggplot(Data_summary_clean_etop_area %>% filter(treatment %in% c('DMSO','etoposide')) %>% group_by(treatment,genotype,line,well,duration) %>% summarize(
  area_well=mean(Total_nucl_area),ctip2_pos_well=mean(ctip2_pos_site),tuj_well=mean(Total_tuj_site),ctip2_neg_well=mean(ctip2_neg_site)),
  aes(x=duration,y=(area_well),color=treatment)) +theme_bw() +
  facet_grid(~treatment) + theme(strip.text.x = element_text(size = 8, colour = "black", angle = 90))+
  ggtitle("Average area of nuclear pIRF3 in \n Primary Cortical Neurons (per site per well)") +
  ylim(0,NA) + ylab("Average area of nuclear pIRF3 in \n Primary Cortical Neurons (per site per well)") + xlab("") + theme(axis.title = element_text(size = 10)) +
  stat_summary(fun.data=mean_se,geom="crossbar",fill="grey50",alpha=0.2) +
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = 'bold')) +
  geom_jitter(height=0,width=0.1,size=2,alpha=0.7) +
  theme(axis.text.x=element_text(angle=30,hjust=1)) +
  scale_color_manual(values=c("grey20","blue2"))

#pIRF3 nucl - area - PER CELL
ggplot(Data_summary_clean_etop_area %>% filter(treatment %in% c('DMSO','etoposide')) %>% group_by(treatment,genotype,line,well,duration) %>% summarize(
  area_well=mean(Total_nucl_area),ctip2_pos_well=mean(ctip2_pos_site),tuj_well=mean(Total_tuj_site),ctip2_neg_well=mean(ctip2_neg_site)),
  aes(x=duration,y=((area_well)/(tuj_well)),color=treatment)) +theme_bw() +
  facet_grid(~treatment) + theme(strip.text.x = element_text(size = 8, colour = "black", angle = 90))+
  ggtitle("Area of nuclear pIRF3 in \n Primary Cortical Neurons (per cell per site per well)") +
  ylim(0,0.4) + ylab("Area of nuclear pIRF3 in \n Primary Cortical Neurons (per cell per site per well)") + xlab("") + theme(axis.title = element_text(size = 10)) +
  stat_summary(fun.data=mean_se,geom="crossbar",fill="grey50",alpha=0.2) +
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = 'bold')) +
  geom_jitter(height=0,width=0.1,size=2,alpha=0.7) +
  theme(axis.text.x=element_text(angle=30,hjust=1)) +
  scale_color_manual(values=c("grey20","blue2"))

# run t-test etoposide/ALL
etop_ALLarea_ttest <- Data_summary_clean_etop_area %>% filter(treatment %in% c('DMSO','etoposide')) %>% filter(duration == "1H") %>% group_by(treatment,genotype,line,well,duration) %>% summarize(
  area_well=mean(Total_nucl_area),tuj_well=mean(Total_tuj_site),y = area_well/tuj_well)

t.test(y ~ treatment, data = etop_ALLarea_ttest)

# save csvs for normalization

Data_summary_clean_etop_area <- Data_summary_clean_etop_area %>% mutate(replicate = "1")
write.csv(Data_summary_clean_etop_area, "pIRF3-area-etoposide.csv")



# dev.off()



### between CTIP2+/- neurons t-tests
## comparing between ctip2+/- neurons
etop_POSarea_ttest$ctip2 = "+"
etop_POSarea_merge1 = etop_POSarea_ttest %>% select(treatment, genotype, line, well, duration, y, ctip2)
etop_NEGarea_ttest$ctip2 = "-"
etop_NEGarea_merge2 = etop_NEGarea_ttest %>% select(treatment, genotype, line, well, duration, y, ctip2)
etop_POSNEGarea_merge = rbind(etop_POSarea_merge1, etop_NEGarea_merge2)

t.test(y ~ ctip2, data = etop_POSNEGarea_merge %>% filter(treatment == "etoposide"))
t.test(y ~ ctip2, data = etop_POSNEGarea_merge %>% filter(treatment == "DMSO"))


## etoposide 
Data_summary_clean_etop_POSarea <- Data_summary_clean_etop_POSarea %>% mutate(replicate = "1",
                                                                              ctip2 = "pos")
write.csv(Data_summary_clean_etop_POSarea, "pIRF3-area-etop-pos.csv")

Data_summary_clean_etop_NEGarea <- Data_summary_clean_etop_NEGarea %>% mutate(replicate = "1",
                                                                              ctip2 = "neg")
write.csv(Data_summary_clean_etop_NEGarea, "pIRF3-area-etop-neg.csv")

