## This script loads, analyzes, and plots RNAseq differential expression

library(readxl) #loads library for reading excel files
library(tidyverse) #loads tidyverse library
library(DESeq2) #loads DESeq2
library(pheatmap)

setwd("C:/Data/Wainger_Lab/RNA/Maornof 2020")
excelsheet = "C:/Data/Wainger_Lab/RNA/Maornof 2020/Maornof_counts.xlsx"

#import dataset
genecounts <- read_excel(excelsheet, sheet =1, col_names = TRUE) #loads data from excel
genecounts <- as.data.frame(genecounts) #turns tibble output from above into data.frame
rownames(genecounts) <- genecounts[,1] #makes first column into rownames
genecounts$Gene_Name <- NULL #removes the Gene_name column (rownames have this info)

#make DESeqDataSet, needs matrix of conditions to organize comparisons
condition <- read_excel(excelsheet, sheet =2, col_names = TRUE) #loads data from excel
condition <-as.data.frame(condition)
rownames(condition) <- condition[,1] #makes first column into rownames
condition$Library <- NULL #removes the Gene_name column (rownames have this info)

dds <- DESeqDataSetFromMatrix(countData = genecounts, colData = condition, design = ~Genotype + Time + Genotype:Time)

dds_lrt <-DESeq(dds, test="LRT", reduced = ~Genotype + Time)
res <- results(dds_lrt)

write.csv(as.data.frame(res),file='C:/Data/Wainger_Lab/RNA/Maornof 2020/TDP43_DE.csv')


genecounts <- read_excel(excelsheet, sheet =3, col_names = TRUE) #loads data from excel
genecounts <- as.data.frame(genecounts) #turns tibble output from above into data.frame
rownames(genecounts) <- genecounts[,1] #makes first column into rownames
genecounts$Gene_Name <- NULL #removes the Gene_name column (rownames have this info)

#make DESeqDataSet, needs matrix of conditions to organize comparisons
condition <- read_excel(excelsheet, sheet =4, col_names = TRUE) #loads data from excel
condition <-as.data.frame(condition)
rownames(condition) <- condition[,1] #makes first column into rownames
condition$Library <- NULL #removes the Gene_name column (rownames have this info)

dds <- DESeqDataSetFromMatrix(countData = genecounts, colData = condition, design = ~Genotype + Time + Genotype:Time)

dds_lrt <-DESeq(dds, test="LRT", reduced = ~Genotype + Time)
res <- results(dds_lrt)

write.csv(as.data.frame(res),file='C:/Data/Wainger_Lab/RNA/Maornof 2020/PR_DE.csv')

## Analysis of NDE for Christine
#PR vs ctrl
excelsheet = "C:/Data/Wainger_Lab/RNA/Maornof 2020/Maornof_AH221020.xlsx"

genecounts <- read_excel(excelsheet, sheet =1, col_names = TRUE) #loads data from excel
genecounts <- as.data.frame(genecounts) #turns tibble output from above into data.frame
rownames(genecounts) <- genecounts[,1] #makes first column into rownames
genecounts$Gene_Name <- NULL #removes the Gene_name column (rownames have this info)

#make DESeqDataSet, needs matrix of conditions to organize comparisons
condition <- read_excel(excelsheet, sheet =2, col_names = TRUE) #loads data from excel
condition <-as.data.frame(condition)
rownames(condition) <- condition[,1] #makes first column into rownames
condition$Library <- NULL #removes the Gene_name column (rownames have this info)

dds <- DESeqDataSetFromMatrix(countData = genecounts, colData = condition, design = ~genotype)

#differential expression
de <- DESeq(dds)
res <- results(de)

write.csv(as.data.frame(res),file='C:/Data/Wainger_Lab/RNA/Maornof 2020/PR_DE2.csv')

#GFP 24hr vs 72hr
excelsheet = "C:/Data/Wainger_Lab/RNA/Maornof 2020/Maornof_AH221020.xlsx"

genecounts <- read_excel(excelsheet, sheet =3, col_names = TRUE) #loads data from excel
genecounts <- as.data.frame(genecounts) #turns tibble output from above into data.frame
rownames(genecounts) <- genecounts[,1] #makes first column into rownames
genecounts$Gene_Name <- NULL #removes the Gene_name column (rownames have this info)

#make DESeqDataSet, needs matrix of conditions to organize comparisons
condition <- read_excel(excelsheet, sheet =4, col_names = TRUE) #loads data from excel
condition <-as.data.frame(condition)
rownames(condition) <- condition[,1] #makes first column into rownames
condition$Library <- NULL #removes the Gene_name column (rownames have this info)

dds <- DESeqDataSetFromMatrix(countData = genecounts, colData = condition, design = ~Time)

#differential expression
de <- DESeq(dds)
res <- results(de)

write.csv(as.data.frame(res),file='C:/Data/Wainger_Lab/RNA/Maornof 2020/GFP_Time_DE.csv')

#GFP vs TDP43
genecounts <- read_excel(excelsheet, sheet =5, col_names = TRUE) #loads data from excel
genecounts <- as.data.frame(genecounts) #turns tibble output from above into data.frame
rownames(genecounts) <- genecounts[,1] #makes first column into rownames
genecounts$Gene_Name <- NULL #removes the Gene_name column (rownames have this info)

#make DESeqDataSet, needs matrix of conditions to organize comparisons
condition <- read_excel(excelsheet, sheet =6, col_names = TRUE) #loads data from excel
condition <-as.data.frame(condition)
rownames(condition) <- condition[,1] #makes first column into rownames
condition$Library <- NULL #removes the Gene_name column (rownames have this info)

dds <- DESeqDataSetFromMatrix(countData = genecounts, colData = condition, design = ~Genotype)

#differential expression
de <- DESeq(dds)
res <- results(de)

write.csv(as.data.frame(res),file='C:/Data/Wainger_Lab/RNA/Maornof 2020/GFPvsTDP43.csv')

#PR plot
genecounts <- read_excel(excelsheet, sheet =9, col_names = TRUE) #loads data from excel
genecounts <- as.data.frame(genecounts) #turns tibble output from above into data.frame
rownames(genecounts) <- genecounts[,1] #makes first column into rownames
genecounts$Gene_Name <- NULL #removes the Gene_name column (rownames have this info)

#make DESeqDataSet, needs matrix of conditions to organize comparisons
condition <- read_excel(excelsheet, sheet =10, col_names = TRUE) #loads data from excel
condition <-as.data.frame(condition)
rownames(condition) <- condition[,1] #makes first column into rownames
condition$Library <- NULL #removes the Gene_name column (rownames have this info)

dds <- DESeqDataSetFromMatrix(countData = genecounts, colData = condition, design = ~Genotype+Time)
rlogcounts = rlog(dds,blind=FALSE) #DESeq2 normalization that log2 transforms and scales for library size
plotPCA(rlogcounts, intgroup=c('Genotype','Time')) #plotting using plotPCA

#TDP43 plot
genecounts <- read_excel(excelsheet, sheet =11, col_names = TRUE) #loads data from excel
genecounts <- as.data.frame(genecounts) #turns tibble output from above into data.frame
rownames(genecounts) <- genecounts[,1] #makes first column into rownames
genecounts$Gene_Name <- NULL #removes the Gene_name column (rownames have this info)

#make DESeqDataSet, needs matrix of conditions to organize comparisons
condition <- read_excel(excelsheet, sheet =12, col_names = TRUE) #loads data from excel
condition <-as.data.frame(condition)
rownames(condition) <- condition[,1] #makes first column into rownames
condition$Library <- NULL #removes the Gene_name column (rownames have this info)

dds <- DESeqDataSetFromMatrix(countData = genecounts, colData = condition, design = ~Genotype+Time)
rlogcounts = rlog(dds,blind=FALSE) #DESeq2 normalization that log2 transforms and scales for library size
plotPCA(rlogcounts, intgroup=c('Genotype','Time')) #plotting using plotPCA
