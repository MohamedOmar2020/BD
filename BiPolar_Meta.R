################################################################
## Mohamed Omar
## 13/09/2019
## Detecting differentially expressed genes between Bipolar disorder patients and normal individuals (Postmortem brain)
#################################################################

rm(list = ls())

setwd("/Users/mohamedomar/Documents/Research/Projects/Bipolar")

library(RMySQL)
library(MetaIntegrator)
library(GEOquery)
library(limma)

#########################################################
## Get the data sets from GEO

#Datasets <- getGEOData(c("GSE120340", "GSE87610", "GSE92538", "GSE78246", "GSE62191", "GSE35977", "GSE12654"))

#save(Datasets, file = "./Data/TrainingData.rda")
#load("./Data/TrainingData.rda")


# Dataset1 <- Datasets$originalData$GSE120340
# Dataset2 <- Datasets$originalData$GSE87610
# Dataset3 <- NewData$originalData$GSE12654
# Dataset4 <- NewData$originalData$GSE12649
# Dataset5 <- Datasets$originalData$GSE62191
# 
# Dataset6 <- Datasets$originalData$GSE92538_GPL10526
# Dataset7 <- Datasets$originalData$GSE92538_GPL17027
# Dataset8 <- NewData$originalData$GSE5388
# Dataset9 <- Datasets$originalData$GSE35977

# TrainingData <- list(Dataset1, Dataset2, Dataset3, Dataset4, Dataset5, Dataset6, Dataset7)
# names(TrainingData) <- c(Dataset1$formattedName, Dataset2$formattedName, Dataset3$formattedName, Dataset4$formattedName, Dataset5$formattedName, Dataset6$formattedName, Dataset7$formattedName)
# 
# TestingData <- list(Dataset8, Dataset9)
# names(TestingData) <- c(Dataset8$formattedName, Dataset9$formattedName)


## Training: GSE120340 - GSE87610 - GSE12654 - GSE12649 - GSE62191 - GSE92538_GPL10526 - GSE92538_GPL17027.
## Testing: GSE5388  - GSE35977
#save(TrainingData, file = "./Data/TrainingData.rda")
#save(TestingData, file = "./Data/TestingData.rda")

# Load the 7 training data sets
load("./Data/TrainingData.rda")

Dataset1 <- TrainingData$GSE120340
Dataset2 <- TrainingData$GSE87610
Dataset3 <- TrainingData$GSE12654
Dataset4 <- TrainingData$GSE12649
Dataset5 <- TrainingData$GSE62191
Dataset6 <- TrainingData$`GSE92538 GPL10526`
Dataset7 <- TrainingData$`GSE92538 GPL17027`

## Get the expression matrix for each data set
Expr1 <- Dataset1$expr
Expr2 <- Dataset2$expr
Expr3 <- Dataset3$expr
Expr4 <- Dataset4$expr
Expr5 <- Dataset5$expr
Expr6 <- Dataset6$expr
Expr7 <- Dataset7$expr


## Get the Phenotype table for each data set
Pheno1 <- Dataset1$pheno
Pheno2 <- Dataset2$pheno
Pheno3 <- Dataset3$pheno
Pheno4 <- Dataset4$pheno
Pheno5 <- Dataset5$pheno
Pheno6 <- Dataset6$pheno
Pheno7 <- Dataset7$pheno


##################################################################
## Expression annotation 
# Annotate Expr1
head(rownames(Expr1))
rownames(Expr1) <- Dataset1$keys
summary(is.na(rownames(Expr1)))
summary(rownames(Expr1) == "")
Expr1 <- Expr1[!is.na(rownames(Expr1)), ]
dim(Expr1)
Dataset1$expr <- Expr1
# The same for the keys (Remove NAs)
Dataset1$keys <- Dataset1$keys[!is.na(Dataset1$keys)]

###########
# Annotate Expr2
head(rownames(Expr2))
rownames(Expr2) <- Dataset2$keys
summary(is.na(rownames(Expr2)))
summary(rownames(Expr2) == "")
Expr2 <- Expr2[!is.na(rownames(Expr2)), ]
dim(Expr2)
Dataset2$expr <- Expr2
# The same for the keys
Dataset2$keys <- Dataset2$keys[!is.na(Dataset2$keys)]

###########
# Annotate Expr3
head(rownames(Expr3))
rownames(Expr3) <- Dataset3$keys
summary(is.na(rownames(Expr3)))
Expr3 <- Expr3[!is.na(rownames(Expr3)), ]
dim(Expr3)                                    # Small number of genes !!
Dataset3$expr <- Expr3
# The same for the keys
Dataset3$keys <- Dataset3$keys[!is.na(Dataset3$keys)]

###########
# Annotate Expr4
head(rownames(Expr4))
rownames(Expr4) <- Dataset4$keys
summary(is.na(rownames(Expr4)))
Expr4 <- Expr4[!is.na(rownames(Expr4)), ]
dim(Expr4)                                 
Dataset4$expr <- Expr4
# The same for the keys
Dataset4$keys <- Dataset4$keys[!is.na(Dataset4$keys)]

###############
# Annotate Expr5
head(rownames(Expr5))
rownames(Expr5) <- Dataset5$keys
summary(is.na(rownames(Expr5)))
Expr5 <- Expr5[!is.na(rownames(Expr5)), ]
dim(Expr5)
Dataset5$expr <- Expr5
# The same for the keys
Dataset5$keys <- Dataset5$keys[!is.na(Dataset5$keys)]

#######################
# Annotate Expr6
head(rownames(Expr6))
rownames(Expr6) <- Dataset6$keys
summary(is.na(rownames(Expr6)))
summary(rownames(Expr6) == "")
Expr6 <- Expr6[!is.na(rownames(Expr6)), ]
dim(Expr6)
Dataset6$expr <- Expr6
# The same for the keys (Remove NAs)
Dataset6$keys <- Dataset6$keys[!is.na(Dataset6$keys)]

#########################
# Annotate Expr7
head(rownames(Expr7))
rownames(Expr7) <- Dataset7$keys
summary(is.na(rownames(Expr7)))
summary(rownames(Expr7) == "")
Expr7 <- Expr7[!is.na(rownames(Expr7)), ]
dim(Expr7)
Dataset7$expr <- Expr7
# The same for the keys (Remove NAs)
Dataset7$keys <- Dataset7$keys[!is.na(Dataset7$keys)]

##################################################################
## Modify the phenotype tables

# Modify Pheno1
# Remove the SCV cases
Pheno1 <- Pheno1[!(Pheno1$`disease state:ch1` == "SCZ"), ]
Pheno1$Disease_Status <- Pheno1$`disease state:ch1`
table(Pheno1$Disease_Status)
Pheno1$Disease_Status[Pheno1$Disease_Status == "BD(-)"] <- "BD"
Pheno1$Disease_Status[Pheno1$Disease_Status == "BD(+)"] <- "BD"
Pheno1$Disease_Status <- as.factor(Pheno1$Disease_Status)
table(Pheno1$Disease_Status)

# Modify the Expr1
Expr1 <- Expr1[, colnames(Expr1) %in% rownames(Pheno1)]
# Check if the sample names are consistent
all(rownames(Pheno1) == colnames(Expr1))
# Replace the old Pheno1 and EXpr1 with the new ones
Dataset1$pheno <- Pheno1
Dataset1$expr <- Expr1

###############################
# Modify Pheno2
# Remove schezo and MDD
Pheno2 <- Pheno2[!(Pheno2$`genotype:ch1` == "Major Depressive Disorder"), ]
Pheno2 <- Pheno2[!(Pheno2$`genotype:ch1` == "Schizophrenia"), ]

Pheno2$Disease_Status <- Pheno2$`genotype:ch1`
Pheno2$Disease_Status[Pheno2$Disease_Status == "Bipolar"] <- "BD"
Pheno2$Disease_Status[Pheno2$Disease_Status == "Unaffected Comparison subject"] <- "control"
Pheno2$Disease_Status <- as.factor(Pheno2$Disease_Status)
table(Pheno2$Disease_Status)

# Modify the Expr2
Expr2 <- Expr2[, colnames(Expr2) %in% rownames(Pheno2)]
# Check if the sample names are consistent
all(rownames(Pheno2) == colnames(Expr2))
# Replace the old Pheno and Expr with the new ones
Dataset2$pheno <- Pheno2
Dataset2$expr <- Expr2

############################
# Modify Pheno3
# Remove Schezo and MDD
Pheno3 <- Pheno3[!(Pheno3$source_name_ch1 == "schizophrenia"), ]
Pheno3 <- Pheno3[!(Pheno3$source_name_ch1 == "depression"), ]

Pheno3$Disease_Status <- as.character(Pheno3$source_name_ch1)
Pheno3$Disease_Status[Pheno3$Disease_Status == "bipolar disorder"] <- "BD"
Pheno3$Disease_Status <- as.factor(Pheno3$Disease_Status)
table(Pheno3$Disease_Status)

# Modify the Expr3
Expr3 <- Expr3[, colnames(Expr3) %in% rownames(Pheno3)]
# Check if the sample names are consistent
all(rownames(Pheno3) == colnames(Expr3))
# Replace the old Pheno and Expr with the new ones
Dataset3$pheno <- Pheno3
Dataset3$expr <- Expr3

#############################
# Modify Pheno4
# Remove Schezo and MDD
Pheno4 <- Pheno4[!(Pheno4$characteristics_ch1.1 == "schizophrenia"), ]

Pheno4$Disease_Status <- as.character(Pheno4$characteristics_ch1.1)
Pheno4$Disease_Status[Pheno4$Disease_Status == "bipolar disorder"] <- "BD"
Pheno4$Disease_Status <- as.factor(Pheno4$Disease_Status)
table(Pheno4$Disease_Status)

# Modify the Expr4
Expr4 <- Expr4[, colnames(Expr4) %in% rownames(Pheno4)]
# Check if the sample names are consistent
all(rownames(Pheno4) == colnames(Expr4))
# Replace the old Pheno and Expr with the new ones
Dataset4$pheno <- Pheno4
Dataset4$expr <- Expr4

##############################
# Modify Pheno5
# Remove Shezo
Pheno5 <- Pheno5[!(Pheno5$`disease state:ch1` == "schizophrenia"), ]

Pheno5$Disease_Status <- Pheno5$`disease state:ch1`
Pheno5$Disease_Status[Pheno5$Disease_Status == "bipolar disorder"] <- "BD"
Pheno5$Disease_Status[Pheno5$Disease_Status == "healthy control"] <- "control"
Pheno5$Disease_Status <- as.factor(Pheno5$Disease_Status)
table(Pheno5$Disease_Status)

# Modify the Expr5
Expr5 <- Expr5[, colnames(Expr5) %in% rownames(Pheno5)]
# Check if the sample names are consistent
all(rownames(Pheno5) == colnames(Expr5))
# Replace the old Pheno and Expr with the new ones
Dataset5$pheno <- Pheno5
Dataset5$expr <- Expr5

##################################
# Modify Pheno6
# Remove the SCV cases
Pheno6 <- Pheno6[!(Pheno6$`diagnosis:ch1` == "Schizophrenia"), ]
Pheno6 <- Pheno6[!(Pheno6$`diagnosis:ch1` == "Major Depressive Disorder"), ]
Pheno6$Disease_Status <- as.character(Pheno6$`diagnosis:ch1`)
table(Pheno6$Disease_Status)
Pheno6$Disease_Status[Pheno6$Disease_Status == "Bipolar Disorder"] <- "BD"
Pheno6$Disease_Status[Pheno6$Disease_Status == "Control"] <- "control"
Pheno6$Disease_Status <- as.factor(Pheno6$Disease_Status)
table(Pheno6$Disease_Status)

# Modify Expr6
Expr6 <- Expr6[, colnames(Expr6) %in% rownames(Pheno6)]
# Check if the sample names are consistent
all(rownames(Pheno6) == colnames(Expr6))
# Replace the old Pheno and EXpr with the new ones
Dataset6$pheno <- Pheno6
Dataset6$expr <- Expr6

#######################################
# Modify Pheno7
# Remove the SCV cases
Pheno7 <- Pheno7[!(Pheno7$`diagnosis:ch1` == "Schizophrenia"), ]
Pheno7 <- Pheno7[!(Pheno7$`diagnosis:ch1` == "Major Depressive Disorder"), ]
Pheno7$Disease_Status <- as.character(Pheno7$`diagnosis:ch1`)
table(Pheno7$Disease_Status)
Pheno7$Disease_Status[Pheno7$Disease_Status == "Bipolar Disorder"] <- "BD"
Pheno7$Disease_Status[Pheno7$Disease_Status == "Control"] <- "control"
Pheno7$Disease_Status <- as.factor(Pheno7$Disease_Status)
table(Pheno7$Disease_Status)

# Modify Expr7
Expr7 <- Expr7[, colnames(Expr7) %in% rownames(Pheno7)]
# Check if the sample names are consistent
all(rownames(Pheno7) == colnames(Expr7))
# Replace the old Pheno and EXpr with the new ones
Dataset7$pheno <- Pheno7
Dataset7$expr <- Expr7


#########################################################################
############################################################################

## Label samples (All samples need to be assigned labels in the $class vector, 1 for ‘disease’ or 0 for ‘control’)
Dataset1 <- classFunction(Dataset1, column = "Disease_Status", diseaseTerms = c("BD"))
Dataset2 <- classFunction(Dataset2, column = "Disease_Status", diseaseTerms = c("BD"))
Dataset3 <- classFunction(Dataset3, column = "Disease_Status", diseaseTerms = c("BD"))
Dataset4 <- classFunction(Dataset4, column = "Disease_Status", diseaseTerms = c("BD"))
Dataset5 <- classFunction(Dataset5, column = "Disease_Status", diseaseTerms = c("BD"))
Dataset6 <- classFunction(Dataset6, column = "Disease_Status", diseaseTerms = c("BD"))
Dataset7 <- classFunction(Dataset7, column = "Disease_Status", diseaseTerms = c("BD"))


############################################################################
############################################################################
#########################################################################
## The metaanalysis

## Creating the meta object
AllDataSets <- list(Dataset1, Dataset2, Dataset3, Dataset4, Dataset5, Dataset6, Dataset7)
names(AllDataSets) <- c(Dataset1$formattedName, Dataset2$formattedName, Dataset3$formattedName, Dataset4$formattedName, Dataset5$formattedName, Dataset6$formattedName, Dataset7$formattedName)

Bipolar_meta <- list()
Bipolar_meta$originalData <- AllDataSets

## Replace keys within each data set
Bipolar_meta$originalData$GSE120340$keys <- rownames(Expr1)
Bipolar_meta$originalData$GSE87610$keys <- rownames(Expr2)
Bipolar_meta$originalData$GSE12654$keys <- rownames(Expr3)
Bipolar_meta$originalData$GSE12649$keys <- rownames(Expr4)
Bipolar_meta$originalData$GSE62191$keys <- rownames(Expr5)
Bipolar_meta$originalData$`GSE92538 GPL10526`$keys <- rownames(Expr6)
Bipolar_meta$originalData$`GSE92538 GPL17027`$keys <- rownames(Expr7)

## Check the meta object before the metaanalysis
checkDataObject(Bipolar_meta, "Meta", "Pre-Analysis") ## If true, Proceed to the meta analysis


## Run the meta analysis
Bipolar_metaanalysis <- runMetaAnalysis(Bipolar_meta, runLeaveOneOutAnalysis = F, maxCores = 3)

## Filter out significant genes from the metaanalysis results (this will be the gene signature that separates Metas from No_Mets)
Bipolar_metaanalysis <- filterGenes(Bipolar_metaanalysis, isLeaveOneOut = F, effectSizeThresh = 0, FDRThresh = 0.15, numberStudiesThresh = 7, heterogeneityPvalThresh = 0.05)

## Assigning a name to the filter
filter <- Bipolar_metaanalysis$filterResults[[1]]

## Summarize filter results
filter_summary <- summarizeFilterResults(metaObject = Bipolar_metaanalysis, getMostRecentFilter(Bipolar_metaanalysis))

## Gene names
PositiveGenes <- filter$posGeneNames
PositiveGenes

NegativeGenes <- filter$negGeneNames
NegativeGenes

# ## Save the tables of positive and negative genes
write.csv(filter_summary$pos, file = "./Objs/Filter_Positive_genes.csv", quote = T)
write.csv(filter_summary$neg, file = "./Objs/Filter_Negative_genes.csv", quote = T)

## Summary ROC curve in all training data sets using the original filter
set.seed(333)
png(filename = "./Figs/SummaryAUC_Training.png", width = 4000, height = 3000, res = 300)
summaryROCPlot(Bipolar_metaanalysis, filter, orderByAUC = FALSE)
dev.off()

## Save the filter
save(filter, filter_summary, PositiveGenes, NegativeGenes, file = "./Objs/FilterResults.rda")

## Modify the gene signature for more accuracy and AUC
# Using forward search 
# New_filter <- forwardSearch(metaObject = Bipolar_metaanalysis, filterObject = filter)
# 
# ## Replace the old filter with the new smaller one
# Bipolar_metaanalysis$filterResults$FDR0.15_es0_nStudies7_looaFALSE_hetero0.05$posGeneNames <- New_filter$posGeneNames
# Bipolar_metaanalysis$filterResults$FDR0.15_es0_nStudies7_looaFALSE_hetero0.05$negGeneNames <- New_filter$negGeneNames
# 
# New_filter <- Bipolar_metaanalysis$filterResults[[1]]
# New_filter_summary <- summarizeFilterResults(metaObject = Bipolar_metaanalysis, getMostRecentFilter(Bipolar_metaanalysis))
# 
# ## Save the tables of positive and negative genes
# write.table(New_filter_summary$pos, file = "./Objs/NewFilter_Positive_genes.csv", quote = T, sep = "\t", col.names = T, row.names = F)
# write.table(New_filter_summary$neg, file = "./Objs/NewFilter_Negative_genes.csv", quote = T, sep = "\t", col.names = T, row.names = F)
# 
# ## Gene names
# PositiveGenes <- New_filter$posGeneNames
# NegativeGenes <- New_filter$negGeneNames
# 
# 
# ## Summary ROC curve in all training data sets using the new filter
# set.seed(333)
# summaryROCPlot(Bipolar_metaanalysis, New_filter)


## Effect size visualization in each study
# For the up-regulated genes
pdf(file = "./Figs/UpGenes_EffectSize.pdf", title = "Effect size of up-regulated genes across the discovery data sets", width = 20, height = 10)
par(mfrow=c(3,5))

forestPlot(metaObject = Bipolar_metaanalysis, geneName = "MYO6", textColor = "black")
forestPlot(metaObject = Bipolar_metaanalysis, geneName = "ERCC8", textColor = "black")
forestPlot(metaObject = Bipolar_metaanalysis, geneName = "SLC12A4", textColor = "black")
forestPlot(metaObject = Bipolar_metaanalysis, geneName = "CA12", textColor = "black")
forestPlot(metaObject = Bipolar_metaanalysis, geneName = "ZNF273", textColor = "black")
forestPlot(metaObject = Bipolar_metaanalysis, geneName = "TP53", textColor = "black")
forestPlot(metaObject = Bipolar_metaanalysis, geneName = "MYRIP", textColor = "black")
forestPlot(metaObject = Bipolar_metaanalysis, geneName = "FGF2", textColor = "black")
forestPlot(metaObject = Bipolar_metaanalysis, geneName = "BCL2L11", textColor = "black")
forestPlot(metaObject = Bipolar_metaanalysis, geneName = "COL11A1", textColor = "black")
forestPlot(metaObject = Bipolar_metaanalysis, geneName = "PRKG1", textColor = "black")
forestPlot(metaObject = Bipolar_metaanalysis, geneName = "RRM1", textColor = "black")
forestPlot(metaObject = Bipolar_metaanalysis, geneName = "PAQR3", textColor = "black")
dev.off()

# For the down-regulated genes
pdf(file = "./Figs/DownGenes_effect_size.pdf", title = "Effect size of Down-regulated genes across the discovery data sets", width = 20, height = 12)
par(mfrow=c(2,5))

forestPlot(metaObject = Bipolar_metaanalysis, geneName = "CD74", textColor = "black")
forestPlot(metaObject = Bipolar_metaanalysis, geneName = "CNP", textColor = "black")
forestPlot(metaObject = Bipolar_metaanalysis, geneName = "HLA-DPB1", textColor = "black")
forestPlot(metaObject = Bipolar_metaanalysis, geneName = "DDX17", textColor = "black")
forestPlot(metaObject = Bipolar_metaanalysis, geneName = "HLA-DPA1", textColor = "black")
forestPlot(metaObject = Bipolar_metaanalysis, geneName = "MRPL33", textColor = "black")
forestPlot(metaObject = Bipolar_metaanalysis, geneName = "CCKBR", textColor = "black")
forestPlot(metaObject = Bipolar_metaanalysis, geneName = "ACAT2", textColor = "black")
forestPlot(metaObject = Bipolar_metaanalysis, geneName = "CD53", textColor = "black")
dev.off()

## Heatmap of the effect sizes of the signature genes
png("./Figs/Heatmap_EffectSizes.png", width = 4000, height = 2000, res = 300)
heatmapPlot(Bipolar_metaanalysis, filterObject = filter)
dev.off()
