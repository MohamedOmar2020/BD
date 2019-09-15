##########################################################################
# Mohamed Omar
# 14/09/2019
## Goal: Validation of the Bipolar gene signature in independent data sets.
##########################################################################

rm(list = ls())

setwd("/Users/mohamedomar/Documents/Research/Projects/Bipolar")

library(MetaIntegrator)
library(GEOquery)
library(pROC)
library(caret)
library(enrichR)
library(ggplot2)
library(xtable)

# Load the 2 testing data sets
load("./Data/TestingData.rda")

Dataset8 <- TestingData$GSE5388
Dataset9 <- TestingData$GSE35977


## Get the expression matrix for each data set
Expr8 <- Dataset8$expr
Expr9 <- Dataset9$expr


## Get the Phenotype table for each data set
Pheno8 <- Dataset8$pheno
Pheno9 <- Dataset9$pheno

##################################################################
## Expression annotation 

# Annotate Expr8
head(rownames(Expr8))
rownames(Expr8) <- Dataset8$keys
summary(is.na(rownames(Expr8)))
summary(rownames(Expr8) == "")
Expr8 <- Expr8[!is.na(rownames(Expr8)), ]
dim(Expr8)
Dataset8$expr <- Expr8
# The same for the keys (Remove NAs)
Dataset8$keys <- Dataset8$keys[!is.na(Dataset8$keys)]

##################################################################
## Expression annotation 
# Annotate Expr9
head(rownames(Expr9))
rownames(Expr9) <- Dataset9$keys
summary(is.na(rownames(Expr9)))
summary(rownames(Expr9) == "")
Expr9 <- Expr9[!is.na(rownames(Expr9)), ]
dim(Expr9)
Dataset9$expr <- Expr9
# The same for the keys (Remove NAs)
Dataset9$keys <- Dataset9$keys[!is.na(Dataset9$keys)]


##################################################################
## Modify the phenotype tables

# Modify Pheno8
Pheno8$Disease_Status <- as.character(Pheno8$`Disease_status:ch1`)
table(Pheno8$Disease_Status)
Pheno8$Disease_Status[Pheno8$Disease_Status == "Bipolar disorder"] <- "BD"
Pheno8$Disease_Status[Pheno8$Disease_Status == "Healthy control"] <- "control"
Pheno8$Disease_Status <- as.factor(Pheno8$Disease_Status)
table(Pheno8$Disease_Status)

# Modify Expr8
Expr8 <- Expr8[, colnames(Expr8) %in% rownames(Pheno8)]
# Check if the sample names are consistent
all(rownames(Pheno8) == colnames(Expr8))
# Replace the old Pheno and EXpr with the new ones
Dataset8$pheno <- Pheno8
Dataset8$expr <- Expr8

#############################################
# Modify Pheno9
# Remove the SCV, depression and NA cases
Pheno9 <- Pheno9[!(Pheno9$`disease status:ch1` == "schizophrenia"), ]
Pheno9 <- Pheno9[!(Pheno9$`disease status:ch1` == "NA"), ]
Pheno9 <- Pheno9[!(Pheno9$`disease status:ch1` == "depression"), ]
Pheno9 <- Pheno9[!(Pheno9$`disease status:ch1` == "bipolar (not bipolar)"), ]
Pheno9$Disease_Status <- Pheno9$`disease status:ch1`
table(Pheno9$Disease_Status)
Pheno9$Disease_Status[Pheno9$Disease_Status == "bipolar"] <- "BD"
Pheno9$Disease_Status[Pheno9$Disease_Status == "unaffected"] <- "control"
Pheno9$Disease_Status <- as.factor(Pheno9$Disease_Status)
table(Pheno9$Disease_Status)

# Modify the Expr9
Expr9 <- Expr9[, colnames(Expr9) %in% rownames(Pheno9)]
# Check if the sample names are consistent
all(rownames(Pheno9) == colnames(Expr9))
# Replace the old Pheno and Expr with the new ones
Dataset9$pheno <- Pheno9
Dataset9$expr <- Expr9



#########################################################################
############################################################################

## Label samples (All samples need to be assigned labels in the $class vector, 1 for ‘disease’ or 0 for ‘control’)
Dataset8 <- classFunction(Dataset8, column = "Disease_Status", diseaseTerms = c("BD"))
Dataset9 <- classFunction(Dataset9, column = "Disease_Status", diseaseTerms = c("BD"))


## Keys
Dataset8$keys <- rownames(Expr8)
Dataset9$keys <- rownames(Expr9)

############################################################################
#### Validation

# Load the filter
load("./Objs/FilterResults.rda")

# ROC curves
png(filename = "./Figs/Val1AUC.png", width = 2000, height = 2000, res = 300)
rocPlot(datasetObject = Dataset8, filterObject = filter)
dev.off()

png(filename = "./Figs/Val2AUC.png", width = 2000, height = 2000, res = 300)
rocPlot(datasetObject = Dataset9, filterObject = filter)
dev.off()

# PRC curves
prcPlot(filterObject = filter, datasetObject = Dataset8)
prcPlot(filterObject = filter, datasetObject = Dataset9)

# Violin plots
violinPlot(filterObject = filter, datasetObject = Dataset8, labelColumn = "Disease_Status")
violinPlot(filterObject = filter, datasetObject = Dataset9, labelColumn = "Disease_Status")


## Calculate ROC data
scoreRs1 <- calculateScore(filterObject = filter, datasetObject = Dataset8)
rocRes1 <- calculateROC(predictions = scoreRs1, labels = Dataset8$class)


scoreRs2 <- calculateScore(filterObject = filter, datasetObject = Dataset9)
rocRes2 <- calculateROC(predictions = scoreRs2, labels = Dataset9$class)

###########################################################
## Calculation of the accuracy, sensitivity and specificity for each testing data set

## Calculate a signature score (Z score) and add it to the phenotype table
Pheno8$score <- calculateScore(filterObject = filter, datasetObject = Dataset8)
Pheno9$score <- calculateScore(filterObject = filter, datasetObject = Dataset9)

## Find the best threshold (for further use for the classification in the Test data)
thr_test1 <- coords(roc(Pheno8$Disease_Status, Pheno8$score, levels = c("control", "BD"),),"best", transpose = TRUE)["threshold"]
thr_test1


thr_test2 <- coords(roc(Pheno9$Disease_Status, Pheno9$score, levels = c("control", "BD"),),"best", transpose = TRUE)["threshold"]
thr_test2


### Get predictions and confusion matrix

# In Dataset8
Test_Predictions1 <- ifelse(Pheno8$score >= thr_test1, "BD", "control")
confusionMatrix(as.factor(Test_Predictions1), Pheno8$Disease_Status, positive = "BD")

# In Dataset9
Test_Predictions2 <- ifelse(Pheno9$score >= thr_test2, "BD", "control")
confusionMatrix(as.factor(Test_Predictions2), Pheno9$Disease_Status, positive = "BD")

###############################################################################
## Gene set enrichment analysis
dbs <- listEnrichrDbs()
dbs <- c("GO_Molecular_Function_2018", "GO_Cellular_Component_2018", "GO_Biological_Process_2018" , "ChEA_2016" ,"KEGG_2019_Human")
Enriched_PositiveGns <- enrichr(genes = PositiveGenes, databases = dbs)
Enriched_NegativeGns <- enrichr(genes = NegativeGenes, databases = dbs)
printEnrich(Enriched_PositiveGns, "PosGnsEnrichment.txt" , sep = "\t", columns = c(1:9))
printEnrich(Enriched_NegativeGns, "NegGnsEnrichment.txt" , sep = "\t", columns = c(1:9))

##############
Pos_GO_BP <- Enriched_PositiveGns["GO_Biological_Process_2018"]
Pos_GO_BP <- Pos_GO_BP$GO_Biological_Process_2018
Pos_GO_BP <- Pos_GO_BP[Pos_GO_BP$P.value <= 0.05, ]

write.csv(Pos_GO_BP, file = "./Objs/Pos_GO_BP.csv", quote = TRUE)
print(xtable(Pos_GO_BP, type = "latex"), file = "Pos_GO_BP")
################

Neg_GO_BP <- Enriched_NegativeGns["GO_Biological_Process_2018"]
Neg_GO_BP <- Neg_GO_BP$GO_Biological_Process_2018
Neg_GO_BP <- Neg_GO_BP[Neg_GO_BP$P.value <= 0.05, ]

write.csv(Neg_GO_BP, file = "./Objs/Neg_GO_BP.csv", quote = TRUE)
print(xtable(Neg_GO_BP, type = "latex"), file = "Neg_GO_BP")

###################

Pos_KEGG <- Enriched_PositiveGns["KEGG_2019_Human"]
Pos_KEGG <- Pos_KEGG$KEGG_2019_Human
Pos_KEGG <- Pos_KEGG[Pos_KEGG$P.value <= 0.05, ]

write.csv(Pos_KEGG, file = "./Objs/Pos_KEGG.csv", quote = TRUE)
print(xtable(Pos_KEGG, type = "latex"), file = "Pos_KEGG")

####################

Neg_KEGG <- Enriched_NegativeGns["KEGG_2019_Human"]
Neg_KEGG <- Neg_KEGG$KEGG_2019_Human
Neg_KEGG <- Neg_KEGG[Neg_KEGG$P.value <= 0.05, ]

write.csv(Neg_KEGG, file = "./Objs/Neg_KEGG.csv", quote = TRUE)
print(xtable(Neg_KEGG, type = "latex"), file = "Neg_KEGG")



################################################################################
################################################################################
