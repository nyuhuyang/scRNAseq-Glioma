############################################
# combine hpca and blueprint_encode
############################################
library(SingleR)
library(genefilter)
source("../R/Seurat_functions.R")
source("../R/SingleR_functions.R")
############################################
# gbm_eset
############################################
library(GSEABase)
library(GSVAdata)
library(SingleR)
data(gbm_VerhaakEtAl)

#------View data------
gbm_eset
head(featureNames(gbm_eset))
gbm_eset$subtype = paste0(gbm_eset$subtype,"_TCGA.GBM")
table(gbm_eset$subtype)
gbm_expr <- exprs(gbm_eset)
dim(gbm_expr)
head(gbm_expr[,1:5])
testMMM(gbm_expr)

#TPM
gbm_expr = TPM(2**(exprs(gbm_eset)), human_lengths)
head(gbm_expr[,1:5])
rownames(gbm_expr) = toupper(rownames(gbm_expr))
gbm_expr = log2(gbm_expr+1)
Median = apply(gbm_expr,2,median)
names(Median[Median < 3])
gbm_expr = gbm_expr[,Median > 3]
dim(gbm_expr)
gbm_expr[gbm_expr<0] = 0
par(mfrow=c(1,1))
boxplot(hpca_data, ylim=c(0,max(gbm_expr))) # slow
boxplot(gbm_expr, ylim=c(0,max(gbm_expr))) # slow
title(main = "Glioblastoma multiforme")

TCGA_gbm <- CreateSinglerReference(name = 'TCGA_GBM',
                                   expr = as.matrix(gbm_expr), 
                                   types = as.character(gbm_eset$subtype[Median > 3]), 
                                   main_types = as.character(gbm_eset$subtype[Median > 3]))
save(TCGA_gbm,file='./data/GeneSets/TCGA_gbm.RData') # it is best to name the object and the file with the same name.


############################################
# combine hpca, blueprint_encode, gbm_eset and IvyGbm
############################################
Iname <- load(file='../SingleR/data/ref_Human.RData');Iname
lname1 <- load(file='./data/GeneSets/TCGA_gbm.RData');lname1
Iname2 <- load(file='./data/GeneSets/Ivy_Glioblastoma.RData');Iname2
testMMM(TCGA_gbm$data)
testMMM(Ivy_Glioblastoma$data)
boxplot(Ivy_Glioblastoma$data) # slow
head(Ivy_Glioblastoma$data[,1:5])
refs_TCGA_IvyGbm <- merge(TCGA_gbm$data,Ivy_Glioblastoma$data,by="row.names",all=FALSE)
dim(refs_TCGA_IvyGbm)
rownames(refs_TCGA_IvyGbm) = refs_TCGA_IvyGbm$Row.names
refs_TCGA_IvyGbm <- refs_TCGA_IvyGbm[-which(colnames(refs_TCGA_IvyGbm)=="Row.names")]
head(refs_TCGA_IvyGbm[,1:5])

Refs_TCGA_IvyGbm <- merge(ref$data,refs_TCGA_IvyGbm,by="row.names",all=FALSE)
dim(Refs_TCGA_IvyGbm)
rownames(Refs_TCGA_IvyGbm) = Refs_TCGA_IvyGbm$Row.names
Refs_TCGA_IvyGbm <- Refs_TCGA_IvyGbm[-which(colnames(Refs_TCGA_IvyGbm)=="Row.names")]
head(Refs_TCGA_IvyGbm[,1:5])

par(mfrow = c(1,1))
boxplot(Refs_TCGA_IvyGbm$data,main="hpca, blueprint+encode, TCGA.BGM and Ivy.Glioblastoma") #slow!
Refs_TCGA_IvyGbm <- CreateSinglerReference(name = 'hpca_blueprint_encode_TCGA_Ivy',
                                   expr = as.matrix(Refs_TCGA_IvyGbm), 
                                   types = c(ref$types,TCGA_gbm$types,Ivy_Glioblastoma$types), 
                                   main_types = c(ref$main_types,TCGA_gbm$main_types,
                                                  Ivy_Glioblastoma$main_types),
                                   de.num = 200,de.main.num = 300)
save(Refs_TCGA_IvyGbm,file='./data/GeneSets/Refs_TCGA_IvyGbm.RData') # it is best to name the object and the file with the same name.

################################################
# extract marker gene list
################################################
In <- load(file='./data/GeneSets/Refs_TCGA_IvyGbm.RData');In # it is best to name the object and the file with the same name.
Refs_TCGA_IvyGbm_main = CreatGenePanelFromSingleR(object = Refs_TCGA_IvyGbm,
                                        main.type = TRUE, Human =TRUE)
write.csv(Refs_TCGA_IvyGbm_main,file = "./output/Refs_TCGA_IvyGbm_main.csv")
Refs_TCGA_IvyGbm_sub = CreatGenePanelFromSingleR(object = Refs_TCGA_IvyGbm,
                                       main.type = FALSE, Human =TRUE)
write.csv(Refs_TCGA_IvyGbm_sub,file = "./output/Refs_TCGA_IvyGbm_sub.csv")
