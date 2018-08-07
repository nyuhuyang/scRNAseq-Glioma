############################################
# prepare Ivy_glioblastoma reference
############################################
library(SingleR)
library(genefilter)
library(plyr)
library(dplyr)
source("../R/Seurat_functions.R")
#====== load Ivy_glioblastoma data  ==========================================
Ivy_gbm <- read.csv(file = "./data/GeneSets/gene_expression_matrix_2014-11-25/fpkm_table.csv",header = T)
dim(Ivy_gbm)
head(Ivy_gbm[,1:5])

# Convert gene id, log transfer
rows_genes <- read.csv(file = "./data/GeneSets/gene_expression_matrix_2014-11-25/rows-genes.csv")
head(rows_genes)
Ivy_gbm$gene_id.rna_well_id <- plyr::mapvalues(x = Ivy_gbm$gene_id.rna_well_id,
                                from = rows_genes$gene_id,
                                to = as.character(rows_genes$gene_symbol))
head(Ivy_gbm[,1:5])
rownames(Ivy_gbm) = Ivy_gbm$gene_id.rna_well_id
Ivy_gbm <- Ivy_gbm[-which(colnames(Ivy_gbm)=="gene_id.rna_well_id")]
head(Ivy_gbm[,1:5])
dim(Ivy_gbm)
boxplot(Ivy_gbm) # slow
testMMM(Ivy_gbm)
Ivy_gbm_log <- log1p(Ivy_gbm)
testMMM(Ivy_gbm_log)
boxplot(Ivy_gbm_log) # slow

# read columns-samples
columns_samples <- read.csv(file = "./data/GeneSets/gene_expression_matrix_2014-11-25/columns-samples.csv",
                            header = T)
columns_samples$rna_well_id = paste0("X",columns_samples$rna_well_id)
head(columns_samples[,1:10])
table(columns_samples$rna_well_id == colnames(Ivy_gbm_log))

# read columns-samples
rna_seq_samples <- read.csv(file = "./data/GeneSets/gene_expression_matrix_2014-11-25/rna_seq_samples_details.csv",
                            header = T)
# duplicated rows with merge by tumor_id, because there are duplicated tumor_id
apply(columns_samples,2,function(x) table(duplicated(x)))
apply(rna_seq_samples,2,function(x) table(duplicated(x)))
colnames(rna_seq_samples)[colnames(rna_seq_samples)=="sample_id"] = "rna_well_id"
rna_seq_samples$rna_well_id = paste0("X",rna_seq_samples$rna_well_id)
columns_rna_seq <- merge(columns_samples,rna_seq_samples, by = "rna_well_id")
dim(columns_rna_seq)
columns_rna_seq = columns_rna_seq[,c("rna_well_id","molecular_subtype")]
columns_rna_seq$molecular_subtype = as.character(columns_rna_seq$molecular_subtype)
#strsplit(columns_rna_seq$molecular_subtype, ", ")
#columns_rna_seq$main_types <- sapply(columns_rna_seq$molecular_subtype,
#                                     function(row) strsplit(row,", ")[[1]][1])
#columns_rna_seq$types <- gsub(", ", "_",columns_rna_seq$molecular_subtype)
types = columns_rna_seq[match(colnames(Ivy_gbm_log),columns_rna_seq$rna_well_id),"molecular_subtype"]
main_types = columns_rna_seq[match(colnames(Ivy_gbm_log),columns_rna_seq$rna_well_id),"molecular_subtype"]
#======== Create reference=====
Ivy_Glioblastoma <- CreateSinglerReference(name = 'Ivy_Glioblastoma',
                              expr = as.matrix(Ivy_gbm_log), 
                              types = paste0(types,"_Ivy.Glioblastoma"), 
                              main_types = paste0(main_types,"_Ivy.Glioblastoma"))

save(Ivy_Glioblastoma,file='./data/GeneSets/Ivy_Glioblastoma.RData') 
