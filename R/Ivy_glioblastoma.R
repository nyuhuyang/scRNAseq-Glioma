############################################
# prepare Ivy_glioblastoma reference
############################################
library(SingleR)
library(genefilter)
library(plyr)
source("../R/Seurat_functions.R")
#====== load Ivy_glioblastoma data  ==========================================
Ivy_gbm <- read.csv(file = "./data/Ref/gene_expression_matrix_2014-11-25/fpkm_table.csv",header = T)
dim(Ivy_gbm)
head(Ivy_gbm[,1:5])

# Convert gene id, log transfer
rows_genes <- read.csv(file = "./data/Ref/gene_expression_matrix_2014-11-25/rows-genes.csv")
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
test(Ivy_gbm)
Ivy_gbm_log <- log1p(Ivy_gbm)
test(Ivy_gbm_log)
boxplot(Ivy_gbm_log) # slow

# read columns-samples
columns_samples <- read.csv(file = "./data/Ref/gene_expression_matrix_2014-11-25/columns-samples.csv",
                            header = T)
columns_samples$rna_well_id = paste0("X",columns_samples$rna_well_id)
head(columns_samples[,1:10])
table(columns_samples$rna_well_id == colnames(Ivy_gbm_log))

#======== Create reference=====
CreateSinglerReference <- function(name,expr,types,main_types){
        name=name
        data = expr # the expression matrix
        types=types # a character list of the types. Samples from the same type should have the same name.
        main_types=main_types # a character list of the main types. 
        ref = list(name=name,data = expr, types=types, main_types=main_types)
        
        # if using the de method, we can predefine the variable genes
        ref$de.genes = CreateVariableGeneSet(expr,types,200)
        ref$de.genes.main = CreateVariableGeneSet(expr,main_types,300)
        
        # if using the sd method, we need to define an sd threshold
        sd = rowSds(expr)
        sd.thres = sort(sd, decreasing = T)[4000] # or any other threshold
        ref$sd.thres = sd.thres
        
        return(ref)
        
}

ref_IvyGbm <- CreateSinglerReference(name = 'Ivy_glioblastoma',
                              expr = as.matrix(Ivy_gbm_log), 
                              types = as.character(columns_samples$structure_abbreviation), 
                              main_types = as.character(columns_samples$structure_abbreviation))

save(ref_IvyGbm,file='./data/Ref/ref_IvyGbm.RData') 
