############################################
# combine hpca and blueprint_encode
############################################
library(SingleR)
library(genefilter)
source("../R/Seurat_functions.R")
#-------skip if ref data was saved-----
data("HPCA")
data("Blueprint_Encode")
head(hpca$types)
head(hpca$main_types)
head(blueprint_encode$types)
dim(hpca$data)
dim(blueprint_encode$data)
head(hpca$data[,1:5])
head(blueprint_encode$data[,1:5])
blueprint_encode$data[is.na(blueprint_encode$data)] = 0
head(colSums(blueprint_encode$data))
boxplot(hpca$data)+title(main="hpca") #slow!
boxplot(blueprint_encode$data)+title(main="blueprint_encode")#slow!
boxplot(blueprint_encode$data[,1:100])#slow!

# remove low quanlity columns
par(mfrow=c(2,1))
hist(colMeans(blueprint_encode$data),breaks=ncol(blueprint_encode$data))
quantile_75 <- apply(blueprint_encode$data,2,function(x) quantile(x,probs =0.75))
hist(quantile_75, breaks=ncol(blueprint_encode$data))
rm_samples <- names(quantile_75)[quantile_75<1]
rm_index <- which(colnames(blueprint_encode$data) %in% rm_samples)
blueprint_encode_rm <- blueprint_encode$data[,-rm_index]
par(mfrow=c(1,1))
boxplot(blueprint_encode_rm)#slow!
title(main="blueprint_encode_rm")

# adjest median
test <- function(x) {
        par(mfrow=c(3,1))
        Max <- apply(x,2,max)
        Median <- apply(x,2,median)
        Min <- apply(x,2,min)
        xlim = c(min(x),max(x))
        hist(Max, xlim = xlim)
        hist(Median, xlim = xlim)
        hist(Min, xlim = xlim)
        par(mfrow=c(1,1))
}

#hpca_median <- mean(apply(hpca$data,2,median))
# make all reads positive-------
hpca_data <- scale(hpca$data) - min(scale(hpca$data)) 
test(hpca$data)
test(hpca_data)
boxplot(hpca_data) #slow!

#blue_encode_median <- mean(apply(blueprint_encode_rm,2,median))
# make all reads positive-------
#blue_encode_data <- scale(blueprint_encode_rm) + blue_encode_median
blue_encode_data <- scale(blueprint_encode_rm) - min(scale(blueprint_encode_rm))
#blue_encode_min <- apply(blue_encode_data,2,min)
#blue_encode_data <- sweep(blue_encode_data,
#                          2, blue_encode_min) # by column
#min(blue_encode_data)
test(scale(blueprint_encode_rm))
test(blue_encode_data)
boxplot(blue_encode_data) #slow!

# merge
hpca_blue_encode <- merge(hpca_data,blue_encode_data
                               ,by="row.names",all=FALSE)
rownames(hpca_blue_encode) = hpca_blue_encode$Row.names
hpca_blue_encode <- hpca_blue_encode[-which(colnames(hpca_blue_encode)=="Row.names")]
head(hpca_blue_encode[,1:5])
dim(hpca_blue_encode)

# scale
boxplot(hpca_blue_encode)#slow!
title(main='hpca_blueprint_encode')
#allmean <- mean(apply(hpca_blue_encode,2,median))  
#hpca_blue_encode <- scale(hpca_blue_encode)
#hpca_blue_encode <- hpca_blue_encode +allmean

name = 'hpca_blueprint_encode'
expr = as.matrix(hpca_blue_encode) # the expression matrix
types = as.character(c(hpca$types,blueprint_encode$types[-rm_index])) # a character list of the types. Samples from the same type should have the same name.
main_types = as.character(c(hpca$main_types,blueprint_encode$main_types[-rm_index])) # a character list of the main types. 

# Fine tune types
FineTune <- function(x){
        x = gsub(" ","_",x)
        x = gsub("B-cells","B_cells",x)
        x = gsub("CD4\\+_T-cells","T_cells",x)
        x = gsub("CD8\\+_T-cells","T_cells",x)
        x = gsub("CD4\\+_Tcm","T_cells:CD4+_central_memory",x)
        x = gsub("CD4\\+_Tem","T_cells:CD4+_effector_memory",x)
        x = gsub("CD8\\+_Tcm","T_cells:CD8+_Central_memory",x)
        x = gsub("CD8\\+_Tem","T_cells:CD8+_effector_memory",x)
        x = gsub("cells","cell",x)
        x = gsub("cell","cells",x)
        x = gsub("Class-switched_memory_B_cells","B_cells:Class-switched_memory",x)
        x = gsub("Macrophages","Macrophage",x)
        x = gsub("Memory_B_cells","B_cells:Memory",x)
        x = gsub("Monocytes","Monocyte",x)
        x = gsub("Neutrophils","Neutrophil",x)
        x = gsub("Smooth_muscle_cells","Smooth_muscle",x)
        return(x)
}
types = FineTune(types)
main_types <- FineTune(main_types)
ref = list(name=name,data = expr, types=types, main_types=main_types)

# if using the de method, we can predefine the variable genes
ref$de.genes = CreateVariableGeneSet(expr,types,200)
ref$de.genes.main = CreateVariableGeneSet(expr,main_types,300)

# if using the sd method, we need to define an sd threshold
sd = rowSds(expr)
sd.thres = sort(sd, decreasing = T)[4000] # or any other threshold
ref$sd.thres = sd.thres

save(ref,file='../SingleR/data/ref.RData') # it is best to name the object and the file with the same name.


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
gbm_eset$subtype = paste0(gbm_eset$subtype,"_Glioblastoma")
table(gbm_eset$subtype)
gbm_expr <- exprs(gbm_eset)
dim(gbm_expr)
head(gbm_expr[,1:5])
test(gbm_expr)

#TPM
gbm_expr = TPM(exp(exprs(gbm_eset)), human_lengths)
head(gbm_expr[,1:5])
rownames(gbm_expr) = toupper(rownames(gbm_expr))
gbm_expr = log1p(gbm_expr)
Median = apply(gbm_expr,2,median)
names(Median[Median < 2])
gbm_expr = gbm_expr[,Median > 2]
dim(gbm_expr)
#allmedian <- mean(colMeans(gbm_expr)) 
# make all reads positive-------
gbm_expr1 = scale(gbm_expr) - min(scale(gbm_expr))
head(gbm_expr1[,1:5])
par(mfrow=c(1,2))
boxplot(gbm_expr) # slow
boxplot(gbm_expr1) # slow
test(gbm_expr1)

############################################
# combine hpca, blueprint_encode and gbm_eset
############################################
Iname <- load(file='../SingleR/data/ref.RData')
Iname
head(ref$data[,1:5])
test(ref$data)
boxplot(ref$data) # slow
Refs_gbm <- merge(ref$data,gbm_expr1,by="row.names",all=FALSE)
dim(Refs_gbm)
rownames(Refs_gbm) = Refs_gbm$Row.names
Refs_gbm <- Refs_gbm[-which(colnames(Refs_gbm)=="Row.names")]
head(Refs_gbm[,1:5])
dim(Refs_gbm)
ncol(ref$data)+ncol(gbm_expr)

# scale
#allmean <- mean(apply(Refs_gbm,2,median))
#Refs_gbm <- scale(Refs_gbm)
#Refs_gbm <- Refs_gbm +allmean

boxplot(Refs_gbm)
test(Refs_gbm)

name = 'hpca_blueprint_encode_gbm'
expr = as.matrix(Refs_gbm) # the expression matrix
types = c(ref$types,as.character(gbm_eset$subtype[Median > 2])) # a character list of the types. Samples from the same type should have the same name.
main_types = c(ref$main_types,
               as.character(gbm_eset$subtype[Median > 2])) # a character list of the main types. 
Refs_gbm = list(name=name,data = expr, types=types, main_types=main_types)

# if using the de method, we can predefine the variable genes
Refs_gbm$de.genes = CreateVariableGeneSet(expr,types,200)
Refs_gbm$de.genes.main = CreateVariableGeneSet(expr,main_types,300)

# if using the sd method, we need to define an sd threshold
sd = rowSds(expr)
sd.thres = sort(sd, decreasing = T)[4000] # or any other threshold
Refs_gbm$sd.thres = sd.thres

save(Refs_gbm,file='./data/Ref/Refs_gbm.RData') # it is best to name the object and the file with the same name.


#---- Skip ------
expr <- Refs_gbm
dim(expr)
head(expr)
tail(colnames(expr),10)
boxplot(expr) # optional
expr = expr - min(expr)
name = "gbm_VerhaakEtAl"
types = gbm_eset$subtype
main_types = as.character(types) # a character list of the main types. 
Refs_gbm_VerhaakEtAl <- list(name=name,data = expr, types=types, main_types=main_types)
# if using the de method, we can predefine the variable genes
Refs_gbm_VerhaakEtAl$de.genes.main = CreateVariableGeneSet(expr,main_types,300)

# if using the sd method, we need to define an sd threshold
sd = rowSds(expr)
sd.thres = sort(sd, decreasing = T)[4000] # or any other threshold
Refs_gbm_VerhaakEtAl$sd.thres = sd.thres

save(Refs_gbm_VerhaakEtAl,file='./data/GeneSets/Refs_gbm_VerhaakEtAl.RData') 
# it is best to name the object and the file with the same name.

#==== brainTxDbSets======
data(brainTxDbSets)
brainTxDbSets
sapply(brainTxDbSets, length)

lapply(brainTxDbSets, head)
gbm_es <- gsva(gbm_eset, brainTxDbSets, mx.diff=FALSE, verbose=FALSE, parallel.sz=1)
