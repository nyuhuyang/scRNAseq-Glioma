# load and split Seurat================
library(Seurat) #v2.3.4 main
library(magrittr)
library(velocyto.R)
library(ggplot2)
source("../R/Seurat_functions.R")
(load(file = "./data/Glioma_Harmony_20181201.Rda"))
Glioma@meta.data$group = gsub("\\-.*","",Glioma@meta.data$orig.ident)
table(Glioma@meta.data$group)
Glioma %<>% SetAllIdent("group")
TSNEPlot(Glioma,do.label = T)
sub_object <- SubsetData(Glioma, ident.use = "PM1005")

sub_object %<>% SetAllIdent("singler2sub")
cell.colors <- sub_object@meta.data$singler2sub.colors
names(cell.colors) = rownames(sub_object@meta.data)
TSNEPlot.1(sub_object,colors.use = cell.colors, do.label = T)
ldat <- read.loom.matrices("data/velocyto/PM1005_merged.loom")

# extract tsne
emb <- sub_object@dr$tsne@cell.embeddings[,1:2]
row.names(emb) <- stringi::stri_replace_last_fixed(row.names(emb),"_",":")
row.names(emb) <- paste0(row.names(emb),"x")
remove(Glioma)
remove(sub_object)
GC()

# load loom file
ldat <- read.loom.matrices("data/velocyto/PM1005_merged.loom")
hist(log10(rowSums(ldat$spliced)+1),col='wheat',xlab='log10[ number of reads + 1]',
     main='number of reads per gene')

emat <- ldat$spliced
nmat <- ldat$unspliced
smat <- ldat$ambiguous
str(emat@Dim);str(nmat@Dim);str(smat@Dim)
# look at the resulting gene set
str(intersect(intersect(rownames(emat),rownames(nmat)),rownames(smat)))

# filter expression matrices based on some minimum max-cluster averages
emat <- filter.genes.by.cluster.expression(emat,cell.colors,min.max.cluster.average = 0.2)
nmat <- filter.genes.by.cluster.expression(nmat,cell.colors,min.max.cluster.average = 0.1)
smat <- filter.genes.by.cluster.expression(smat,cell.colors,min.max.cluster.average = 0.02)
str(emat@Dim);str(nmat@Dim);str(smat@Dim)
# look at the resulting gene set
str(intersect(intersect(rownames(emat),rownames(nmat)),rownames(smat)))