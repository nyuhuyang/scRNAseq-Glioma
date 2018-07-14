########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################

library(Seurat)
library(dplyr)
source("./R/Seurat_functions.R")
########################################################################
#
#  1 Seurat Alignment 
# 
# ######################################################################
#======1.1 Setup the Seurat objects =========================
# Load the dataset

# setup Seurat objects since both count matrices have already filtered
# cells, we do no additional filtering here

# rename all "_" into "-" in sample names 
Glioma_raw <- list()
Glioma_Seurat <- list()
samples <- c("PM1258-2D","PM1258-GLICO","PM1258-TO")
conditions <- c("2D","GLICO","TO")

for(i in 1:length(samples)){
    Glioma_raw[[i]] <- Read10X(data.dir = paste0("./data/",
                                                     samples[i],"/outs/filtered_gene_bc_matrices/hg19/"))
    colnames(Glioma_raw[[i]]) <- paste0(samples[i],"_",colnames(Glioma_raw[[i]]))
    Glioma_Seurat[[i]] <- CreateSeuratObject(Glioma_raw[[i]],
                                           min.cells = 3,
                                           min.genes = 200,
                                           names.delim = "_")
    Glioma_Seurat[[i]]@meta.data$conditions <- conditions[i]
}
Glioma_Seurat <- lapply(Glioma_Seurat, FilterCells, 
                            subset.names = "nGene", 
                            low.thresholds = 200, 
                            high.thresholds = Inf)
# Calculate median UMI per cell
Glioma_raw_data <- lapply(Glioma_Seurat,function(x) as.matrix(x = x@raw.data))
lapply(Glioma_raw_data,function(x) mean(colSums(x)))
lapply(Glioma_raw_data,function(x) median(colSums(x)))

par(mfrow = c(1, length(samples)), cex = 1.5)
for(i in 1:length(samples)) {
        boxplot(colSums(Glioma_raw_data[[i]]))
        title(samples[i])
}

Glioma_Seurat <- lapply(Glioma_Seurat, NormalizeData)
Glioma_Seurat <- lapply(Glioma_Seurat, ScaleData)
Glioma_Seurat <- lapply(Glioma_Seurat, FindVariableGenes, do.plot = FALSE)

# we will take the union of the top 1k variable genes in each dataset for
# alignment note that we use 1k genes in the manuscript examples, you can
# try this here with negligible changes to the overall results
genes.use <- lapply(Glioma_Seurat, function(x) head(rownames(x@hvg.info), 1000))
genes.use <- unique(unlist(genes.use))
for(i in 1:length(samples)){
        genes.use <- intersect(genes.use, rownames(Glioma_Seurat[[i]]@scale.data))
}
length(genes.use) # 1/10 of total sample size


#======1.2 Perform a canonical correlation analysis (CCA) =========================
# run a canonical correlation analysis to identify common sources
# of variation between the two datasets.
remove(Glioma_raw)
GC()
Glioma <- RunMultiCCA(object.list = Glioma_Seurat, 
                    genes.use = genes.use,
                    niter = 25, num.ccs = 50,
                    standardize =TRUE)
save(Glioma, file = "./data/Glioma_alignment.Rda")

# CCA plot CC1 versus CC2 and look at a violin plot
p1 <- DimPlot(object = Glioma, reduction.use = "cca", group.by = "conditions", 
              pt.size = 0.5, do.return = TRUE)
p2 <- VlnPlot(object = Glioma, features.plot = "CC1", group.by = "conditions", 
              do.return = TRUE)
plot_grid(p1, p2)

p3 <- MetageneBicorPlot(Glioma, grouping.var = "conditions", dims.eval = 1:50, 
                        display.progress = TRUE, smooth = TRUE) # run on cluster
p3 + geom_smooth(method = 'loess')

PrintDim(object = Glioma, reduction.type = "cca", dims.print = 1:2, genes.print = 10)

DimHeatmap(object = Glioma, reduction.type = "cca", cells.use = 500, dim.use = c(1:9), 
           do.balanced = TRUE)
DimHeatmap(object = Glioma, reduction.type = "cca", cells.use = 500, dim.use = c(10:18), 
           do.balanced = TRUE)
#======1.3 QC =========================
Glioma <- CalcVarExpRatio(object = Glioma, reduction.type = "pca",
                              grouping.var = "conditions", dims.use = 1:20)
Glioma <- SubsetData(Glioma, subset.name = "var.ratio.pca",accept.low = 0.5) #8052 out of 8103

mito.genes <- grep(pattern = "^MT-", x = rownames(x = Glioma@data), value = TRUE)
percent.mito <- Matrix::colSums(Glioma@raw.data[mito.genes, ])/Matrix::colSums(Glioma@raw.data)
Glioma <- AddMetaData(object = Glioma, metadata = percent.mito, col.name = "percent.mito")

#Now we can run a single integrated analysis on all cells!
p1 <- GenePlot.1(object = Glioma, gene1 = "nUMI", gene2 = "percent.mito",col.use=gg_color_hue(3))
p2 <- GenePlot.1(object = Glioma, gene1 = "nUMI", gene2 = "nGene",col.use=gg_color_hue(3))
plot_grid(p1,p2)
g1 <- VlnPlot(object = Glioma, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 1)

Glioma <- FilterCells(object = Glioma, subset.names = c("nGene", "percent.mito"), 
                          low.thresholds = c(800, -Inf), high.thresholds = c(6000, 0.1))

p3 <- GenePlot.1(object = Glioma, gene1 = "nUMI", gene2 = "percent.mito",col.use=gg_color_hue(3))
p4 <- GenePlot.1(object = Glioma, gene1 = "nUMI", gene2 = "nGene",col.use=gg_color_hue(3))
g2 <- VlnPlot(object = Glioma, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 1)
plot_grid(p1, p3, p2, p4, ncol=2)
plot_grid(g1, g2, ncol=2)

Glioma <- ScaleData(object = Glioma, genes.use = genes.use, display.progress = FALSE, 
                    vars.to.regress = "percent.mito")
#======1.4 align seurat objects =========================
#Now we align the CCA subspaces, which returns a new dimensional reduction called cca.aligned
set.seed(42)
Glioma <- AlignSubspace(object = Glioma, reduction.type = "cca", grouping.var = "conditions", 
                            dims.align = 1:20)
#Now we can run a single integrated analysis on all cells!

Glioma <- FindClusters(object = Glioma, reduction.type = "cca.aligned", dims.use = 1:20, 
                           resolution = 3, force.recalc = T, save.SNN = TRUE)

Glioma <- RunTSNE(object = Glioma, reduction.use = "cca.aligned", dims.use = 1:20, 
                      do.fast = TRUE)

p1 <- TSNEPlot(Glioma, do.return = T, pt.size = 1, no.legend = TRUE,
               do.label = TRUE,group.by = "orig.ident",label.size = 8) + 
        ggtitle("Samples") + 
        theme(plot.title = element_text(hjust = 0.5))
p2 <- TSNEPlot(Glioma, do.return = T, pt.size = 1, no.legend = TRUE,
               do.label = TRUE, group.by = "ident",label.size = 8) + 
        ggtitle("Cluster ID") + 
        theme(plot.title = element_text(hjust = 0.5))
#png('./output/TSNESplot_alignment.png')
plot_grid(p1, p2)
#dev.off()
TSNEPlot(object = Glioma,do.label = TRUE, group.by = "ident", 
         do.return = TRUE, no.legend = TRUE,
         pt.size = 1,label.size = 8 )+
        ggtitle("tSNEplot for all cell clusters")+
        theme(text = element_text(size=20),     #larger text including legend title							
              plot.title = element_text(hjust = 0.5)) #title in middle
save(Glioma, file = "./data/Glioma_alignment.Rda")
#======1.5 PCA =========================
Glioma <- FindVariableGenes(object = Glioma, mean.function = ExpMean, dispersion.function = LogVMR, 
                         do.plot = FALSE)
hv.genes <- head(rownames(Glioma@hvg.info), 1000)
Glioma <- RunPCA(object = Glioma, pc.genes = hv.genes, pcs.compute = 100, do.print = TRUE, 
              pcs.print = 1:5, genes.print = 5)
PCElbowPlot(object = Glioma, num.pc = 100)
PCHeatmap(Glioma, pc.use = c(1:3, 70:75), cells.use = 500, do.balanced = TRUE)
Glioma <- StashIdent(object = Glioma, save.name = "cca_3.0")
Glioma <- FindClusters(object = Glioma, reduction.type = "pca", dims.use = 1:75, resolution = 3, 
                    save.SNN = TRUE, n.start = 10,force.recalc=T, nn.eps = 0.5, print.output = FALSE)
Glioma <- RunTSNE(object = Glioma, reduction.use = "pca", dims.use = 1:75, tsne.method = "FIt-SNE", 
               nthreads = 4, reduction.name = "FItSNE", reduction.key = "FItSNE_", 
               fast_tsne_path = "/Users/yah2014/src/FIt-SNE/bin/fast_tsne", 
               max_iter = 2000)
library(cowplot)
p1 <- DimPlot(object = Glioma, reduction.use = "FItSNE", no.legend = TRUE,
              do.return = TRUE,vector.friendly = F, pt.size = 1,
              do.label = TRUE,label.size = 8, group.by = "conditions") + 
        ggtitle("Samples") + 
        theme(plot.title = element_text(hjust = 0.5))

p2 <- DimPlot(object = Glioma, reduction.use = "FItSNE", no.legend = TRUE,
              do.return = TRUE,vector.friendly = F, pt.size = 1,
              do.label = TRUE,label.size = 8, group.by = "ident") + 
        ggtitle("Cluster ID") + 
        theme(plot.title = element_text(hjust = 0.5))

plot_grid(p1, p2)
Glioma <- StashIdent(object = Glioma, save.name = "pca_3.0")
save(Glioma, file = "./data/Glioma_alignment.Rda")
