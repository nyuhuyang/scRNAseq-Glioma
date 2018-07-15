library(SingleR)
library(Seurat)
library(reshape2)
library(pheatmap)
library(kableExtra)

source("./R/Seurat_functions.R")
#====== 3.1 Create Singler Object  ==========================================
lnames = load(file = "./data/Glioma_alignment.Rda")
lnames
#cca
Glioma <- SetAllIdent(object = Glioma, id = "cca_3.0")
DimPlot(object = Glioma, reduction.use = "tsne", no.legend = TRUE,
        do.return = TRUE,vector.friendly = F, pt.size = 1,
        do.label = TRUE,label.size = 8, group.by = "ident") + 
        ggtitle("Cluster ID") + 
        theme(plot.title = element_text(hjust = 0.5))
singler = CreateSinglerObject(as.matrix(Glioma@data), annot = Glioma@ident,
                                 project.name=Glioma@project.name,
                              min.genes = 500,technology = "10X", species = "Human", 
                              citation = "",ref.list = list(), normalize.gene.length = F,
                              variable.genes = "de",
                              fine.tune = T, do.signatures = F, clusters = NULL) # run in cluster
singler$seurat = Glioma # (optional)
singler$meta.data$orig.ident = Glioma@meta.data$orig.ident # the original identities, if not supplied in 'annot'
singler$meta.data$xy = Glioma@dr$tsne@cell.embeddings # the tSNE coordinates
singler$meta.data$clusters = Glioma@ident # the Seurat clusters (if 'clusters' not provided)
save(singler,file="./data/singler_Glioma.RData")
# pca
Glioma <- SetAllIdent(object = Glioma, id = "pca_3.0")
DimPlot(object = Glioma, reduction.use = "FItSNE", no.legend = TRUE,
        do.return = TRUE,vector.friendly = F, pt.size = 1,
        do.label = TRUE,label.size = 8, group.by = "ident") + 
        ggtitle("Cluster ID") + 
        theme(plot.title = element_text(hjust = 0.5))
singler$meta.data$xy = Glioma@dr$FItSNE@cell.embeddings # the tSNE coordinates
singler$meta.data$clusters = Glioma@ident # the Seurat clusters (if 'clusters' not provided)
save(singler,file="./data/singler_Glioma.RData")
#====== 3.2 SingleR specifications ==========================================
# Step 1: Spearman coefficient
lnames = load(file = "./data/singler_Glioma.RData")
lnames
DimPlot(object = singler$seurat, reduction.use = "FItSNE", no.legend = TRUE,
        do.return = TRUE,vector.friendly = F, pt.size = 1,
        do.label = TRUE,label.size = 8, group.by = "ident") + 
        ggtitle("Cluster ID") + 
        theme(plot.title = element_text(hjust = 0.5))

SingleR.DrawScatter(sc_data = singler$seurat@data,cell_id = 10, 
                    ref = immgen, sample_id = 232)

# Step 2: Multiple correlation coefficients per cell types are aggregated 
# to provide a single value per cell type per single-cell. 
# In the examples below we use the 80% percentile of correlation values.
# for visualization purposes we only present a subset of cell types (defined in labels.use)
out = SingleR.DrawBoxPlot(sc_data = singler$seurat@data,cell_id = 100, 
                          ref = hpca,main_types = T,
                          labels.use=c('Astrocyte','iPS_cells','Neuroepithelial_cell','Neurons',
                                       'Tissue_stem_cells','Embryonic_stem_cells','Fibroblasts','Endothelial cells'))
print(out$plot)
SingleR.DrawHeatmap(singler$singler[[1]]$SingleR.single.main, top.n = Inf,
                    clusters = singler$meta.data$orig.ident)
#Or by all cell types (showing the top 50 cell types):
SingleR.DrawHeatmap(singler$singler[[1]]$SingleR.single, top.n = 50,
                    clusters = singler$meta.data$orig.ident)
#Next, we can use the fine-tuned labels to color the t-SNE plot:
out = SingleR.PlotTsne(singler$singler[[1]]$SingleR.single,
                       singler$meta.data$xy,do.label=F,
                       do.letters = F,labels = singler$singler[[1]]$SingleR.single$labels,
                       label.size = 4, dot.size = 3,do.legend = TRUE)
out$p

SplitSingleR.PlotTsne(singler = singler, split.by = "conditions",do.label=F,do.legend = T)
output <- SplitSingleR.PlotTsne(singler = singler, split.by = "conditions",
                              return.plots=T,do.label=T,do.legend = T)
output[[1]]
output[[2]]
output[[2]]
#Finally, we can also view the labeling as a table compared to the original identities:

table(singler$singler[[1]]$SingleR.single$labels, singler$meta.data$orig.ident) %>%
        kable %>%  kable_styling()
