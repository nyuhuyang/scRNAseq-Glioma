library(SingleR)
library(Seurat)
library(reshape2)
library(pheatmap)
library(kableExtra)

source("../R/Seurat_functions.R")
#====== 3.1 Create Singler Object  ==========================================
lnames = load(file = "./data/Glioma_alignment.Rda")
lnames
lname = load(file='./data/GeneSets/Refs_TCGA_IvyGbm.RData') 
lname
length(Refs_TCGA_IvyGbm$types)
length(unique(Refs_TCGA_IvyGbm$types))
length(unique(Refs_TCGA_IvyGbm$main_types))
#pca
DimPlot(object = Glioma, reduction.use = "tsne", no.legend = TRUE,
        do.return = TRUE,vector.friendly = F, pt.size = 1,
        do.label = TRUE,label.size = 8, group.by = "ident") + 
        ggtitle("Cluster ID") + 
        theme(plot.title = element_text(hjust = 0.5))
singler = CreateSinglerObject(as.matrix(Glioma@data), annot = Glioma@ident,
                                 project.name=Glioma@project.name,
                              do.main.types = T,
                              min.genes = 500,technology = "10X", species = "Human", 
                              citation = "",ref.list = list(Refs_TCGA_IvyGbm), 
                              normalize.gene.length = F,
                              variable.genes = "de",
                              fine.tune = F, do.signatures = F, clusters = NULL) # run in cluster

singler$meta.data$orig.ident = Glioma@meta.data$orig.ident # the original identities, if not supplied in 'annot'
singler$meta.data$xy = Glioma@dr$tsne@cell.embeddings # the tSNE coordinates
singler$meta.data$clusters = Glioma@ident # the Seurat clusters (if 'clusters' not provided)
save(singler,file="./data/singler_Verh_IvyGbm.RData")
#save(singler,file="./data/singler_Refs_TCGA_IvyGbm1.RData")
#save(singler,file="./data/singler_Glioma.RData")

#====== 3.2 SingleR specifications ==========================================
# Step 1: Spearman coefficient
lnames = load(file = "./data/singler_Verh_IvyGbm.RData")
lnames
singler$seurat = Glioma # (optional)

TSNEPlot(object = Glioma,do.label = TRUE, group.by = "orig.ident", 
         do.return = TRUE, no.legend = TRUE,
         pt.size = 1,label.size = 8 )+
        ggtitle("tSNEplot for all samples")+
        theme(text = element_text(size=20),     #larger text including legend title							
              plot.title = element_text(hjust = 0.5)) #title in middle

SingleR.DrawScatter(sc_data = singler$seurat@data,cell_id = 10, 
                    ref = immgen, sample_id = 232)

# Step 2: Multiple correlation coefficients per cell types are aggregated 
# to provide a single value per cell type per single-cell. 
# In the examples below we use the 80% percentile of correlation values.
# for visualization purposes we only present a subset of cell types (defined in labels.use)
out = SingleR.DrawBoxPlot(sc_data = singler$seurat@data,cell_id = 10, 
                          ref = Refs_TCGA_IvyGbm,main_types = T,
                          labels.use=c('Astrocyte','Fibroblasts',
                                       'Myocytes','Mesangial_cells','Neuroepithelial_cells:ESC-derived'))
print(out$plot)
SingleR.DrawHeatmap(singler$singler[[1]]$SingleR.single.main, top.n = Inf,
                    clusters = singler$meta.data$orig.ident)
#Or by all cell types (showing the top 50 cell types):
SingleR.DrawHeatmap(singler$singler[[1]]$SingleR.single, top.n = 50,
                    clusters = singler$meta.data$orig.ident)
SingleR.DrawHeatmap(singler$singler[[1]]$SingleR.single.main,top.n = 60,
                    normalize = F,clusters = singler$meta.data$orig.ident)
#Next, we can use the fine-tuned labels to color the t-SNE plot:
# types-------
out = SingleR.PlotTsne.1(singler$singler[[1]]$SingleR.single,
                         singler$meta.data$xy,do.label=T,
                         do.letters = F,labels = singler$singler[[1]]$SingleR.single$labels,
                         label.size = 5, dot.size = 1,do.legend = F,alpha = 1,
                         label.repel = T,force=2)
out$p+  ggtitle("Supervised cell labeling by HPCA, Blueprint, Encode, TCGA.GBM and Ivy Glioblastoma")+#ggplot title
        theme(text = element_text(size=20),     #larger text including legend title
              plot.title = element_text(hjust = 0.5,size = 17, face = "bold")) #title in middle
# main types-------
out = SingleR.PlotTsne.1(singler$singler[[1]]$SingleR.single.main,
                         singler$meta.data$xy,do.label=T,
                         do.letters = F,labels = singler$singler[[1]]$SingleR.single.main$labels,
                         label.size = 5, dot.size = 1,do.legend = F,alpha = 1,
                         label.repel = T,force=2)
out$p+  ggtitle("Supervised cell labeling by HPCA, Blueprint, Encode, TCGA.GBM and Ivy.Glioblastoma")+#ggplot title
        theme(text = element_text(size=20),     #larger text including legend title
              plot.title = element_text(hjust = 0.5,size = 15, face = "bold")) #title in middle

SplitSingleR.PlotTsne(singler = singler, split.by = "orig.ident",do.label=T,
                      do.letters = T,do.legend = FALSE,force=2)
output <- SplitSingleR.PlotTsne(singler = singler, split.by = "conditions",
                                return.plots=T,do.label=T,do.legend = F,alpha = 0.5,
                                label.repel = F, force=2)
output[[1]]
output[[2]]
#Finally, we can also view the labeling as a table compared to the original identities:

kable(table(singler$singler[[1]]$SingleR.single.main$labels,
            singler$meta.data$orig.ident)) %>%
        kable_styling()
kable(table(singler$meta.data$orig.ident,
            singler$singler[[1]]$SingleR.single.main$labels)) %>%
        kable_styling()
kable(table(singler$meta.data$orig.ident,singler$seurat@ident)) %>%
        kable_styling()

cells_prop <- as.data.frame(table(singler$singler[[1]]$SingleR.single$labels,
                                  singler$meta.data$orig.ident))
cells_prop <- dcast(cells_prop,Var2~Var1)
rownames(cells_prop) = cells_prop$Var2
cells_prop <- cells_prop[,-1]
total  <- colSums(cells_prop)
Cells_prop <- cells_prop/total
Cells_prop %>% kable %>% kable_styling()

All_cells <- as.data.frame(table(Glioma@meta.data$orig.ident))
glioblastoma$total <- All_cells$Freq
glioblastoma$percentage <- glioblastoma$Freq/glioblastoma$total
glioblastoma
