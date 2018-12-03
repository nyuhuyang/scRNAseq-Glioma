########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################

suppressPackageStartupMessages({
        library(Seurat)
        library(magrittr)
        library(harmony)
        library(dplyr)
        library(kableExtra)
        source("../R/Seurat_functions.R")
})
path <- paste0("./output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
########################################################################
#
#  1 harmony Alignment 
# 
# ######################################################################
#======1.1 read sample file =========================
# Load the mouse.eyes dataset
# setup Seurat objects since both count matrices have already filtered
# cells, we do no additional filtering here
df_samples <- readxl::read_excel("doc/181002_Single_cell_sample list.xlsx")
sample_n = which(df_samples$tests %in% paste0("test",1:5))
table(df_samples$tests)
df_samples[sample_n,] %>% kable() %>% kable_styling()
samples <- df_samples$samples[sample_n]
projects <- df_samples$projects[sample_n]
conditions <- df_samples$conditions[sample_n]
tests <- df_samples$tests[sample_n]

Glioma_raw <- list()
Glioma_Seurat <- list()
for(i in 1:length(samples)){
    Glioma_raw[[i]] <- Read10X(data.dir = paste0("./data/",
                            samples[i],"/outs/filtered_gene_bc_matrices/hg19/"))
    colnames(Glioma_raw[[i]]) <- paste0(samples[i],"_",colnames(Glioma_raw[[i]]))
    Glioma_Seurat[[i]] <- CreateSeuratObject(Glioma_raw[[i]],
                                        min.cells = 3,
                                        min.genes = 200,
                                        project = projects[i],
                                        names.delim = "_")
    Glioma_Seurat[[i]]@meta.data$conditions <- conditions[i]
}
#======1.1.2 QC before merge =========================
cell.number <- sapply(Glioma_Seurat, function(x) length(x@cell.names))
QC_list <- lapply(Glioma_Seurat, function(x) as.matrix(x = x@raw.data))
median.nUMI <- sapply(QC_list, function(x) median(colSums(x)))
median.nGene <- sapply(QC_list, function(x) median(apply(x,2,function(y) sum(length(y[y>0])))))

min.nUMI <- sapply(QC_list, function(x) min(colSums(x)))
min.nGene <- sapply(QC_list, function(x) min(apply(x,2,function(y) sum(length(y[y>0])))))

QC.list <- cbind(df_samples,cell.number, median.nUMI,median.nGene,min.nUMI,min.nGene,
                row.names = samples)
write.csv(QC.list, paste0(path,"QC_list.csv"))
QC.list %>% kable() %>% kable_styling()
remove(QC_list);GC()

#========1.1.3 merge ===================================
Glioma <- Reduce(function(x, y) MergeSeurat(x, y, do.normalize = F), Glioma_Seurat)
remove(Glioma_raw,Glioma_Seurat);GC()
mito.genes <- grep(pattern = "^MT-", x = rownames(x = Glioma@data), value = TRUE)
percent.mito <- Matrix::colSums(Glioma@raw.data[mito.genes, ])/Matrix::colSums(Glioma@raw.data)
Glioma <- AddMetaData(object = Glioma, metadata = percent.mito, col.name = "percent.mito")

Glioma@ident = factor(Glioma@ident,levels = samples)

g1 <- VlnPlot(object = Glioma, features.plot = c("nGene", "nUMI", "percent.mito"),
            nCol = 1,point.size.use = 0.2,
            x.lab.rot = T, do.return = T,return.plotlist =T)

save(g1, file = "./data/g1_18_20181201.Rda")

#======1.2 load  SingleCellExperiment =========================
(load(file = "./data/sce_18_20181201.Rda"))
names(sce_list)
Glioma_Seurat <- lapply(sce_list, as.seurat) %>%
        lapply(NormalizeData) %>%
        #lapply(ScaleData) %>%
        lapply(FindVariableGenes, do.plot = FALSE)

for(i in 1:length(samples)){
        Glioma_Seurat[[i]]@meta.data$conditions <- conditions[i]
        Glioma_Seurat[[i]]@meta.data$tests <- tests[i]
}
# we will take the union of the top 1k variable genes in each dataset for alignment
genes.use <- Glioma_Seurat %>% lapply(function(object) head(rownames(object@hvg.info), 1000)) %>%
                unlist %>% unique
length(genes.use)

#========1.3 merge ===================================
Glioma <- Reduce(function(x, y) MergeSeurat(x, y, do.normalize = F), Glioma_Seurat)
Glioma@var.genes = genes.use
remove(sce_list,Glioma_Seurat);GC()

#Glioma@meta.data$orig.ident <- gsub("Pt-MD", "MD", Glioma@meta.data$orig.ident)
#Glioma = SetAllIdent(Glioma, id = "orig.ident")
#======1.4 mito, QC, filteration =========================
mito.genes <- grep(pattern = "^MT-", x = rownames(x = Glioma@data), value = TRUE)
percent.mito <- Matrix::colSums(Glioma@raw.data[mito.genes, ])/Matrix::colSums(Glioma@raw.data)
Glioma <- AddMetaData(object = Glioma, metadata = percent.mito, col.name = "percent.mito")

(load(file = "./data/g1_18_20181201.Rda"))

Glioma <- FilterCells(object = Glioma, subset.names = c("nGene","nUMI","percent.mito"),
                   low.thresholds = c(1000,3000, -Inf), 
                   high.thresholds = c(Inf,Inf, 0.15))

par(mfrow = c(1, 2))
GenePlot(object = Glioma, gene1 = "nUMI", gene2 = "percent.mito")
GenePlot(object = Glioma, gene1 = "nUMI", gene2 = "nGene")

Glioma@ident = factor(Glioma@ident,levels = samples)

g2 <- VlnPlot(object = Glioma, features.plot = c("nGene", "nUMI", "percent.mito"), 
              nCol = 1,point.size.use = 0.2,
              x.lab.rot = T, do.return = T,return.plotlist =T)
save(g2,file = "./data/g2_18_20181201.Rda")
jpeg(paste0(path,"/S1_nGene.jpeg"), units="in", width=10, height=7,res=600)
print(plot_grid(g1[[1]]+ggtitle("nGene in raw data")+ 
                    scale_y_log10(limits = c(200,10000)),#+ylim(c(0,1000)),
                g2[[1]]+ggtitle("nGene after filteration")+ 
                    scale_y_log10(limits = c(200,10000))))
dev.off()
jpeg(paste0(path,"/S1_nUMI.jpeg"), units="in", width=10, height=7,res=600)
print(plot_grid(g1[[2]]+ggtitle("nUMI in raw data")+ 
                    scale_y_log10(limits = c(1500,100000)),#+ylim(c(0,1000)),
                g2[[2]]+ggtitle("nUMI after filteration")+ 
                    scale_y_log10(limits = c(1500,100000))))
dev.off()
jpeg(paste0(path,"/S1_mito.jpeg"), units="in", width=10, height=7,res=600)
print(plot_grid(g1[[3]]+ggtitle("mito % in raw data")+ 
                    ylim(c(0,0.5)),
                g2[[3]]+ggtitle("mito % after filteration")+ 
                    ylim(c(0,0.5))))
dev.off()
# After removing unwanted cells from the dataset, the next step is to normalize the data.
Glioma <- NormalizeData(object = Glioma, normalization.method = "LogNormalize", 
                     scale.factor = 10000)
Glioma <- FindVariableGenes(object = Glioma, mean.function = ExpMean, 
                         dispersion.function = LogVMR, do.plot = FALSE, 
                         x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
length(Glioma@var.genes)

#======1.5 Add Cell-cycle score =========================
# Read in a list of cell cycle markers, from Tirosh et al, 2015
cc.genes <- readLines(con = "../R/seurat_resources/regev_lab_cell_cycle_genes.txt")
s.genes <- HumanGenes(Glioma,cc.genes[1:43])
g2m.genes <- HumanGenes(Glioma,cc.genes[44:97])
Glioma <- CellCycleScoring(object = Glioma, s.genes = s.genes, g2m.genes = g2m.genes, 
                        set.ident = TRUE)
RidgePlot(object = Glioma, features.plot = HumanGenes(Glioma,c("CCND1","CDK4","CCND2","CDK6","CCND3","RB1")), 
          nCol = 2)
Glioma@meta.data$CC.Difference <- Glioma@meta.data$S.Score - Glioma@meta.data$G2M.Score
Glioma@meta.data$S.Score = Glioma@meta.data$S.Score - min(Glioma@meta.data$S.Score)
Glioma@meta.data$G2M.Score = Glioma@meta.data$G2M.Score - min(Glioma@meta.data$G2M.Score)
head(x = Glioma@meta.data)
save(Glioma, file = "./data/Glioma_Harmony_20181201.Rda")
#======1.5 1st run of pca-tsne  =========================
Glioma %<>%  ScaleData %<>%
         RunPCA(pc.genes = Glioma@var.genes, pcs.compute = 100, do.print = F)
#Glioma %<>% RunICA(ic.genes = Glioma@var.genes, ics.compute = 50, print.results = F)
jpeg(paste0(path,"/S1_DimElbowPlot_pca.jpeg"), units="in", width=10, height=7,res=600)
DimElbowPlot(Glioma, reduction.type = "pca", dims.plot = 100)
dev.off()
system.time(Glioma %<>% RunHarmony("orig.ident", 
                                theta = 2, plot_convergence = TRUE,
                                nclust = 50, max.iter.cluster = 100))

Glioma@ident %<>% factor(levels = samples)
p1 <- DimPlot(object = Glioma, reduction.use = "harmony", pt.size = .1, group.by = "orig.ident", do.return = T)
p2 <- VlnPlot(object = Glioma, features.plot = "Harmony1", group.by = "orig.ident", do.return = TRUE,
              x.lab.rot = T)
jpeg(paste0(path,"/S1_Harmony_vplot.jpeg"), units="in", width=10, height=7,res=600)
plot_grid(p1,p2)
dev.off()

jpeg(paste0(path,"/S1_Harmony_DimHeatmap.jpeg"), units="in", width=10, height=7,res=600)
DimHeatmap(object = Glioma, reduction.type = "harmony", cells.use = 500, dim.use = 1:6, do.balanced = TRUE)
dev.off()

#========1.6 Seurat tSNE Functions for Integrated Analysis Using Harmony Results=======
system.time({
        Glioma %<>% RunTSNE(reduction.use = "harmony", dims.use = 1:20, do.fast = T)
        Glioma %<>% FindClusters(reduction.type = "harmony", resolution = 0.6, dims.use = 1:20, 
                              force.recalc = TRUE, print.output = FALSE)
})

p3 <- TSNEPlot(Glioma, do.return = T, pt.size = 0.5, group.by = "orig.ident")
p4 <- TSNEPlot(Glioma, do.label = T, do.return = T, pt.size = 0.5)
jpeg(paste0(path,"/S1_Harmony_TSNEPlot.jpeg"), units="in", width=10, height=7,res=600)
plot_grid(p3, p4)
dev.off()

g_Harmony <- TSNEPlot.1(object = Glioma, do.label = F, group.by = "ident", 
                 do.return = TRUE, no.legend = T, 
                 colors.use = ExtractMetaColor(Glioma),
                 pt.size = 1,label.size = 6 )+
        ggtitle("Tsne plot of all cell types")+
        theme(text = element_text(size=15),							
              plot.title = element_text(hjust = 0.5,size = 18, face = "bold")) 

jpeg(paste0(path,"/TSNEplot-Harmony.jpeg"), units="in", width=10, height=7,res=600)
print(g_Harmony)
dev.off()

Glioma.scale.data <- Glioma@scale.data
save(Glioma.scale.data, file = "./data/Glioma.scale.data_Harmony_20181201.Rda")
remove(Glioma.scale.data);GC()

Glioma@scale.data =NULL
save(Glioma, file = "./data/Glioma_Harmony_20181201.Rda")
GC()



jpeg(paste0(path,"/TSNEplot-alignment~.jpeg"), units="in", width=10, height=7,res=600)
do.call(plot_grid, list(g_CCA,g_MNN,g_Harmony))
dev.off()
