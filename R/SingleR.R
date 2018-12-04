library(SingleR)
library(Seurat)
library(reshape2)
library(pheatmap)
library(kableExtra)
source("../R/Seurat_functions.R")
source("../R/SingleR_functions.R")
path <- paste0("./output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
#====== 3.1 Create Singler Object  ==========================================
(load(file = "./data/Glioma_Harmony_20181201.Rda"))
(load(file='./data/GeneSets/ref_GBMLGG_RSEM.RData'))
length(ref$types)
length(unique(ref$types))
length(unique(ref$main_types))
#pca
DimPlot(object = Glioma, reduction.use = "tsne", no.legend = TRUE,
        do.return = TRUE,vector.friendly = F, pt.size = 1,
        do.label = TRUE,label.size = 8, group.by = "ident") + 
        ggtitle("Cluster ID") + 
        theme(plot.title = element_text(hjust = 0.5))
singler = CreateSinglerObject(Glioma@data, annot = Glioma@ident,
                              project.name=Glioma@project.name,
                              min.genes = 500,technology = "10X", species = "Human", 
                              citation = "",ref.list = list(ref), 
                              normalize.gene.length = F,
                              variable.genes = "de",
                              fine.tune = T, do.signatures = F, clusters = NULL) # run in cluster
# if singler didn't find all cell labels
length(singler$singler[[1]]$SingleR.single$labels) == ncol(Glioma@data)
singler$meta.data$orig.ident = Glioma@meta.data$orig.ident # the original identities, if not supplied in 'annot'
singler$meta.data$xy = Glioma@dr$tsne@cell.embeddings # the tSNE coordinates
singler$meta.data$clusters = Glioma@ident # the Seurat clusters (if 'clusters' not provided)
save(singler,file="./output/singler_Glioma_20181203.RData")

#====== 3.2 SingleR specifications ==========================================
# Step 1: Spearman coefficient
(load(file="./output/singler_Glioma_20181203.RData"))

##############################
# add singleR label to Seurat
###############################
singlerDF = data.frame("singler1sub"=singler$singler[[1]]$SingleR.single$labels,
                       "singler1main"=singler$singler[[1]]$SingleR.single.main$labels,
                       "singler2sub"=singler$singler[[2]]$SingleR.single$labels,
                       "singler2main"=singler$singler[[2]]$SingleR.single.main$labels,
                       "kang" = FineTune(singler$other, main.type = FALSE),
                       row.names = rownames(singler$singler[[1]]$SingleR.single$labels))

table(rownames(singlerDF) %in% Glioma@cell.names)

#knowDF = data.frame("cell.names"= Glioma@cell.names)
#ident.DF = full_join(singlerDF,knowDF, by="cell.names")
#ident.DF<- apply(ident.DF,2,as.character)
#rownames(ident.DF) = ident.DF[,"cell.names"]
#ident.DF = ident.DF[,-which(colnames(ident.DF) == "cell.names")]
apply(singlerDF,2,function(x) length(unique(x)))
#ident.DF[is.na(ident.DF)] <- "unknown"
Glioma <- AddMetaData(object = Glioma,
                   metadata = singlerDF)
Glioma <- SetAllIdent(object = Glioma, id = "singler1sub")

##############################
# check the spearman correlation
###############################
#Or by all cell types (showing the top 50 cell types):
jpeg(paste0(path,"/DrawHeatmap_sub1.jpeg"), units="in", width=10, height=7,
     res=600)
print(SingleR.DrawHeatmap(singler$singler[[1]]$SingleR.single, top.n = 8,normalize = F))
dev.off()
jpeg(paste0(path,"/DrawHeatmap_sub2.jpeg"), units="in", width=10, height=7,
     res=600)
print(SingleR.DrawHeatmap(singler$singler[[2]]$SingleR.single,top.n = 50,normalize = F))
dev.off()

#Finally, we can also view the labeling as a table compared to the original identities:

kable(table(singler$meta.data$orig.ident,
            singler$singler[[1]]$SingleR.single$labels)) %>%
        kable_styling()
Glioma@meta.data$singler1sub %>% table() %>% kable() %>% kable_styling()
Glioma@meta.data$singler2sub %>% table() %>% kable() %>% kable_styling()

##############################
# process color scheme
##############################
singler_colors <- readxl::read_excel("./doc/singler.colors.xlsx")
singler_colors1 = as.vector(singler_colors$singler.color1[!is.na(singler_colors$singler.color1)])
singler_colors2 = as.vector(singler_colors$singler.color2[!is.na(singler_colors$singler.color2)])
singler_colors1[duplicated(singler_colors1)];singler_colors2[duplicated(singler_colors2)]
length(singler_colors1);length(singler_colors2)
apply(Glioma@meta.data[,c("singler1sub","singler1main","singler2sub","singler2main")],
      2,function(x) length(unique(x)))
Glioma@meta.data[,c("singler2sub")] %>% table() %>% kable() %>% kable_styling()
Glioma <- AddMetaColor(object = Glioma, label= "singler1sub", colors = singler_colors1[1:8])
Glioma <- AddMetaColor(object = Glioma, label= "singler2sub", colors = singler_colors2[1:21])
Glioma <- SetAllIdent(object = Glioma, id = "singler1sub")
TSNEPlot.1(Glioma, colors.use = ExtractMetaColor(Glioma),no.legend = F)

##############################
# draw tsne plot
##############################
p3 <- TSNEPlot.1(object = Glioma, do.label = T, group.by = "ident", 
                 do.return = TRUE, no.legend = T, 
                 colors.use = ExtractMetaColor(Glioma),
                 pt.size = 1,label.size = 5,force = 2)+
        ggtitle("Supervised cell type labeling by TCGA GBM and LGG")+
        theme(text = element_text(size=10),							
              plot.title = element_text(hjust = 0.5,size = 18, face = "bold")) 

jpeg(paste0(path,"PlotTsne_sub1.jpeg"), units="in", width=10, height=7,
     res=600)
print(p3)
dev.off()

save(Glioma,file="./data/Glioma_Harmony_20181201.Rda")
##############################
# subset Seurat
###############################
table(Glioma@meta.data$orig.ident)
table(Glioma@ident)

df_samples <- readxl::read_excel("doc/181002_Single_cell_sample list.xlsx")
tests <- paste0("test",c(1:5))
for(test in tests){
        sample_n = which(df_samples$tests %in% test)
        samples <- unique(df_samples$samples[sample_n])
        print(samples)
        
        cell.use <- rownames(Glioma@meta.data)[Glioma@meta.data$orig.ident %in% samples]
        subset.Glioma <- SubsetData(Glioma, cells.use = cell.use)
        g <- SplitTSNEPlot(subset.Glioma,group.by = "ident",split.by = "orig.ident",
                           no.legend = T,do.label =T,label.size=3,
                           return.plots =T, label.repel = T,force=2)
        jpeg(paste0(path,test,"_TSNEPlot.jpeg"), units="in", width=10, height=7,
             res=600)
        print(do.call(plot_grid, g))
        dev.off()
}
##############################
# draw barchart
##############################
for(test in tests){
      sample_n = which(df_samples$tests %in% test)
      samples <- unique(df_samples$samples[sample_n])
      print(samples)
      
      cell.use <- rownames(Glioma@meta.data)[Glioma@meta.data$orig.ident %in% samples]
      subset.Glioma <- SubsetData(Glioma, cells.use = cell.use)
      g <- SplitBarchart(subset.Glioma,group.by = "ident",split.by = "orig.ident",
                         no.legend = T,label.size=3,do.print =T,
                         return.plots =T)
      jpeg(paste0(path,test,"_SplitBarChart.jpeg"), units="in", width=10, height=7,res=600)
      print(do.call(cowplot::plot_grid, c(g, align = "hv")))
      dev.off()
}

# Split Seurat by certein criteria and make tsne plot
#' @param object Seurat object
#' @param split.by the criteria to split, can be gene name, or any variable in meta.data
#' @param select.plots output order, default to NULL. If want to change,use c(2,1) for example
#' @param return.data TRUE/FASLE, return splited ojbect or not.
#' @export p ggplot object from barchart
#' @example SplitBarchart(Glioma, group.by = "ident",split.by = "orig.ident")
SplitBarchart <- function(object, split.by = "orig.ident",select.plots = NULL, 
                          do.return = TRUE, do.print = FALSE,
                          do.label = T, group.by = "ident", no.legend = TRUE,
                          pt.size = 1,label.size = 5,size=20,... ){
  
  
      object.subsets <- SplitSeurat(object = object, split.by = split.by)
      levels <- object.subsets[[length(object.subsets)]]
      
      p <- list()
      if(is.null(select.plots)) select.plots <- 1:length(levels)
      for(i in 1:length(select.plots)){
        p[[i]] <- ggplot(object.subsets[[select.plots[i]]]@meta.data, 
                         aes(x = singler1sub, fill=singler1sub))+
          geom_bar()+ ggtitle(levels[select.plots[i]])+
          ylab("Cell Number")+
          theme(text = element_text(size = size),
                legend.position="none",
                plot.title = element_text(hjust = 0.5),
                axis.title.x=element_blank(),
                axis.text.x = element_text(angle = 90, hjust = 1))+
          scale_y_log10()+
          scale_fill_manual(values=ExtractMetaColor(object.subsets[[select.plots[i]]]))
      }
      p <- p[lapply(p,length)>0] # remove NULL element
      if(do.print) {
        path <- paste0("./output/",gsub("-","",Sys.Date()),"/")
        if(!dir.exists(path)) dir.create(path, recursive = T)
        jpeg(paste0(path,"SplitBarChart_",deparse(substitute(object)),".jpeg"), units="in", width=10, height=7,res=600)
        print(do.call(cowplot::plot_grid, c(p, align = "hv")))
        dev.off()
      }
      if(do.return) return(p) 
}
