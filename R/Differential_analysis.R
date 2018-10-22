########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################

library(Seurat)
library(dplyr)
source("../R/Seurat_functions.R")

#3.1  Compare DE across all major cell types==================
#We would need the data for all clusters, as well the subclusters.
#detect changes in gene expression between young and aged, 
#in the different cell types and subtypes. 
#It will also be interesting to check if there is some subtype enriched in young compared to aged or viceversa. 

# answer 08/09 request
lnames = load(file = "./data/Glioma_alignment.Rda")
lnames
table(Glioma@ident)
ident.vector <- as.factor(Glioma@meta.data$orig.ident)
names(x = ident.vector) <- names(Glioma@ident)
Glioma@ident = ident.vector
TSNEPlot(object = Glioma,do.label = T, group.by = "ident", 
         do.return = TRUE, no.legend = T,
         pt.size = 1,label.size = 8 )+
        ggtitle("Glioma")+
        theme(text = element_text(size=20),     #larger text including legend title							
              plot.title = element_text(hjust = 0.5,size = 25, face = "bold")) 
Glioma.markers <- FindAllMarkers.UMI(Glioma)
write.csv(Glioma.markers,"./output/Glioma_markers.csv")
top <- Glioma.markers %>% group_by(cluster) %>% top_n(30, avg_logFC)
# setting slim.col.label to TRUE will print just the cluster IDS instead of
# every cell name
DoHeatmap.1(object = Glioma, genes.use = top$gene, 
            slim.col.label = TRUE, remove.key = T,cex.row = 8,
            group.cex = 15,
            rotate.key = T,group.label.rot = F)+
        ggtitle("Expression heatmap of top 30 differential expression genes in each group")+
        theme(text = element_text(size=20),
              plot.title = element_text(hjust = 0.5))


SingleFeaturePlot.1(object = Glioma, "FGF13",threshold = 0.1)
