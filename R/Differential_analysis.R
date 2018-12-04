########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################

library(Seurat)
library(dplyr)
source("../R/Seurat_functions.R")
path <- paste0("./output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)

(load(file = "./data/Glioma_Harmony_20181201.Rda"))

# heatmap.2 =================
markers <- HumanGenes(Glioma,c("TP53","ATRX","CIC","FUBP1","NOTCH1",
                               "NIPBL","PIK3CA","PIK3R1","PDGFRA",
                               "CDKN2A","CDKN2B","CCND2","CDK4","RB1",
                               "MDM4","MDM2","EGFR","NF1","BRAF",
                               "PIK3R1","FGF13","PTEN"),unique = T)

y = Glioma@data[markers,]
## Column clustering (adjust here distance/linkage methods to what you need!)
hc <- fastcluster::hclust.vector(dist(1-cor(as.matrix(y), method="spearman")),
                                 method="single") #Run in cluster
cc = GBMLGG@meta.data[hc$labels,"subtype.colors"] %>% as.vector

jpeg(paste0(path,"/Heatmap2_TCGA_GBMLGG.jpeg"), units="in", width=10, height=7,res=600)
heatmap.2(as.matrix(y),
          Colv = as.dendrogram(hc),Rowv= FALSE,
          ColSideColors = cc, trace ="none",labCol = FALSE,dendrogram = "column",
          key.xlab = "scale log RSEM",
          cexRow = 0.01,
          margins = c(5,5),
          breaks = seq(-3,3,length.out = 101),
          #scale = "row",
          col = bluered,
          main = paste("TCGA GBM and LGG bulk RNA-seq"))
par(lend = 1)           # square line ends for the color legend
legend(-0.02, 0.8,       # location of the legend on the heatmap plot
       legend = GBMLGG@meta.data[hc$labels,"subtype"] %>% as.vector %>% unique, # category labels
       col = GBMLGG@meta.data[hc$labels,"subtype.colors"] %>% as.vector %>% unique,  # color key
       lty= 1,             # line style
       lwd = 10            # line width
)
dev.off()


# Doheatmap =================
markers <- HumanGenes(Glioma,c("TP53","ATRX","CIC","FUBP1","NOTCH1",
                               "NIPBL","PIK3CA","PIK3R1","PDGFRA",
                               "CDKN2A","CDKN2B","CCND2","CDK4","RB1",
                               "MDM4","MDM2","EGFR","NF1","BRAF",
                               "PIK3R1","FGF13","PTEN"),unique = T)
df_samples <- readxl::read_excel("doc/181002_Single_cell_sample list.xlsx")

jpeg(paste0(path,"Doheatmap.jpeg"), units="in", width=10, height=7,res=600)
DoHeatmap(object = Glioma, genes.use = markers, group.label.rot =T,
          cex.row = 8, remove.key =T, use.scaled =FALSE,
          title = paste("Expression heatmap of marker genes in Glioma"))
dev.off()











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
