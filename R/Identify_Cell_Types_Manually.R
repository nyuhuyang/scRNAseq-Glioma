library(Seurat)
library(dplyr)
source("../R/Seurat_functions.R")

(load(file = "./data/Glioma_Harmony_20181201.Rda"))
##############################
# geom_density
##############################

markers <- HumanGenes(Glioma,c("EGFR","FGF13","TP53","CCND2","ATRX"))
tests <- paste0("test",1:5)
for(test in tests){
    sample_n = which(df_samples$tests %in% test)
    samples <- unique(df_samples$samples[sample_n])
    print(samples)
    
    cell.use <- rownames(Glioma@meta.data)[Glioma@meta.data$orig.ident %in% samples]
    subset.Glioma <- SubsetData(Glioma, cells.use = cell.use)
    
    g <- split(rownames(subset.Glioma@meta.data), 
               subset.Glioma@meta.data[,"orig.ident"]) %>% 
        lapply(function(cells_use) {
    single.Glioma <- SubsetData(subset.Glioma, cells.use = cells_use)
    sample <- unique(single.Glioma@meta.data$orig.ident)
    data.use <- single.Glioma@data[markers,] %>% as.matrix %>% t %>% as.data.frame %>%
            tidyr::gather(key = markers, value = ave.expr)
    ggplot(data.use, aes(x = ave.expr, color = markers)) + 
            geom_density(size = 1) +
            scale_y_sqrt() + ylim(0, 1)+
            xlab("log nUMI")+
            ggtitle(sample)+
            theme(text = element_text(size=15),
                  #legend.position="none", 
                  legend.position=c(0.3,0.85) ,
                  plot.title = element_text(hjust = 0.5,size = 15, face = "bold"))
    })
    jpeg(paste0(path,"density_",test,".jpeg"), units="in", width=10, height=7,res=600)
    print(do.call(plot_grid,c(g,nrow = 1))+ #   plot_grid(g[[1]],g[[5]],g[[4]],g[[3]],g[[2]])
              ggtitle("Density plot for markers in Glioma")+
              theme(text = element_text(size=15),							
                    plot.title = element_text(hjust = 0.5,size = 15, face = "bold")))
    dev.off()
}


##############################
# SingleFeaturePlot
##############################
df_samples <- readxl::read_excel("doc/181002_Single_cell_sample list.xlsx")
tests <- paste0("test",c(1:5))
markers <- HumanGenes(Glioma,c("TP53","EGFR","CDKN2A","FGF13","PTEN"))
for(test in tests){
    sample_n = which(df_samples$tests %in% test)
    samples <- unique(df_samples$samples[sample_n])
    print(samples)
    
    cell.use <- rownames(Glioma@meta.data)[Glioma@meta.data$orig.ident %in% samples]
    subset.Glioma <- SubsetData(Glioma, cells.use = cell.use)
    SplitSingleFeaturePlot(subset.Glioma,group.by = "ident",split.by = "orig.ident",
                       no.legend = T,label.size=3,do.print =T,markers = markers,
                       threshold = 0.1)
}

#====== 2.1 identify phenotype for each cluster  ==========================================
lnames = load(file = "./data/Glioma_alignment.Rda")
lnames

featureplot <- function(x,object = Glioma,...){
    p <- FeaturePlot(object = object, 
                     reduction.use = "tsne",
                     features.plot = x, min.cutoff = NA, 
                     cols.use = c("lightgrey","blue"), pt.size = 0.5,...)
    return(p)
}
# Cell Type Marker gene database
Housekeeping <- HumanGenes(Glioma,c("RNR2","Rpl4","Actb","Gnas",
                                    "Tubb","Kras","Calb1"))
# Stem cell=======
Embryonic_SCs <- HumanGenes(Glioma,c("ALPL","ALPP","ALPI","ALPPL2","TNFRSF8","TDGF1",
                                   "PODXL","NR6A1","POU5F1","B3GALNT1",
                                   "FUT1","FUT2","KIT","TERT","TERC","TEP1",
                                   "NANOG","SOX2","Dppa5","Klf4","zfp42","Spp1"))
Ectoderm <- HumanGenes(Glioma,c("NES","NCAM1","PAX6","VIM"))# VIM
Mesoderm <- HumanGenes(Glioma,c("BMP4","TBXT"))
Endoderm <- HumanGenes(Glioma,c("AFP","GATA4","HNF4A"))
Pluripotent_Stem_Cells <- unique(c(Embryonic_SCs,Ectoderm,Mesoderm,Endoderm))

mouse_lin <- HumanGenes(Glioma,c("CD2","CD3G","CD3D","CD4","CD5","CD8A","KLRB1","PTPRC","Ly76","Ly6g"))
human_lin <- HumanGenes(Glioma,c("CD3G","CD3D","CD14","FCGR3A","CD19","MS4A1","NCAM1"))
lineage_Cells <- unique(c(mouse_lin,human_lin))

MSC <- HumanGenes(Glioma,c("BMPR2","BMPR1A","BMPR1B","CD34","ATXN1","CD44","THY1",
                         mouse_lin,human_lin),unique = T)
HSC <- HumanGenes(Glioma,c("CD34","CD38","KIT","ATXN1","THY1",
                         mouse_lin,human_lin),unique = T)
Collagen <- HumanGenes(Glioma,c("COL2A1","COL1A1","COL1A2"),unique =T)


# Blood Vessel=====
Endothelium <- HumanGenes(Glioma,c("Cdh5","Pecam1","Flt1","Plvap","Kdr","ptprb",
                                 "Vwf","EMCN","Car4","VEGFA"))
Smooth_muscle_cells <- HumanGenes(Glioma,c("Acta2","MYH11"))
# Bone ==
Osteoblast <- HumanGenes(Glioma,c("ALPL","SPP1","IBSP","BGLAP"))

# Fat
Adipocytes <- HumanGenes(Glioma,c("TRIP4","ASCC1","SLC36A2","P2RX5","MYF5","UCP1"))
Fat_stem_cell <- HumanGenes(Glioma,c("FABP3","FABP1","FABP2","FABP7","FABP4",
                                   "FABP5","FABP12","FABP6","SLC27A1","SLC27A2",
                                   "SLC27A3","SLC27A4","SLC27A5","SLC27A6"))
# Gender
Y_chromosome <- HumanGenes(Glioma,c("UTY","DDX3Y"))
X_chromosome <- HumanGenes(Glioma,c("Xist"))

# Liver
Hepatocyte <- HumanGenes(Glioma,c("ALB","ITGB1"))

# Nervous System
Neuron2A <- HumanGenes(Glioma,c("Tmod2","Skil","Slc30a1","Erbb2ip","PCDHA5",
                                "Vgf","Gabrb3"))

RPE <- HumanGenes(Glioma,c("Rpe65","Rlbp1"))

S1Pyr <-HumanGenes(Glioma,c("Tbr1","Rasgrf2","RASGRF1","Pvrl3","Cux2","Rorb","Plcxd2",
                            "Thsd7a","Kcnk2","Cplx3","Sulf2","Foxp2","Syt6","Rprm",
                            "Nr4a2","Synpr","Pcp4"))

# Bone Marrow and Blood ===
Mesenchymal <- HumanGenes(Glioma,c("Pdgfrb","Vim","Has2","Dcn"))

#--Hematopoietic----
Hematopoietic <- HumanGenes(Glioma,c("PTPRC","LAPTM5","SRGN"))
#------Myeloid----
megakaryocytes <-  HumanGenes(Glioma,c("PPBP","GNG11"))
erythrocyte <-  HumanGenes(Glioma,c("HBA2","HBB"))
MastCells <- HumanGenes(Glioma,c("Cma1","Mcpt4","Tpsb2","Cpa3"))
Neutrophil <- HumanGenes(Glioma,c("ADAM8","MSMO1","FUT4","FCGR3A","CEACAM8"))
CD14_Monocytes <-  HumanGenes(Glioma,c("CD14","LYZ","S100A9","CCL2"))
CD16_Monocytes <- HumanGenes(Glioma,c("FCGR3A","MS4A7","VMO1"))
Monocytes <- unique(c(CD14_Monocytes,CD16_Monocytes))
Macrophages <- HumanGenes(Glioma,c("LYZ","CD68","MARCO","EMR1","ITGB2"))
DendriticCells <- HumanGenes(Glioma,c("Itgax","GPR183","CST3","HLA-DQA1","FCER1A","TSPAN13",
                                   "IL3RA","IGJ"))
Platelet <- HumanGenes(Glioma,c("VWF","CD36","ITGB3","GP9","GP6","GP1BA","SELP","PDGFRB")) # CD36 is not specific
Myeloid_all <-  HumanGenes(Glioma,c(megakaryocytes,erythrocyte,MastCells,
                             CD14_Monocytes,CD16_Monocytes,Macrophages,DendriticCells),unique=T)
#------Lymphoid----
Lymphoid <- HumanGenes(Glioma,c("Cd19","CD79A","MS4A1",
                             "GNLY","Ncr1","CCL5","KLRD1","NKG7"))
# T cell
T_Cell <- HumanGenes(Glioma,c("CD2","CD3G","CD3D","CD4","CD8A","IL2RA","FOXP3",
                           "IL7R","SELL","IL2RG","GIMAP5"))

Treg <- HumanGenes(Glioma,c("FOXP3","CD4","IL2RA","CTLA4","PDCD1","ENTPD1","CD38",
                         "ICOS","TNFSF9","TNFRSF9"))
CD4_Naive_T <- HumanGenes(Glioma,c("CD4","IL7R","GIMAP5","SELL","IL2RG"))
Natural_killer_T <- HumanGenes(Glioma,c("NKG7","CCL5","NCAM1","FCGR3A","Ncr1","KLRD1"))

T_Cell_all <- unique(c(T_Cell,Treg,CD4_Naive_T,Natural_killer_T))
# B cell
B_Cell <-HumanGenes(Glioma,c("CD19","MS4A1","CD79A","CD40","CD22","FCER2","HLA-DRB1",
                          "CXCR4","SOX11","CD5","PAX5","CD27","IL4R"))

B_StemCell <- HumanGenes(Glioma,c("SPN","CD20"))
Pre_Pro_B <- HumanGenes(Glioma,c("CD34","MME","CD38"))
Pro_B <- HumanGenes(Glioma,c("MME","CD19","SPN","CD38","CD24","IL7","IL3RA"))
Pre_B <- HumanGenes(Glioma,c("MME","CD19","MS4A1","CD24","CD38","IL7","IL3RA","IL4R"))
Immature_B <- HumanGenes(Glioma,c("MME","CD19","MS4A1","CR2","CD40","CD24","CD38","IL4R"))
Transitional_B <- HumanGenes(Glioma,c("CD19","MS4A1","CD5","CR2","CD24","CD38"))
Marginal_zone_B <- HumanGenes(Glioma,c("CD1C","CD19","MS4A1","CR2","CD27"))
Regulatory_B <- HumanGenes(Glioma,c("CD1D","CD5","CD19","CR2","CD24"))
Follicular_B <- HumanGenes(Glioma,c("CD19","MS4A1","CR2","CD22","FCER2","CD24",
                                 "HLA-DRB1","HLA-DQB1","HLA-DRA","HLA-DQA1"))
Activated_B <- HumanGenes(Glioma,c("CD27","CD19","MS4A1","IL2RA","TNFRSF8","CD69","CD80","CD86","FLT3"))
Germinal_center_B <- HumanGenes(Glioma,c("MME","CD19","MS4A1","FCER2","CD27","CD38","TNFRSF17"))
Plasma_blast <- HumanGenes(Glioma,c("CD19","CD38","CD27","TNFRSF17","HLA-DRB1"))
Plasma_cell_long_lived <- HumanGenes(Glioma,c("CXCR4","CD27","CD38","CD138","CD269"))
Memory_B <- HumanGenes(Glioma,c("CD19","MS4A1","CD40","CD27","CXCR4","CXCR5","ACKR3"))

B_Cell_all <- unique(c(B_Cell,B_StemCell,Pre_Pro_B,Pro_B,Pre_B,Immature_B,Transitional_B,
                       Marginal_zone_B,Regulatory_B,Follicular_B,Activated_B,Germinal_center_B,
                       Plasma_blast,Plasma_cell_long_lived,Memory_B))
Lymphoid_all <- unique(c(Lymphoid,T_Cell_all,B_Cell_all))
Hematopoietic_all <- unique(c(Hematopoietic,Myeloid_all,Lymphoid_all))

# other
Melanocytes <- HumanGenes(Glioma,c("Pmel","Mlana"))

Myelinating_Schwann_cells <- HumanGenes(Glioma,c("MBP","MPZ"))
Pericytes <- HumanGenes(Glioma,c("Pdgfrb","Cspg4","Anpep","Rgs5",
                              "Myh11","Mylk","Des","Vtn","Ifitm1"))
Stromal_fibroblasts <- HumanGenes(Glioma,c("DCN","COL6A1","TIMP3","PDGFRA"))
Neurons <- HumanGenes(Glioma,c("Ihh","Gli1", "Ptch1", "Hhip"))
CellCycle <- HumanGenes(Glioma,c("CCND1","CCND2","CCND3","CDK4","CDK6","PCNA","SOX11",
                               "RB1","E2F1","TK1","CCNA2","MKI67","CDK1"))
Fibroblast <- HumanGenes(Glioma,c("FGF1","FGF9","SFRP1"))
Epithelium <- HumanGenes(Glioma,c("Epcam","KRT19","KRT5",
                                  "MUC1","SCGB3A2","SCGB1A1","SCGB3A1","SFTPB","FOXJ1","Rpe65",
                                  "Rlbp1","Msln","Upk3b","Lrrn4"))

# "Single-cell RNA-seq supports a developmental hierarchy in human oligodendroglioma" nature20123
IDH_oligodendroglioma <- HumanGenes(Glioma,c("OLIG2", "OMG"))
IDH_astrocytoma <- HumanGenes(Glioma,c("APOE", "ALDOC", "SOX9", "GFAP"))
stemness <- HumanGenes(Glioma,c("SOX4","SOX11", "SOX2","CCND2"))
cell_proliferation <- HumanGenes(Glioma,c("MKI67"))
glioma_CSCs <- HumanGenes(Glioma,c("NFIB", "ASCL1", "CHD7", "CD24", "BOC", "TCF4"))
nature20123 <- c(IDH_oligodendroglioma,IDH_astrocytoma,stemness,cell_proliferation, glioma_CSCs)

SingleFeaturePlotSave(object = Glioma, feature.plot = nature20123,mkdir=T)

# "Single-cell RNA-seq highlights intratumoral heterogeneity in primary glioblastoma" science.1254257
MGH30 <- HumanGenes(Glioma,c("EGFR", "PDGFRA", "PDGFA","FGFR1", "FGF1", "NOTCH2", "JAG1"))
quiescence <- HumanGenes(Glioma,c("HES1", "TSC22D1", "KDM5B","NFIB"))
stem_cell <- HumanGenes(Glioma,c("PROM1"))
tumor_propagation <- HumanGenes(Glioma,c("POU3F2"))
neural_stem_cells_self_renewal <- HumanGenes(Glioma,c("NFIA"))
science.1254257 <- c(MGH30,quiescence,stem_cell,tumor_propagation,neural_stem_cells_self_renewal)

SingleFeaturePlotSave(object = Glioma, feature.plot = science.1254257,mkdir=T)

# "Decoupling genetics, lineages, and microenvironment in IDH-mutant gliomas by single-cell RNA-seq" science.aai8478
IDH_O <- HumanGenes(Glioma,c("OLIG2","TERT"))
IDH_A <- HumanGenes(Glioma,c("TP53","ATRX","GFAP"))
Microglia_macrophage <- HumanGenes(Glioma,c("CD14", "AIF1", "CSF1R"))
Oligodendrocytes <-  HumanGenes(Glioma,c("MBP", "MOBP", "PLLP", "CLDN11"))
neuro_developmental_TF= HumanGenes(Glioma,c("SOX4","SOX11", "TCF4"))
Microglia_like <- HumanGenes(Glioma,c("CX3CR1","P2RY12","P2RY13","SELPLG"))
Macrophage_like<- HumanGenes(Glioma,c("CD14","CD163","IFITM2","IFITM3","TAGLN2",
                                      "F13A1","TGFBI","IFNGR1"))

science.aai8478 <- unique(c(IDH_O,IDH_A,Microglia_macrophage,Oligodendrocytes,
                            neuro_developmental_TF,Microglia_like,Macrophage_like))
SingleFeaturePlotSave(object = Glioma, feature.plot = science.aai8478)
#############


Refs_TCGA_IvyGbm_main = read.csv("./output/Refs_TCGA_IvyGbm_main.csv",header = T,row.names = 1)

SearchMarker <- function(df, marker){
        result = apply(df, 2, function(x) which(grepl(marker, x)))
        result_df = data.frame(sort(unlist(result)))
        colnames(result_df) = marker
        return(result_df)
        }

SearchAllMarkers <- function(df, markers){
        results <- lapply(markers,function(x) SearchMarker(df,x))
        names(results) <- markers
        temp_df = results[[1]]
        for(i in 2:length(results)) {
                temp_df = merge(temp_df,results[[i]], by="row.names",all=TRUE)
                rownames(temp_df) = temp_df$Row.names
                temp_df = temp_df[,-1]
        }
        removeNA(temp_df)
        return(temp_df)
}
removeNA <- function(df){
        AllNA = apply(df,2,function(x) length(which(!is.na(x))))
        df = df[,AllNA !=0]
        return(df)
}

publised_markers <- HumanGenes(Glioma,sort(c(nature20123,science.1254257,
                                             science.aai8478)),
                               unique = T)
Summary <- SearchAllMarkers(df = Refs_TCGA_IvyGbm_main, 
                            markers = publised_markers)
IDH_astrocytoma = removeNA(Summary["Astrocytes",])
IDH_astrocytoma = colnames(IDH_astrocytoma)


SingleFeaturePlotSave(object = Glioma, feature.plot = IDH_astrocytoma)
SingleFeaturePlotSave(object = Glioma, feature.plot = publised_markers,mkdir=F,
                      folder.name = "ALL",threshold= 0.5)

SingleFeaturePlotSave(object = Glioma, feature.plot = "OLIG2", mkdir=F,
                      folder.name = "ALL",threshold= 1.5)
SingleFeaturePlotSave(object = Glioma, feature.plot = "SOX2", mkdir=F,
                      folder.name = "ALL",threshold= 2.0)
SingleFeaturePlotSave(object = Glioma, feature.plot = "ALDOC", mkdir=F,
                      folder.name = "ALL",threshold= 1.5)
SingleFeaturePlotSave(object = Glioma, feature.plot = "ATRX", mkdir=F,
                      folder.name = "ALL",threshold= 1.5)
SingleFeaturePlotSave(object = Glioma, feature.plot = "VIM", mkdir=F,
                      folder.name = "ALL",threshold= 5.)
SingleFeaturePlotSave(object = Glioma, feature.plot = stemness, mkdir=F,
                      folder.name = "ALL",threshold= 1.5)

Neurons = HumanGenes(Glioma,brainTxDbSets$neuronal_up)
Neurons
SingleFeaturePlotSave(object = Glioma, feature.plot = Neurons,mkdir=T,
                      threshold= 0.1)
SingleFeaturePlotSave(object = Glioma, feature.plot = "NEUROD6", mkdir=F,
                      folder.name = "Neurons",threshold= .1)
Astroglia = HumanGenes(Glioma,brainTxDbSets$astroglia_up)
SingleFeaturePlotSave(object = Glioma, feature.plot = Astroglia,mkdir=T,
                      threshold= 0.1)
SingleFeaturePlotSave(object = Glioma, feature.plot = "IGFBP3",mkdir=F,
                      folder.name = "Astroglia", threshold= 1.0)

Mesenchymal = HumanGenes(Glioma,Refs_TCGA_IvyGbm_main$Fibroblasts[1:100])
SingleFeaturePlotSave(object = Glioma, feature.plot = Mesenchymal,mkdir=T,
                      threshold= 0.1)
quiescence <- HumanGenes(Glioma,c("HES1", "TSC22D1", "KDM5B","NFIB"))
SingleFeaturePlotSave(object = Glioma, feature.plot = quiescence,mkdir=T,
                      threshold= 1.5)

Neurons = HumanGenes(Glioma,c("NEFL","GABRA1","SYT1","SLC12A5"))
SingleFeaturePlotSave(object = Glioma, feature.plot = "SYT1", mkdir=F,
                      folder.name = "Neurons",threshold= .1)
SingleFeaturePlotSave(object = Glioma, feature.plot = "SLC12A5", mkdir=F,
                      folder.name = "Neurons",threshold= .1)
