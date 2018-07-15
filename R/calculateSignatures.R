calculateSignatures.1 <- function (sc_data=sc.data.gl, species = "Human", signatures = NULL, n.break = 1000) 
{
        sc_data = as.matrix(sc_data)
        data("signatures")
        if (is.null(signatures)) {
                if (species == "Human") {
                        egc = human.egc
                }
                else if (species == "Mouse") {
                        egc = mouse.egc
                }
                rownames(sc_data) = tolower(rownames(sc_data)) # re-write
        }
        else {
                egc = signatures
                if(species == "Human") rownames(sc_data) = toupper(rownames(sc_data)) # re-write
                if(species == "Mouse") rownames(sc_data) = Hmisc::capitalize(tolower(rownames(sc_data))) # re-write
        }
        numClusters = min(parallel::detectCores(all.tests = FALSE, logical = TRUE) - 
                                  1, 4)
        scores = matrix(NA, length(egc), ncol(sc_data))
        wind = seq(1, ncol(sc_data), by = n.break)
        print(paste("Using sets of", n.break, "cells. Running", 
                    length(wind), "times."))
        for (i in wind) {
                last = min(ncol(sc_data), i + n.break - 1)
                a = GSVA::gsva(sc_data[, i:last], egc, method = "ssgsea", 
                               ssgsea.norm = F, parallel.sz = numClusters, parallel.type = "FORK")
                scores[, i:last] = a
        }
        mmin = matrixStats::rowMins(scores)
        mmax = matrixStats::rowMaxs(scores)
        scores = scores/(mmax - mmin)
        rownames(scores) = rownames(a)
        output = data.frame(t(scores))
        if (is.null(signatures) && species == "Human") {
                output$Cell_Cycle = rowMeans(t(scores[c("G1S", "G2M"), 
                                                      ]))
                output[, c("G1S", "G2M")] = c()
        }
        else {
                colnames(output) = c("G1/S", "G2/M", "M", "M/G1", "S")
        }
        output
}


library(GSEABase)
library(GSVAdata)
library(Biobase)
source("http://bioconductor.org/biocLite.R")
biocLite("genefilter")
library(genefilter)
library(limma)
library(RColorBrewer)
library(GSVA)
data(gbm_VerhaakEtAl)
gbm_eset
head(featureNames(gbm_eset))
table(gbm_eset$subtype)

data(brainTxDbSets)
sapply(brainTxDbSets, length)

lapply(brainTxDbSets, head)
gbm_es <- gsva(gbm_eset, brainTxDbSets, mx.diff=FALSE, verbose=FALSE, parallel.sz=1)
