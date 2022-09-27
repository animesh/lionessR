#setup####
libLocal <- "/home/ash022/scripts/Singh/rLib"
install.packages("infotheo", lib = libLocal)
install.packages("BiocManager", lib = libLocal)
BiocManager::install("GenomicRanges", lib = libLocal)
BiocManager::install("MatrixGenerics", lib = libLocal)
BiocManager::install("matrixStats", lib = libLocal)
BiocManager::install("SummarizedExperiment", lib = libLocal)
BiocManager::install("lionessR", lib = libLocal)
BiocManager::install("igraph", lib = libLocal)
library("infotheo", lib = libLocal)
library("matrixStats", lib = libLocal)
library("GenomicRanges", lib = libLocal)
library("MatrixGenerics", lib = libLocal)
library("SummarizedExperiment", lib = libLocal)
library("lionessR", lib = libLocal)
library("igraph", lib = libLocal)
options(nwarnings = 1000000)
summary(warnings(1000000))
#data####
M <- read.csv("/home/ash022/scripts/Singh/dataTmm.csv", header = T, row.names = 1)
dat <- as.matrix(M)
summary(dat)
#MI####
netFunMI <- function(x, ...) {
    mutinformation(discretize(t(x), disc = "equalwidth", nbins = 100), method = "emp")
}
#LIONESS####
cormat <- lioness(dat, netFunMI)
head(cormat)
z <- cormat@assays@data[[1]]
rownames(z) <- rownames(cormat)
head(z)
toptable_edges <- t(matrix(unlist(c(strsplit(row.names(z), "_"))), 2))
head(toptable_edges)
#sample####
chkGene <- "LAGE3"
for (patients in 1:ncol(z)) {
    print(paste0("patient", patients, chkGene))
    z1 <- cbind(toptable_edges, data.frame(z[, patients]))
    z2 <- z1[abs(z1[, 3]) > 1, ]
    g <- graph.data.frame(z2, directed = F)
    dg <- degree(g)
    print(summary(dg))
    dg_max <- sort.int(dg, decreasing = T, index.return = FALSE)
    dg_df <- as.matrix(dg_max)
    print(dg_df[rownames(dg_df) == chkGene])
    print(head(dg_df))
    save(g, file = paste0("patient", patients, ".MI.RData"))
    write.csv(dg_df, file = paste0("patient", patients, ".MI.csv"))
}
summary(warnings(1000000))
savehistory(file = "setup.Rhistory")
#end####