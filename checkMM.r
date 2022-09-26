summary(warnings(1000000))
libLocal <- "/home/ash022/scripts/Singh/rLib"
BiocManager::install("GenomicRanges", lib = libLocal)
BiocManager::install("MatrixGenerics", lib = libLocal)
BiocManager::install("matrixStats", lib = libLocal)
BiocManager::install("SummarizedExperiment", lib = libLocal)
BiocManager::install("lionessR", lib = libLocal)
BiocManager::install("igraph", lib = libLocal)
install.packages("infotheo", lib = libLocal)
library("infotheo", lib = libLocal)
library("matrixStats", lib = libLocal)
library("GenomicRanges", lib = libLocal)
library("MatrixGenerics", lib = libLocal)
library("SummarizedExperiment", lib = libLocal)
library("lionessR", lib = libLocal)
library("igraph", lib = libLocal)
M <- read.csv("/home/ash022/scripts/Singh/dataTmm.csv", header = T, row.names = 1)
cvar <- apply(as.array(as.matrix(M)), 1, sd)
nsel <- nrow(M)
dat <- cbind(cvar, M)
dat <- dat[order(dat[, 1], decreasing = T), ]
dat <- dat[1:nsel, -1]
dat <- as.matrix(dat)
rN <- rownames(dat)
cN <- colnames(dat)
datDisc <- discretize(t(dat))
rownames(datDisc) <- cN
colnames(datDisc) <- rN
summary(datDisc)
testMI = mutinformation(datDisc, method = "emp")
head(testMI)
library(lionessR) # , help, pos = 2, lib.loc = NULL)
netFunMI <- function(x, ...) {
    mutinformation(datDisc, method = "emp")
}
dim(datDisc)
cormat <- lioness(as.matrix(dat), netFunMI)
summary(warnings(1000000))
head(cormat)
cormat1 <- cormat
cormat_ani <- cormat1@assays@data[[1]]
rownames(cormat_ani) <- rownames(cormat1)
z <- cormat_ani
head(z)
toptable_edges <- t(matrix(unlist(c(strsplit(row.names(z), "_"))), 2))
head(toptable_edges)
chkGene <- "LAGE3"
for (patients in 1:ncol(z)) {
    print(paste0("patient", patients, chkGene))
    z1 <- cbind(toptable_edges, data.frame(z[, patients]))
    z2 <- z1[abs(z1[, 3]) > 1, ]
    g <- graph.data.frame(z2, directed = F)
    dg <- degree(g)
    dg_max <- sort.int(dg, decreasing = T, index.return = FALSE)
    dg_df <- as.matrix(dg_max)
    print(dg_df[rownames(dg_df) == chkGene])
    save(g, file = paste0("patient", patients, ".MI.RData"))
    write.csv(dg_df, file = paste0("patient", patients, ".MI.csv"))
}
summary(warnings(1000000))
savehistory(file = "setup.Rhistory")
