summary(warnings(1000000))
BiocManager::install("GenomicRanges", lib = "/home/ash022/scripts/Singh/rLib")
BiocManager::install("MatrixGenerics", lib = "/home/ash022/scripts/Singh/rLib")
BiocManager::install("matrixStats", lib = "/home/ash022/scripts/Singh/rLib")
BiocManager::install("SummarizedExperiment", lib = "/home/ash022/scripts/Singh/rLib")
BiocManager::install("lionessR", lib = "/home/ash022/scripts/Singh/rLib")
BiocManager::install("igraph", lib = "/home/ash022/scripts/Singh/rLib")
library("matrixStats", lib = "/home/ash022/scripts/Singh/rLib")
library("GenomicRanges", lib = "/home/ash022/scripts/Singh/rLib")
library("MatrixGenerics", lib = "/home/ash022/scripts/Singh/rLib")
library("SummarizedExperiment", lib = "/home/ash022/scripts/Singh/rLib")
library("lionessR", lib = "/home/ash022/scripts/Singh/rLib")
library("igraph", lib = "/home/ash022/scripts/Singh/rLib")
M <- read.csv("/home/ash022/scripts/Singh/dataTmm.csv", header = T, row.names = 1)
cvar <- apply(as.array(as.matrix(M)), 1, sd)
nsel <- nrow(M)
dat <- cbind(cvar, M)
dat <- dat[order(dat[, 1], decreasing = T), ]
dat <- dat[1:nsel, -1]
dat <- as.matrix(dat)
netFunS <- function(x, ...) {
    stats::cor(t(x), method = "spearman")
}
cormat1 <- lioness(dat, netFunS)
cormat_ani <- cormat1@assays@data[[1]]
rownames(cormat_ani) <- rownames(cormat1)
z <- cormat_ani
toptable_edges <- t(matrix(unlist(c(strsplit(row.names(z), "_"))), 2))
# toptable_edges <- t(unlist(c(strsplit(row.names(z), "_"))))
for (patients in 1:ncol(z)) {
    print(paste0("patient", patients, ".RData"))
    z1 <- cbind(toptable_edges, data.frame(z[, patients]))
    z2 <- z1[abs(z1[, 3]) > 1, ]
    g <- graph.data.frame(z2, directed = F)
    dg <- degree(g)
    dg_max <- sort.int(dg, decreasing = T, index.return = FALSE)
    dg_df <- as.matrix(dg_max)
    print(dg_df[rownames(dg_df) == "LAGE3"])
    save(g, file = paste0("patient", patients, ".spearman.RData"))
    write.csv(dg_df, file = paste0("patient", patients, ".spearman.csv"))
}
summary(warnings(1000000))
savehistory(file = "setup.Rhistory")