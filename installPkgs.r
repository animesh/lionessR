dirLocal <- getwd()
dir.create(file.path(dirLocal, "rLib"), showWarnings = FALSE)
libLocal <- paste(dirLocal,"rLib",sep="/")
reposLib="http://cran.us.r-project.org"
install.packages(c("infotheo","BiocManager","jsonlite"), lib = libLocal, repos = reposLib)
BiocManager::install(c("BiocGenerics","GenomeInfoDb","Biobase","S4Vectors","IRanges", "GenomicRanges","MatrixGenerics", "matrixStats", "SummarizedExperiment","reshape2","limma","lionessR","igraph"),  lib = libLocal)
