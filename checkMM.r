
#https://cloud.r-project.org/bin/linux/ubuntu/
#sudo vim /etc/apt/sources.list and remove cran-* lines
#sudo apt update -qq
#sudo apt install --no-install-recommends software-properties-common dirmngr
# add the signing key (by Michael Rutter) for these repos
# To verify key, run gpg --show-keys /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc 
# Fingerprint: E298A3A825C0D65DFD57CBB651716619E084DAB9
#wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | sudo tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc
# add the R 4.0 repo from CRAN -- adjust 'focal' to 'groovy' or 'bionic' as needed
#sudo add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/"
#sudo apt update 
#sudo apt upgrade 
#sudo apt install --no-install-recommends r-base
#install.packages("BiocManager")
#install.packages("jsonlite")
#BiocManager::install("lionessR")
M<-read.csv("dataTmm.csv",header=T,row.names=1)
cvar <- apply(as.array(as.matrix(M)), 1, sd)
# https://github.com/orgs/community/discussions/26316 for plot , unblock cookies
hist(cvar)
dat <- cbind(cvar, M)
dat <- dat[order(dat[,1], decreasing=T),]
dat <- dat[cvar>1, -1]
dat <- as.matrix(dat)
dim(dat)
hist(dat)
#library(lionessR)#, help, pos = 2, lib.loc = NULL)
cormat <- lionessR::lioness(as.matrix(dat))
dim(cormat)
cData <- cormat@assays@data[[1]]
rownames(cData) <- rownames(cormat)
dim(cData)
head(cData)
hist(cData)
cData<-cData[abs(cData[,3])>1,]
hist(cData)
# https://github.com/mararie/lionessR/blob/master/vignettes/lionessR.Rmd
edgeCdata <- t(matrix(unlist(c(strsplit(row.names(cData), "_"))), 2))
#edgeCdata <- unlist(c(strsplit(row.names(cData), "_")))
tail(edgeCdata)
z <- cbind(edgeCdata, cData)
#sudo apt install gfortran
#sudo apt install libblas-dev
#sudo apt install liblapack-dev
#BiocManager::install("igraph")
g <- igraph::graph.data.frame(z, directed = F)
plot(g)
dg<-igraph::degree(g)
sortDG<-sort.int(dg,decreasing=T,index.return=F)
max(sortDG)
min(sortDG)
plot(sortDG)
head(sortDG,25)
tail(sortDG, 25)
