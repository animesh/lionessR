### LIONESS ##

#############################################
### Pearson correlation between MM and ML ###
#############################################


# remove loaded object form current R session
rm(list=ls(all=T))
ls()

## set path for working directory
setwd("C:/Users/asi073/Desktop/Singh/")
path  <- getwd()
setwd(path)

#setup####
#sudo apt install gfortran libblas-dev liblapack-dev
options(nwarnings = 1000000)
#summary(warnings(1000000))
set.seed(2042)
#dirLocal <- getwd()


library("BiocGenerics")
library("S4Vectors")
library("IRanges")
library("GenomeInfoDb")
library("Biobase")
library("matrixStats")
library("GenomicRanges")
library("MatrixGenerics")
library("SummarizedExperiment")
library("lionessR")
library("igraph")
library("limma")
library("reshape2")

#data####
M <- read.csv("dataTmmS42.csv", header = T, row.names = 1)
#select####
dat <- as.matrix(M)
dim(dat)
hist(dat)
summary(dat)
rN <- rownames(dat)
cN <- colnames(dat)
cormat<-readRDS("cfunctionxnnstatscortxmethodpearsonn.cormat.RDS")
head(cormat)
#sampleN####
## label sample group names
sampleN <- c('MGUS1','MGUS2','MGUS3','MGUS4','MGUS5','MGUS6','MGUS7','MGUS8','MGUS9',
             'Ml1','Ml2','Ml3','Ml4','Ml5','Ml6','MM7','MM8','MM9','MM10','MM11','MM12',
             'MM13','MM14','MM15','MM16','MM17','MM18','MM19','MM20','MM21','MM22','MM23',
             'MM24','MM25','MM26','MM27','MM28','MM29','MM30','MM31','MM32','MM33','MM34',
             'MM35','MM36','MM37')
grpSample <- c('MGUS','MGUS','MGUS','MGUS','MGUS','MGUS','MGUS','MGUS','MGUS',
               'Ml','Ml','Ml','Ml','Ml','Ml','MM','MM','MM','MM','MM','MM','MM',
               'MM','MM','MM','MM','MM','MM','MM','MM','MM','MM','MM','MM','MM',
               'MM','MM','MM','MM','MM','MM','MM','MM','MM','MM','MM')
## add matching sample names in trait object
trait<-data.frame(cbind(sampleN,grpSample))
head(trait)
dim(M)
colnames(M) <- sampleN
#[1] 3956   46
head(M)

#sampleMM####
datMM<-data.frame(M[,trait$grpSample == "MM"])
summary(datMM)
dim(datMM)
nsel = nrow(M)
#nsel <- 500
cvar <- apply(as.array(as.matrix((datMM))), 1, sd)
sum(cvar==0)
length(cvar)
hist(cvar)
sum(datMM==0)
netyes <- cor(t(datMM))
dim(netyes)
hist(netyes)

#sampleML####
datMG<-M[,trait$grpSamp == "Ml"]
summary(datMG)
cvar <- apply(as.array(as.matrix((datMG))), 1, sd)
sum(cvar==0)
hist(cvar)
datMG <- as.matrix(datMG)
summary(datMG)
dim(datMG)
netno  <- cor(t(datMG))
dim(netno)
hist(netno)

netdiff <- netyes-netno
hist(netdiff)
dim(netdiff)
zeroDiff=netdiff==0
sum(zeroDiff)

## convert these adjacency matrices to edgelists
cormat2 <- rep(1:nsel, each=nsel)
cormat1 <- rep(1:nsel,nsel)
el      <- cbind(cormat1, cormat2, c(netdiff))

hist(el[,3])
melted  <- melt(upper.tri(netdiff))
summary(melted)
melted  <- melted[which(melted$value),]
values  <- netdiff[which(upper.tri(netdiff))]
melted  <- cbind(melted[,1:2], values)
summary(melted)

genes      <- row.names(netdiff)
melted[,1] <- genes[melted[,1]]
melted[,2] <- genes[melted[,2]]


row.names(melted) <- paste(melted[,1], melted[,2], sep="_")
tosub <- melted
hist(tosub[,3])
length(tosub[,3])
corThr<-1
tosel <- row.names(tosub[which(abs(tosub[,3])>corThr),])
length(tosel)
corsub <- cormat[which(row.names(cormat)  %in% tosel) ,]
CC <- corsub@assays@data[[1]]
rownames(CC) <- rownames(corsub)
head(CC)
hist(CC)
dim(CC)
group <- factor(trait[,2])
design <- model.matrix(~0+group)
cont.matrix <- makeContrasts(yesvsno = (groupMM - groupMl), levels = design)
fit <- lmFit(CC, design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2e <- eBayes(fit2)
toptable <- topTable(fit2e, number=nrow(CC), adjust.method="fdr")
head(toptable,10)
hist(toptable$logFC)
hist(toptable$P.Value)
hist(toptable$adj.P.Val)
hist(toptable$t)
summary(toptable$t)
toptable_edges <- t(matrix(unlist(c(strsplit(row.names(toptable), "_"))),2))
N <- nrow(CC)

N<-100
z <- cbind(toptable_edges[1:N,], as.data.frame(toptable$logFC[1:N]))
summary(z)
hist(z[,3])
#z <- z[abs(z[,3]) > 1.5, ]
g <- graph.data.frame(z, directed=F)
g <- delete.vertices(simplify(g), degree(g)==0 )
g <- simplify(g, remove.multiple = T, remove.loops = T)
g <- delete.vertices(g, which(E(g)$comm.freq == 4))


cl <- components(g)
ccG<-lapply(seq_along(cl$csize)[cl$csize > 2], function(x) V(g)$name[cl$membership %in% x])
length(ccG)
ccG

#ccG1<-ccG[[1]]

for (i in 1:length(ccG)){
  #print(ccG[[i]])
  maxlen <- max(c(ccG[[i]]))   
  genelist<- unlist(ccG[[i]])
  #print(genelist)
  print(strsplit(genelist, "[.]"))
}


#f <- function(s) strsplit(s, ".")
#tt<-sapply(ccG, f)

x <- sapply(ccG, function(ccG) strsplit(ccG, "-")[[1]], USE.NAMES=FALSE)

#figGeneList <- read.table("EIF2.txt", header = F)
geneppr <- read.table("geneListppr.txt", header = F)

##---------------------------------------------------------

#Exploratory analysis

colnames(geneppr) <- "name"
#length(unique(geneppr$name))
#geneppr$name %in% rownames(M)

#strsplit(rownames(M),'.')

b <- gsub("\\..*","",V(g)$name)
which(geneppr$name %in% b)
geneppr$name[geneppr$name %in% b]

#length(unique(figGeneList$name)) == nrow(figGeneList)
#n_occur <- data.frame(table(geneppr$name))
#n_occur[n_occur$Freq > 1,]
#figGeneList[figGeneList$name %in% n_occur$Var1[n_occur$Freq > 1],]

##---------------------------------------------------------------------

## Color the edges using limma T-statistics

labelOffset <- strwidth(names(V(g)),units = 'in')
vertex.label.dist = labelOffset*2.6
#E(g)$weight <- abs(z[,3])
tt <- unlist(ccG)
#tt <- tt[34:length(ccG)]
#g <- delete_vertices(g, tt)
#plot(g)
E(g)$weight<- as.numeric(z[,3])
E(g)$color[E(g)$weight<0] <- "blue"
E(g)$color[E(g)$weight>0] <- "red"

topgeneslist <- unique(c(toptable_edges[1:N,]))
fit <- lmFit(M, design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2e <- eBayes(fit2)
topDE <- topTable(fit2e, number=nrow(M), adjust="fdr")
hist(topDE$logFC)
topDE <- topDE[which(row.names(topDE) %in% topgeneslist),]
topgenesDE <- cbind(row.names(topDE), topDE$t)

##----------------------------------------------------------------------

#We color nodes based on the t-statistic from the LIMMA analysis:

# add t-statistic to network nodes
nodeorder <- cbind(V(g)$name, 1:length(V(g)))
nodes <- merge(nodeorder, topgenesDE, by.x=1, by.y=1)
nodes <- nodes[order(as.numeric(as.character(nodes[,2]))),]
nodes[,3] <- as.numeric(as.character(nodes[,3]))
nodes <- nodes[,-2]
V(g)$weight <- nodes[,2]

# make a color palette
mypalette4 <- colorRampPalette(c("blue","white","white","red"), space="Lab")(256)
breaks2a <- seq(min(V(g)$weight), 0, length.out=128)
breaks2b <- seq(0.00001, max(V(g)$weight)+0.1,length.out=128)
breaks4 <- c(breaks2a,breaks2b)

# select bins for colors
bincol <- rep(NA, length(V(g)))
for(i in 1:length(V(g))){
  bincol[i] <- min(which(breaks4>V(g)$weight[i]))
}
bincol <- mypalette4[bincol]

# add colors to nodes
V(g)$color <- bincol

par(mar=c(0,0,0,0))
layout <- layout.reingold.tilford(g, circular=T)
plot(g, vertex.label.cex=0.5, directed=F,
     vertex.size=3, 
     vertex.label.color = "black",
     vertex.label.font=2, vertex.label.dist = labelOffset*1.6,
     edge.width=E(g)$width,
     #edge.width=10*(abs(as.numeric(z[,3]))-0.7),
     #vertex.color=V(g)$color,
     layout= layout)
#layout=layout_with_graphopt)
#layout = layout_as_tree)

hitGeneMM_ML <- geneppr$name[geneppr$name %in% b]
mylist <- hitGeneMM_ML 
write.table(mylist,"hitGeneMM_ML.txt",sep=",")
write.table(b,"GeneMM_ML.txt", sep = ",")
