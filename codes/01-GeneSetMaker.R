# Library Imports; Make sure all packages are installed
libList <- c("CHRONOS","pathview","KEGGgraph","Rgraphviz","igraph",
             "org.Mm.eg.db","tm","stringi","stringr","dplyr","sna",
             "RBGL","tidyr","AnnotationDbi","org.Hs.eg.db","annotate",
             "biomaRt","Hmisc","broom","xtable")


for(i in libList) {
    if(!require(i)){
        install.packages(i)
    }
    if(!require(i)){
        BiocManager::install(i)
    }
}

library(CHRONOS)
library(pathview)
library(KEGGgraph)
library(Rgraphviz)
library(igraph)
library(org.Mm.eg.db)
library(tm)
library(stringi)
library(stringr)
library(dplyr)
library(sna)
library(RBGL)
library(tidyr)
library("AnnotationDbi")
library("org.Hs.eg.db")
library(annotate)
library(biomaRt)
library(Hmisc)
library(broom)
library(xtable)

if(!dir.exists("PATHWAYLIST/")){
    dir.create("PATHWAYLIST")
    a <- CHRONOS::downloadKEGGPathwayList("hsa")
    pathview::download.kegg(pathway.id = a$Id, kegg.dir = "./PATHWAYLIST/")
    paths_list  <- as_data_frame(a)
    write.csv(paths_list, file = "data/pathlist.txt", row.names = F)
}


a <- read.table("data/pathlist.txt", sep = ",", header = T,
                colClasses = c("character","character"))
# Retreiving Pathway graphs and information for KGML files
paths.address  <- paste("PATHWAYLIST/hsa",a$Id,".xml",sep = "")

graphs.list    <- sapply(paths.address, function(X)(
    KEGGgraph::parseKGML2Graph(X, expandGenes=TRUE)))




pathway.objs       <- sapply(paths.address, KEGGgraph::parseKGML)
pathway.titles     <- sapply(pathway.objs, getTitle)
pathway.titles     <- unname(pathway.titles)
names(graphs.list) <- paste0("Pathway.KEGG_",pathway.titles)

graphs.hsap    <-  sapply(graphs.list, function(X) ({
    temp <- translateKEGGID2GeneID(nodes(X))
    nodes(X) <- unname(temp)
    return(X)
}))

graphs.nodes  <-  sapply(graphs.list, function(X) ({
    temp <- translateKEGGID2GeneID(nodes(X))
    return(temp)
}))


saveRDS(graphs.nodes,"KEGG_PCxN/KEGG_PCxN.RDS")

