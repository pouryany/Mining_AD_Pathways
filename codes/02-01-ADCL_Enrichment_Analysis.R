# Enrichment analysis of AD pathways and AD genelists. 


rm(list=ls())
library(dplyr)
library(parallel)
library(clusterProfiler)
library(org.Hs.eg.db)




path.thresh  <- 9
genesets.dir <- "data/KEGG_PCxN.RDS"
ADsets.dir  <- "data/ADGeneSets.RDS"  
num.cors     <- detectCores()/2
#input.dir    <-"Data/temp_pourya/"
output.dir   <- "data"

if(!dir.exists(output.dir)){
    dir.create(output.dir)
}




    
    
    
    kegg.paths   <- readRDS(genesets.dir)
    kegg.sel     <- sapply(kegg.paths,length)
    kegg.paths   <- kegg.paths[kegg.sel > path.thresh]
    kegg.ref     <- Reduce(union,kegg.paths)
    AD.sets2     <- readRDS(ADsets.dir)
    AD.sets3     <- lapply(AD.sets2, function(X){
        X[X %in% kegg.ref]})

    
    



    
    ##
    ###
    sel.vec      <- sapply(AD.sets3,length)
    AD.sets3    <- AD.sets3[sel.vec > path.thresh]
    
    
    options(stringsAsFactors = F)
    iterator     <- (merge(names(AD.sets3),names(kegg.paths)))
    names(iterator)
    
    #iterator <- as.character(iterator)
    
    iterator     <- iterator %>% mutate_all(.,as.character)
    all          <- length(kegg.ref)
    
    
    
    
    
    enrichs   <- mclapply(1:nrow(iterator), function(Y){
        
        X <- iterator[Y,]
        q <- length(intersect(unlist(kegg.paths[X[[2]]]), 
                              unlist(AD.sets3[X[[1]]])))  
        
        m <- length(unlist(AD.sets3[X[[1]]]))
        n <- all - m
        k <- length(unlist(kegg.paths[X[[2]]]))
        pval <- phyper(q-1,m,n,k,lower.tail = F,log.p = F)
        return(c("pval" = pval,
                 "Intersect"    = q, 
                 "ADset_Size"  = m,
                 "not_ADset"   = n, 
                 "pathway_Size" = k,
                 "ratio_in"     = q/(k-q),
                 "ratio_out"    = (m-q)/n,
                 "ratio_ratios" = (q/(k-q))/((m-q)/n) ))
    }, mc.cores = num.cors )
    
    
    temp      <- do.call(rbind,enrichs)
    temp      <- as.data.frame(temp)
    
    
    iterator  <- cbind(iterator,temp)
    
    
    iterator  <- iterator %>% group_by(.,x) %>% 
                 mutate(., adj_pval_listWise = p.adjust(pval, method = "fdr"))
    
    
    save.dir <- paste0(output.dir,"/ADCL_ENRICHMENT.RDS")
    saveRDS(iterator,save.dir)
    write.csv(iterator,"data/ADCL_EnrichmentAnalysis.csv")
    
    
    #saveRDS(iterator,"Data/temp_pourya/PATHWAY_KEGG_ADSET_ENRICHMENT.rds")
    
    AD.sets3$open_targets

    
    str        <- read.csv("data/ADCL.csv",header = F)
    str        <- str$V1
    names(str) <- "ADCL"
    
    str        <- bitr(unlist(str),fromType = "SYMBOL",
                       toType = "ENTREZID",OrgDb = "org.Hs.eg.db")
    
    
    str <- str[str$ENTREZID %in% AD.sets3$ADCL,]
    write.csv(str,"outData/ADCL_geneTable.csv")
    
    
  
    
    #open.targets <- read.csv("data/AD_genes_Open_Targets.txt", sep = "\t")
    open.targets0<- read.csv("data/my_AD_list.csv")
    open.targets <- read.csv("data/my_AD_list.csv")
    open.targets <- open.targets[open.targets$association_score.overall >= 0.1,]
    open.targets <- open.targets$target.gene_info.symbol
    open.targets <- bitr(unlist(open.targets),fromType = "SYMBOL",
                         toType = "ENTREZID",OrgDb = "org.Hs.eg.db")
    
    open.targets <- open.targets[open.targets$ENTREZID %in% AD.sets3$open_targets,]
    
    
    open.targets0 <-open.targets0[open.targets0$target.gene_info.symbol %in% open.targets$SYMBOL,]
    nrow(open.targets0)
    write.csv(open.targets0,"outData/openTargets_geneTable.csv")
    



