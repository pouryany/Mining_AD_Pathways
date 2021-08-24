rm(list = ls())
library(dplyr)

ADlists = c("ADCL","ADGeneList","open_targets")
for (i in ADlists) {
    enrich.table <- readRDS("data/ADCL_ENRICHMENT.RDS")
    
    enrich.table <- enrich.table %>% filter(.,x ==i)
    path.dir <- read.csv("data/Overall_rank_minusDUPandADDs.txt", sep = "\t") 
    path.dir$Pathway  <- paste0("Pathway.KEGG_",path.dir$Pathway)
    
    t <- enrich.table %>% filter(.,y %in%path.dir$Pathway ) %>%
        mutate(., adj_p = p.adjust(pval, method = "fdr")) 
    t <- t %>% mutate(.,enrich_rank = rank(pval), 
                      is_rich = ifelse(adj_p < 0.05,T,F))
    t <- dplyr::left_join(path.dir,t,by = c("Pathway"= "y"))
    t <- t[!is.na(t$x),]
    
    
    
    rich_paths <- t[t$is_rich,]$Overall_rank
    rest_paths <- t[!t$is_rich,]$Overall_rank
    
    wilcox.test(rich_paths,rest_paths,alternative = "less")
    
    t <- t[,c("Pathway","Overall_rank","pval","Intersect","ADset_Size",
              "pathway_Size","adj_p","enrich_rank")   ]   
    t <- t[order(t$adj_p),]
    write.csv(t,paste0("outData/2020_12_01_",i,"_Enrichment.csv"))
    
}



