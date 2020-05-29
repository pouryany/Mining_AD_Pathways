rm(list = ls())
library(dplyr)

ADlists = c("ADCL","open_targets")
for (i in ADlists) {
    enrich.table <- readRDS("data/ADCL_ENRICHMENT.RDS")
    
    enrich.table <- enrich.table %>% filter(.,x ==i)
    path.dir <- read.csv("new_data/Overall_rank_minusDUPandADDs.txt", sep = "\t") 
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
    write.csv(t,paste0("outData/",i,"_Enrichment.csv"))
    
}



# colnames(t)
# 
# t <- t[t$adj_p < 0.05,]
# 
# library(ggplot2)
# 
# ggplot(t, aes(Overall_rank, -log(adj_p))) +
#     geom_point() +      
#     geom_smooth(method = "lm", se = FALSE) +
#     labs( x = "Overall rank according to text mining",
#           y = "Enrichment score",
#           title = "KEGG Pathways with significant enrichment\nwithin AD genes from GWAS")+
#     theme_classic(base_size = 20)
# 
# cor(t$Overall_rank,-log2(t$adj_p),method = "spearman")
# 
# 
# 

