# processing the AD curated list into a list
# The table can be found in PCxN publication Tables S.4 and S.5
# PCxN publication here: https://doi.org/10.1371/journal.pcbi.1006042
library(org.Hs.eg.db)
library(clusterProfiler)


str        <- readr::read_file("data/ADCL.txt")
str        <- gsub("\\s","",str)
str        <- (strsplit(str,split = ","))
names(str) <- "ADCL"

str        <- bitr(unlist(str),fromType = "SYMBOL",
                   toType = "ENTREZID",OrgDb = "org.Hs.eg.db")

saveRDS(str,"data/ADCL.RDS")


list2   <- readxl::read_excel("data/journal.pcbi.1006042.s009.xlsx")

saveRDS(list2,"data/ADGeneList.RDS")

#open.targets <- read.csv("data/AD_genes_Open_Targets.txt", sep = "\t")
open.targets <- read.csv("data/my_AD_list.csv")
open.targets <- open.targets[open.targets$association_score.overall >= 0.1,]
open.targets <- open.targets$target.gene_info.symbol
open.targets <- bitr(unlist(open.targets),fromType = "SYMBOL",
                   toType = "ENTREZID",OrgDb = "org.Hs.eg.db")
write.csv(open.targets,"outData/open_targets_genes.csv")

list3   <- list("ADCL" = str$ENTREZID, "ADGeneList" = list2$ENTREZID, 
                "open_targets" = open.targets$ENTREZID)
saveRDS(list3,"data/ADGeneSets.RDS")


# Listing genes in pathways

genesets.dir <- "data/KEGG_PCxN.RDS"
kegg.paths   <- readRDS(genesets.dir)

kegg.genes   <- sapply(kegg.paths, function(X){ 
                        gene.df <- bitr(X, fromType ="ENTREZID",toType = "SYMBOL",
                                                     OrgDb = org.Hs.eg.db)
                        return(paste(gene.df$SYMBOL,collapse =", "))
                      })

path.tab    <- data.frame("Pathway" = names(kegg.genes), "Genes" = kegg.genes )
rownames(path.tab ) <- NULL
write.table(path.tab,"data/PathGeneTab.tsv", sep = "\t")

