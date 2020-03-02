library(GSA)

temp <- GSA::GSA.read.gmt(filename = "GenesetUpdate/c2.cp.v7.0.entrez.gmt")


temp2 <- temp$genesets 
names(temp2) <- paste0("Pathway.",temp$geneset.names)
names(temp2)
temp$geneset.descriptions
meta <- cbind(names(temp2),temp$geneset.descriptions)
saveRDS(temp2,"MSigDBV7_Canonical.RDS")
saveRDS(meta,"PathwaysMeta.RDS")
