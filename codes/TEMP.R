rm(list = ls())
print("Mapping to PCxN")
options(StringAsFactors = F)

library(RColorBrewer)
library(igraph)
library("GGally")
library(intergraph)
library(intergraph)
library(dplyr)

output.dir  <- paste0("Output")
enrich.fdr  <- 0.05
edge.fdr    <- 0.05
cor.thresh  <- 0

# Loading PCxN
    pcxn.dir   <- paste0("data/improved_PCxN_KEGG_PCxN.RDS")
    PCxN       <- readRDS(pcxn.dir)

    
    enrich.map.dir   <- paste0("data/KEGG_KEGG_ENRICHMENT.RDS")
    enrich.map       <- readRDS(enrich.map.dir)
    
    part.A <- dplyr::inner_join(enrich.map, PCxN, 
                                by = c("x" = "Pathway.A", "y" = "Pathway.B"))
    part.B <- dplyr::inner_join(enrich.map,
                                PCxN, by = c("x" = "Pathway.B", "y" = "Pathway.A"))
    
    part.AB <- rbind(part.A,part.B)    
    
    part.AB <- part.AB %>% dplyr::select(., x, y, adj_pval_listWise,PathCor,p.Adjust)
    
    part.AB <-  part.AB %>% mutate(.,enriched = ifelse(adj_pval_listWise < enrich.fdr,T,F),
                         correlated = ifelse(p.Adjust < edge.fdr, T, F),
                         enrich_weight = -log2(adj_pval_listWise))
    
    
    
    
    # Phenylalanine, tyrosine and tryptophan biosynthesis = 5 
    # The above pathway will be excluded from the top list
    
    path.dir   <- read.csv("new_data/Annotated_with_KEGG_Class.txt", sep = "\t") 
    path.cls   <- path.dir$Pclass
    shrt.name  <- as.character(path.dir$Short_name)
    path.rank  <- path.dir$Overall_rank
    path.smpl  <- path.dir$Class_Simple
    
    names(shrt.name) <- path.dir$Pathway
    names(path.rank) <- path.dir$Pathway
    names(path.smpl) <- path.dir$Pathway
    
      
    path.dir2  <- path.dir$Pathway
    San.paths  <- path.dir2
    San.paths  <- paste0("Pathway.KEGG_",San.paths)
    all.paths  <- union(unique(part.AB$x),unique(part.AB$y))
    San.paths  <- San.paths[which(((San.paths) %in% all.paths))]
    San.paths  <- San.paths[1:30]
    San.paths  <- gsub("Pathway.KEGG_","", San.paths)
    
    
    
    san.names        <- path.rank[San.paths]
    #names(san.names) <- San.paths
    
    
    
    path.cls   <- path.smpl[San.paths]
    path.cls   <- as.factor(as.character(path.cls))
    names(path.cls) <- San.paths
    
    
    
    net <- part.AB
    net <- igraph::graph_from_data_frame(net, directed = F)
    
    
    V(net)$name <- gsub("Pathway.KEGG_","", V(net)$name)
    
    sub.mods <- induced_subgraph(net,(San.paths))
    
    
    # remove links that are neither enriched nor correlated. 
    sub.mods <- sub.mods - E(sub.mods)[E(sub.mods)$enriched == F & E(sub.mods)$correlated ==F]
    
    sub.rich <- sub.mods - E(sub.mods)[E(sub.mods)$enriched == F]
    E(sub.rich)$width <- (abs(E(sub.rich)$enrich_weight) /70) + 1
    
    
    sub.cor  <- sub.mods - E(sub.mods)[E(sub.mods)$correlated == F]
    E(sub.cor)$width <- (abs(E(sub.cor)$PathCor) *10) + 0.5
    
    
    set.seed(1)
    layers   <- layout_with_graphopt(sub.rich, charge = 0.05)
    # Setting cholestrol metabolism in a proper position
    
    layers[11,1] <- layers[11,1] - 1200
    layers[2,1] <- layers[2,1] + 100
    layers[7,2] <- layers[7,2] + 50
    
    library(RColorBrewer)
    getPalette = colorRampPalette(brewer.pal(9, "Set1"))
    colourCount = length(unique(path.cls))
    getPalette(colourCount)[path.cls]
    
    pdf("Images/myfig.pdf", height = 24, width = 10)
    par(mfrow=c(3,1))
    
    plot(sub.rich,
         vertex.color = getPalette(colourCount)[path.cls[V(sub.rich)$name]],
         layout = layers,
         vertex.size = 12,
         vertex.label= san.names[V(sub.rich)$name], 
         vertex.label.font=1,
         vertex.label.family = "sans",
         main = "Gene overlap between KEGG pathways in \n top 10% AD-pathways via text-mining")
    plot(sub.cor,
         vertex.color = getPalette(colourCount)[path.cls[V(sub.rich)$name]], 
         layout = layers,
         vertex.label= shrt.name[V(sub.rich)$name], 
         vertex.size = 12,
         vertex.label.font=1,
         vertex.label.cex = 1.4 ,
         vertex.label.family = "sans",
         main = "Canonical correlation (based on gene expression) between KEGG pathways in \n top 10% AD-pathways via text-mining")
    
    plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
    legend(x=0, y=0.6, as.character(unique(path.cls)), pch=21,
           col="#777777", pt.bg=getPalette(colourCount)[unique(path.cls)],
           pt.cex=4, cex=1.9, bty="n", ncol=1)
    
    dev.off()
    

    

