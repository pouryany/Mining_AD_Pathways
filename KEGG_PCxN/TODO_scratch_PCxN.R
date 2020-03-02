
# Problems to be solved. 

# Throws error with only one set selected



####### estimate00.R

fname = paste0("../data/", geneset_file, sep = "")


######### PCxN_estimate01.R



        rels_char <- "7,11"
        geneset_file <- "sok_test_L1000.RDS"
        pcor_choice <- "1"
        output_folder <- "output_pour"
        id <- "1"


        for(j in gse){
            # subset expression ranks
            seriesn <- tissue_seriesn[j]
            series <- tissue_series[j]
            
            ind_series <- colnames(tissue_rnk) %in% tissue_meta$sample[tissue_meta$series == seriesn]
            exprs_rnk <- tissue_rnk[ ,ind_series ]
            cat(
                "+==========================+\n",
                "Tissue:",tissue_select,"\n",
                "Series:",series,"\n",
                "+==========================+\n\n\n"
            )
            
            # ---- Get pathway summaries for disjoint gene sets ----
            summary_dis_list <-list()
            
            for (gs in 1:length(gs_lst)) {
                nm <- names(gs_lst[gs])
                tmp <- list(GetSummary(dat=exprs_rnk,gs=gs_lst[[gs]],colMeans))
                summary_dis_list[[nm]] <- tmp
            }
            
            # get experiment-level estimates (parallel loop for pathway pairs)
            pb = txtProgressBar(min=0,max=number_of_pathways,style=3,initial=0)
            cat("\n")
            res = mclapply(input,ProcessElement,mc.cores=as.numeric(4))
            # Remove NULLs
            res[sapply(res, is.null)] <- NULL
            close(pb)
            
            # save experiment-level estimates
            fname = paste0(make.names(tissue_select),"_",series)
            saveRDS(res, paste0("../",output_folder,"/",fname,"_cpad_pathcor.RDS"))

            accepted_pairs <- length(res)
            actual_pairs <- data.frame("actual_number_of_pairs" = accepted_pairs)
            saveRDS(actual_pairs, paste0("../",output_folder,"/pairs.RDS"))

        }





        # Pre-existing
        r_geneset_file <- geneset_file
        r_cores <- 4
        r_output_folder <- output_folder
        r_pcor <- pcor_choice
        r_rels <- rels
        r_pathways <- number_of_pathways


#######





# PCxN_estimate02



            # ==== INPUTS ====
        geneset_file <- "sok_test_L1000.RDS"

        output_folder <- "output_pour"
        id <- 1




        # ==== INPUTS ====
geneset_file <- "sok_test_L1000.RDS"
output_folder <- "output_pour"







#ion
#In addition: Warning message:
#In gzfile(file, "rb") :
#  cannot open compressed file '../output_pour/mean_pcor2_barcode/pbmc_GSE7753_cpad_pathcor.RDS', probable reason 'No such file or directory'
# 


# directory with gene expression background
barcode_dir = "../data/HGU133plus2/"

# ==== Functions ====
# adjust p-values and correlation estimates,
# otherwise the functions to combine the p-values
# cannot handle values near 0 and values near 1
AdjustPmat = function(p_mat,eps=1E-16){
  res = t(apply(p_mat,1,function(pval){
    pval[pval <= eps] = eps
    pval[pval >= 1-eps] = 1 - eps
    return(pval)
  }))
  return(res)
}

# ==== Barcode Annotation ====
# Sample annotation for the gene expression background
tissue_annot <- readRDS( "../data/Barcode3.tissue.RDS" )

# ==== GSE Series ====
# GSE series per tissue
tmp <- subset(tissue_annot,select=c(tissue,series))
gse_lst <- lapply(split(tmp, tmp$tissue),function(x){table(x$series)})
# order by number of samples
gse_lst <- gse_lst[order(sapply(gse_lst,sum),decreasing=T)]
rm(tmp)
# filter series with at least 5 samples
res <- lapply(gse_lst,function(x){x[which(x >= 10)]})
res <- res[lapply(res,length)>0]
# order by number of samples
res <- res[order(sapply(res,sum),decreasing=T)]


fname = c()
for(tissue_select in names(res)){
    # tissue meta-data
    tissue_meta <- subset(tissue_annot[tissue_annot$tissue == tissue_select,],select=c(sample,series))
    # tissue GSE series
    tissue_series <- names(res[[ tissue_select ]])
    tissue_series <- gsub(";","_",tissue_series)
    # RDS file name
    fname = c(fname, paste0(make.names(tissue_select),"_",tissue_series))
}

# ==== Pathway Annotation ====
# Filtered gene set annotation
gs_lst = readRDS(paste("../data/",geneset_file, sep=""))

# indices for pathway pairs
npairs <- readRDS(paste0("../",output_folder,"/mean_pcor2_barcode/res/pairs.RDS"))
number_of_pathways = npairs[[1]]
# split pathway pairs in chunks
pairs_chunks <- split(1:number_of_pathways, ceiling(1:number_of_pathways/100000))

# ==== Pathway Names =====
# load a set of experiment-level estimates
my_rds = paste0("../",output_folder,"/mean_pcor2_barcode/",fname[13],"_cpad_pathcor.RDS")
tML = readRDS(my_rds)

info <- data.frame("parts" = length(pairs_chunks))

saveRDS(info, paste0("../",output_folder,"/mean_pcor2_barcode/res/parts.RDS"))

# Start looping through parts here
#for (cp in 1:length(pairs_chunks)) {
for (cp in id:id){
    print("In loop")
    print(cp)

    myLst = tML[ pairs_chunks[[ cp ]] ]
    print(length(myLst))

    # get names for pathway pairs
    t0 <- sapply(myLst,function(x){ c(x[["Pathway.A"]],x[["Pathway.B"]]) })

    t1 <- t(t0)
    print(dim( t1 ))
    #print(t1[1,])
    #print(t1[,1])
    res = as.data.frame( t1 )
    colnames(res) = c("Pathway.A","Pathway.B")

    # ==== Overlap Coefficient =====
    # get overlap coefficients for pathway pairs
    res$Overlap.Coefficient = unlist( sapply(myLst,function(x){x[[ "Overlap.Coeff" ]]}), use.names = F )

    # ===== Correlation Estimates ====
    GetCorEstimate = function(ic){
        # load experiment-level estimates
        my_rds = paste0("../",output_folder,"/mean_pcor2_barcode/",fname[ic],"_cpad_pathcor.RDS")
        myLst = readRDS(my_rds)[ pairs_chunks[[ cp ]] ]
        # extract correlation estimates
        tmp = unlist( sapply(myLst, function(x){x[[ "estimate" ]]}), use.names = F )
        setTxtProgressBar(pb,ic)
        return( tmp )
    }

    pb = txtProgressBar(min=0,max=length(fname),style=3,initial=0,char="-")
    cat("\n")
    cor_estimates = mclapply( seq_along(fname), GetCorEstimate, mc.cores = as.numeric("2") )
    close(pb)

    cor_estimates = Reduce(f = cbind, x = cor_estimates)
    colnames(cor_estimates) = fname

    # number of samples per experiment
    n_vec = readRDS( paste0( "../",output_folder,"/mean_pcor2_barcode/res/n_vec.RDS" ) )

    # weighted average for the correlation estimates
    n_mult = n_vec/sum(n_vec)
    res$PathCor = c( cor_estimates%*%n_mult )

    # ==== P-Values ====
    GetPvals = function(ic){
        # load experiment-level estimates
        my_rds = paste0("../",output_folder,"/mean_pcor2_barcode/",fname[ic],"_cpad_pathcor.RDS")
        myLst = readRDS(my_rds)[ pairs_chunks[[ cp ]] ]
        # extract p-value
        tmp = unlist( sapply(myLst, function(x){x[[ "p.value" ]]}), use.names = F )
        setTxtProgressBar(pb,ic)
        return( tmp )
    }

    pb = txtProgressBar(min=0,max=length(fname),style=3,initial=0,char="-")
    cat("\n")
    pvals = mclapply( seq_along(fname), GetPvals, mc.cores = as.numeric("2") )
    close(pb)

    pvals = Reduce(f = cbind, x = pvals)
    colnames(pvals) = fname

    # adjust p-values, otherwise the functions to combine the p-values
    # cannot handle values near 0 and values near 1
    pvals = AdjustPmat(pvals)

    # Liptak's Method to combine p-values
    CombinePval = function(ic){
        if( any( is.na( pvals[ic,] ) ) ){
            tmp = NA
        }else{
            tmp = c(sumz(p=pvals[ic,], weights=n_vec)$p)
        }
        setTxtProgressBar(pb,ic)
        return(tmp)
    }

    pb = txtProgressBar(min=0,max=length(fname),style=3,initial=0,char="-")
    cat("\n")
    combined_pvals = mclapply( 1:nrow(pvals), CombinePval, mc.cores = as.numeric("2") )
    close(pb)

    res$p.value = unlist( combined_pvals )

    # save results in data frame
    saveRDS(res, paste0("../",output_folder,"/mean_pcor2_barcode/res/pcxn_mean_pcor2_barcode_part",cp,".RDS"))
}

# Finish looping through parts

r_parts <- length(pairs_chunks)
r_output_folder <- output_folder
report <- matrix(nrow = 0, ncol = 2)
report <- rbind(report,c("Parts                                 ",r_parts))

temp2 <- paste("../",r_output_folder,"/report.txt",sep = "")
write.table(report, file = temp2, sep = " : ", dec = ".", append = TRUE,
            col.names = FALSE,row.names = FALSE,quote = FALSE)

rm(list = ls())
gc()

