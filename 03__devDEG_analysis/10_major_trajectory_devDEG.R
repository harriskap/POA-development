suppressPackageStartupMessages({
  library("edgeR")
  library("stringr")
  library("limma")
  library("RColorBrewer")
  library("Glimma")
  library("dplyr")
  library("clusterProfiler")
  library("org.Mm.eg.db")
})

dt_file_list <- list.files( "./data/limma_voom", full.names=TRUE)
obs_files <- dt_file_list[grepl("obs-cts", dt_file_list)] 
blk_files <- dt_file_list[grepl("pseudo-bulk", dt_file_list)] 

length( obs_files)==length( blk_files)

for (itr in 1:length(obs_files)){
    # get files for each cell type
    # itr =1 
    obs_fl_itr <- obs_files[itr]
    blk_fl_itr <- blk_files[itr]
    counts <- read.delim( blk_fl_itr, row.names=1, sep=',')
    obs    <- read.delim( obs_fl_itr, row.names=1, sep=',')
    
    obs_file_name  <- basename(file.path(obs_fl_itr))
    # was dt_nm
    cell_type <- str_split(obs_file_name, "_obs")[[1]][1]
    dt_nm <- cell_type
    print(dt_nm)
    
    save_path <- "./data/limma_voom"
    dt_nm1 <- paste0("glimma_", dt_nm)
    res_dir <- file.path(save_path, dt_nm)
    if (dir.exists(res_dir)){
       next
    }
    
    # create edgeR DGE object
    x <- DGEList(counts)
    rownames(obs) <- make.names(rownames(obs))      
    group <- as.factor(obs$age)
    x$samples$group <- group
    # Rename all levels, by name
    # levels(group) <- list('Fetal'='Fetal','Neonatal'='Neonatal','Infancy'='Infancy','Childhood'='Childhood',
    #                    'Adolescence'='Adolescence','Adult'='Adult')
  
    # add all wanted features to x$samples
    #feats_oi <- c( 'chem', 'Sex', 'Library.Prep.Lot')
    #for (feat_itr in feats_oi){
    #  fact <- as.factor( obs[feat_itr])
    #  x$samples[feat_itr] <- obs[feat_itr]
    #}
    #x$samples[,feats_oi] <- mutate_all( x$samples[,feats_oi], as.character)     
    # control for lot 2 only, all others show no effects
    #prep_lot <- x$samples$Library.Prep.Lot
    #x$samples$Library.Prep.Lot[prep_lot!=2] <- 999
    keep.exprs <- filterByExpr(x, group=group)
    x <- x[keep.exprs,, keep.lib.sizes=FALSE]
    # split cellranger gene ENSEMBL and SYMBOLS 
    #comb_geneids <- rownames(x)
    #split_nms <- sapply( comb_geneids, function(x) {strsplit( x, "--", fixed=TRUE)[[1]]}[c(1,2)])
    # re-format and label name types
    #geneid <- t(split_nms)
    #geneid <- t(comb_geneids)
    #colnames(geneid) <- c("ENSEMBL","SYMBOL")
    #row.names(geneid) <- 1:nrow(geneid)
    #rownames(x) <- geneid[,2]
    #x$genes <- geneid[,1]
    x$genes <- rownames(x)
    # TMM normalize
    x <- calcNormFactors(x, method ="TMM") # "TMM"
    lcpm <- cpm(x, log=TRUE)
    # create design matrix
    # chem <- as.factor( x$samples$chem)
    # sex  <- as.factor( x$samples$Sex)
    # lot  <- as.factor( x$samples$Library.Prep.Lot)
    # nlvls = c( nlevels( chem), nlevels( sex), nlevels( lot))
    #factors <- c( "+ chem","+ sex","+ lot")
    #des_form = as.formula( paste("~0 + group", paste( factors[nlvls>1], collapse="")))
    des_form <- as.formula("~0 + group")
    design <- model.matrix( des_form)
  
    colnames(design) <- gsub("group", "", colnames(design))
    colnames(design) <- make.names(colnames( design))
    ### need to be sure to only make comparison of stages present in data
    # pull stage_ids in design matrix
    #contr_nms <- colnames( design)
    # pull intersect with stage_order
    #inter_nms <- intersect( make.names( obs$stage_id), contr_nms)
    # create contrast vector for comparisons across stages
    num_nms <- length(colnames(design))
    inter_nms <- colnames(design)
    # empty vector to hold contrast comparisons
    comp_vec <- c()
    # append comparisons
    for (itr in 1:num_nms){
        itr_nms2 <- inter_nms[-c(1:itr)]
        itr_nms1 <- rep(c(inter_nms[itr]), times=length(itr_nms2))
        comp_vec <- append( comp_vec, paste( itr_nms1, itr_nms2, sep='-'))
    }
    contr.matrix <- makeContrasts( contrasts=comp_vec, levels=colnames( design))
    v <- voom(x, design, plot=FALSE)
    vfit <- lmFit(v, design)
    vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
    efit <- eBayes(vfit)
    tfit <- treat( vfit, lfc=1)
    dt <- decideTests(efit)
  
    iii <- 1
    for (itr in comp_vec){
      html_dir <- paste0("data/limma_voom/", dt_nm)
      glMDPlot(efit, coef=iii, status=dt, main=colnames(efit)[iii],
             side.main="SYMBOL", counts=lcpm, groups=group, launch=FALSE, folder=html_dir, html=itr) # html="MD-Plot"
      iii <- iii + 1
    }
    dt_nm <- paste0("glimma_", dt_nm)
    dir.create(file.path(save_path, dt_nm), showWarnings = FALSE)
    logTMM_fil_nm <- paste(save_path, dt_nm, 'logTMM_cts.csv', sep="/")
    write.csv(lcpm, logTMM_fil_nm)
    results_fl_nm <- paste(save_path, dt_nm, 'results_file.txt', sep="/")
    write.fit(efit, dt, file = results_fl_nm)
    DEG_fl_nm <- paste(save_path, dt_nm, 'DEGlist.RDS', sep="/")
    saveRDS(x, DEG_fl_nm)
}