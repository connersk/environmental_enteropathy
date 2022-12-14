require(Seurat)
require(parallel)
require(tidyverse)

#'
#' returns output directory of figs/ srobjs/ endclusts/ 
#' srobj with base level clusters as column
GenTieredClusters <- function(srobj, cluster_assay = "SCT", cells = NULL, tier=0, clustN = 0,
                              PreProcess_fun = PreProcess_sctransform,
                              BaseCondition_fun = BaseCondition_default,
                              ChooseOptimalClustering_fun = ChooseOptimalClustering_default,
                              saveSROBJdir=NULL, figdir=NULL, SaveEndNamesDir=NULL, SaveEndFileName=NULL) {
  ######################################################################################################
  #' make sure output directories exist
  ######################################################################################################
  if (!is.null(saveSROBJdir)) dir.create(saveSROBJdir, showWarnings = F, recursive = T)
  if (!is.null(figdir)) dir.create(figdir, showWarnings = F, recursive = T)
  if (!is.null(SaveEndNamesDir) & tier==0) dir.create(SaveEndNamesDir, showWarnings = F, recursive = T)
  SaveEndFileName <- paste0(SaveEndFileName, sprintf("tier%s_clust%s", tier, clustN))
  
  ######################################################################################################
  #' make sure QC metadata exists
  ######################################################################################################
  if (is.null(srobj@meta.data$nFeature_RNA) & tier == 0) {srobj$nFeature_RNA <- Matrix::colSums(srobj@assays$RNA@counts > 0)}
  if (is.null(srobj@meta.data$nCount_RNA) & tier == 0) {srobj$nCount_RNA <- Matrix::colSums(srobj@assays$RNA@counts)}
  if (is.null(srobj@meta.data$percent.mt) & tier == 0) {srobj$percent.mt <- PercentageFeatureSet(srobj, pattern = "^MT-")}
  
  ######################################################################################################
  #' basic processing
  ###################################################################################################### 
  #' subset to provided cells
  if (is.null(cells)) {
    cells <- colnames(srobj)
  }
  working_srobj <- subset(srobj, cells=cells)
  working_srobj@misc$tier <- tier
  
  #' keep track of posistion in tree for logging
  message(paste0("Starting tier: ", tier, ", cluster: ", clustN, ", with ", ncol(working_srobj), " cells" ))
  
  #' preprocess
  working_srobj <- PreProcess_fun(working_srobj, figdir=figdir)
  
  #' get optimum cluster res and clusters
  res <- ChooseOptimalClustering_fun(working_srobj,
                                     assay=cluster_assay,
                                     PreProcess_fun = PreProcess_fun,
                                     downsample_num = 7500,
                                     figdir=figdir)
  working_srobj <- FindClusters(working_srobj,
                                assay = cluster_assay,
                                resolution = res,
                                k.param=ceiling(0.5*sqrt(ncol(working_srobj))))
  ######################################################################################################
  #' base conditions
  ######################################################################################################
  # 0 = not reached, 1 = too few cells, 2 = only one cluster, 3 = 2 indistinguishable clusters
  reached_base_condition <- BaseCondition_fun(working_srobj)
  
  ######################################################################################################
  #' logging & plotting
  ###################################################################################################### 
  cell.num.per.clust <- table(Idents(working_srobj))
  message("Cell counts per cluster")
  print(cell.num.per.clust)
  
  if (!is.null(figdir)) {
    
    #' save basic cell counts
    cat(sprintf("Srobj: with %s genes & %s cells", nrow(working_srobj), ncol(working_srobj)),
        file = sprintf("%s/basic_counts.txt", figdir), sep="\n")
    cat(sprintf("\n Choose %s PCs \n", working_srobj@misc$nPCs), file = sprintf("%s/basic_counts.txt", figdir), append = T)
    cat("\nTable of num cells per cluster:\n", file = sprintf("%s/basic_counts.txt", figdir), append = T)
    write.table(cell.num.per.clust, file = sprintf("%s/basic_counts.txt", figdir), row.names = F, append = T)  
    
    #' violin plot of nFeature, nCount, & percent mito
    old.idents <- Idents(working_srobj) 
    Idents(working_srobj) <- rep("all.data", ncol(working_srobj))
    
    VlnPlot(working_srobj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = -1)
    ggsave(sprintf("%s/QC_vln_plot.pdf", figdir))
    
    #' basic QC scatter plots
    FeatureScatter(working_srobj, feature1 = "nCount_RNA", feature2 = "percent.mt")
    ggsave(sprintf("%s/QC_scatter_plot_nCount_pctMito.pdf", figdir))
    FeatureScatter(working_srobj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
    ggsave(sprintf("%s/QC_scatter_plot_nCount_nFeat.pdf", figdir))
    
    #' set idents back to clusters
    Idents(working_srobj) <- old.idents
    
    #' basic umap of clusters (using umap-learn, as only it is allowed on premade graphs)
    working_srobj <- RunUMAP(working_srobj, graph=paste0(cluster_assay, "_snn"), umap.method="umap-learn")
    DimPlot(working_srobj, reduction = "umap", label = T)
    ggsave(sprintf("%s/cluster_umap.pdf", figdir))
    
    if (!is.null(working_srobj@meta.data$orig.ident)) {
      DimPlot(working_srobj, group.by = "orig.ident", reduction = "umap", label = F)
      ggsave(sprintf("%s/orig_ident_umap.pdf", figdir))
    }
    #   DimPlot(working_srobj, group.by = "smp_id", reduction = "umap", label = T)
    #   DimPlot(working_srobj, group.by = "inflammation", reduction = "umap", label = T)
    
    markers <- FindAllMarkers(working_srobj, assay = "RNA", only.pos = T, max.cells.per.ident = 2000)
    if ( !(nrow(markers) < 1) ) {
      top10markers <- markers %>% group_by(cluster) %>% filter(p_val_adj < 0.05) %>%  top_n(10, avg_logFC)
      write.table(top10markers, sprintf("%s/top10markers.tsv", figdir), sep="\t", row.names = F)
      write.table(markers, sprintf("%s/allmarkers.tsv", figdir), sep="\t", row.names = F)
      
      #' output heatmap of markers
      if (nrow(top10markers) < 5) {top10markers <- markers %>% group_by(cluster) %>% top_n(10, avg_logFC)}
      tryCatch({ 
        working_srobj <- ScaleData(working_srobj, features = top10markers$gene)
        DoHeatmap(working_srobj, features = top10markers$gene, raster = F)
        ggsave(sprintf("%s/top10markers_heatmap.pdf", figdir))
     }, error = function(e) {message("Failed to produce heatmap"); print(paste("HEATMAP ERRORED ON:  ", e))})
    }
  }
  if (!is.null(saveSROBJdir)) {
    if (reached_base_condition == 3) # merge final two clusters
    {
      Idents(working_srobj) <- 0
      working_srobj$two_not_quite_clusters <- working_srobj$seurat_clusters
      working_srobj$seurat_clusters <- 0
    }
    message(paste("saving srobj to", saveSROBJdir))
    saveRDS(working_srobj, sprintf("%s/subset_srobj.rds", saveSROBJdir))
  }


  ######################################################################################################
  #' setup and recurse
  ######################################################################################################
  
  #' how to stop 
  if (reached_base_condition > 0) {
    if (!is.null(SaveEndNamesDir)) {
      write.table(as.matrix(colnames(working_srobj)),
                  sprintf("%s/%s.tsv", SaveEndNamesDir, SaveEndFileName),
                  sep="\t", row.names = F, col.names = F)
    }
    return(working_srobj)
  }
  
  
  #' split all clusters into separate srobjs
  message("continuing to recurse")
  subsets <- SplitSrobjOnIdents(working_srobj, paste0("tier", tier))
  print(subsets)
  #' we haven't reached a base condition continue
  return(lapply(seq_along(subsets), 
                function(i) {
                  if (!is.null(saveSROBJdir)) {
                    saveSROBJdir_nextTier <- sprintf("%s/%s", saveSROBJdir, paste0("tier", tier+1, "_cluster", i-1))
                  }
                  if (!is.null(figdir)) {
                    figdir_nextTier <- sprintf("%s/%s", figdir, paste0("tier", tier+1, "_cluster", i-1))
                  }
                  GenTieredClusters(srobj, 
                                    cells=colnames(subsets[[i]]),
                                    tier=tier+1, clustN = i-1,
                                    saveSROBJdir = saveSROBJdir_nextTier,
                                    figdir = figdir_nextTier,
                                    SaveEndNamesDir = SaveEndNamesDir,
                                    SaveEndFileName = SaveEndFileName
                                    )
        }))
}


PreProcess <- function(srobj, ChoosePCs_fun = ChoosePCs_default, ChooseHVG_fun = ChooseHVG_default, figdir=NULL) {
  #' seurat v3 processing
  srobj <- NormalizeData(object = srobj)
  srobj <- ChooseHVG_fun(srobj)
  srobj <- ScaleData(object = srobj, features=VariableFeatures(srobj))
  srobj <- RunPCA(object = srobj, npcs = min(50, round(ncol(srobj)/2)), verbose=F)
  nPCs  <- ChoosePCs_fun(srobj, figdir=figdir)
  srobj@misc$nPCs <- nPCs
  message(paste0("chose ", length(nPCs), "PCs, added choices to srobj@misc$nPCs "))
  srobj <- FindNeighbors(object = srobj, dims=nPCs)
  return(srobj)
}

PreProcess_sctransform <- function(srobj, ChoosePCs_fun = ChoosePCs_default, figdir=NULL) {
  #' seurat v3 processing
  srobj <- SCTransform(object = srobj)
  srobj <- NormalizeData(srobj, assay = "RNA")
  srobj <- RunPCA(object = srobj, npcs = min(50, round(ncol(srobj)/2)), verbose=F)
  nPCs  <- ChoosePCs_fun(srobj, figdir=figdir)
  srobj@misc$nPCs <- nPCs
  message(paste0("chose ", length(nPCs), "PCs, added choices to srobj@misc$nPCs "))
  srobj <- FindNeighbors(object = srobj, dims=nPCs)
  return(srobj)
}

ChoosePCs_default <- function(srobj, improved_diff_quantile=.85, significance=0.01, figdir=NULL) {
  #' if num cells big use elbow plot if small use jackstraw, 
  #'   if none sig use first 2 pcs
  if (ncol(srobj) > 500) {
    eigValues <- (srobj@reductions$pca@stdev)^2  ## EigenValues
    varExplained <- eigValues / sum(eigValues)
    nPCs <- 1:max(which(diff(varExplained) < quantile(diff(varExplained), 1-improved_diff_quantile)) + 1)
    if (!is.null(figdir)) {
      ElbowPlot(srobj, ndims=ncol(srobj@reductions$pca@cell.embeddings)) + 
        geom_vline(xintercept  = length(nPCs), color="red")
      ggsave(sprintf("%s/ChoosenPCs.pdf", figdir))
    }
  } 
  else {
    suppressWarnings({srobj <- JackStraw(srobj, dims=50)})
    srobj <- ScoreJackStraw(srobj, dims=1:ncol(srobj@reductions$pca@cell.embeddings))
    nPCs <- which(JS(srobj$pca)@overall.p.values[,"Score"] < significance)
    if (!is.null(figdir)) {
      JackStrawPlot(srobj, dims = 1:ncol(srobj@reductions$pca@cell.embeddings)) + NoLegend()
      ggsave(sprintf("%s/ChoosenPCs.pdf", figdir))
    }
  }
  if (length(nPCs) < 2) {
    nPCs <- 1:2
  }
  return(nPCs)
}

ChoosePCs_JackStraw <- function(srobj, threshold=0.01, figdir=NULL) {
  suppressWarnings({srobj <- JackStraw(srobj, dims=50)})
  srobj <- ScoreJackStraw(srobj, dims=1:ncol(srobj@reductions$pca@cell.embeddings))
  nPCs <- which(JS(srobj$pca)@overall.p.values[,"Score"] < threshold)
  if (!is.null(figdir)) {
    JackStrawPlot(srobj, dims = 1:ncol(srobj@reductions$pca@cell.embeddings)) + NoLegend()
    ggsave(sprintf("%s/ChoosenPCs.pdf", figdir))
  }
  if (length(nPCs) < 2) {
    nPCs <- 1:2
  }
  return(nPCs)
}

ChooseHVG_default <- function(srobj) {
  FindVariableFeatures(srobj, nfeatures=2000, selection.method = "disp")
}

ChooseOptimalClustering_default <- function(srobj, downsample_num = Inf,
                                            PreProcess_fun = PreProcess, ...) {
  #' get resolution from paramsweep of heavily downsampled dataset
  
  # downsample to 10% of data
  # downsample input srobj?
  if ((downsample_num < Inf) & (downsample_num <  ncol(srobj))) {
    cells <- sample(colnames(srobj), downsample_num)
    srobj <- subset(srobj, cells=cells)
    srobj <- PreProcess_fun(srobj)
  }
  srobj <- ChooseClusterResolutionDownsample(srobj, ...)
  return(srobj@misc$resolution.choice)

}

EnoughDiffExp <- function(srobj, ident1, ident2, nUpregulated = 10, nDownregulated = 25) {
  
  markers = FindMarkers(srobj, ident.1 = ident1, ident.2 = ident2, assay = "RNA", logfc.threshold = 0.25, min.pct = 0.2, test.use = "wilcox")
  sig.markers = markers %>% filter(p_val_adj < 0.01)
  sig.markers = sig.markers %>% filter(avg_logFC > log(2))  
  return((nrow(sig.markers) > nUpregulated))
  
}

BaseCondition_default <- function(srobj, min_cluster_size=100, max_tiers=10) {
  
  #' if min cells or max tiers: stop
  if ( (ncol(srobj) < min_cluster_size) | (srobj@misc$tier > max_tiers) ) {
    message("found end state with min number of cells")
    return(1)
  }
  
  #' if only one cluster return
  if ( !(length(unique(Idents(srobj))) > 1) ) {
    message("found end state with one cluster")
    return(2)
  }
  
  #' if 2 clusts & different with diffexp ? continue : else stop
  else if ( (length(unique(Idents(srobj))) == 2) & 
       !(EnoughDiffExp(srobj, levels(Idents(srobj))[1], levels(Idents(srobj))[2])) ) {
    message("found end state with two indistinguishable clusters")
    return(3)
  }

  return(0)
}

SplitSrobjOnIdents <- function(srobj, proj_nm_suffix) {
  numbers_only <- function(x) !grepl("\\D", x)
  if (any(!numbers_only(levels(srobj)))) {
    idts <- sort(levels(Idents(srobj)))
  } else {
    idts <- sort(as.numeric(levels(Idents(srobj))))
  }
  lapply(idts,
         function(idt) {tmp <- subset(srobj, idents=idt); 
           tmp@project.name=paste0(tmp@project.name, idt, proj_nm_suffix, sep="_");
           return(tmp)})
}

ChooseClusterResolutionDownsample <- function(
  input.srobj, assay="RNA", n.pcs = input.srobj@misc$nPCs, sample.name =  format(Sys.time(), "%a-%b-%d-%X-%Y-%Z"),
  res.low = 0.2, res.high=0.8, res.n = 40, bias = "under", figdir=NULL) {#I modified the resolutions here - ck feb21
  
  ######## step 1: save the input seurat object as a new temporary object, 
  ########         dont want to overwrite or change the original one with all of the parameter scans
  
  srobj.tmp = input.srobj 
  # in case there have been other things calculated in the metadata, just cut down to simplify/avoid errors
  srobj.tmp@meta.data = srobj.tmp@meta.data[,c("nCount_RNA","nFeature_RNA")] # should just be the nUMI and nGene	
  
  ######## step 2: calculate the FindClusters over a large range of resolutions
  print("Performing parameter scan over multiple resolutions...")
  
  set.res = round(exp(seq(log(res.low), log(res.high), length.out=res.n)), digits=3)
  srobj.tmp = FindClusters(srobj.tmp, assay=assay, dims.use = n.pcs, k.param=ceiling(0.5*sqrt(ncol(srobj.tmp))),
                           resolution=set.res[1], save.SNN=T, plot.SNN=F,
                           force.recalc=TRUE, verbose=FALSE)
  
  res.that.work <- rep(T, length(set.res))
  for(i in 2:length(set.res)){
    tryCatch({
    srobj.tmp = FindClusters(
      srobj.tmp, assay=assay, resolution=set.res[i], reuse.SNN=TRUE, plot.SNN=FALSE, verbose=FALSE, force.recalc=FALSE)
    }, error = function(e) {res.that.work[i] <<- F; message(paste("errored on", set.res[i], "Flagging as doesn't work"))})
    print(paste("          ", round(100*i/length(set.res)), "% done with parameter scan", sep=""))
  }
  set.res <- set.res[res.that.work]

  
  ######## step 3: output plot of how the resolution changes the number of clusters you get
  n.clusters = vector(mode="numeric", length=length(set.res))
  names(n.clusters) = set.res
  for(i in 1:length(n.clusters)){
    n.clusters[i] = length(table(as.vector(srobj.tmp@meta.data[,paste0(assay, "_snn_res.", names(n.clusters)[i])])))
  }
  
  ######## step 4: calculate the silhouette width for each resolution
  print("Computing a silhouette width for each cell, for each resolution...")
  require(cluster)
  
  dist.temp = cor(t(srobj.tmp@reductions$pca@cell.embeddings[,n.pcs]), method="pearson")
  # Downsample to 5k
  random.cells.choose = if ( nrow(dist.temp) > 5000 ) sample(1:nrow(dist.temp), 5000) else 1:nrow(dist.temp)
  dist.temp.downsample = dist.temp[random.cells.choose, random.cells.choose]
  sil.all.matrix = matrix(data=NA, nrow=nrow(dist.temp.downsample), ncol=0)
  
  for(i in 1:length(set.res)){
    clusters.temp = as.numeric(as.vector(
      srobj.tmp@meta.data[random.cells.choose,paste0(assay, "_snn_res.", set.res[i])]))
    if(length(table(clusters.temp))>1){
      sil.out = silhouette(clusters.temp, as.dist(1-as.matrix(dist.temp.downsample)))
      sil.all.matrix = cbind(sil.all.matrix, sil.out[,3])
    }
    if(length(table(clusters.temp))==1){
      sil.all.matrix = cbind(sil.all.matrix, rep(0, length(clusters.temp)))
    }
    print(paste("          ", round(100*i/length(set.res)), "% done with silhouette calculation", sep=""))
    
  }
  
  
  ######## step 5: calculate summary metric to compare the silhouette distributions,
  ########         average has worked well so far... could get fancier
  
  print("Identifying a best resolution to maximize silhouette width")
  sil.average = setNames(colMeans(sil.all.matrix), set.res)
  sil.medians <- setNames(apply(sil.all.matrix, 2, median), set.res)
  
  ######## step 6: automate choosing resolution that maximizes the silhouette 
  hist.out = hist(sil.average, length(sil.average)/1.2,  plot=FALSE)
  
  #  take the ones that fall into the top bin, 
  #  and the max OR MIN of those  ******* can change this to under vs over cluster
  if(bias=="over"){
    resolution.choice = as.numeric(max(
      names(sil.average[which(sil.average>hist.out$breaks[length(hist.out$breaks)-1])])))
  }
  if(bias=="under"){
    resolution.choice = as.numeric(min(
      names(sil.average[which(sil.average>hist.out$breaks[length(hist.out$breaks)-1])])))
  }
  
  # get the silhouette of the best resolution: 
  silhouette.best = as.numeric(sil.average[paste(resolution.choice)])
  
  print(paste("Best Resolution Choice: ", resolution.choice, ", with average silhouette score of: ",
              round(silhouette.best, digits=3), ", giving ", as.numeric(n.clusters[paste(resolution.choice)]),
              " clusters", sep=""))
  
  
  ######### step 7: output plot and data 
  
  if (!is.null(figdir)) {
    
    print(paste0("Ouptutting summary statistics and returning seurat object... ",
                 "This will create a pdf in your output directory,",
                 " and will return your input seurat object ammended with the best choice",
                 " for clusters (found as Best.Clusters in the meta.data matrix, and set to your new ident)..."))
    
    pdf(paste(figdir, "/", sample.name, ".pdf", sep=""),
        width=10, height=4, useDingbats=FALSE)
    par(mfrow=c(1,3))
    # Resolution vs # of Clusters
    plot(set.res, n.clusters, col="black", pch=19,
         type="p", xlab="Resolution", ylab="# Clusters",
         main="Resolution vs. # Clusters")
    
    # Resolution vs Average Silhouette
    plot(set.res, sil.average, col="black", pch=19,
         type="p", xlab="Resolution", ylab="Average Silhouette",
         main="Resolution vs. Average Silhouette")
    points(set.res, sil.medians, col="red", pch=15)
    abline(h=hist.out$breaks[length(hist.out$breaks)-1], col="firebrick3", lty=2)
    abline(v=resolution.choice, col="dodgerblue2", lty=2)
    
    # N Clusters vs Average Silhouette
    plot(n.clusters, sil.average, col="black", pch=19,
         type="p", xlab="# Clusters", ylab="Average Silhouette",
         main="# Clusters vs. Average Silhouette")
    points(n.clusters, sil.medians, col="red", pch=15)
    abline(h=hist.out$breaks[length(hist.out$breaks)-1], col="firebrick3", lty=2)
    abline(v=as.numeric(n.clusters[paste(resolution.choice)]), col="dodgerblue2", lty=2)
    dev.off()
  }
  
  ######## step 8: return the original seurat object, with the metadata containing a 
  ########         concatenated vector with the clusters defined by the best choice here,
  ########         as well as the ident set to this new vector
  
  Best.Clusters = srobj.tmp@meta.data[,paste0(assay, "_snn_res.", resolution.choice)]
  
  input.srobj$Best.Clusters = Best.Clusters
  Idents(input.srobj) = input.srobj$Best.Clusters
  input.srobj@misc <- list("resolution.choice" = resolution.choice)
  return(input.srobj)
}



#######################################################################################################
#' Adopted from
#' https://raw.githubusercontent.com/broadinstitute/BipolarCell2016/master/class.R
#' 
#' Ported to work with Seurat v3 by Benjamin Doran Oct 2019
#' Added columns for FDR and posFrac to output by default
#' 
#######################################################################################################


markers.binom <- function(object, clust.1,clust.2=NULL,effect.size=log(2), TPM.mat=NULL, Count.mat=NULL) {
  genes.use=rownames(object@assays$RNA@data)
  clust.use=Idents(object)
  cells.1=names(clust.use[which(clust.use%in%clust.1)])
  
  if (is.null(clust.2)) {
    clust.2="rest"
    cells.2=names(clust.use)
    cells.2=cells.2[!(cells.2%in%cells.1)]
  } else {
    cells.2=names(clust.use[which(clust.use%in%clust.2)])
  }
  
  Count.mat = object@assays$RNA@counts
  if (is.null(TPM.mat)) TPM.mat = exp(object@assays$RNA@data[, c(cells.1, cells.2)])-1
  if (is.null(Count.mat)) Count.mat = object@assays$RNA@counts[genes.use, c(cells.1, cells.2)]
  result=binomcount.test(object, cells.1, cells.2, effect.size, TPM.mat, Count.mat)
  
  if (nrow(result) < 1) {
    return(setNames(data.frame(matrix(ncol = 7, nrow = 0)),
                    c("log.effect", "pval", "fdr", "posFrac.1", "posFrac.2", "nTrans_1", "nTrans_2")))
  }
  
  result[,"posFrac.1"] <- posFrac.1 <- apply(object@assays$RNA@data[rownames(result),cells.1, drop=F],1,function(x) round(sum(x > 0)/length(x),2))
  result[,"posFrac.2"] <- posFrac.2 <- apply(object@assays$RNA@data[rownames(result),cells.2, drop=F],1,function(x) round(sum(x > 0)/length(x),2))
  
  if (clust.2=="rest"){
    genes.include = posFrac.1 >= 0.1
  } else{
    genes.include = (posFrac.1 >= 0.1) | (posFrac.2 >= 0.1)
  }
  
  result = result[genes.include,]
  result = result[order(abs(result$log.effect), decreasing=TRUE),]
  
  if (nrow(result) < 1) {
    return(setNames(data.frame(matrix(ncol = 7, nrow = 0)),
                    c("log.effect", "pval", "fdr", "posFrac.1", "posFrac.2", "nTrans.1", "nTrans.2")))
  }
  
  #Mean number of transcripts per cell
  if (!is.null(attr(object@assays$RNA,"counts"))){
    nTrans.1 = apply(object@assays$RNA@counts[rownames(result), cells.1, drop=F], 1, function(x) round(mean(x),3))
    nTrans.2 = apply(object@assays$RNA@counts[rownames(result), cells.2, drop=F], 1, function(x) round(mean(x),3))
    result[,"nTrans.1"] = nTrans.1
    result[,"nTrans.2"] = nTrans.2
  }
  
  return(result)
} 


binomcount.test <- function(object, cells.1, cells.2, effect.size, TPM.mat, Count.mat) {
  
  x=TPM.mat
  y=Count.mat
  
  #Test for enrichments in cluster #1
  m = apply(x[, cells.2], 1, function(x) sum(x>0)) #Number of cells expressing marker in cluster #2
  m1 = m; m1[m==0]=1; # Regularization. Add a pseudo count of 1 to unexpressed markers to avoid false positives
  n = apply(x[, cells.1], 1, function(x) sum(x>0)) #Number of cells expressing marker in cluster #1
  #Find the probability of finding n or more cells +ve for marker in cluster 1 given the fraction in cluster 2
  pv1 = pbinom(n, length(cells.1), m1/length(cells.2), lower.tail = FALSE) + dbinom(n, length(cells.1), m1/length(cells.2))
  
  log_fold_express = log(n*length(cells.2)/(m*length(cells.1))) #log proportion of expressing cells
  d1 <- data.frame(log.effect=log_fold_express,pval=pv1)
  d1 <- subset(d1, log.effect >= effect.size)
  d1 <- d1[order(d1$pval,decreasing=FALSE),]
  
  #Enrichments in cells.2
  n1 = n; n1[n==0]=1; # Regularization.  Add a pseudo count of 1 to unexpressed markers to avoid false positives
  #Find the probability of finding m or more cells +ve for marker in cluster 2 given the fraction in cluster 1
  pv2 = pbinom(m, length(cells.2), n1/length(cells.1), lower.tail=FALSE) + dbinom(m, length(cells.2), n1/length(cells.1))
  d2 <- data.frame(log.effect=log_fold_express,pval=pv2)
  d2 <- subset(d2, log.effect <= -effect.size)
  d2 <- d2[order(d2$pval,decreasing=FALSE),]
  
  d = rbind(d1, d2);
  d = d[order(d$pval, decreasing=FALSE),]
  d$fdr <- p.adjust(d$pval, method="fdr")
  return(d)
} 

