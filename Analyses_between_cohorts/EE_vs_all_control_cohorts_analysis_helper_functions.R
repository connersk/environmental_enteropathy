# Helper functions for comparing EE to all three control cohorts, broken down by analysis type


####################################################################################################
#################################      Differential expression     #################################
####################################################################################################


# For a given cell type, these functions

  # read in DE results from the following analyses
    # EE vs. all three control cohorts
    # EE vs. EoE cohort
    # EE. vs. Resection cohort
    # EE vs. Durban cohort

  # Returns genes that are 
    # significantly upregulated in EE vs. all three cohorts
    # display the same direction of regulation when comparing EE vs. each cohort individually

combine_markers = function(markers_all,markers_eoe,markers_wp,markers_durban,sig_threshold=0.05){
  bool_for_sig_vs_outgroups = markers_all$p_val_adj < sig_threshold
  marker_genes = markers_all$gene
  bools_for_direction = as.vector(sapply(marker_genes,function(x){check_marker(x,markers_eoe,markers_wp,markers_durban,1.0)}))
  final_genes = marker_genes[bool_for_sig_vs_outgroups & bools_for_direction]
  return(markers_all %>% filter(gene %in% final_genes))
}


check_marker = function(marker,markers_eoe,markers_wp,markers_durban,sig_threshold=1.0){
  in_eoe = marker %in% markers_eoe$gene[markers_eoe$p_val_adj <= sig_threshold]
  in_wp = marker %in% markers_wp$gene[markers_wp$p_val_adj <= sig_threshold]
  in_durban = marker %in% markers_durban$gene[markers_durban$p_val_adj <= sig_threshold]
  if(in_eoe & in_wp & in_durban){
    if(sign(markers_eoe$avg_logFC[markers_eoe$gene==marker])==sign(markers_wp$avg_logFC[markers_wp$gene==marker])&sign(markers_eoe$avg_logFC[markers_eoe$gene==marker])==sign(markers_durban$avg_logFC[markers_durban$gene==marker])){
      return(TRUE)
    }else{
      return(FALSE)
    }
  }else{
    return(FALSE)
  }
}

# Reading in the marker genes from the different sources
read_markers = function(path,all_three_prefix,eoe_prefix,wp_prefix,durban_prefix){
  markers.vs.all = read.csv(paste0(all_three_prefix,path,"all.csv"),row.names=1)
  markers.vs.all$cluster = "EE"
  markers.vs.all$gene = rownames(markers.vs.all)
  
  markers.vs.eoe = read.csv(paste0(eoe_prefix,path,"EoE.csv"),row.names = 1)
  markers.vs.eoe$cluster = "EE"
  markers.vs.eoe$gene = rownames(markers.vs.eoe)
 
  
  markers.vs.wp = read.csv(paste0(wp_prefix,path,"Resection.csv"),row.names = 1)
  markers.vs.wp$cluster = "EE"
  markers.vs.wp$gene = rownames(markers.vs.wp)
  
  markers.vs.durban = read.csv(paste0(durban_prefix,path,"Durban.csv"),row.names = 1)
  markers.vs.durban$cluster = "EE"
  markers.vs.durban$gene = rownames(markers.vs.durban)
  
  combined_markers = combine_markers(markers.vs.all,markers.vs.eoe,markers.vs.wp,markers.vs.durban)
  return(combined_markers)
}

####################################################################################################
#################################          Cell composition        #################################
####################################################################################################

# Reads in all_meta (@meta.data slot from Seurat Object)
# Returns the counts of each cell type
generate_all_counts = function(all_meta){
  all_meta$orig.ident = as.factor(all_meta$orig.ident)
  all_counts <- all_meta %>%
    group_by(cell_types, orig.ident,study,.drop=FALSE) %>%
    tally()
  temp <- rep(0,nrow(all_counts))
  for(i in 1:length(temp)){
    temp[i] <- all_counts$n[i]/(sum(all_counts$n[all_counts$orig.ident==all_counts$orig.ident[i]]))
  }
  all_counts$percent_of_sample <- temp
  temp <- rep(0,nrow(all_counts))
  for(i in 1:length(temp)){
    temp[i] <- all_counts$n[i]/(sum(all_counts$n[all_counts$cell_types==all_counts$cell_types[i]]))
  }
  all_counts$proportion_of_cell_type <- temp
  orig_score <- all_meta  %>% group_by(orig.ident,study) %>% tally()
  for(i in 1:nrow(all_counts)){
    if(is.na(all_counts$study[i])){
      all_counts$study[i] <- orig_score$study[orig_score$orig.ident==all_counts$orig.ident[i]]
    }
  }
  all_counts$orig.ident = as.character(all_counts$orig.ident)
  all_counts$cell_types = as.character(all_counts$cell_types)
  return(all_counts)
}

# generates the cell counts from each sample
generate_counts_by_sample = function(all_counts){
  origs <- unique(all_counts$orig.ident)
  cells <- unique(all_counts$cell_types)
  
  sample_all_counts <- data.frame(matrix(ncol = length(cells), nrow = 0))
  colnames(sample_all_counts) <- cells
  
  for(orig in origs){
    temp <- rep(0, length(cells))
    for(i in 1:length(cells)){        
      possible_cells <- all_counts$cell_types[all_counts$orig.ident==orig]
      if(cells[i] %in% possible_cells){
        temp[i] <- all_counts$n[(all_counts$orig.ident==orig & all_counts$cell_types==cells[i])]
      }
    }
    sample_all_counts <- rbind(sample_all_counts, temp)
  }
  rownames(sample_all_counts) <- origs
  colnames(sample_all_counts) <- cells
  return(sample_all_counts)
}

# generates the percent of each sample that each cell type constitutes
generate_percent_by_sample = function(all_counts){
  origs <- unique(all_counts$orig.ident)
  cells <- unique(all_counts$cell_types)
  
  sample_all_counts <- data.frame(matrix(ncol = length(cells), nrow = 0))
  colnames(sample_all_counts) <- cells
  
  for(orig in origs){
    temp <- rep(0, length(cells))
    for(i in 1:length(cells)){        
      possible_cells <- all_counts$cell_types[all_counts$orig.ident==orig]
      if(cells[i] %in% possible_cells){
        temp[i] <- all_counts$percent_of_sample[(all_counts$orig.ident==orig & all_counts$cell_types==cells[i])]
      }
    }
    sample_all_counts <- rbind(sample_all_counts, temp)
  }
  rownames(sample_all_counts) <- origs
  colnames(sample_all_counts) <- cells
  return(sample_all_counts)
}
generate_sample_meta = function(all_meta){
  sample_unique_pos <- c()
  samples <- c()
  for(i in 1:nrow(all_meta)){
    if(!(all_meta$orig.ident[i] %in% samples)){
      sample_unique_pos <- c(sample_unique_pos, i)
      samples <- c(samples, all_meta$orig.ident[i])
    }
  }
  
  sample_meta <- all_meta[sample_unique_pos,]
  rownames(sample_meta) <- 1:nrow(sample_meta)
  return(sample_meta)
}

# Runs Fisher's exact test to find cell populations with significantly different relative abundances
  # between EE and outgroup
run_fischer_test = function(s_obj,ident.group=Idents(s_obj),column_name = "study",group_name="EE",cell_types){
  pvals = rep(0,length(cell_types))
  signs = rep(0,length(cell_types))
  for(i in 1:length(cell_types)){
    clust_not_group = sum(ident.group==cell_types[i] &  s_obj@meta.data[,column_name] !=group_name)
    clust_group = sum(ident.group==cell_types[i] & s_obj@meta.data[,column_name]==group_name)
    not_clust_not_group = sum(ident.group!=cell_types[i] & s_obj@meta.data[,column_name]!=group_name)
    not_clust_group = sum(ident.group!=cell_types[i] & s_obj@meta.data[,column_name]==group_name)
    
    percent_not_group = clust_not_group/(clust_not_group +not_clust_not_group)
    percent_group = clust_group/(clust_group +not_clust_group)
    
    m = matrix(c(clust_not_group,clust_group,not_clust_not_group,not_clust_group),nrow=2)
    pvals[i] = fisher.test(m)$p.value
    signs[i] = sign(percent_group - percent_not_group)
  }
  p_adj = p.adjust(pvals,method="BH")
  return(data.frame(clusters=cell_types,adj_pval=p_adj,sign=signs))
}

# Final function for running composition testing
  # uses a leave one out approach, where Fisher's exact test is rerun
  # with each data point removed and the maximum pvalue from these leave-one-out tests
  # is taken as the final p value
fisher_composition_test <- function(s_obj,unique_types,cell_type_column_name,group_column_name,in_group_name,threshold=0.05){
  Idents(s_obj) <-cell_type_column_name
  s_obj$cell_types = s_obj@meta.data[,cell_type_column_name]
  all_counts = generate_all_counts(s_obj@meta.data)
  sample_all_counts = generate_counts_by_sample(all_counts)
  
  fisher_results = run_fischer_test(s_obj,column_name=group_column_name,group_name=in_group_name,cell_types = unique_types)
  fisher_results_sig = fisher_results %>% filter(adj_pval <= threshold)
  
  pvals = matrix(rep(1,nrow(sample_all_counts)*length(unique_types)),
                 nrow=nrow(sample_all_counts),
                 ncol=length(unique(unique_types)))
  colnames(pvals) <- fisher_results$clusters
  for(i in 1:length(unique(s_obj$orig.ident))){
    temp = subset(s_obj, orig.ident!=unique(s_obj$orig.ident)[i])
    pval_table = run_fischer_test(temp,column_name=group_column_name,group_name=in_group_name,cell_types = unique_types)
    
    for(j in 1:nrow(pval_table)){
      pvals[i,colnames(pvals)==pval_table$clusters[j]] = pval_table$adj_pval[j]
    }
    
  }
  
  max_pvals = apply(pvals,2,max)
  signs = fisher_results$sign
  signs = signs[order(max_pvals)]
  max_pvals = max_pvals[order(max_pvals)]
  max_table = data.frame(max_adj_pvals=max_pvals,sign=signs)
  max_table$cell_type = rownames(max_table)
  return(max_table)
  
}


####################################################################################################
##########################              Progeny & Dorothea             #############################
####################################################################################################

# Summarizes progeny scores for each cell type in Seurat object
get_progeny_for_cell_type = function(sobj){
  
  CellsDiseases <- data.frame(Cell = names(Idents(sobj)), 
                              cellType = as.character(sobj$study),
                              stringsAsFactors = FALSE)
  ## We transform Progeny scores into a data frame to better handling the results
  progeny_scores_df <- 
    as.data.frame(t(GetAssayData(sobj, slot = "scale.data", 
                                 assay = "progeny"))) %>%
    rownames_to_column("Cell") %>%
    gather(Pathway, Activity, -Cell) 
  
  ## We match Progeny scores with the cell clusters.
  progeny_scores_df <- inner_join(progeny_scores_df, CellsDiseases)
  ## We summarize the Progeny scores by cellpopulation
  summarized_progeny_scores <- progeny_scores_df %>% 
    group_by(Pathway, cellType) %>%
    summarise(avg = mean(Activity), std = sd(Activity))
  cohenD = c()
  for(i in 1:nrow(summarized_progeny_scores)){
    if(i %% 2 ==0){
      cohenD = c(cohenD, (summarized_progeny_scores$avg[i-1]-summarized_progeny_scores$avg[i])/(sqrt((summarized_progeny_scores$std[i-1]^2 + summarized_progeny_scores$std[i]^2)/2)))
    }
  }
  return(tibble(Pathway=unique(summarized_progeny_scores$Pathway),cohensD=cohenD))
}  

# Finds the differenes in progeny score between EE and outgroup
get_sobj_progeny_two_way = function(sobj,meta){
  
  # initialize matrix for pvalues
  meta="study"
  pvals = matrix(1, nrow=14,ncol=3)
  temp = GetAssayData(sobj,assay="progeny",slot="data")
  rownames(pvals) = rownames(temp)
  colnames(pvals) = c("pval","p_val_adj","sign") #1 if mean higher in EE, -1 if mean lower in EE
  
  for(i in 1:14){ # loop thru pathways
    pvals[i,1] = wilcox.test(temp[i,sobj[[meta]][,1]=="EE"],temp[i,sobj[[meta]][,1]!="EE"])$p.value
    if(mean(temp[i,sobj[[meta]][,1]=="EE"])>mean(temp[i,sobj[[meta]][,1]!="EE"])){
      pvals[i,3] = 1
    }else{
      pvals[i,3] = -1
    }
  }
  pvals[,2] = p.adjust(pvals[,1],method="BH")
  return(as.data.frame(pvals))
  
}

# Combines progeny scores in the same manner as DE functions to find significant progeny results

combine_markers_progeny = function(markers_all,markers_eoe,markers_wp,markers_durban,sig_threshold=0.05){
  sig_all = rownames(markers_all)[markers_all$p_val_adj < sig_threshold]
  marker_bools = as.vector(sapply(sig_all,function(x){check_marker_progeny(x,markers_eoe,markers_wp,markers_durban,1)}))
  sig_all = sig_all[marker_bools]
  return(markers_all[sig_all,])
}

check_marker_progeny = function(marker,markers_eoe,markers_wp,markers_durban,sig_threshold=1){
  in_eoe = marker %in% rownames(markers_eoe)[markers_eoe$p_val_adj < sig_threshold]
  in_wp = marker %in% rownames(markers_wp)[markers_wp$p_val_adj < sig_threshold]
  in_durban = marker %in% rownames(markers_durban)[markers_durban$p_val_adj < sig_threshold]
  if(in_eoe & in_wp & in_durban){
    if(sign(markers_eoe[marker,]$sign)==sign(markers_wp[marker,]$sign) & sign(markers_eoe[marker,]$sign)==sign(markers_durban[marker,]$sign)){
      return(TRUE)
    }else{
      return(FALSE)
    }
  }else{
    return(FALSE)
  }
}


# Summarizes Dorothea scores for each cell type in Seurat object
get_dorothea_for_cell_type = function(sobj){
  
  CellsDiseases <- data.frame(Cell = names(Idents(sobj)), 
                              cellType = as.character(sobj$study),
                              stringsAsFactors = FALSE)
  ## We transform dorothea scores into a data frame to better handling the results
  dorothea_scores_df <- 
    as.data.frame(t(GetAssayData(sobj, slot = "scale.data", 
                                 assay = "dorothea"))) %>%
    rownames_to_column("Cell") %>%
    gather(TF, Activity, -Cell) 
  
  ## We match dorothea scores with the cell clusters.
  dorothea_scores_df <- inner_join(dorothea_scores_df, CellsDiseases)
  ## We summarize the dorothea scores by cellpopulation
  summarized_dorothea_scores <- dorothea_scores_df %>% 
    group_by(TF, cellType) %>%
    # group_by(TF) %>%
    summarise(avg = mean(Activity), std = sd(Activity))
  cohenD = c()
  for(i in 1:nrow(summarized_dorothea_scores)){
    if(i %% 2 ==0){
      cohenD = c(cohenD, (summarized_dorothea_scores$avg[i-1]-summarized_dorothea_scores$avg[i])/(sqrt((summarized_dorothea_scores$std[i-1]^2 + summarized_dorothea_scores$std[i]^2)/2)))
    }
  }
  return(tibble(TF=unique(summarized_dorothea_scores$TF),cohensD=cohenD))
}  

# Finds the differenes in dorothea score between EE and outgroup
get_sobj_dorothea_two_way = function(sobj,meta){
  
  # initialize matrix for pvalues
  meta="study"
  pvals = matrix(1, nrow=nrow(GetAssayData(sobj,assay="dorothea",slot="data")),ncol=3)
  temp = GetAssayData(sobj,assay="dorothea",slot="data")
  rownames(pvals) = rownames(temp)
  colnames(pvals) = c("pval","p_val_adj","sign") #1 if mean higher in EE, -1 if mean lower in EE
  
  for(i in 1:nrow(GetAssayData(sobj,assay="dorothea",slot="data"))){ # loop thru TFs
    pvals[i,1] = wilcox.test(temp[i,sobj[[meta]][,1]=="EE"],temp[i,sobj[[meta]][,1]!="EE"])$p.value
    if(mean(temp[i,sobj[[meta]][,1]=="EE"])>mean(temp[i,sobj[[meta]][,1]!="EE"])){
      pvals[i,3] = 1
    }else{
      pvals[i,3] = -1
    }
  }
  pvals[,2] = p.adjust(pvals[,1],method="BH")
  return(as.data.frame(pvals))
  
}

# Combines Dorothea scores in the same manner as DE functions to find significant progeny results


combine_markers_dorothea = function(markers_all,markers_eoe,markers_wp,markers_durban,sig_threshold=0.05){
  sig_all = rownames(markers_all)[markers_all$p_val_adj < sig_threshold]
  marker_bools = as.vector(sapply(sig_all,function(x){check_marker_dorothea(x,markers_eoe,markers_wp,markers_durban,1)}))
  sig_all = sig_all[marker_bools]
  return(markers_all[sig_all,])
}

check_marker_dorothea = function(marker,markers_eoe,markers_wp,markers_durban,sig_threshold=1){
  in_eoe = marker %in% rownames(markers_eoe)[markers_eoe$p_val_adj < sig_threshold]
  in_wp = marker %in% rownames(markers_wp)[markers_wp$p_val_adj < sig_threshold]
  in_durban = marker %in% rownames(markers_durban)[markers_durban$p_val_adj < sig_threshold]
  if(in_eoe & in_wp & in_durban){
    if(sign(markers_eoe[marker,]$sign)==sign(markers_wp[marker,]$sign) & sign(markers_eoe[marker,]$sign)==sign(markers_durban[marker,]$sign)){
      return(TRUE)
    }else{
      return(FALSE)
    }
  }else{
    return(FALSE)
  }
}




