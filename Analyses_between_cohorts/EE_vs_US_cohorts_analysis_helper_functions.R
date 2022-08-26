# Helper functions for comparing EE to both U.S. cohorts, broken down by analysis type


####################################################################################################
#################################      Differential expression     #################################
####################################################################################################


# For a given cell type, these functions

# read in DE results from the following analyses
# EE vs. both U.S. cohorts
# EE vs. EoE cohort
# EE. vs. Resection cohort

# Returns genes that are 
# significantly upregulated in EE vs. both U.S. cohorts
# display the same direction of regulation when comparing EE vs. each cohort individually

combine_markers = function(markers_all,markers_eoe,markers_wp,sig_threshold=0.05){
  bool_for_sig_vs_outgroups = markers_all$p_val_adj < sig_threshold
  
  marker_genes = markers_all$gene
  bools_for_direction = as.vector(sapply(marker_genes,function(x){check_marker(x,markers_eoe,markers_wp,1.0)}))
  final_genes = marker_genes[bool_for_sig_vs_outgroups & bools_for_direction]
  return(markers_all %>% filter(gene %in% final_genes))
}

check_marker = function(marker,markers_eoe,markers_wp,sig_threshold=1.0){
  in_eoe = marker %in% markers_eoe$gene[markers_eoe$p_val_adj <= sig_threshold]
  in_wp = marker %in% markers_wp$gene[markers_wp$p_val_adj <= sig_threshold]
  if(in_eoe & in_wp){
    if(sign(markers_eoe$avg_logFC[markers_eoe$gene==marker])==sign(markers_wp$avg_logFC[markers_wp$gene==marker])){
      return(TRUE)
    }else{
      return(FALSE)
    }
  }else{
    return(FALSE)
  }
}


read_markers = function(path){
  markers.vs.all = read.csv(paste0("both/",path,"both.csv"),row.names = 1)
  markers.vs.all$cluster = "EE"
  markers.vs.all$gene = rownames(markers.vs.all)
  markers.vs.eoe = read.csv(paste0("eoe/",path,"EoE.csv"),row.names = 1)
  markers.vs.eoe$cluster = "EE"
  markers.vs.eoe$gene = rownames(markers.vs.eoe)
  markers.vs.wp = read.csv(paste0("resection/",path,"Resection.csv"),row.names = 1)
  markers.vs.wp$cluster = "EE"
  markers.vs.wp$gene = rownames(markers.vs.wp)
  combined_markers = combine_markers(markers.vs.all,markers.vs.eoe,markers.vs.wp)
  return(combined_markers)
}

read_markers_durban = function(path){
  print("here")
  print(path)
  markers.vs.all = read.csv(paste0("~/Dropbox (MIT)/Zambia/durban/DE/SA_vs_both_US/",path,"all.csv"),row.names = 1)
  markers.vs.all$cluster = "Durban"
  markers.vs.all$gene = rownames(markers.vs.all)
  markers.vs.eoe = read.csv(paste0("~/Dropbox (MIT)/Zambia/reseq_analysis/fixed_de/eoe/",path,"EoE.csv"),row.names = 1)
  markers.vs.eoe$cluster = "Durban"
  markers.vs.eoe$gene = rownames(markers.vs.eoe)
  markers.vs.wp = read.csv(paste0("~/Dropbox (MIT)/Zambia/reseq_analysis/fixed_de/resection/",path,"Resection.csv"),row.names = 1)
  markers.vs.wp$cluster = "Durban"
  markers.vs.wp$gene = rownames(markers.vs.wp)
  combined_markers = combine_markers(markers.vs.all,markers.vs.eoe,markers.vs.wp)
  return(combined_markers)
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
get_all_way_progeny_progeny = function(sobj,meta){
  sobj_test = subset(sobj, study %in% c("EoE","EE"))
  results_eoe = get_sobj_progeny_two_way(sobj_test,"study")
  sobj_west = subset(sobj,study %in% c("EE","Uninvolved cancer resection"))
  results_whipple = get_sobj_progeny_two_way(sobj_west,"study")
  sobj$study[sobj$study=="Uninvolved cancer resection"] = "not EE"
  sobj$study[sobj$study=="EoE"] = "not EE"
  results_all =get_sobj_progeny_two_way(sobj,"study")
  results_all = cbind(results_all,results_eoe,results_whipple)
  results_all = results_all[,-c(1,3,5)]
  colnames(results_all) = c("adj_pval_all","adj_pval_eoe","adj_pval_whipple")
  return(results_all[order(rowSums(results_all),decreasing=F),])
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


combine_markers_progeny = function(markers_all,markers_eoe,markers_wp,sig_threshold=0.05){
  sig_all = rownames(markers_all)[markers_all$p_val_adj < sig_threshold]
  marker_bools = as.vector(sapply(sig_all,function(x){check_marker_progeny(x,markers_eoe,markers_wp,1)}))
  sig_all = sig_all[marker_bools]
  return(markers_all[sig_all,])
}

check_marker_progeny = function(marker,markers_eoe,markers_wp,sig_threshold=1){
  in_eoe = marker %in% rownames(markers_eoe)[markers_eoe$p_val_adj < sig_threshold]
  in_wp = marker %in% rownames(markers_wp)[markers_wp$p_val_adj < sig_threshold]
  if(in_eoe & in_wp){
    if(sign(markers_eoe[marker,]$sign)==sign(markers_wp[marker,]$sign)){
      return(TRUE)
    }else{
      return(FALSE)
    }
  }else{
    return(FALSE)
  }
}

# Combines progeny scores in the same manner as DE functions to find significant progeny results

combine_markers_dorothea = function(markers_all,markers_eoe,markers_wp,sig_threshold=0.05){
  sig_all = rownames(markers_all)[markers_all$p_val_adj < sig_threshold]
  marker_bools = as.vector(sapply(sig_all,function(x){check_marker_dorothea(x,markers_eoe,markers_wp,1)}))
  sig_all = sig_all[marker_bools]
  return(markers_all[sig_all,])
}

check_marker_dorothea = function(marker,markers_eoe,markers_wp,sig_threshold=1){
  in_eoe = marker %in% rownames(markers_eoe)[markers_eoe$p_val_adj < sig_threshold]
  in_wp = marker %in% rownames(markers_wp)[markers_wp$p_val_adj < sig_threshold]
  if(in_eoe & in_wp){
    if(sign(markers_eoe[marker,]$sign)==sign(markers_wp[marker,]$sign)){
      return(TRUE)
    }else{
      return(FALSE)
    }
  }else{
    return(FALSE)
  }
}

# Finds the differenes in dorothea score between EE and outgroup

get_sobj_dorothea_two_way = function(sobj,meta){
  
  # initialize matrix for pvalues
  meta="study"
  pvals = matrix(1, nrow=292,ncol=3)
  temp = GetAssayData(sobj,assay="dorothea",slot="data")
  rownames(pvals) = rownames(temp)
  colnames(pvals) = c("pval","p_val_adj","sign") #1 if mean higher in EE, -1 if mean lower in EE
  
  for(i in 1:292){ # loop thru pathways
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

# Summarizes Dorothea scores for each cell type in Seurat object
get_dorothea_for_cell_type = function(sobj){
  viper_scores_df <- GetAssayData(sobj, slot = "scale.data",  assay = "dorothea") %>% data.frame() %>% t()
  
  CellsClusters <- data.frame(cell = names(Idents(sobj)), 
                              study = as.character(sobj$study),
                              stringsAsFactors = FALSE)
  
  viper_scores_clusters <- viper_scores_df  %>% data.frame() %>% rownames_to_column("cell") %>% gather(tf, activity, -cell) %>% inner_join(CellsClusters)
  summarized_viper_scores <- viper_scores_clusters %>% 
    group_by(tf, study) %>%
    summarise(avg = mean(activity),
              std = sd(activity))
  
  cohenD = c()
  for(i in 1:nrow(summarized_viper_scores)){
    if(i %% 2 ==0){
      cohenD = c(cohenD, (summarized_viper_scores$avg[i-1]-summarized_viper_scores$avg[i])/(sqrt((summarized_viper_scores$std[i-1]^2 + summarized_viper_scores$std[i]^2)/2)))
    }
  }
  
  return(tibble(tf=unique(summarized_viper_scores$tf),cohensD=cohenD))

}




