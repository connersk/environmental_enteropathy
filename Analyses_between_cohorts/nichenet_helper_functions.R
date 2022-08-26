
get_cell_subset_name = function(de_string){
  x = strsplit(de_string,"tier")[[1]]
  substring(strsplit(x[length(x)],"combined.csv")[[1]][1],2)
}

get_tier_for_sobj = function(de_string){
  x = strsplit(de_string,"_")[[1]]
  substring(x[length(x)],1,5)
}

run_nichenet = function(sobj,receiver,sender_celltypes,cell_type_markers){
  
  # adding extra filter to make sure that the DE response is strict
  combined_markers = cell_type_markers %>% filter(cluster==receiver) %>% filter(p_val_adj < 0.05)
  senders = sender_celltypes
  expressed_genes_receiver = get_expressed_genes(receiver, sobj, pct = 0.10)
  background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]
  
  list_expressed_genes_sender = senders %>% unique() %>% lapply(get_expressed_genes, sobj, 0.10) # lapply to get the expressed genes of every sender cell type separately here
  expressed_genes_sender = list_expressed_genes_sender %>% unlist() %>% unique()
  
  # A
  
  geneset_oi = combined_markers %>% filter(p_val_adj <= 0.05 & abs(avg_logFC) >= 0.25) %>% pull(gene)
  geneset_oi = geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]
  
  ligands = lr_network %>% pull(from) %>% unique()
  receptors = lr_network %>% pull(to) %>% unique()
  
  expressed_ligands = intersect(ligands,expressed_genes_sender)
  expressed_receptors = intersect(receptors,expressed_genes_receiver)
  
  potential_ligands = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% pull(from) %>% unique()
  
  ligand_activities = predict_ligand_activities(geneset = as.character(geneset_oi), background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)
  ligand_activities = ligand_activities %>% arrange(-pearson) %>% mutate(rank = rank(desc(pearson)))
  best_upstream_ligands = ligand_activities %>% top_n(20, pearson) %>% arrange(-pearson) %>% pull(test_ligand) %>% unique()
  
  active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 200) %>% bind_rows() %>% drop_na()
  
  active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix, cutoff = 0.33)
  
  order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev() %>% make.names()
  order_targets = active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links)) %>% make.names()
  rownames(active_ligand_target_links) = rownames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23
  colnames(active_ligand_target_links) = colnames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23
  
  vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()
  p_ligand_target_network = vis_ligand_target %>% make_heatmap_ggplot("Prioritized ligands","Predicted target genes", color = "purple",legend_position = "top", x_axis_position = "top",legend_title = "Regulatory potential")  + theme(axis.text.x = element_text(face = "italic")) + scale_fill_gradient2(low = "whitesmoke",  high = "purple", breaks = c(0,0.0045,0.0090))
  
  lr_network_top = lr_network %>% filter(from %in% best_upstream_ligands & to %in% expressed_receptors) %>% distinct(from,to)
  best_upstream_receptors = lr_network_top %>% pull(to) %>% unique()
  
  lr_network_top_df_large = weighted_networks_lr %>% filter(from %in% best_upstream_ligands & to %in% best_upstream_receptors)
  
  lr_network_top_df = lr_network_top_df_large %>% spread("from","weight",fill = 0)
  lr_network_top_matrix = lr_network_top_df %>% select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df$to)
  
  dist_receptors = dist(lr_network_top_matrix, method = "binary")
  hclust_receptors = hclust(dist_receptors, method = "ward.D2")
  order_receptors = hclust_receptors$labels[hclust_receptors$order]
  
  dist_ligands = dist(lr_network_top_matrix %>% t(), method = "binary")
  hclust_ligands = hclust(dist_ligands, method = "ward.D2")
  order_ligands_receptor = hclust_ligands$labels[hclust_ligands$order]
  
  order_receptors = order_receptors %>% intersect(rownames(lr_network_top_matrix))
  order_ligands_receptor = order_ligands_receptor %>% intersect(colnames(lr_network_top_matrix))
  
  vis_ligand_receptor_network = lr_network_top_matrix[order_receptors, order_ligands_receptor]
  rownames(vis_ligand_receptor_network) = order_receptors %>% make.names()
  colnames(vis_ligand_receptor_network) = order_ligands_receptor %>% make.names()
  p_ligand_receptor_network = vis_ligand_receptor_network %>% t() %>% make_heatmap_ggplot("Ligands","Receptors", color = "mediumvioletred", x_axis_position = "top",legend_title = "Prior interaction potential")
  
  
  lr_network_strict = lr_network %>% filter(database != "ppi_prediction_go" & database != "ppi_prediction")
  ligands_bona_fide = lr_network_strict %>% pull(from) %>% unique()
  receptors_bona_fide = lr_network_strict %>% pull(to) %>% unique()
  
  lr_network_top_df_large_strict = lr_network_top_df_large %>% distinct(from,to) %>% inner_join(lr_network_strict, by = c("from","to")) %>% distinct(from,to)
  lr_network_top_df_large_strict = lr_network_top_df_large_strict %>% inner_join(lr_network_top_df_large, by = c("from","to"))
  
  lr_network_top_df_strict = lr_network_top_df_large_strict %>% spread("from","weight",fill = 0)
  lr_network_top_matrix_strict = lr_network_top_df_strict %>% select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df_strict$to)
  
  dist_receptors = dist(lr_network_top_matrix_strict, method = "binary")
  hclust_receptors = hclust(dist_receptors, method = "ward.D2")
  order_receptors = hclust_receptors$labels[hclust_receptors$order]
  
  dist_ligands = dist(lr_network_top_matrix_strict %>% t(), method = "binary")
  hclust_ligands = hclust(dist_ligands, method = "ward.D2")
  order_ligands_receptor = hclust_ligands$labels[hclust_ligands$order]
  
  order_receptors = order_receptors %>% intersect(rownames(lr_network_top_matrix_strict))
  order_ligands_receptor = order_ligands_receptor %>% intersect(colnames(lr_network_top_matrix_strict))
  
  vis_ligand_receptor_network_strict = lr_network_top_matrix_strict[order_receptors, order_ligands_receptor]
  rownames(vis_ligand_receptor_network_strict) = order_receptors %>% make.names()
  colnames(vis_ligand_receptor_network_strict) = order_ligands_receptor %>% make.names()
  p_ligand_receptor_network_strict = vis_ligand_receptor_network_strict %>% t() %>% make_heatmap_ggplot("Ligands","Receptors", color = "mediumvioletred", x_axis_position = "top",legend_title = "Prior interaction potential\n(bona fide)")
  
  # ligand activity heatmap
  ligand_pearson_matrix = ligand_activities %>% select(pearson) %>% as.matrix() %>% magrittr::set_rownames(ligand_activities$test_ligand)
  rownames(ligand_pearson_matrix) = rownames(ligand_pearson_matrix) %>% make.names()
  colnames(ligand_pearson_matrix) = colnames(ligand_pearson_matrix) %>% make.names()
  
  vis_ligand_pearson = ligand_pearson_matrix[order_ligands, ] %>% as.matrix(ncol = 1) %>% magrittr::set_colnames("Pearson")
  p_ligand_pearson = vis_ligand_pearson %>% make_heatmap_ggplot("Prioritized ligands","Ligand activity", color = "darkorange",legend_position = "top", x_axis_position = "top", legend_title = "Pearson correlation coefficient\ntarget gene prediction ability)") + theme(legend.text = element_text(size = 9))
  
  # ligand expression Seurat dotplot
  order_ligands_adapted = order_ligands
  order_ligands_adapted[order_ligands_adapted == "H2.M3"] = "H2-M3" # cf required use of make.names for heatmap visualization | this is not necessary if these ligands are not in the list of prioritized ligands!
  order_ligands_adapted[order_ligands_adapted == "H2.T23"] = "H2-T23" # cf required use of make.names for heatmap visualization | this is not necessary if these ligands are not in the list of prioritized ligands!
  #     temp = subset(sobj,tier2 %in% senders)
  #     rotated_dotplot = DotPlot(temp, features = order_ligands_adapted %>% rev(), cols = "RdYlBu") + coord_flip() + theme(legend.text = element_text(size = 10), legend.title = element_text(size = 12)) # flip of coordinates necessary because we want to show ligands in the rows when combining all plots
  
  DE_table_all = cell_type_markers
  DE_table_all$gene = as.character(DE_table_all$gene)
  
  # Combine ligand activities with DE information
  ligand_activities_de = ligand_activities %>% select(test_ligand, pearson) %>% rename(ligand = test_ligand) %>% left_join(DE_table_all %>% rename(ligand = gene))
  print(head(ligand_activities_de))
  ligand_activities_de[is.na(ligand_activities_de)] = 0
  ligand_activities_de = ligand_activities_de %>% filter(avg_logFC!=0)
  
  # # # make LFC heatmap
  lfc_matrix = ligand_activities_de  %>% select(-pearson,-p_val,-pct.1,-pct.2,-p_val_adj)
  lfc_matrix  = lfc_matrix %>% spread(cluster,avg_logFC)
  lfc_matrix [is.na(lfc_matrix )] = 0
  ligands = lfc_matrix$ligand
  lfc_matrix = lfc_matrix %>% select(-ligand) %>% as.matrix
  rownames(lfc_matrix) = ligands
  # rownames(lfc_matrix) = rownames(lfc_matrix) %>% make.names()
  
  
  order_ligands = order_ligands[order_ligands %in% rownames(lfc_matrix)]
  vis_ligand_lfc = lfc_matrix[order_ligands,]
  # rownames(vis_ligand_lfc) = paste0(rownames(vis_ligand_lfc)," ",vis_ligand_lfc[,6])
  # temp = rownames(vis_ligand_lfc)
  # vis_ligand_lfc = apply(vis_ligand_lfc[,-6],2,as.numeric)
  # rownames(vis_ligand_lfc) = temp
  # vis_ligand_lfc = vis_ligand_lfc[vis_ligand_lfc[,1]!=0,]
  
  colnames(vis_ligand_lfc) = vis_ligand_lfc %>% colnames() %>% make.names()
  p_ligand_lfc = vis_ligand_lfc %>% make_threecolor_heatmap_ggplot("Prioritized ligands","LFC in Sender", low_color = "midnightblue",mid_color = "white", mid = median(vis_ligand_lfc), high_color = "red",legend_position = "top", x_axis_position = "top", legend_title = "LFC") + theme(axis.text.y = element_text(face = "italic"))
  
  return(list(p_ligand_pearson,p_ligand_target_network,order_ligands_adapted,p_ligand_lfc))
}

run_nichenet_data = function(sobj,receiver,sender_celltypes,cell_type_markers){
  
  combined_markers = cell_type_markers %>% filter(cluster==receiver) %>% filter(p_val_adj < 0.05)
  senders = sender_celltypes
  expressed_genes_receiver = get_expressed_genes(receiver, sobj, pct = 0.10)
  background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]
  
  list_expressed_genes_sender = senders %>% unique() %>% lapply(get_expressed_genes, sobj, 0.10) # lapply to get the expressed genes of every sender cell type separately here
  expressed_genes_sender = list_expressed_genes_sender %>% unlist() %>% unique()
  
  # goi
  
  geneset_oi = combined_markers %>% filter(p_val_adj <= 0.05 & abs(avg_logFC) >= 0.25) %>% pull(gene)
  geneset_oi = geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]
  
  ligands = lr_network %>% pull(from) %>% unique()
  receptors = lr_network %>% pull(to) %>% unique()
  
  expressed_ligands = intersect(ligands,expressed_genes_sender)
  expressed_receptors = intersect(receptors,expressed_genes_receiver)
  
  potential_ligands = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% pull(from) %>% unique()
  
  ligand_activities = predict_ligand_activities(geneset = as.character(geneset_oi), background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)
  ligand_activities = ligand_activities %>% arrange(-pearson) %>% mutate(rank = rank(desc(pearson)))
  best_upstream_ligands = ligand_activities %>% top_n(20, pearson) %>% arrange(-pearson) %>% pull(test_ligand) %>% unique()
  
  active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 200) %>% bind_rows() %>% drop_na()
  
  active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix, cutoff = 0.33)
  
  order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev() %>% make.names()
  order_targets = active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links)) %>% make.names()
  rownames(active_ligand_target_links) = rownames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23
  colnames(active_ligand_target_links) = colnames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23
  
  vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()
  p_ligand_target_network = vis_ligand_target %>% make_heatmap_ggplot("Prioritized ligands","Predicted target genes", color = "purple",legend_position = "top", x_axis_position = "top",legend_title = "Regulatory potential")  + theme(axis.text.x = element_text(face = "italic")) + scale_fill_gradient2(low = "whitesmoke",  high = "purple", breaks = c(0,0.0045,0.0090))
  
  lr_network_top = lr_network %>% filter(from %in% best_upstream_ligands & to %in% expressed_receptors) %>% distinct(from,to)
  best_upstream_receptors = lr_network_top %>% pull(to) %>% unique()
  
  lr_network_top_df_large = weighted_networks_lr %>% filter(from %in% best_upstream_ligands & to %in% best_upstream_receptors)
  
  lr_network_top_df = lr_network_top_df_large %>% spread("from","weight",fill = 0)
  lr_network_top_matrix = lr_network_top_df %>% select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df$to)
  
  dist_receptors = dist(lr_network_top_matrix, method = "binary")
  hclust_receptors = hclust(dist_receptors, method = "ward.D2")
  order_receptors = hclust_receptors$labels[hclust_receptors$order]
  
  dist_ligands = dist(lr_network_top_matrix %>% t(), method = "binary")
  hclust_ligands = hclust(dist_ligands, method = "ward.D2")
  order_ligands_receptor = hclust_ligands$labels[hclust_ligands$order]
  
  order_receptors = order_receptors %>% intersect(rownames(lr_network_top_matrix))
  order_ligands_receptor = order_ligands_receptor %>% intersect(colnames(lr_network_top_matrix))
  
  vis_ligand_receptor_network = lr_network_top_matrix[order_receptors, order_ligands_receptor]
  rownames(vis_ligand_receptor_network) = order_receptors %>% make.names()
  colnames(vis_ligand_receptor_network) = order_ligands_receptor %>% make.names()
  p_ligand_receptor_network = vis_ligand_receptor_network %>% t() %>% make_heatmap_ggplot("Ligands","Receptors", color = "mediumvioletred", x_axis_position = "top",legend_title = "Prior interaction potential")
  
  
  lr_network_strict = lr_network %>% filter(database != "ppi_prediction_go" & database != "ppi_prediction")
  ligands_bona_fide = lr_network_strict %>% pull(from) %>% unique()
  receptors_bona_fide = lr_network_strict %>% pull(to) %>% unique()
  
  lr_network_top_df_large_strict = lr_network_top_df_large %>% distinct(from,to) %>% inner_join(lr_network_strict, by = c("from","to")) %>% distinct(from,to)
  lr_network_top_df_large_strict = lr_network_top_df_large_strict %>% inner_join(lr_network_top_df_large, by = c("from","to"))
  
  lr_network_top_df_strict = lr_network_top_df_large_strict %>% spread("from","weight",fill = 0)
  lr_network_top_matrix_strict = lr_network_top_df_strict %>% select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df_strict$to)
  
  dist_receptors = dist(lr_network_top_matrix_strict, method = "binary")
  hclust_receptors = hclust(dist_receptors, method = "ward.D2")
  order_receptors = hclust_receptors$labels[hclust_receptors$order]
  
  dist_ligands = dist(lr_network_top_matrix_strict %>% t(), method = "binary")
  hclust_ligands = hclust(dist_ligands, method = "ward.D2")
  order_ligands_receptor = hclust_ligands$labels[hclust_ligands$order]
  
  order_receptors = order_receptors %>% intersect(rownames(lr_network_top_matrix_strict))
  order_ligands_receptor = order_ligands_receptor %>% intersect(colnames(lr_network_top_matrix_strict))
  
  vis_ligand_receptor_network_strict = lr_network_top_matrix_strict[order_receptors, order_ligands_receptor]
  rownames(vis_ligand_receptor_network_strict) = order_receptors %>% make.names()
  colnames(vis_ligand_receptor_network_strict) = order_ligands_receptor %>% make.names()
  p_ligand_receptor_network_strict = vis_ligand_receptor_network_strict %>% t() %>% make_heatmap_ggplot("Ligands","Receptors", color = "mediumvioletred", x_axis_position = "top",legend_title = "Prior interaction potential\n(bona fide)")
  
  # ligand activity heatmap
  ligand_pearson_matrix = ligand_activities %>% select(pearson) %>% as.matrix() %>% magrittr::set_rownames(ligand_activities$test_ligand)
  rownames(ligand_pearson_matrix) = rownames(ligand_pearson_matrix) %>% make.names()
  colnames(ligand_pearson_matrix) = colnames(ligand_pearson_matrix) %>% make.names()
  
  vis_ligand_pearson = ligand_pearson_matrix[order_ligands, ] %>% as.matrix(ncol = 1) %>% magrittr::set_colnames("Pearson")
  p_ligand_pearson = vis_ligand_pearson %>% make_heatmap_ggplot("Prioritized ligands","Ligand activity", color = "darkorange",legend_position = "top", x_axis_position = "top", legend_title = "Pearson correlation coefficient\ntarget gene prediction ability)") + theme(legend.text = element_text(size = 9))
  
  # ligand expression Seurat dotplot
  order_ligands_adapted = order_ligands
  order_ligands_adapted[order_ligands_adapted == "H2.M3"] = "H2-M3" # cf required use of make.names for heatmap visualization | this is not necessary if these ligands are not in the list of prioritized ligands!
  order_ligands_adapted[order_ligands_adapted == "H2.T23"] = "H2-T23" # cf required use of make.names for heatmap visualization | this is not necessary if these ligands are not in the list of prioritized ligands!
  #     temp = subset(sobj,tier2 %in% senders)
  #     rotated_dotplot = DotPlot(temp, features = order_ligands_adapted %>% rev(), cols = "RdYlBu") + coord_flip() + theme(legend.text = element_text(size = 10), legend.title = element_text(size = 12)) # flip of coordinates necessary because we want to show ligands in the rows when combining all plots
  
  DE_table_all = cell_type_markers
  DE_table_all$gene = as.character(DE_table_all$gene)
  
  # Combine ligand activities with DE information
  ligand_activities_de = ligand_activities %>% select(test_ligand, pearson) %>% rename(ligand = test_ligand) %>% left_join(DE_table_all %>% rename(ligand = gene))
  print(head(ligand_activities_de))
  ligand_activities_de[is.na(ligand_activities_de)] = 0
  ligand_activities_de = ligand_activities_de %>% filter(avg_logFC!=0) %>% filter(!(is.na(avg_logFC)))
  
  # # # make LFC heatmap
  lfc_matrix = ligand_activities_de  %>% select(-pearson,-p_val,-pct.1,-pct.2,-p_val_adj)
  lfc_matrix  = lfc_matrix %>% spread(cluster,avg_logFC)
  lfc_matrix [is.na(lfc_matrix )] = 0
  ligands = lfc_matrix$ligand
  lfc_matrix = lfc_matrix %>% select(-ligand) %>% as.matrix
  rownames(lfc_matrix) = ligands
  # rownames(lfc_matrix) = rownames(lfc_matrix) %>% make.names()
  
  
  order_ligands = order_ligands[order_ligands %in% rownames(lfc_matrix)]
  vis_ligand_lfc = lfc_matrix[order_ligands,]
  # rownames(vis_ligand_lfc) = paste0(rownames(vis_ligand_lfc)," ",vis_ligand_lfc[,6])
  # temp = rownames(vis_ligand_lfc)
  # vis_ligand_lfc = apply(vis_ligand_lfc[,-6],2,as.numeric)
  # rownames(vis_ligand_lfc) = temp
  # vis_ligand_lfc = vis_ligand_lfc[vis_ligand_lfc[,1]!=0,]
  
  colnames(vis_ligand_lfc) = vis_ligand_lfc %>% colnames() %>% make.names()
  p_ligand_lfc = vis_ligand_lfc %>% make_threecolor_heatmap_ggplot("Prioritized ligands","LFC in Sender", low_color = "midnightblue",mid_color = "white", mid = median(vis_ligand_lfc), high_color = "red",legend_position = "top", x_axis_position = "top", legend_title = "LFC") + theme(axis.text.y = element_text(face = "italic"))
  
  return(list(vis_ligand_pearson,vis_ligand_target,order_ligands_adapted,vis_ligand_lfc))
}

