# show colors -------
library(scales)
#show_col(col)

# Read_files ---------
func_read_files<-function(sample.name,datadir){
  seg = read.table(file = dir(datadir, full.names=T, pattern="uber.*.seg.txt"), header = T, sep = "\t")
  bin = read.table(file = dir(datadir, full.names=T, pattern="uber.*.bin.txt"), header = T, sep = "\t") 
  ratio = read.table(file = dir(datadir, full.names=T, pattern="uber.*.ratio.txt"), header = T, sep = "\t")
  seg = seg %>% dplyr::filter(!chrom %in% 24)
  bin = bin %>% dplyr::filter(!chrom %in% 24)
  ratio = ratio %>% dplyr::filter(!chrom %in% 24)
  seg_cndata = data.frame(seg[,-c(1:3)])
  bin_cndata = data.frame(bin[,-c(1:3)])
  ratio_cndata = data.frame(ratio[,-c(1:3)])
  chrom_chrompos_abspos = bin[1:3]
  files = list(seg,bin,ratio,seg_cndata,bin_cndata,ratio_cndata,chrom_chrompos_abspos); 
  names(files) = c("seg","bin","ratio","seg_cndata","bin_cndata","ratio_cndata","chrom_chrompos_abspos")
  message("import seg, bin, ratio matrics")
  return(files)
}

# get_groups ------

get_timepoint <- function(x, ifrowname = T) {
  #x = as.data.frame(t(files$seg_cndata))
  Primary = grep(rownames(x),pattern = "P|_P_|primary",ignore.case = T) #P|
  Recurrence = grep(rownames(x),pattern = "RE|recurrence",ignore.case = T) #R|_R_|
  allcluster = c("Primary", "Recurrence")
  Group = rep(NA, nrow(x))
  Group[Primary] <- "Primary"
  Group[Recurrence] <- "Recurrence"
  xcluster = unique(Group)
  levels = allcluster[allcluster %in% xcluster]
  Groups <- data.frame(Group); rownames(Groups) <- rownames(x); table(Groups); colnames(Groups)<-"G"
  Groups[which(is.na(Groups)),] <- "Primary"
  Groups$G = factor(Groups$G, levels = levels)
  Groups_2<-factor(Groups[,1],levels = levels)
  if(ifrowname == T ){return(Groups)}
  if(ifrowname == F ){return(Groups_2)}
}

Conv_ratio2int <- function(x, prim_ploidy, rec_ploidy) {
  
  # get the df's
  x_df = x$seg_cndata %>% tibble::as_tibble()
  y_df = x$well %>% tibble::as_tibble()
  z_df = x$Groups %>% tibble::as_tibble()
  
  x_df_tr = t(x_df) %>%  tibble::as_tibble()
  x_df_tr$cells = as.character(colnames(x_df))
  y_df$cells = as.character(y_df$cells)
  z_df$cells = rownames(x$Groups)
  #table(colnames(x_df) == y_df$cells)
  
  # merge the df
  xy_df = dplyr::left_join(x_df_tr, y_df)
  xy_df = dplyr::left_join(xy_df, z_df)
  
  # convert to long
  xy_df_long = reshape2::melt(xy_df, id.vars = c("well", "cells", "G"),
                              variable.name = "seg", value.name = "seg_rat")
  
  # convert to int CN
  xy_df_long = xy_df_long %>% dplyr::mutate(seg_rat_int = dplyr::case_when((well == "a_well") & (G == "Primary") ~ seg_rat*prim_ploidy,
                                                                           (well %in% "a_well") & (G %in% "Recurrence") ~ seg_rat*rec_ploidy,
                                                                           (well %in% c("d_well", "pop_well") ~ seg_rat*2)))
  
  # convert to int CN, all cells
  xy_df_long = xy_df_long %>% dplyr::mutate(seg_rat_int_all = dplyr::case_when(G == "Primary" ~ seg_rat*prim_ploidy,
                                                                               G == "Recurrence" ~ seg_rat*rec_ploidy))
  
  xy_df_long$seg_rat_int = round(xy_df_long$seg_rat_int)
  xy_df_long$seg_rat_int_all = round(xy_df_long$seg_rat_int_all)
  
  # convert back to wide
  xy_df_int = reshape2::recast(xy_df_long, cells ~ seg, measure.var = "seg_rat_int")
  xy_df_int1 = xy_df_int %>% dplyr::select(-cells)
  rownames(xy_df_int1) = xy_df_int$cells
  
  xy_df_int_all = reshape2::recast(xy_df_long, cells ~ seg, measure.var = "seg_rat_int_all")
  xy_df_int_all1 = xy_df_int_all %>% dplyr::select(-cells)
  rownames(xy_df_int_all1) = xy_df_int_all$cells
  
  #Groups_2<-factor(Groups[,1],levels = levels)
  {return(list(int_col = t(xy_df_int1), int_all_col = t(xy_df_int_all1)))}
}


Conv_ratio2int_pcf <- function(x, prim_ploidy, rec_ploidy) {
  
  # get the df's
  x_df = x$pcf_seg_exp %>% tibble::as_tibble()
  y_df = x$well %>% tibble::as_tibble()
  z_df = x$Groups %>% tibble::as_tibble()
  
  #x_df_tr = t(x_df) %>%  tibble::as_tibble()
  x_df$cells = as.character(rownames(x$pcf_seg_exp))
  y_df$cells = as.character(y_df$cells)
  z_df$cells = rownames(x$Groups)
  #table(colnames(x_df) == y_df$cells)
  
  # merge the df
  xy_df = dplyr::left_join(x_df, y_df)
  xy_df = dplyr::left_join(xy_df, z_df)
  
  # convert to long
  xy_df_long = reshape2::melt(xy_df, id.vars = c("well", "cells", "G"),
                              variable.name = "seg", value.name = "seg_rat")
  
  # convert to int CN
  xy_df_long = xy_df_long %>% dplyr::mutate(seg_rat_int = dplyr::case_when((well == "a_well") & (G == "Primary") ~ seg_rat*prim_ploidy,
                                                                           (well %in% "a_well") & (G %in% "Recurrence") ~ seg_rat*rec_ploidy,
                                                                           (well %in% c("d_well", "pop_well") ~ seg_rat*2)))
  
  # convert to int CN, all cells
  xy_df_long = xy_df_long %>% dplyr::mutate(seg_rat_int_all = dplyr::case_when(G == "Primary" ~ seg_rat*prim_ploidy,
                                                                               G == "Recurrence" ~ seg_rat*rec_ploidy))
  
  xy_df_long$seg_rat_int = round(xy_df_long$seg_rat_int)
  xy_df_long$seg_rat_int_all = round(xy_df_long$seg_rat_int_all)
  
  # convert back to wide
  xy_df_int = reshape2::recast(xy_df_long, cells ~ seg, measure.var = "seg_rat_int")
  xy_df_int1 = xy_df_int %>% dplyr::select(-cells)
  rownames(xy_df_int1) = xy_df_int$cells
  
  xy_df_int_all = reshape2::recast(xy_df_long, cells ~ seg, measure.var = "seg_rat_int_all")
  xy_df_int_all1 = xy_df_int_all %>% dplyr::select(-cells)
  rownames(xy_df_int_all1) = xy_df_int_all$cells
  
  #Groups_2<-factor(Groups[,1],levels = levels)
  {return(list(int_col = t(xy_df_int1), int_all_col = t(xy_df_int_all1)))}
}

# func_extract_well -------

func_extract_well <- function(x, return_df = TRUE) {
  
  a_well = grep(rownames(x),pattern = "*A_",ignore.case = T)
  d_well = grep(rownames(x),pattern = "*D_",ignore.case = T)
  pop_well = grep(rownames(x),pattern = "*POP_",ignore.case = T)
  allwells = c("a_well", "d_well", "pop_well")
  Well = rep(NA, nrow(x))
  Well[a_well] <- "a_well"
  Well[d_well] <- "d_well"
  Well[pop_well] <- "pop_well"
  xcluster = unique(Well)
  levels = allwells[allwells %in% xcluster]
  Well <- data.frame(Well); rownames(Well) <- rownames(x); table(Well); colnames(Well)<-"W"
  Well[which(is.na(Well)),] <- "a_well"
  Well$W = factor(Well$W, levels = levels)
  Well_df = data.frame(cells = rownames(x), well = Well$W)
  
  
  if(return_df == TRUE ){return(Well_df)}
  if(return_df == FALSE ){return(Well)}
  
}

# function to make heatmap ------
fun_htmap <- function(data, chr_lengths,  show_row_names, Groups, ifCN = F){
  if(ifCN == T){
    title = "Copy Number"
    breaks = c(0,1,2,4,6)
    dat_mat = as.data.frame.matrix(data$seg_cndata_int_all) %>% dplyr::mutate(across(.fns = ~case_when(. > 8 ~ "8", TRUE ~ as.character(.))))
    dat_mat = as.matrix(t((dat_mat)))
    col_vec = structure(c("#181C43", "#2641A2", "#247ABA", "#73A8BD", "white","#E4CAC3", "#CF8971", "#B8462D", "#850E28"), #, "#3C0912"
                        names = c(0:8)
    )
    
  }else{
    breaks = c(-2,0,2)
    dat_mat = as.matrix(t(log2(data$seg_cndata)))
    title = "Log2 (Ratio)"
    col_vec = circlize::colorRamp2(breaks =breaks, 
                                   c("dodgerblue3", "white", "firebrick3"))
  }
  
  data$G$cells = rownames(data$G)
  row_anno_df = dplyr::left_join(data$well, data$G, by = "cells")
  
  col_grp = c("Primary" = "maroon", "Recurrence" = "#009999")
  col_peak = c("a_well" = "grey", "d_well" = "seagreen", "pop_well" = "black")
  
  ha_row=rowAnnotation(df = row_anno_df %>% dplyr::select(G, well), 
                       col = list(G = col_grp, well = col_peak), 
                       show_annotation_name = T)
  chr_lengths = readRDS(file = "/seg_chrom_chrompos_abspos.rds")
  chr_lengths = chr_lengths %>% 
    dplyr::select(abspos, chrom, chrompos) %>% 
    dplyr::group_by(chrom) %>%
    dplyr::summarize(n = dplyr::n()) %>%
    pull(n)
  chr_lengths = chr_lengths[-24]
  chr_binary = rep(c(2,1), 11);chr_binary = c(chr_binary, 2)
  chr1 = data.frame(chr = rep.int(x = chr_binary, times = chr_lengths))
  # getting lengths for chr numbers annotation
  chr_rl_c = c(1, cumsum(chr_lengths))
  # creating a data frame to calculate rowMeans
  chr_df =  data.frame(a = chr_rl_c[1:length(chr_rl_c)-1],b= chr_rl_c[2:length(chr_rl_c)])
  chr_l_means = round(rowMeans(chr_df))
  chrom.names = c(1:22,"X")
  # creating the vector for chr number annotations
  v = vector(mode = "character",length = cumsum(chr_lengths)[length(chr_lengths)])
  v[chr_l_means] = chrom.names
  v[is.na(v)] = ""
  chr_bar = HeatmapAnnotation(chr_text = anno_text(v, gp = gpar(fontsize = 8)),
                              df = chr1,
                              show_legend = FALSE,
                              which = "column",
                              col = list(chr = c("1" = "grey88", "2" = "black")))
  
  
  
  
  ht_all = Heatmap(dat_mat,
                   cluster_columns = FALSE,
                   border = TRUE,
                   cluster_rows = FALSE,
                   show_row_names = show_row_names,
                   row_names_gp = gpar(fontsize = 5),
                   row_names_side = "left",
                   show_column_names = FALSE,
                   use_raster = TRUE,
                   top_annotation = chr_bar,
                   left_annotation = ha_row,
                   col = col_vec,
                   heatmap_legend_param = list(title = title, title_gp = gpar(fontsize = 15, fontface = "bold"), 
                                               labels_gp = gpar(fontsize = 12)))
  pdf(file = paste0("CN_",ifCN, "_nonfilter_heatmap0.pdf"), width = 10, height = 8)
  hm = draw(ht_all, heatmap_legend_side = "right")
  print(hm)
  dev.off()
}  

# func_rowAnno ---------
func_rowAnno <- function(clusterT = NA, clusterN = NA, timepoint_info, ifCF=F, ifCL =F){
  if(ifCF == T){
    #set color
    Groups_set <- c("Primary","Recurrence")
    col<-c("maroon", "#009999")
    names(col) <- Groups_set
  }
  else{
    #set color
    Groups_set <- c("Primary", "Recurrence")
    col<-c("maroon", "#009999")
    col = as.vector(col)
    names(col) <- Groups_set
  }
  
  if(ifCL == T){
    
    n_colors = length(unique(clusterT))
    row_col = c(pinks, escher, ernst)[1:n_colors]
    names(row_col) = 1:n_colors
    
    heat_row_col <- rowAnnotation(
      Cl = c(clusterT),
      Tx = timepoint_info,
      col = list(Cl = row_col,
                 Tx= col), 
      annotation_legend_param = list(Cl = list(title_gp = gpar(fontsize = 12, fontface = "bold"), 
                                               labels_gp = gpar(fontsize = 12)),
                                     Tx = list(title_gp = gpar(fontsize = 12, fontface = "bold"), 
                                               labels_gp = gpar(fontsize = 12))))
  }
  else{
    
    
    heat_row_col <- rowAnnotation(
      Tx = timepoint_info,
      col = list(Tx= col), 
      annotation_legend_param = list(
        Tx = list(title_gp = gpar(fontsize = 12,fontface = "bold"), 
                  labels_gp = gpar(fontsize = 12))))
  }
  
  return(heat_row_col)
}

# Heatmap_hclust ----------
library(circlize)
library(ComplexHeatmap)
func_htmap_clus <- function(sample, FC,
                           data, chr_lengths, row_dend, row_split, ha_row, ifCN = F,
                           show_row_names = FALSE){
  if(ifCN == T){
    breaks = c(0,1,3)
  }else{breaks = c(-1.5,0,1.5)}
  chr_lengths = chr_lengths %>% 
    dplyr::select(abspos, chrom, chrompos) %>% 
    dplyr::group_by(chrom) %>%
    dplyr::summarize(n = dplyr::n()) %>%
    pull(n)
  chr_lengths1 = chr_lengths[-24]
  chr_binary = rep(c(2,1), 11);chr_binary = c(chr_binary, 2)
  chr = data.frame(chr = rep.int(x = chr_binary, times = chr_lengths1))
  # getting lengths for chr numbers annotation
  chr_rl_c = c(1, cumsum(chr_lengths))
  # creating a data frame to calculate rowMeans
  chr_df =  data.frame(a = chr_rl_c[1:length(chr_rl_c)-1],b= chr_rl_c[2:length(chr_rl_c)])
  chr_l_means = round(rowMeans(chr_df))
  chrom.names = c(1:22,"X")
  # creating the vector for chr number annotations
  v = vector(mode = "character",length = cumsum(chr_lengths1)[length(chr_lengths1)])
  v[chr_l_means] = chrom.names
  v[is.na(v)] = ""
  chr_bar = HeatmapAnnotation(chr_text = anno_text(v, gp = gpar(fontsize = 8)),
                              df = chr,
                              show_legend = FALSE,
                              which = "column",
                              col = list(chr = c("1" = "grey88", "2" = "black")))
  col_fun = circlize::colorRamp2(breaks = breaks, c("dodgerblue3", "white", "firebrick3"))
  lgd = Legend(col_fun = col_fun, title = "logRatio")
  if (show_row_names == FALSE){ ht = Heatmap(as.matrix(log2(data)),
                                             name = "logRatio",
                                             row_title = "single cells",
                                             column_title = "genomic coordinates",
                                             column_title_side = "bottom",
                                             cluster_rows = row_dend,
                                             row_title_gp = gpar(col = c("hotpink1", "royalblue", "darkgreen","purple","pink","yellow","orange"), font = 1:3),
                                             row_split = row_split,
                                             circlize::colorRamp2(breaks = breaks, 
                                                                  c("dodgerblue3", "white", "firebrick3")),
                                             top_annotation = chr_bar, 
                                             left_annotation = ha_row,
                                             show_row_dend = TRUE, row_dend_width=unit(40,"mm"), cluster_columns = FALSE, 
                                             show_row_names = F, show_column_names = FALSE, use_raster = TRUE)}
  else { ht = Heatmap(as.matrix(log2(data)),
                      cluster_rows = row_dend,
                      name = "logRatio",
                      row_title = "single cells",
                      column_title = "genomic coordinates",
                      column_title_side = "bottom",
                      row_title_gp = gpar(col = c("hotpink1", "royalblue", "darkgreen","purple","pink","yellow","orange"), font = 1:3),
                      row_split = row_split,
                      circlize::colorRamp2(breaks = breaks, 
                                           c("dodgerblue3", "white", "firebrick3")),
                      top_annotation = chr_bar, 
                      left_annotation = ha_row,
                      show_row_dend = TRUE, row_dend_width=unit(40,"mm"), cluster_columns = FALSE, 
                      show_row_names = T, show_column_names = FALSE, use_raster = TRUE)}
  
  pdf(file = paste0("CF_NULL_CN_",ifCN, "_filter_heatmap1.pdf"), width = 10, height = 8)
  hm = draw(ht, heatmap_legend_side = "right")
  print(hm)
  dev.off()
  
}

# filter for outlier cells --------
filter_density<-function(object, resolution = 0.8, k = 5){
  seg = object$seg_cndata
  dst <- cor(seg)
  df_dist <- apply(as.matrix(dst), 1, function(x) {
    mean(sort(x, decreasing = T)[2:(k + 1)])
  }) %>% tibble::enframe(name = "sample", value = "cor")
  
  df_dist <- df_dist %>% dplyr::mutate(filtered = dplyr::case_when(cor >= resolution ~ "kept", cor < resolution ~ "removed"))
  object$df_dist = df_dist
  data<-as.data.frame(t(object$seg_cndata));dim(data);rownames(data)
  d <- amap::Dist(data, method = "manhattan", nbproc = 50)
  row_dend = hclust(d, method = "ward.D2");row_dend = color_branches(row_dend, k = 10)
  
  object$well$cells = as.character(object$well$cells)
  object$well$well = as.character(object$well$well)
  object$well = tibble(object$well)
  df_dist2 = dplyr::left_join(df_dist, object$well, by = c("sample" = "cells"))
  
  col_fil = c("kept" = "#a8e6cf", "removed" = "firebrick3")
  col_peak = c("a_well" = "grey", "d_well" = "seagreen", "pop_well" = "black")
  
  heat_row_col =  ComplexHeatmap::rowAnnotation(df = df_dist2 %>% dplyr::select(filtered, well) %>% as.data.frame(), 
                                                col = list(filtered = col_fil, well = col_peak), 
                                                show_annotation_name = T)
  
  object$seg_cndata_fil = object$seg_cndata[,which(df_dist$filtered == 'kept')]
  return(object)
}

func_umaps<-function(x, colour, title, shape, group){
  if(shape == F){
    x %>% 
      rownames_to_column(var = "Cell") %>%
      ggplot() + 
      geom_point(aes(x = umap1, 
                     y = umap2,
                     colour = colour
                     #shape = "21"
      ),
      size = 0.2,#alpha = 0.5,stroke = 0.5,
      show.legend = T) + 
      # geom_text(aes(x = umap1, 
      #               y = umap2, label=colour), size=3) +
      guides(color = guide_legend(title = title, override.aes = list(size=5))) + 
      theme(legend.position = "right",
            legend.background = element_blank(),
            legend.title = element_text(size=15),
            legend.text = element_text(size=10),
            legend.key = element_blank(),
            panel.grid = element_blank(),
            panel.border = element_rect(colour = "black", fill=NA, size=1),
            axis.line = element_blank(),
            panel.background = element_blank(),
            axis.text = element_text(size = 20),
            axis.ticks = element_blank(),
            axis.title=element_text(size=20))
  }
  else{
    x %>% 
      rownames_to_column(var = "Cell") %>%
      ggplot() +
      geom_point(aes(x = umap1, 
                     y = umap2,
                     colour = colour,
                     shape = group),
                 size = 1,
                 show.legend = T) +
      guides(color = guide_legend(title = title, override.aes = list(size=5))) + 
      theme(legend.position = "right",
            legend.background = element_blank(),
            legend.title = element_text(size=15),
            legend.text = element_text(size=10),
            legend.key = element_blank(),
            panel.grid = element_blank(),
            panel.border = element_rect(colour = "black", fill=NA, size=1),
            axis.line = element_blank(),
            panel.background = element_blank(),
            axis.text = element_text(size = 20),
            axis.ticks = element_blank(),
            axis.title=element_text(size=20))
  }
}

func_htmap <- function(ht.data, split, heat_row_col,
                      #top_annotation = readRDS(file = "/volumes/lab/users/tapsi/projects/ecis/analysis/20200610_dcis_cnv/king01_v3/chr_bar.rds"),
                      chr_top_annotation = chr_bar,
                      gene_annotation = NULL,
                      heatmap_legend_side = "right", 
                      ifconcensus = F, save_pdf = TRUE, peak = T, ifgene_anno = F){
 
  breaks = c(-1.5,0,1.5)
  dat_mat = as.matrix(t(log2(ht.data)))
  title = "Log2 (Ratio)"
  col_vec = circlize::colorRamp2(breaks = breaks, 
                                 c("dodgerblue3", "white", "firebrick3"))
  
    ht = Heatmap(t(dat_mat),
                 # cluster_rows = function(x){
                 #   fastcluster::hclust(amap::Dist(x, method = "manhattan"), method = "ward.D2")
                 # },
                 cluster_columns = FALSE,
                 border = TRUE,
                 cluster_rows = FALSE,
                 name = title,
                 row_title = "single cells",
                 column_title = "genomic coordinates",
                 column_title_side = "bottom",
                 show_row_names = F,
                 row_names_gp = gpar(fontsize = 8),
                 row_names_side = "right",
                 show_column_names = F,
                 split = split,
                 show_parent_dend_line = FALSE,
                 use_raster = TRUE,
                 top_annotation = chr_top_annotation,
                 bottom_annotation = gene_annotation, 
                 left_annotation = heat_row_col,
                 col = col_vec,
                 heatmap_legend_param = list(title = title, title_gp = gpar(fontsize = 12, fontface = "bold"), 
                                             labels_gp = gpar(fontsize = 12)))
 
  
  hm = draw(ht, heatmap_legend_side = heatmap_legend_side)
    pdf(file = paste0("clusteredhtmap_paper1.pdf"), width = 8, height = 7)
    print(hm)
    dev.off()
    
    return(hm)
 
  
  
  
}

# consensus heatmaps
func_consensus <- function(ht.data, split, heat_row_col,
                              chr_top_annotation = chr_bar,
                               anno_bar_df = df_grp_count_prop,
                               gene_annotation = NULL, 
                               heatmap_legend_side = "right", 
                               gene_anno = TRUE){
  breaks = c(-1.5,0,1.5)
  dat_mat = as.matrix(t(log2(ht.data)))
  title = "Log2 (Ratio)"
  col_vec = circlize::colorRamp2(breaks = breaks, 
                                 c("dodgerblue3", "white", "firebrick3"))
  col_tp = c("Primary" = "maroon", "Recurrence" =  "#009999")
  
  if(bar_anno == T){
    ha = HeatmapAnnotation('Cell counts' = anno_barplot(anno_bar_df, which = "row", gp=gpar(fill=col_tp,border = NA, lty = "blank")), border = F, which = "row", 
                           width = unit(1.2, "cm"),  
                           #size = unit(0.95, "mm"), 
                           #bar_width = 0.85,
                           show_legend = T)
    legend_list = list(Legend(labels = c(names(col_tp)), title = "Timepoint", legend_gp =gpar(fill=col_tp)))
  
      ht = Heatmap(dat_mat,
                   cluster_columns = FALSE,
                   border = TRUE,
                   cluster_rows = FALSE,
                   name = title,
                   row_title = "Subclones",
                   column_title = "genomic coordinates",
                   column_title_side = "bottom",
                   show_row_names = F,
                   row_names_gp = gpar(fontsize = 8),
                   row_names_side = "right",
                   show_column_names = F,
                   show_parent_dend_line = FALSE,
                   use_raster = TRUE,
                   bottom_annotation  = gene_annotation,
                   top_annotation = chr_top_annotation,
                   left_annotation = heat_row_col,
                   right_annotation = ha,
                   col = col_vec,
                   heatmap_legend_param = list(title = title, title_gp = gpar(fontsize = 12, fontface = "bold"), 
                                               labels_gp = gpar(fontsize = 12)))
    
  }
  
    hm = draw(ht, heatmap_legend_side = heatmap_legend_side, annotation_legend_list = legend_list)
    pdf(file = paste0("consensus_clusteredhtmap.pdf"), width = 10, height = 4.5)
    print(hm)
    dev.off()
    return(hm)
  
  
}

gene_anno <- function(gene_list, remove_Y = TRUE) {
  bins_in_cna_pipeline <- read.delim("bins_in_cna_pipeline_bands.txt", sep = "\t",
                                     stringsAsFactors = FALSE, 
                                     header = T)
  library(GenomicRanges)
  library(Homo.sapiens)
  hg_genes_grange <- genes(Homo.sapiens, columns = "SYMBOL")
  hg_genes_grange_df <- as.data.frame(hg_genes_grange) 
  hg_genes_grange_df$SYMBOL <- unlist(unname(as.vector(hg_genes_grange_df$SYMBOL)))
  gene_anno_df <- hg_genes_grange_df[which(hg_genes_grange_df$SYMBOL %in% gene_list),]
  names(gene_anno_df)[6] <- "gene"
  
  if (remove_Y == TRUE){
    bins_in_cna_pipeline <- bins_in_cna_pipeline %>% dplyr::filter(chr != "chrY")  
  }
  
  genes.annotation = character(length = nrow(bins_in_cna_pipeline))
  
  for (i in 1:nrow(bins_in_cna_pipeline)) {
    for (j in 1:nrow(gene_anno_df)) {
      if ((findInterval(gene_anno_df$start[j],c(bins_in_cna_pipeline$start[i], bins_in_cna_pipeline$end[i])) == 1) && (bins_in_cna_pipeline$chr[i] == gene_anno_df$seqnames[j])){
        genes.annotation[i] <- as.character(gene_anno_df$gene[j])
      }
    }
  }
  
  genes.annotation[genes.annotation == ""] <- NA
  
  labels <- genes.annotation[!is.na(genes.annotation)]
  positions <- which(!is.na(genes.annotation))
  padding = unit.c(unit(2, "cm"), unit(2, "cm"),
                   unit(c(2, 1), "cm"))
  
  gene_anno_return <- HeatmapAnnotation(link = anno_mark(at = positions, side = "bottom", labels = labels, labels_gp = gpar(fontsize = 10, col = "black")), height = unit(2, "cm"), which = "column")
  
  return(gene_anno_return)
  
} 


library("ape")
func_make_tree <- function(x,add.edge = 100, dismethod = "manhattan", tree_method = "NJ"){
 
  clusters_consensus = x$clusters_consensus
  clusters_consensus = clusters_consensus
  colnames(clusters_consensus) = c(1:ncol(clusters_consensus))
  outgroupclone = x$outgroupclone
  cluster_col = x$cluster_col
  seg_data <- clusters_consensus %>% as.data.frame()
  normal1 = data.frame(normal1 = c(rep(1, dim(clusters_consensus)[1])))
  seg_data1 = cbind(seg_data, normal1)
  colnames(seg_data1) = c(1:length(seg_data),"normal")
  dist_mat <- amap::Dist(t(seg_data1), method = dismethod, nbproc = 25)
  
  
  tree <- ape::nj(dist_mat)

  # define the root
  tree_r <- ape::root(tree, outgroup = "normal", resolve.root = TRUE) # outgroupclone
  tree_r <- ape::drop.tip(tree_r, "normal")
  tree_r$edge.length = tree_r$edge.length+add.edge
  tree_r$root.edge = add.edge
  #plot(tree_r)
  
  list_samples <- split(rownames(t(clusters_consensus)), rownames(t(clusters_consensus)))
  names(list_samples) = c(1:length(list_samples))
  tree_r <- ggtree::groupOTU(tree_r, list_samples); x$tree_r = tree_r
  
  treeplt<-ggtree::ggtree(tree_r,
                          ladderize = FALSE,
                          size = .6, root.position = 0) + 
    ggtree::geom_tippoint(aes(color = group), size = 6) +
    ggtree::geom_tiplab(size=5, color="black")+
    scale_colour_manual(values = c(cluster_col)) + 
    theme(legend.position = "right") +
    labs(color = "Clones") + geom_rootedge()
  
  print(treeplt)
  message("saving tree plot...")
  x$treeplt = treeplt
  return(x)
}
