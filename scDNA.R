#' Scripts to process scDNA data
#' These snipetts were created with the help of the following contributors: Tapsi Kumar, Darlan Minussi, Yiyun lin, Yuheui Zhao from Navinlab.

# load directories ----
# ****************** --------------------
rm(list = ls()) 
args = commandArgs(trailingOnly=TRUE)
sample.name =  "samp_1" # without Y chr
save_folder = "samp_1_v1"
wkdir = paste0("path/",save_folder,"/");dir.create(wkdir,recursive = T)
setwd(wkdir)

datadir = paste0("/path_segmented_file/")
datadir2 = paste0("/path_segmented_file/")
`%!in%` = Negate(`%in%`)
Sys.setenv(RETICULATE_PYTHON = "/path/miniconda3/bin/python")
library(reticulate)
reticulate::use_python("/path/miniconda3/bin/python", required = T)
reticulate::py_discover_config(required_module = "umap")
reticulate::py_discover_config(required_module = "leidenalg")
library(leiden)
source("/path/utils.R")
library(tidyverse)
library(Rphenograph)
library(dplyr)
library(dbscan)
library(dendextend)
library(ComplexHeatmap)
library(ggplot2)
library(cowplot)
require(RColorBrewer)
library(fastcluster)
library(ape)
library(ggtree)
#install.packages("devtools")
#devtools::install_github("johannesbjork/LaCroixColoR")
library("LaCroixColoR")

# define colors -------------------
ernst = c("#91323A", "#3A4960", "#6D7345", "#554540", "#D7C969")
escher = c("#C1395E", "#AEC17B", "#E07B42", "#89A7C2", "#F0CA50")
colr = c("#000000", "#C70E7B", "#A6E000", "#1BB6AF", "#824915", "#F25146", "#00258C", "#E6B43A", "#4E8857", ernst)#, ""
col_grp = c("Primary" = "darkorchid1", "Recurrence" = "darkorange1")
col_peak = c("a_well" = "gold1", "d_well" = "seagreen", "pop_well" = "black")


# read files --------------------
# ****************** --------------------
files = func_read_files(sample.name, datadir)
dim(files$seg_cndata) # 12167  2032
dim(files$bin_cndata)
dim(files$ratio_cndata)
files$sample.name = sample.name

# annotate data -------
# ****************** --------------------
# Get the groups ----
Groups = get_timepoint(t(files$seg_cndata));table(Groups) #get meta data groups #Primary Recurrence 226        946 
files$Groups = Groups

# Get the wells did not work----
Well = func_extract_well(t(files$seg_cndata), return_df = TRUE);table(Well) #get the wells info
files$well = Well

# convert to integer CN ----
# Get avg FACS ploidy 
ratio_pri = 85051.37/43002.36
ratio_rec = 83540.24/42276.78
ploidy_pri = 2*ratio_pri
ploidy_rec = 2*ratio_rec

# Get the cn 
intger_cn = Conv_ratio2int(x = files$seg_cndata, prim_ploidy = ploidy_pri, rec_ploidy = ploidy_rec) 
files$seg_cndata_int = intger_cn

# quick check the data 
# make log ratio plot
fun_htmap(data=files$seg_cndata, chr_lengths = files$chrom_chrompos_abspos, 
              show_row_names = F, Groups = Groups, ifCF = T, ifCN = F)


# Filter cells -----
# ****************** --------------------
# clean up bad daata (use k correlation rank)
files = filter_density(object = files, resolution = 0.85, k=8) #9

# hclust run heatmap for kept data
data<-as.data.frame(t(files$seg_cndata_fil));dim(data);rownames(data);files$data = data
data_int<-as.data.frame(t(files$seg_cndata_int))
data_int <- subset(data_int, rownames(data_int) %in% rownames(data));files$data_int = data_int

d <- amap::Dist(data, method = "euclidean", nbproc = 50)
row_dend = hclust(d, method = "ward.D2"); row_dend = color_branches(row_dend, k = 2, col = c("red", "blue", "green",
                                                                                             "yellow", "grey", "orange"))
cl = data.frame(cutree(row_dend,4));table(cl); colnames(cl) = "group" 

# Get the timepoint ----
index = get_timepoint(x=data, ifrowname = T);table(index) 
files$Groups_filt = index

#make the group named with the color
heat_row_col = func_rowAnno(timepoint_info = get_timepoint(x=data, ifrowname = F), clusterT = T, clusterN = T)

# make the heatmap
c = func_htmap_clus(data=data, chr_lengths = files$chrom_chrompos_abspos, 
                   row_dend = row_dend, ha_row = heat_row_col, row_split = 6)



# extract tumor cells --------
## extract only tumor cells based on hclust heatmap

data_tumor = data[which(cl$group %in% c("1", "2",  "4")),] # this varies for each sample
files$data_tumor = data_tumor; dim(data_tumor) 

data_tumor_int = data_int[which(cl$group %in% c("1", "2",  "4")),] 
files$data_tumor_int = data_tumor_int; dim(data_tumor_int) #761 12167

# select tumor cells ----
tum_names = rownames(files$data_tumor)
tum_names = tum_names[!tum_names %in% pop_wells]
select_tumor = c(names(files$ratio)[1:2], tum_names);length(select_tumor) #763

files_tumor_ratio = files$ratio[, select_tumor]
files_tumor_seg = files$seg[, select_tumor] 
files_tumor_bin = files$bin[, select_tumor] 

# run multiPCF to make the segments
library("copynumber")
winsorized.logratios = winsorize(files_tumor_ratio)
multi.seg <- multipcf(data = winsorized.logratios, gamma = 10, verbose=TRUE)

files_tumor_multiseg2 <- multi.seg[,6:dim(multi.seg)[2]] 
files_tumor_multiseg2 <- t(files_tumor_multiseg2)

# return the multiseg file to 12167 rows, not in log ratio
pop_expand<-function(popseg_file) {
  popseg_times <- popseg_file$n.probes
  popseg_clean <- dplyr::select(popseg_file, -c("chrom", "arm", "start.pos", "end.pos", "n.probes"))
  popseg_t <- as.data.frame(t(popseg_clean))
  popseg_long <- as.data.frame(apply(popseg_t, 1, function(m) { rep.int(m, popseg_times) }))
  popseg_long_t <- as.data.frame(t(popseg_long))
  return(popseg_long_t)
}

files_tumor_multiseg1 <- pop_expand(multi.seg)
files$seg_exp = files_tumor_multiseg1; dim(files$seg_exp) # 761 12167
files$seg = files_tumor_multiseg2;dim(files$seg) # 761 98

# umap ----- 
# ****************** --------------------
set.seed(31)
nnbs = 35
umap.data = data_tumor
umap.data.int = data_tumor_int
timepoint_info = get_timepoint(x=data_tumor, ifrowname = F)
files$umap_grp = timepoint_info
files$dtumap = data.frame(uwot::umap(umap.data,
                                     metric = "manhattan",
                                     min_dist = 0.2,
                                     n_neighbors = nnbs))
files$dtumap = files$dtumap %>% 
  dplyr::rename("umap1" = "X1",
                "umap2" = "X2") ;rownames(files$dtumap) = rownames(umap.data)

p1<-func_umaps(x=files$dtumap , colour = (get_timepoint(x=data.frame(umap.data), ifrowname=T))$G,
           title = paste0("Timepoint"), shape = F) + scale_colour_manual(values=c("lightblue", "lightpink"));p1
library(cowplot)
save_plot(filename = paste0(wkdir,"umap1.pdf"), plot = p1, base_height = 4, base_width = 5)


# findcluster ------
# ****************** --------------------
# here i replaced minor with major
cluster<-findcluster(files$dtumap, k_minor = 30, k_major = 35); files$cluster = cluster

cluster1 <- findcluster_dbscan(x = files$dtumap, k_minor = NULL, k_major = NULL); 
cluster1$subclones = factor(cluster1$subclones, levels = c("1","2","3","4","5", "6", "7", "8", "9", "10", "11")) # define levels for each sample

files$clusterdb = cluster1
p3 <- func_umaps(x=files$dtumap, colour = as.factor(cluster$major ), 
             title = paste0("subclones"), shape = F) +
  labs(x ="umap1",
       y ="umap2")+ 
  scale_color_manual(values = c(pinks, escher, ernst));p3

p3db <- func_umaps(x=files$dtumap, colour = as.factor(cluster1$subclones), 
               title = paste0("subclones"), shape = F) +
  labs(x ="umap1",
       y ="umap2")+ 
  scale_color_manual(values = c(pinks, escher, ernst ));p3db

cluster = cluster %>% mutate(minor_grp = minor)
cluster = cluster %>% mutate(major_grp = major)
cluster = cluster %>% mutate(subclones = cluster1$subclones)
# manually added later
tree_tips_order
cluster$subclones2 = factor(cluster$subclones, levels = tree_tips_order) %>% as.numeric()

cluster$subclones1 = factor(cluster$subclones, levels = seq(1:length(unique(cluster$subclones))))
cluster$subclones = factor(cluster$subclones2, levels = seq(1:length(unique(cluster$subclones2))))

cluster_or <- cluster[order(cluster$major_grp),];files$cluster_or = cluster_or
p3minor <- func_umaps(x=files$dtumap, colour = as.factor(cluster$minor_grp ), 
                  title = paste0("minor_clones"), shape = F) +
  labs(x ="umap1",
       y ="umap2")+ 
  scale_color_manual(values = c(pinks, escher, ernst))
p3 <- func_umaps(x=files$dtumap, colour = as.factor(cluster$major_grp ), 
             title = paste0("major_clones"), shape = F) +
  labs(x ="umap1",
       y ="umap2")+ 
  scale_color_manual(values = c(pinks,escher, ernst))
p3_subclones <- func_umaps(x=files$dtumap, colour = as.factor(cluster$subclones1), 
                       title = paste0("Subclones"), shape = F) +
  labs(x ="umap1",
       y ="umap2")+ 
  scale_color_manual(values = c(pinks, escher, ernst))


p4 = plot_grid(p1,p3minor, p3, p3_subclones, rel_widths = c(5,5, 5,5),rel_heights = 4,ncol =2)
pdf(file = paste0(wkdir,sample.name,"_UMAP1.pdf"), width = 10, height = 7)
p4
dev.off()

p1<-func_umaps(x=files$dtumap , colour = (get_timepoint(x=data.frame(umap.data), ifrowname=T))$G,
           title = paste0("Timepoint"), shape = F) + scale_colour_manual(values=c("lightblue", "lightpink"))

fig_final = plot_grid(p1, p3_subclones, rel_widths = c(5,5),rel_heights = 4, ncol =2)
pdf(file = paste0(wkdir,sample.name,"_UMAP1_final_figure.pdf"), width = 10, height = 3.5)
fig_final
dev.off()

# clustered heatmap ----------
files$outgroupclone = length(unique(cluster_or$subclones1))
ht.data1 = umap.data.int[cluster_or$cell,]
clusterT = cluster_or$subclones1; files$clusterT =clusterT
ht.cluster = clusterT
index = get_timepoint(x=ht.data1, ifrowname=F)
heat_row_col = func_rowAnno(timepoint_info = index, clusterT = clusterT, clusterN = NA ,ifCF=T, ifCL = T)
func_htmap(ht.data = ht.data1, split = ht.cluster, heat_row_col = heat_row_col,  ifconcensus = F,ifCN =T, save_pdf = TRUE)

# make consensus ---------
# ****************** --------------------
cluster_col = colr[1:filescbs$outgroupclone]
names(cluster_col) = c(1:filescbs$outgroupclone)

filescbs$cluster_col= cluster_col 

df_grp_count = tibble(Timepoint = row_anno_df$Timepoint, cell = rownames(ht.datacbs), consensus_cluster = ht.clustercbs)
df_grp_count_prop = table(df_grp_count$consensus_cluster, df_grp_count$Timepoint)
df_grp_count_prop_fill = prop.table(df_grp_count_prop, margin = 1)

# ratio
clusters_list = base::split(ht.datacbs, ht.clustercbs)
clusters_consensus = purrr::map_dfr(clusters_list, function(x) apply(x, 2, median))
clusters_consensus = as.data.frame(t(clusters_consensus)); dim(clusters_consensus); 
filescbs$clusters_consensus = clusters_consensus 
clusters_consensus[is.na(clusters_consensus)] = 1
clusters_consensus1 = t(clusters_consensus)

# set colors of sublcones
n_colors = ncol(clusters_consensus)
row_col = colr[1:n_colors]
names(row_col) = 1:n_colors

consensus_row_anno = rowAnnotation(Subclones = as.factor(unique(sort(ht.clustercbs))), #ht.clustercbs
                                   col = list(Subclones = row_col),
                                   annotation_legend_param = list(
                                     Subclones = list(title_gp = gpar(fontsize = 12,fontface = "bold"), 
                                                      labels_gp = gpar(fontsize = 12))))

func_consensus(ht.data = clusters_consensus, split = rownames(clusters_consensus), heat_row_col=consensus_row_anno,
               chr_top_annotation = chr_bar,
               gene_annotation = NULL, ifCN = F,
               anno_bar_df = df_grp_count_prop_fill,
               heatmap_legend_side = "right", ifhclust = FALSE,
               ifconcensus = T, save_pdf = T, peak = FALSE, bar_anno = TRUE)

# gene annotation -------
files1cbs = filescbs
genes_list1 = c("FHIT", "AKT3",   "EGFR",  "FGFR1", "MYC", 
                "CDKN2A", "GATA3",   "PGR", "CDK4",  "BRCA2", "RB1", "TP53", "BRCA1", "CCNE1","PPM1D","NCOA3", "AURKA", "ERBB2")


gene_ht_anno = gene_anno(genes_list1)
# consensus heatmaps
func_consensus(ht.data = clusters_consensus, split = rownames(clusters_consensus), heat_row_col=consensus_row_anno,
               chr_top_annotation = chr_bar,
               gene_annotation = gene_ht_anno, ifCN = F,
               anno_bar_df = df_grp_count_prop_fill,
               heatmap_legend_side = "right", ifhclust = FALSE,
               ifconcensus = T, save_pdf = T, peak = FALSE, bar_anno = TRUE)


# make phylogenetic tree by consensus ratio -------
# ****************** --------------------

files = func_make_tree(files, dismethod = "manhattan", add.edge = 100, tree_method = "NJ")

# NJ tree
pdf(file = paste0(wkdir,sample.name, dismethod,"_consensus_NJ.pdf"), width = 5, height = 5)
files$treeplt
dev.off()






# get order from tree ----
is_tip <- tree_r$edge[,2] <= length(tree_r$tip.label)
ordered_tips_index <- tree_r$edge[is_tip, 2]
tree_tips_order <- tree_r$tip.label[ordered_tips_index] %>% rev()

#files2pcf = dataReorderpcf(files)
clusterT = files2pcf$cl_orpcf$subclones
ht.cluster3 = c(clusterT) #,clusterN

rownames(files2pcf$cluster_or_pcf) = files2pcf$cluster_or_pcf$cell
index = get_timepoint(x=files2pcf$cluster_or_pcf, ifrowname=T)
index$cell = rownames(index)
timepoints=levels(index$G)


# cal fraction
files2pcf$Groups2 = data_frame(Timepoint = files2pcf$Groups$G, cell = rownames(files2pcf$Groups))
fq_tumor <- data.frame(left_join(files2pcf$cluster_or_pcf, files2pcf$Groups2, by = "cell"))

fq_tumor_pre <- fq_tumor[which(fq_tumor$Timepoint == "Primary"),];
fq_tumor_post <- fq_tumor[which(fq_tumor$Timepoint == "Recurrence"),]  

frac_pre<-data.frame(table(fq_tumor_pre$subclones)) #subclones3
frac_pre$Freq = frac_pre$Freq*100/nrow(fq_tumor_pre);frac_pre$Var1 = as.numeric(as.character(frac_pre$Var1))
frac_pre$Timepoint <- "Primary"
colnames(frac_pre)<-c("clone_id","clonal_prev","timepoint")

frac_post<-data.frame(table(fq_tumor_post$subclones))
frac_post$Freq = frac_post$Freq*100/nrow(fq_tumor_post);frac_post$Var1 = as.numeric(as.character(frac_post$Var1))
frac_post$Timepoint  = "Recurrence"
colnames(frac_post)<-c("clone_id","clonal_prev","timepoint")

clonal_prev <- rbind(frac_pre,frac_post) #this is the summary we need for clonal prevalence
clonal_prev$clone_id = as.character(clonal_prev$clone_id)
clonal_prev = tibble(clonal_prev)
clonal_prev$timepoint = as.character(clonal_prev$timepoint)

# get tree structure
library(trqwe)
tree_edges <- tree_r$edge %>% set_colnames(c("source","target")) %>% data.frame()
# get color for clone
tip_colours	<- cbind(names(files2pcf$cluster_col),files2pcf$cluster_col)
# create color for fake nodes
shadesOfGrey <- colorRampPalette(c("grey65", "grey87"))
alltree = c(tree_edges$source,tree_edges$target)
alltree = unique(alltree)
fake.node = alltree[alltree %!in% tip_colours[,1]]
fn_colours = cbind(fake.node, shadesOfGrey(length(fake.node)))
clone_colours <- rbind(tip_colours,fn_colours) %>% data.frame() %>% set_colnames(c("clone_id","colour"))
clone_colours$colour = as.character(clone_colours$colour)

Fishplot1 <-timescape(clonal_prev, tree_edges, clone_colours, mutations = "NA",  
                      xaxis_title = "Time Point", yaxis_title = "Clonal Prevalence",
                      phylogeny_title = "Clonal Phylogeny", alpha = 0.5,
                      genotype_position = "stack", perturbations = "NA", sort = F,
                      show_warnings = TRUE, width = 800, height = 300)
Fishplot1
saveWidget(Fishplot1, file = "tree_timescape_njtree.html", background = "transparent")

# ****************** DONE ****************** -------
