
#' dimred_plot
#'
#' @description
#' dimensional reduction plot
#'
#' @param df : a data frame including embedding and variables to plot 
#' @param ed.names: a length 2 tuple of embedding names 
#' @param group.by: variable used to group cells
#' @param pt.size: point size, default is 0.3
#' @param centers: a data frame specifying center location of each group
#' @param legend.size: legend label size, default is 3
#' @param legend.ncol: number of columns of legend, default is 1
#' @param label: whether to show labels, default is TRUE
#' @param palette: color theme jama is used or not
#' @param label.size: label font size
#' @return ggplot2 object
#' 
#' 
dimred_plot = function(
  df, ed.names, group.by, pt.size=0.3, centers=NULL, legend.size=3, legend.ncol=1, label=T, palette=F, label.size=4) 
{
  eb.u=ceiling(max(df[,ed.names[1]], df[,ed.names[2]])); eb.l=floor(min(df[,ed.names[1]], df[,ed.names[2]]))
  p = ggplot(df, aes_string(x=ed.names[1], y=ed.names[2])) + geom_point(aes_string(color=group.by), size=pt.size)
  p = p + coord_fixed() + xlim(eb.l, eb.u)+ylim(eb.l,eb.u)
  p = p + theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
  p = p + guides(color = guide_legend(override.aes = list(size = 3), ncol=legend.ncol)) 
  # p = p + labs(color = "Clusters")
  if (label)
  p = p + annotate("text", x=centers$tsne1.c, y=centers$tsne2.c, label=names(table(centers[,1])), check_overlap=T, size=label.size)
  if (palette)
  p = p + scale_color_jama()
  return(p)
}


#' feature_plot
#'
#' @description
#' Colors single cells on a dimensional reduction plot according to expression of a specific gene
#'
#' @param df : a data frame including embedding and variables to plot
#' @param gname : a gene's names to plot 
#' @param ed.names: a length 2 tuple of embedding names 
#' @param pt.size: point size, default is 0.3
#' @param cutoff: max cutoff, default is 100 or no cutoff
#' @return ggplot2 object
#' 
#' 
feature_plot = function(df, gname, ed.names, pt.size=0.4, cutoff=100){
  # cut color spectral based on log counts after excluding extreme big values
  ind=df[,gname]<=quantile(df2[,gname],cutoff/100)
  colors <- colorRampPalette(c("gray75", "purple3"))(100)
  plotcol <- colors[cut(log1p(df[ind, gname]), breaks=100)]
  
  eb.u=ceiling(max(df[,ed.names[1]], df[,ed.names[2]])); eb.l=floor(min(df[,ed.names[1]], df[,ed.names[2]]))
  
  p = ggplot(df[ind,], aes_string(x=ed.names[1], y=ed.names[2])) + geom_point(color=plotcol, size=pt.size) 
  p = p + coord_fixed() + xlim(eb.l,eb.u) + ylim(eb.l,eb.u)  
  p = p + theme_bw() +theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())+ggtitle(gname)
  p = p + scale_color_gradient(low="gray90", high="purple3")
  return (p)
  # ggsave(paste0("/Users/kangbowei/Dropbox/SingleCell/figures/figS3/A_new/",gname,".png"))
}


#' stacked_vln_plot
#'
#' @description
#' stacked violin plot of marker gene for each cluster
#'
#' @param obj: Seurat object
#' @param features: tuple of gene marker names 
#' @param pt.size: point size, default is 0
#' @param plot.margin: plot margin in ggplot2
#' @return plot
#' 
#'
stacked_vln_plot<- function(obj, features, pt.size = 0, plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"), ...) {
  plot_list <- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x, group.by='idnew', col=color_seurat$cols[ordering],...))
  plot_list[[1]]<- plot_list[[1]] + theme(axis.text.x=element_text(), axis.ticks.x = element_line()) + scale_x_discrete(position = "top")
  # change the y-axis tick to only max value 
  ymaxs<- purrr::map_dbl(plot_list, extract_max)
  plot_list<- purrr::map2(plot_list, ymaxs, function(x,y) x + scale_y_continuous(breaks = c(y)) + expand_limits(y = y))
  p <- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
  return (p)
}


#' helper function used in stacked_vln_plot
#' plot.margin to adjust the white space between each plot. pass any arguments to VlnPlot in Seurat
modify_vlnplot<- function(obj, feature, pt.size = 0, plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"), ...) {
  p <- VlnPlot(obj, features = feature, pt.size = pt.size, ... ) + xlab("") + ylab(feature) + ggtitle("") 
  p <- p + theme(legend.position = "none", axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.y = element_text(size = rel(1), angle = 0), axis.text.y = element_blank(), plot.margin = plot.margin ) 
  p <- p + scale_x_discrete(position = "top")
  return(p)
}

#' helper function used in stacked_vln_plot
#' extract the max value of the y axis
extract_max<- function(p){
  ymax<- max(ggplot_build(p)$layout$panel_scales_y[[1]]$range$range)
  return(ceiling(ymax))
}



#' pie_plot
#'
#' @description
#' pie chart showing distribution by days
#'
#' @param df: data frame including column "day"
#' @return pie chart
#' 
#'
pie_plot = function(df) {
  df_pie = data.frame(table(df$day))
  colnames(df_pie) = c('Day','Freq')
  df_pie <- df_pie %>% arrange(desc(Day)) %>% mutate(prop = Freq / sum(df_pie$Freq) *100) 
  df_pie <- df_pie %>% mutate(ypos = cumsum(prop)- 0.5*prop )
  p = ggplot(df_pie, aes(x="", y=Freq, fill=Day)) + geom_bar(stat="identity", width=1, color="white") 
  p = p + coord_polar("y", start=0) + theme_void() + scale_fill_manual(values=color_stages$cols) 
  p = p + geom_text(aes(label = paste0(Freq,', ',round(prop),'%')), size=3, position = position_stack(vjust = 0.5))
  return (p)
}


#' bar_plot
#'
#' @description
#' bar chart showing distribution by days and clusters
#'
#' @param df: data frame including column "day" and "cluster"
#' @param position: "stack" or "fill". default is stack
#' @return bar chart
#' 
#'
bar_plot = function(df, position='stack') {
  df_bar = df[,c('cluster', 'day')] %>% group_by(cluster,day) %>% summarise(freq = n())
  p = ggplot(df_bar, aes(fill=day, y=freq, x=cluster)) 
  p = p + geom_bar(stat="identity", position=position) + labs(fill='Day')
  p = p + theme(panel.background=element_blank(),axis.title = element_blank())+ scale_fill_manual(values=color_stages$cols)
  return (p)
  # ggsave('/Users/kangbowei/Dropbox/SingleCell/figures/figS1/B_new/bar_days_2.png')
}


#' two_way_summ_tbl
#'
#' @description
#' a summary table to show distribution across cells by cluster and day
#'
#' @param df : a data frame including "cluster" and "day" 
#' @param subset: a subset of cluster labels
#' @return a list of tables: counts and percentage
#' 
#' 
two_way_summ_tbl = function(df, subset) {
  tbl1 = table(df$cluster,df$day)
  tbl1 = tbl1[subset+1,]
  tbl1 = cbind(tbl1, apply(tbl1,1,sum))
  tbl1 = rbind(tbl1, apply(tbl1,2,sum))
  return (list(tbl1, round(t(apply(tbl1,1,function(t) t/sum(t[-7]))),3)))
}



#' line_plot
#'
#' @description
#' lines plot showing gene expression changes along developmental stages
#'
#' @param dat: Seurat object with stages named "orig.ident"
#' @param gene_list: list of genes to plot
#' @return line chart
#' 
#'
line_plot = function(dat, gene_list) {
  df = data.frame(day=dat@meta.data$orig.ident)
  for (i in 1:length(gene_list)) df = cbind(df, dat@assays$RNA@data[gene_list[i],])
  colnames(df)[-1]=gene_list
  df = df %>% gather('gene','expr',2:(length(gene_list)+1)) %>% group_by(gene,day) %>% summarise(mean_expr=mean(expr))
  p = ggplot(df, aes(x=day, y=mean_expr, group=gene, col=gene))
  p = p + geom_line(aes(linetype=gene))+ geom_point() + theme_bw()
  return (p)
}


#' vln_plot_compare
#'
#' @description
#' violin plot comparison (E12 vs E14)
#'
#' @param dat: Seurat object with stages named "orig.ident"
#' @param gname: gene name
#' @return violin chart
#' 
#'
vln_plot_compare = function(dat, gname) {
  df = data.frame(dat@assays$RNA@data[gname,], day=dat@meta.data$orig.ident)
  colnames(df)[1] = 'gene'
  p = ggplot(df.tmp,aes(x=day, y=gene)) + geom_violin(trim=T,aes(colour=day,fill=day),scale='area') 
  p = p + geom_jitter(shape=16, position=position_jitter(0.2),size=0.9, color='gray30', alpha=0.5)
  p = p + theme_bw() + ggtitle(gname)+ theme(axis.title = element_blank(),panel.grid=element_blank(),legend.position="bottom")
  p = p + scale_fill_manual(values=color_stages$cols[2:3]) + scale_color_manual(values=color_stages$cols[2:3])
  return (p)
}


#' find_marker_12vs14
#'
#' @description
#' find marker genes E12 vs E14 for a specific cluster
#' require that the log fold change to be greater than 0.6
#'
#' @param dat: Seurat object 
#' @param cc: cluster to consider
#' @return table of marker genes ordering by avr_logFC
#' 
#'
find_marker_12vs14 = function(dat, cc) {
  dat.sub = subset(dat, subset=X19clusters==cc)
  cells.e12=data.frame(dat.sub[['orig.ident']],names=colnames(dat.sub)) %>% dplyr::filter(orig.ident=='E12') 
  cells.e14=data.frame(dat.sub[['orig.ident']],names=colnames(dat.sub)) %>% dplyr::filter(orig.ident=='E14') 
  dat.sub@active.ident = factor(dat.sub@meta.data$orig.ident)
  E12vs14.markers <- FindMarkers(dat.sub, only.pos=F,  ident.1=cells.e12$names, ident.2=cells.e14$names, min.pct = 0.25)
  E12vs14.markers.select=E12vs14.markers %>% dplyr::filter(abs(avg_logFC)>0.6) %>% arrange(avg_logFC)
  return(E12vs14.markers.select)
}


#' find_marker_1vs3
#'
#' @description
#' find marker genes c1 vs c3 (different cell cycle) for a specific cluster
#' require that the log fold change to be greater than 1
#'
#' @param dat: Seurat object 
#' @param day: stage to consider
#' @return table of marker genes ordering by avr_logFC
#' 
#'
find_marker_1vs3 = function(dat, day) {
  dat.sub = subset(dat, subset=orig.ident==day)
  cells.c1=data.frame(dat.sub[['orig.ident']],names=colnames(dat.sub)) %>% dplyr::filter(seurat_clusters==1) 
  cells.c3=data.frame(dat.sub[['orig.ident']],names=colnames(dat.sub)) %>% dplyr::filter(seurat_clusters==3) 
  dat.sub@active.ident = factor(dat.sub@meta.data$seurat_clusters)
  c1vsc3.markers <- FindMarkers(dat.sub, only.pos=F,  ident.1=cells.c1$names, ident.2=cells.c3$names, min.pct = 0.25)
  c1vsc3.markers.select = c1vsc3.markers %>% dplyr::filter(abs(avg_logFC)>1) %>% arrange(avg_logFC)
  return(c1vsc3.markers.select)
}


#' get_Eomes_plus_subset
#'
#' @description
#' get Eomes+ subset 
#'
#' @param dat: Seurat object for full raw data 
#' @return Seurat object for the Eomes+ subset (after QC, PCA, clustering, tsne)
#' 
#'
get_Eomes_plus_subset = function(dat) {
  # sub-clustering c0+4 to get Eomes+ subset
  dat.c04 = subset(dat, subset=X19clusters%in%c(0,4))
  dat.c04 <- NormalizeData(dat.c04, normalization.method = "LogNormalize", scale.factor = 10000)
  dat.c04 <- NormalizeData(dat.c04)
  dat.c04 <- FindVariableFeatures(dat.c04, selection.method = "vst", nfeatures = 2000)
  dat.c04 <- ScaleData(dat.c04, features = rownames(dat.c04))
  dat.c04 <- RunPCA(dat.c04, features = VariableFeatures(object = dat.c04))
  dat.c04 <- FindNeighbors(dat.c04, dims = 1:10)
  dat.c04 <- FindClusters(dat.c04, resolution = 2)
  dat.c04 <- RunTSNE(dat.c04, dims = 1:10)
  
  # select subset by comparing dimplot by new clusters and Eomes feature plot
  # DimPlot(dat.c04, reduction = "tsne",label=T)
  # FeaturePlot(dat.c04,'Eomes',reduction = 'tsne')
  cells.sub1 = colnames(subset(dat.c04, subset=seurat_clusters%in%c(0,5,4,9,12,11,16,15,14))) # 1,255 cells

  # sub-clustering c13 to get Eomes+ subset
  dat.c13 = subset(dat, X19clusters==13)
  dat.c13 <- NormalizeData(dat.c13, normalization.method = "LogNormalize", scale.factor = 10000)
  dat.c13 <- FindVariableFeatures(dat.c13, selection.method = "vst", nfeatures = 2000)
  all.genes <- rownames(dat.c13)
  dat.c13 <- ScaleData(dat.c13, features = all.genes)
  dat.c13 <- RunPCA(dat.c13, features = VariableFeatures(object = dat.c13))
  dat.c13 <- FindNeighbors(dat.c13, dims = 1:20)
  dat.c13 <- FindClusters(dat.c13, resolution =0.8)
  dat.c13 <- RunTSNE(dat.c13, dims = 1:10)
  
  # select subset by comparing dimplot by new clusters and Eomes feature plot
  # DimPlot(dat.c13, reduction = "tsne",label=T)
  # FeaturePlot(dat.c13, 'Eomes', reduction='tsne')
  cells.sub2 = colnames(subset(dat.c13, subset=seurat_clusters == 1)) # 74 cells
  
  # c4 + 7 + 0(part) + 13(part)
  dat[['cc04713_flag']] = ifelse(colnames(dat)%in%c(cells.sub1, cells.sub2) | dat@meta.data$X19clusters%in%c(4,7), 0, 1)
  dat.c04713 = subset(dat, subset=cc04713_flag==0) # 1,849 cells
  dat.c04713 <- NormalizeData(dat.c04713, normalization.method = "LogNormalize", scale.factor = 10000)
  dat.c04713 <- FindVariableFeatures(dat.c04713, selection.method = "vst", nfeatures = 2000)
  dat.c04713 <- ScaleData(dat.c04713, features = rownames(dat.c04713))
  dat.c04713 <- RunPCA(dat.c04713, features = VariableFeatures(object = dat.c04713))
  dat.c04713 <- FindNeighbors(dat.c04713, dims = 1:20)
  dat.c04713 <- FindClusters(dat.c04713, resolution =.8)
  dat.c04713 <- RunTSNE(dat.c04713, dims = 1:10)
  
  return(dat.c04713)
}


#' bubble_plot
#'
#' @description
#' bubble plot for a list of selected genes 
#'
#' @param df: data frame including columns of "cluster" and "day"
#' @param dat: Seurat object including scaled gene expressions 
#' @param genes: list of genes to plot
#' @param clusters: clusters to plot
#' @param lab.display: cluster labels to display (from bottom to top)
#' @return ggplot2 object 
#' 
#'
bubble_plot = function(df, dat, genes, clusters, lab.display) {
  for (i in 1:length(genes)){
    df[, genes[i]] = dat@assays$RNA@scale.data[genes[i],]
  }
  df = df %>% dplyr::filter(cluster %in% clusters)
  df = df %>% gather('gene','value',2:(length(genes)+1)) 
  df = df %>% group_by(gene, cluster) %>% dplyr::summarise(expr=mean(value))
  df = df %>% mutate(cluster=factor(cluster, labels=paste0('c',clusters))) %>% mutate(cluster=factor(cluster, levels=lab.display))
  df = df %>% mutate(gene=factor(gene, levels=genes))
  p = ggplot(df) + geom_point(aes(x=gene,y=cluster,size=expr), color='steelblue')
  p = p + theme_bw() + theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, size=12), axis.text.y=element_text(size=12), legend.position="top", legend.text=element_text(size=12)) 
  p = p + ylab('') +xlab("") + scale_size(range = c(-1, 12))
  return(p)
}


#' scatter_plot
#'
#' @description
#' scatter plot to show gene co-expressions
#'
#' @param df: data frame including columns of "cluster", "day", "cols" and selected gene expression
#' @param genes: list of genes to compare
#' @param clus: clusters to plot
#' @return Seurat object for the Eomes+ subset (after QC, PCA, clustering, tsne)
#' 
#'
scatter_plot = function(df, clus, genes) {
  df.sub = df %>% dplyr::filter(cluster==clus)
  Eomes.sub = df.sub$Eomes
  idx = which(names(df.sub) %in% genes)
  df.sub = df.sub %>% select(c(1:2, idx, ncol(df.sub))) 
  df.sub = df.sub %>% gather('gene', 'expr', 3:(ncol(df.sub)-1))
  df.sub$Eomes = rep(Eomes.sub, length(idx))
  p = ggplot(df.sub) + geom_point(aes(x=Eomes, y=expr, col=day)) + facet_grid(cols=vars(gene)) 
  p = p + theme_bw() + theme(panel.grid=element_blank(), legend.position='top') + scale_color_manual(values=col_vec)
  p = p + ggtitle(paste0('c', clus)) 
  return(p)
}


#' cc_scores
#'
#' @description
#' calculate cell-cycle scores
#' reference:
#' https://satijalab.org/seurat/v3.1/cell_cycle_vignette.html 
#' https://hbctraining.github.io/scRNA-seq/lessons/cell_cycle_scoring.html
#'
#' @param dat: Seurat object with ccg regressed out
#' @param ref: reference file cell-cycle genes
#' @return Seurat object with 3 new meta.data columns added: S.Score, G2M.Score, Phase
#' 
#'
cc_scores = function(dat, ref){
  dat@active.assay = 'RNA'
  # sum(is.na(dat2[['pdt']]$pdt)) # 451
  dat[['nan_pdt']] = is.na(dat[['pdt']]$pdt)
  dat = subset(dat, subset=nan_pdt==FALSE) # 9115 cells
  dat[['SCT']] <- NULL
  
  # reference genes
  cell_cycle_genes = read.csv(ref, header = T)
  # Connect to AnnotationHub
  ah <- AnnotationHub()
  # Access the Ensembl database for organism
  ahDb <- query(ah, pattern = c("Homo sapiens", "EnsDb"), ignore.case = TRUE)
  # Acquire the latest annotation files
  id <- ahDb %>% mcols() %>% rownames() %>% tail(n = 1)
  # Download the appropriate Ensembldb database
  edb <- ah[[id]]
  # Extract gene-level information from database
  annotations <- genes(edb, return.type = "data.frame")
  # Select annotations of interest
  annotations <- annotations %>% dplyr::select(gene_id, gene_name, seq_name, gene_biotype, description)
  # Get gene names for Ensembl IDs for each gene
  cell_cycle_markers <- dplyr::left_join(cell_cycle_genes, annotations, by = c("geneID" = "gene_id"))
  # Acquire the S and G2M phase genes
  s_genes <- cell_cycle_markers %>% dplyr::filter(phase == "S") %>% pull("gene_name")
  g2m_genes <- cell_cycle_markers %>% dplyr::filter(phase == "G2/M") %>% pull("gene_name")
  
  # Perform cell cycle scoring
  dat <- CellCycleScoring(dat, g2m.features = g2m_genes, s.features = s_genes)
  # dat[[]] %>% dplyr::select(S.Score, G2M.Score, Phase)
  return (dat)
}



#' heat_maps
#'
#' @description
#' plot heat map for selected genes
#' 
#'
#' @param df: data frame including 4 columns: cell_id, rank_pdt, S.Score, G2M.Score
#' @param dat: Seurat object including scaled gene expressions in Assay RNA
#' @param genes: genes to plot
#' @return two ggplot2 objects: gene expression vs pseudotime, clustering vs pseudotime
#' 
#'
heat_maps = function(df, dat, genes){
  # top panel: gene expression vs pseudotime
  # scale cell_cycle scores 
  mmax = apply(df[,3:4], 2, max) 
  mmin = apply(df[,3:4], 2, min)
  df = df %>% mutate(S.Score=(S.Score-mmin[1])/(mmax[1]-mmin[1]), G2M.Score=(G2M.Score-mmin[2])/(mmax[2]-mmin[2]))
  colnames(df)[3:4] = c("orig_S.Score","orig_G2M.Score")
  
  # add gene expressions and scale
  for (i in 1: length(genes)) {
    df = cbind(df, dat@assays$RNA@data[genes[i],]/max(dat@assays$RNA@data[genes[i],]))
  }
  colnames(df)[-c(1:4)] = paste(genes, "orig", sep="_")
  df = df %>% arrange(rank_pdt)
  
  # smoothing
  for (i in 1:(2+length(genes))) {
    df = cbind(df, ma(df[,2+i]))
  }
  
  colnames(df)[(length(df)-length(genes)-1):length(df)] = c('S.Score', 'G2M.Score', genes)
  df = df %>% gather('gene','scale_expr', (length(df)-length(genes)-1):length(df)) %>% dplyr::filter(!is.na(scale_expr))
  df$gene = factor(df$gene, levels=c(genes, 'S.Score', 'G2M.Score'))
  
  # plot
  p1 = ggplot(df, aes(rank_pdt, gene, fill= scale_expr)) + geom_tile() 
  p1 = p1 + theme_void() + theme(axis.text.x=element_text(face="plain", color="black", size=7), axis.text.y=element_text(face="plain", color="black", size=7), legend.title=element_text(size=7),legend.text=element_text(size=7))
  p1 = p1 + scale_fill_viridis()
  # +scale_fill_gradient(low='steelblue',high='yellow')
  
  # bottom panel: clustering vs pseudotime
  df.oc = data.frame(cell_id=colnames(dat), rank_pdt=dat@meta.data$rank_pdt)
  df.oc$c7  = dat@meta.data$X19clusters == 7
  df.oc$c4  = dat@meta.data$X19clusters == 4
  df.oc$c0  = dat@meta.data$X19clusters == 0
  df.oc$c13 = dat@meta.data$X19clusters ==13
  df.oc = df.oc %>% arrange(rank_pdt) %>% gather('cluster','flag', 3:6) 
  df.oc$cluster = factor(df.oc$cluster, levels=c('c7','c4','c0','c13'))
  p2 = ggplot(df.oc, aes(rank_pdt, cluster, fill= flag)) + geom_tile() 
  p2 = p2 + theme_void() + theme(axis.text.x=element_text(face="plain", color="black", size=7), axis.text.y=element_text(face="plain", color="black", size=7), legend.title=element_text(size=7),legend.text=element_text(size=7))
  p2 = p2 + scale_fill_viridis(discrete = T, option='E')
  
  return(list(p1, p2))
}

#' helper function for heat_maps
#' smoothing by two-sided moving average
#' @param x: 1-dim times series to work on
#' @param n: window width. Default is 9
#' @return 1-dim times series (may have smaller length than original series)
#' 
ma <- function(x, n = 9){
  as.numeric(stats::filter(x, rep(1 / n, n), sides = 2))
}

# END #
