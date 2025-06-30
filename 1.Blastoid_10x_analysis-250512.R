# Title: Blastoid Blastocyst Analysis —— Blastocyst 10x data analysis
# Author: Gaozy
# Time: 2025-04-30

# comput172-r_env-R
# 0. Basic settings ----
setwd("/home/qcao02/gaozy/blastocyst/")

.libPaths(c("/usr/lib64/R/library",
            "/usr/share/R/library",
            "/home/lushi02/.local/share/r-miniconda/envs/r-reticulate",
            "/home/lushi02/anaconda3/lib/R/library"))

# 1. library ----
library(Matrix, lib.loc = "/usr/lib64/R/library") 
library(Seurat, lib.loc = "/usr/lib64/R/library") 
library(ggplot2)
library(ggpubr)
library(ComplexHeatmap)
library(RColorBrewer)
library(ggcorrplot, lib.loc = "/home/qcao02/R/x86_64-conda-linux-gnu-library/4.4")
library(Hmisc)
library(dplyr, lib.loc = "/usr/lib64/R/library")
library(reshape2)
library(ggrepel)
library(patchwork)  # For combining plots

if (dir.exists("250512") == F) dir.create('250512')
setwd("250512")

on.exit({
  while(dev.cur() > 1) dev.off()
})

# 2. 分析需求1：t-SNE 分析 in house blastoid-10X 和public blastocyst 10X ----
## 1) load data ----
merge_data <- readRDS("../250430/raw_merge_data_GSM4026212.rds")
integrated_obj <- readRDS("../250430/raw_integrated_data_GSM4026212.rds")

## 2) PCA - inte ----
DefaultAssay(integrated_obj) <- "integrated"
Idents(merge_data) <- "cell_info"

RhpcBLASctl::blas_set_num_threads(8)
# Variable features
### b. pca ----
integrated_obj <- ScaleData(integrated_obj, features = VariableFeatures(integrated_obj), do.center = TRUE) %>%
  RunPCA(npcs = 30, features = VariableFeatures(integrated_obj))

## 3) UMAP & tSNE - inte ----
library(clustree)
library(ggplot2)

### a. compare diff resolutions ----
integrated_obj <- FindNeighbors(integrated_obj, dims = 1:30)

### b. select a resolution ----
integrated_obj <- FindClusters(integrated_obj, resolution = 1.2)
integrated_obj <- RunTSNE(integrated_obj, dims = 1:30 )
integrated_obj <- RunUMAP(integrated_obj, dims = 1:30 )

p1 <- DimPlot(integrated_obj, reduction = "umap", group.by = "cell_info", label = TRUE, repel = TRUE, pt.size = 2, cols = c("#1f77b4", "#ff7f0e"))+ guides(color = guide_legend(ncol = 1))
p2 <- DimPlot(integrated_obj, reduction = "umap", group.by = "integrated_snn_res.1.2", label = TRUE, repel = TRUE, pt.size = 2)+ guides(color = guide_legend(ncol = 1))

pdf( "umap_inte.pdf" , width=12, height=5, onefile = F)
print( p1 | p2)
dev.off()

p3 <- DimPlot(integrated_obj, reduction = "tsne", group.by = "cell_info", label = TRUE, repel = TRUE, pt.size = 2, cols = c("#1f77b4", "#ff7f0e"))+ guides(color = guide_legend(ncol = 1))
p4 <- DimPlot(integrated_obj, reduction = "tsne", group.by = "integrated_snn_res.1.2", label = TRUE, repel = TRUE, pt.size = 2)+ guides(color = guide_legend(ncol = 1))

pdf( "tsne_inte.pdf" , width=12, height=5, onefile = F)
print( p3 | p4)
dev.off()

## 4) PCA - merge ----
Idents(merge_data) <- "cell_info"

RhpcBLASctl::blas_set_num_threads(8)
# Variable features
### b. pca ----
merge_data <- FindVariableFeatures(merge_data, nfeatures = 2000)

merge_data <- ScaleData(merge_data, features = VariableFeatures(merge_data), do.center = TRUE) %>%
  RunPCA(npcs = 30, features = VariableFeatures(merge_data))

## 5) UMAP & tSNE - merge ----
### a. compare diff resolutions ----
merge_data <- FindNeighbors(merge_data, dims = 1:30)

### b. select a resolution ----
merge_data <- FindClusters(merge_data, resolution = 0.6)
merge_data <- RunTSNE(merge_data, dims = 1:30 )
merge_data <- RunUMAP(merge_data, dims = 1:30 )

p1 <- DimPlot(merge_data, reduction = "umap", group.by = "cell_info", label = TRUE, repel = TRUE, pt.size = 2, cols = c("#1f77b4", "#ff7f0e"))+ guides(color = guide_legend(ncol = 1))
p2 <- DimPlot(merge_data, reduction = "umap", group.by = "RNA_snn_res.0.6", label = TRUE, repel = TRUE, pt.size = 2)+ guides(color = guide_legend(ncol = 1))

pdf( "umap_merge.pdf" , width=12, height=5, onefile = F)
print( p1 | p2)
dev.off()

p3 <- DimPlot(merge_data, reduction = "tsne", group.by = "cell_info", label = TRUE, repel = TRUE, pt.size = 2, cols = c("#1f77b4", "#ff7f0e"))+ guides(color = guide_legend(ncol = 1))
p4 <- DimPlot(merge_data, reduction = "tsne", group.by = "RNA_snn_res.0.6", label = TRUE, repel = TRUE, pt.size = 2)+ guides(color = guide_legend(ncol = 1))

pdf( "tsne_merge.pdf" , width=12, height=5, onefile = F)
print( p3 | p4)
dev.off()

# 2. 分析需求1.4: 分群注释出TE，PrE和ICM/EPI三群 ----
markers <- list(
  EPI = c("Pou5f1", "Nanog", "Sox2", "Tdgf1", "Klf4", "Klf2", "Klf15"),
  TE = c("Cdx2", "Gata3", "Tfap2c", "Tead4", "Elf5", "Ascl2", "Id2", "Eomes", "Krt18", "Krt8", "Krt7", "Gata2"),
  PrE = c("Sox17", "Gata4", "Gata6", "Sox7", "Pdgfra", "Hnf4a")
)

## 1) anno - plot ----
create_violin <- function(gene, object, slot_data) {
  assay_to_use <- ifelse(gene %in% rownames(object@assays$integrated), "integrated", "RNA")
  
  data_to_plot <- data.frame(
    expression = GetAssayData(object, assay = assay_to_use, slot = slot_data)[gene, ],
    group = object$seurat_clusters
  )
  
  ggplot(data_to_plot, aes(x = group, y = expression, fill = group)) +
    geom_violin(scale = "width", trim = TRUE) +
    #scale_fill_manual(values = c("#1f77b4", "#ff7f0e", "#2ca02c")) +
    labs(title = gene, y = "Expression") +
    theme_void() +
    theme(
      panel.background = element_blank(),
      plot.background = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
      axis.title = element_blank(),
      #axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      plot.title = element_text(size = 8, hjust = 0.5, face = "bold"),
      legend.position = "none",
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 9))
}

### a. vlnplot - inte ----
available_genes <- intersect(c(unlist(markers), "T", "Mixl1", "Mesp1"), rownames(integrated_obj@assays$RNA))
plot_list <- lapply(available_genes, function(g) create_violin(g, integrated_obj, "data"))

combined_plot <- wrap_plots(plot_list, ncol = 5) +
  plot_layout(guides = "collect") & 
  theme(legend.position = "right")

ggsave("anno_marker_expression_violins_inte.pdf", 
       combined_plot, 
       width = 5 * 4, 
       height = ceiling(length(c(unlist(markers), "T", "Mixl1", "Mesp1"))/5)*1.5,
       units = "in")

### b. vlnplot - merge ----
available_genes <- intersect(c(unlist(markers), "T", "Mixl1", "Mesp1"), rownames(merge_data@assays$RNA@scale.data))
plot_list <- lapply(available_genes, function(g) create_violin(g, merge_data, "data"))

combined_plot <- wrap_plots(plot_list, ncol = 5) +
  plot_layout(guides = "collect") & 
  theme(legend.position = "right")

ggsave("anno_marker_expression_violins_merge.pdf", 
       combined_plot, 
       width = 5 * 4, 
       height = ceiling(length(c(unlist(markers), "T", "Mixl1", "Mesp1"))/5)*1.5,
       units = "in")

### c. heatmap - inte ----
exp_data <- GetAssayData(integrated_obj, assay = "RNA", slot = "data")
exp_data <- exp_data[intersect(rownames(exp_data), c(unlist(markers), "T", "Mixl1", "Mesp1")), ]
colnames(exp_data) <- integrated_obj$seurat_clusters

meta <- integrated_obj@meta.data
exp_scale = exp_data
#exp_scale = t(scale(t(exp_data), center = T, scale = T))

cluster_levels <- sort(unique(meta$seurat_clusters))
cluster_colors <- setNames(
  brewer.pal(length(cluster_levels), "Set3")[1:length(cluster_levels)],
  cluster_levels
)

# Create column annotation
column_ha <- HeatmapAnnotation(
  Cluster = meta$seurat_clusters,
  col = list(Cluster = cluster_colors),
  show_annotation_name = TRUE
)

pdf('marker genes_heatmap_scale.pdf', width = 12, height = 5, onefile = F)
Heatmap(as.matrix(exp_scale),
        column_split = meta$seurat_clusters,
        top_annotation = column_ha,
        show_row_names = T,
        show_column_names = F,  # Show column names
        show_column_dend = T,  # Show column dendrogram
        column_title = "Marker genes expression heatmap",
        use_raster = F,
        col = colorRampPalette(rev(brewer.pal(9,"Spectral")))(100),
        heatmap_legend_param = list(title=''))
dev.off()

### d. heatmap - merge ----
exp_data <- GetAssayData(merge_data, assay = "RNA", slot = "scale.data")
exp_data <- exp_data[intersect(rownames(exp_data), c(unlist(markers), "T", "Mixl1", "Mesp1")), ]
colnames(exp_data) <- merge_data$seurat_clusters

meta <- merge_data@meta.data
exp_scale = exp_data
#exp_scale = t(scale(t(exp_data)))

cluster_levels <- sort(unique(meta$seurat_clusters))
cluster_colors <- setNames(
  brewer.pal(length(cluster_levels), "Set3")[1:length(cluster_levels)],
  cluster_levels
)

# Create column annotation
column_ha <- HeatmapAnnotation(
  Cluster = meta$seurat_clusters,
  col = list(Cluster = cluster_colors),
  show_annotation_name = TRUE
)

pdf('marker genes_heatmap_scale_merge.pdf', width = 12, height = 5, onefile = F)
Heatmap(as.matrix(exp_scale),
        column_split = meta$seurat_clusters,
        top_annotation = column_ha,
        show_row_names = T,
        show_column_names = F,  # Show column names
        show_column_dend = T,  # Show column dendrogram
        column_title = "Marker genes expression heatmap",
        use_raster = F,
        col = colorRampPalette(rev(brewer.pal(9,"Spectral")))(100),
        heatmap_legend_param = list(title=''))
dev.off()

## 2) anno - merge ----
merge_data <- AddModuleScore( merge_data, features = markers, name = c("EPI_score", "TE_score", "PrE_score")#, ctrl = 100  # 控制基因数量
)
scores <- merge_data@meta.data[, c("EPI_score1", "TE_score2", "PrE_score3")]
merge_data$cell_type <- gsub("_score|1|2|3","",colnames(scores)[max.col(scores)])

p1 <- DimPlot(merge_data, reduction = "umap", group.by = "cell_info", label = TRUE, repel = TRUE, pt.size = 2, cols = c("#1f77b4", "#ff7f0e"))+ guides(color = guide_legend(ncol = 1))
p2 <- DimPlot(merge_data, reduction = "umap", group.by = "cell_type", label = TRUE, repel = TRUE, pt.size = 2)+ guides(color = guide_legend(ncol = 1))
p3 <- DimPlot(merge_data, reduction = "umap", group.by = "RNA_snn_res.0.6", label = TRUE, repel = TRUE, pt.size = 2)+ guides(color = guide_legend(ncol = 1))
p4 <- DimPlot(merge_data, reduction = "tsne", group.by = "cell_info", label = TRUE, repel = TRUE, pt.size = 2, cols = c("#1f77b4", "#ff7f0e"))+ guides(color = guide_legend(ncol = 1))
p5 <- DimPlot(merge_data, reduction = "tsne", group.by = "cell_type", label = TRUE, repel = TRUE, pt.size = 2)+ guides(color = guide_legend(ncol = 1))
p6 <- DimPlot(merge_data, reduction = "tsne", group.by = "RNA_snn_res.0.6", label = TRUE, repel = TRUE, pt.size = 2)+ guides(color = guide_legend(ncol = 1))

pdf( "cell_type_merge.pdf" , width=18, height=10, onefile = F)
print( (p1 | p2 | p3) / (p4 | p5 | p6))
dev.off()

saveRDS(merge_data, "merge_data_with_anno.rds")

## 3) anno - inte ----
integrated_obj$cell_type <- ifelse(integrated_obj$integrated_snn_res.1.2 %in% c(0,1,2,5,8,10,11), "EPI",
                                   ifelse(integrated_obj$integrated_snn_res.1.2 %in% c(3,4,7,9), "TE", "PrE"))

p1 <- DimPlot(integrated_obj, reduction = "umap", group.by = "cell_info", label = TRUE, repel = TRUE, pt.size = 2, cols = c("#1f77b4", "#ff7f0e"))+ guides(color = guide_legend(ncol = 1))
p2 <- DimPlot(integrated_obj, reduction = "umap", group.by = "cell_type", label = TRUE, repel = TRUE, pt.size = 2)+ guides(color = guide_legend(ncol = 1))
p3 <- DimPlot(integrated_obj, reduction = "umap", group.by = "integrated_snn_res.1.2", label = TRUE, repel = TRUE, pt.size = 2)+ guides(color = guide_legend(ncol = 1))
p4 <- DimPlot(integrated_obj, reduction = "tsne", group.by = "cell_info", label = TRUE, repel = TRUE, pt.size = 2, cols = c("#1f77b4", "#ff7f0e"))+ guides(color = guide_legend(ncol = 1))
p5 <- DimPlot(integrated_obj, reduction = "tsne", group.by = "cell_type", label = TRUE, repel = TRUE, pt.size = 2)+ guides(color = guide_legend(ncol = 1))
p6 <- DimPlot(integrated_obj, reduction = "tsne", group.by = "integrated_snn_res.1.2", label = TRUE, repel = TRUE, pt.size = 2)+ guides(color = guide_legend(ncol = 1))

pdf( "cell_type_inte.pdf" , width=18, height=10, onefile = F)
print( (p1 | p2 | p3) / (p4 | p5 | p6))
dev.off()

saveRDS(integrated_obj, "inte_data_with_anno.rds")

# 3. 分析需求3：Violin plot of  TE marker/EPI marker and prE marker expression of in house blastoid-10X 和public blastocyst 10X ----
epi_markers <- c("Pou5f1", "Nanog", "Klf4", "Sox2", "Klf2", "Tdgf1", "Klf15")
pre_markers <- c("Sox17", "Gata4", "Gata6", "Sox7", "Pdgfra", "Hnf4a")
te_markers  <- c("Cdx2", "Gata3", "Tfap2c", "Tead4", "Elf5", "Ascl2", "Id2", 
                 "Eomes", "Krt18", "Krt8", "Krt7", "Gata2")

create_violin <- function(gene, object, slot_data) {
  assay_to_use <- ifelse(gene %in% rownames(object@assays$integrated), "integrated", "RNA")
  
  data_to_plot <- data.frame(
    expression = GetAssayData(object, assay = assay_to_use, slot = slot_data)[gene, ],
    group = object$cell_type
  )
  
  ggplot(data_to_plot, aes(x = group, y = expression, fill = group)) +
    geom_violin(scale = "width", trim = TRUE) +
    scale_fill_manual(values = c("#1f77b4", "#ff7f0e", "#11999e")) +
    labs(title = gene, y = "Expression") +
    theme_void() +
    theme(
      panel.background = element_blank(),
      plot.background = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      plot.title = element_text(size = 8, hjust = 0.5, face = "bold"),
      legend.position = "none",
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 9))
}

plot_save_vln <- function(input_data_expr, input_data, label, title){
  available_genes <- intersect(all_genes, rownames(input_data_expr))
  plot_list <- lapply(available_genes, function(g) create_violin(g, input_data, "data"))
  
  combined_plot <- wrap_plots(plot_list, ncol = 5) +
    plot_layout(guides = "collect") & 
    theme(legend.position = "right")
  
  combined_plot <- combined_plot +
    plot_annotation(
      title = title,
      theme = list(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
    )
  
  ggsave(paste0("marker_expression_violins_",label,".pdf"), 
         combined_plot, 
         width = 5 * 1.5, 
         height = ceiling(length(all_genes)/5)*0.8,
         units = "in")
}

all_genes <- c(epi_markers, pre_markers, te_markers)

## 1) marker expr vln - inte ----
plotdata <- subset(integrated_obj, cell_info == "In-house Blastoid")
plot_save_vln(input_data_expr = plotdata@assays$integrated, input_data = plotdata, label = "inte_blastoid", title = "In-house Blastoid")

plotdata <- subset(integrated_obj, cell_info == "Public Blastocyst")
plot_save_vln(input_data_expr = plotdata@assays$integrated, input_data = plotdata, label = "inte_blastocyst", title = "Public Blastocyst")

## 2) marker expr vln - merge ----
plotdata <- subset(merge_data, cell_info == "In-house Blastoid")
plot_save_vln(input_data_expr = plotdata@assays$RNA@scale.data, input_data = plotdata, label = "merge_blastoid", title = "In-house Blastoid")

plotdata <- subset(merge_data, cell_info == "Public Blastocyst")
plot_save_vln(input_data_expr = plotdata@assays$RNA@scale.data, input_data = plotdata, label = "merge_blastocyst", title = "Public Blastocyst")

# 6. 分析需求5：分析in house blastoid和 public EPSC-blastoid (10X，GSM4026211)，我们的in house blastoid不表达除 EPI/TE/prE 三个lineage以外的中间态细胞（这群中间态细胞表达mesoderm marker: T / Mixl1 / Mesp1） ----
## 1) blastoid 10x ----
inhouse_filtered <- readRDS("../250430/inhouse_blastoid_filtered.rds")

## 2) GSM4026211 ----
data <- Read10X(data.dir = "../rawdata/GSE135701/GSM4026211/")
assay_data <- CreateAssayObject(counts = data,
                                min.cells = 3, min.features = 200)
GSM4026211 <- CreateSeuratObject(counts = assay_data, project = "GSM4026211", assay = "RNA",
                                 min.cells = 3, min.features = 200)
GSM4026211$cell_info <- "EPSC-public blastoid"

## 3) GSM4026212 ----
data <- Read10X(data.dir = "../rawdata/GSE135701/GSM4026212/")

GSM4026212 <- CreateSeuratObject(counts = data, project = "GSM4026212", assay = "RNA",
                                 min.cells = 3, min.features = 82)

GSM4026212$cell_info <- "Public Blastocyst"

## 3) merge data ----
### a. merge data ----
merge_data <- merge(inhouse_filtered, GSM4026211, add.cell.ids = c("Inhouse", "Pub")) %>% 
  merge(., GSM4026212)
saveRDS(merge_data, "raw_merge_data_GSM4026211_GSM4026212.rds")

### b. integrate data ----
inhouse_filtered <- NormalizeData(inhouse_filtered) %>% FindVariableFeatures( selection.method = "vst", nfeatures = 5000) %>% 
  ScaleData( do.center = FALSE)
GSM4026211 <- NormalizeData(GSM4026211) %>% FindVariableFeatures( selection.method = "vst", nfeatures = 5000) %>% 
  ScaleData( do.center = FALSE)
GSM4026212 <- NormalizeData(GSM4026212) %>% FindVariableFeatures( selection.method = "vst", nfeatures = 5000) %>% 
  ScaleData( do.center = FALSE)

anchors <- FindIntegrationAnchors(object.list = list(inhouse_filtered, GSM4026211, GSM4026212), dims = 1:30)

integrated_obj <- IntegrateData(anchorset = anchors, dims = 1:30)
saveRDS(integrated_obj, "raw_integrated_data_GSM4026211_GSM4026212.rds")

## 4) mesoderm marker: T / Mixl1 / Mesp1 ----
library(ComplexHeatmap)
library(circlize)

mesoderm <- c("T", "Mixl1", "Mesp1")
create_violin <- function(gene, object, slot_data) {
  assay_to_use <- ifelse(gene %in% rownames(object@assays$integrated), "integrated", "RNA")
  
  data_to_plot <- data.frame(
    expression = GetAssayData(object, assay = assay_to_use, slot = slot_data)[gene, ],
    group = object$cell_info
  )
  
  ggplot(data_to_plot, aes(x = group, y = expression, fill = group)) +
    geom_violin(scale = "width", trim = TRUE) +
    scale_fill_manual(values = c("#1f77b4", "#ff7f0e", "#2ca02c")) +
    labs(title = gene, y = "Expression") +
    theme_void() +
    theme(
      panel.background = element_blank(),
      plot.background = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
      axis.title = element_blank(),
      #axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      plot.title = element_text(size = 8, hjust = 0.5, face = "bold"),
      legend.position = "none",
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 9))
}

## 1) marker expr vln - inte ----
available_genes <- intersect(mesoderm, rownames(integrated_obj@assays$RNA))
plot_list <- lapply(available_genes, function(g) create_violin(g, integrated_obj, "data"))

combined_plot <- wrap_plots(plot_list, ncol = 3) +
  plot_layout(guides = "collect") & 
  theme(legend.position = "right")

ggsave("mesoderm_marker_expression_violins_inte.pdf", 
       combined_plot, 
       width = 3 * 2, 
       height = ceiling(length(mesoderm)/3)*1.5,
       units = "in")

## 2) marker expr vln - merge ----
available_genes <- intersect(mesoderm, rownames(merge_data@assays$RNA@scale.data))
plot_list <- lapply(mesoderm, function(g) create_violin(g, merge_data, "data"))

combined_plot <- wrap_plots(plot_list, ncol = 3) +
  plot_layout(guides = "collect") & 
  theme(legend.position = "right")

ggsave("mesoderm_marker_expression_violins_merge.pdf", 
       combined_plot, 
       width = 3 * 2, 
       height = ceiling(length(mesoderm)/3)*1.5,
       units = "in")

# 7. 分析需求6: blatocytst和blastoid  prE/TE/ICM 的DEGs火山图 ----
volcanoFUN <- function(dataset=NULL,
                       title=NULL,
                       sampleoutpath=NULL,
                       sample=NULL,
                       use_adj.pvalue=TRUE,
                       cut_off_geneLabel=10,
                       cut_off_logFC=NULL,
                       cut_off_adj.pvalue=NULL,
                       cut_off_pvalue=0.05,
                       plotCol=c("#BC3C28","grey","#0072B5"),
                       labelUp=NULL,labelDown=NULL,
                       annoData.x=2,
                       annoData.y=1.5,
                       w=4,h=4.6){
  if (use_adj.pvalue) {
    top_upregulated <- subset(dataset,p_val_adj < cut_off_adj.pvalue) %>% arrange(desc(avg_log2FC)) %>%
      head(cut_off_geneLabel) %>% rownames
    top_downregulated <- subset(dataset,p_val_adj < cut_off_adj.pvalue) %>% arrange(avg_log2FC) %>% 
      head(cut_off_geneLabel) %>% rownames
  } else {
    top_upregulated <- subset(dataset,p_val < cut_off_pvalue) %>% arrange(desc(avg_log2FC)) %>%
      head(cut_off_geneLabel) %>% rownames
    top_downregulated <- subset(dataset,p_val < cut_off_pvalue) %>% arrange(avg_log2FC) %>% 
      head(cut_off_geneLabel) %>% rownames
  }
  top_genes <- c(top_upregulated, top_downregulated)
  
  dataset$label = ifelse((rownames(dataset) %in% top_genes),
                         as.character(rownames(dataset)),"")
  
  library(ggrepel)
  cut_off_pvalue = cut_off_pvalue #0.05
  cut_off_adj.pvalue = cut_off_adj.pvalue
  cut_off_logFC = cut_off_logFC
  dataset = dataset[dataset$p_val_adj!=1,]
  dataset$gene = rownames(dataset)
  if(cut_off_adj.pvalue <= 0.05){
    dataset$change = ifelse(dataset$p_val < cut_off_pvalue & dataset$p_val_adj < cut_off_adj.pvalue & 
                              abs(dataset$avg_log2FC) >= cut_off_logFC, 
                            ifelse(dataset$avg_log2FC> cut_off_logFC ,labelUp,labelDown),
                            'NoSig') %>% factor(levels = c(labelUp,'NoSig',labelDown))
  } else {
    dataset$change = ifelse(dataset$p_val < cut_off_pvalue & 
                              abs(dataset$avg_log2FC) >= cut_off_logFC, 
                            ifelse(dataset$avg_log2FC> cut_off_logFC ,labelUp,labelDown),
                            'NoSig') %>% factor(levels = c(labelUp,'NoSig',labelDown))
  }
  
  annoData <- data.frame(x=c(-annoData.x,annoData.x),y=c(annoData.y,annoData.y),
                         label=c(length(which(dataset$change == labelDown)),
                                 length(which(dataset$change == labelUp))))
  
  pdf(paste0(sampleoutpath,"/",sample,"Volcanoplot_",title,".pdf"),width = w,height = h)
  
  if(cut_off_adj.pvalue<=0.05){
    p = ggplot(
      dataset,aes(x = avg_log2FC,
                  y = -log10(p_val_adj),
                  colour=change)) +
      geom_point(alpha=1, size=0.5) +
      scale_color_manual(values=plotCol)+
      
      geom_text_repel(aes(label=label), 
                      segment.color = "black",
                      show.legend = FALSE,
                      box.padding=unit(0.15, "lines"), 
                      point.padding=unit(0.5, "lines"), 
                      #color="white",
                      max.overlaps = 10000000000)+
      geom_text(data = annoData, aes(x = x, y = y, label = label), 
                color = "black", size = 8, angle = 0, fontface = "bold")+ 
      theme_bw()+
      
      geom_vline(xintercept=c(-cut_off_logFC,cut_off_logFC),lty=4,col="grey",lwd=0.6) +
      geom_hline(yintercept = -log10(cut_off_adj.pvalue),lty=4,col="grey",lwd=0.6) +
      
      theme_classic(base_size = 15) + 
      theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
            legend.text.align = 0,
            legend.title = element_blank(),
            panel.grid = element_blank(),
            axis.text = element_text(size = 13),
            axis.title = element_text(size=15),
            legend.position = "top",
            legend.text = element_text(size = 13)) +
      guides(colour = guide_legend(override.aes = list(size=3),reverse = T))+ #可使得legend中的圆点变大
      
      labs(x=expression("Log"["2"]*"(Fold change)"),
           y=expression("-Log"["10"]*"(adjust P vaule)"),
           title = title)
  } else {
    p = ggplot(
      dataset,aes(x = avg_log2FC,
                  y = -log10(p_val),
                  colour=change)) +
      geom_point(alpha=1, size=0.5) +
      scale_color_manual(values=plotCol)+
      
      geom_text_repel(aes(label=label), 
                      segment.color = "black",
                      show.legend = FALSE,
                      box.padding=unit(0.15, "lines"), 
                      point.padding=unit(0.5, "lines"), 
                      #color="white",
                      max.overlaps = 10000000000)+
      geom_text(data = annoData, aes(x = x, y = y, label = label), 
                color = "black", size = 8, angle = 0, fontface = "bold")+ 
      theme_bw()+
      
      geom_vline(xintercept=c(-cut_off_logFC,cut_off_logFC),lty=4,col="grey",lwd=0.6) +
      geom_hline(yintercept = -log10(cut_off_pvalue),lty=4,col="grey",lwd=0.6) +
      
      theme_classic(base_size = 15) + 
      theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
            legend.text.align = 0,
            legend.title = element_blank(),
            panel.grid = element_blank(),
            axis.text = element_text(size = 13),
            axis.title = element_text(size=15),
            legend.position = "top",
            legend.text = element_text(size = 13)) +
      guides(colour = guide_legend(override.aes = list(size=3),reverse = T))+ #可使得legend中的圆点变大
      
      labs(x=expression("Log"["2"]*"(Fold change)"),
           y=expression("-Log"["10"]*"(P vaule)"),
           title = title)
  }
  print(p)
  dev.off()
}

## 1) DEG - inte ----
Idents(integrated_obj) <- "cell_info"

cell_types <- c("EPI", "PrE", "TE")

deg_list <- list()
for (ct in cell_types) {
  cells_in_type <- WhichCells(integrated_obj, expression = cell_type == ct)
  
  subset_obj <- subset(integrated_obj, cells = cells_in_type)
  
  Idents(subset_obj) <- "cell_info"
  
  deg <- FindMarkers(
    subset_obj,
    ident.1 = "In-house Blastoid",
    ident.2 = "Public Blastocyst",
    logfc.threshold = 0.25,
    min.pct = 0.1,
    test.use = "wilcox"
  )
  
  volcanoFUN(dataset = deg,
             title = ct,
             sample = "inte_", sampleoutpath = ".", 
             cut_off_geneLabel = 5,
             use_adj.pvalue = FALSE,
             cut_off_logFC = 0.25, cut_off_adj.pvalue = 0.5, cut_off_pvalue = 0.05,
             plotCol = c("#BC3C28", "#0072B5"),
             labelUp = "In-house Blastoid",
             labelDown = "Public Blastocyst",
             annoData.x=2,
             annoData.y=4,
             w=8, h=6)
  
  deg$cell_type <- ct
  deg$gene <- rownames(deg)
  
  deg_list[[ct]] <- deg
}

deg_all <- bind_rows(deg_list)
write.table(deg_all, "DEG_all_inte.txt", quote = F, row.names = T)

## 2) DEG - merge ----
Idents(merge_data) <- "cell_info"

cell_types <- c("EPI", "PrE", "TE")

deg_list <- list()
for (ct in cell_types) {
  cells_in_type <- WhichCells(merge_data, expression = cell_type == ct)
  
  subset_obj <- subset(merge_data, cells = cells_in_type)
  
  Idents(subset_obj) <- "cell_info"
  
  deg <- FindMarkers(
    subset_obj,
    ident.1 = "In-house Blastoid",
    ident.2 = "Public Blastocyst",
    logfc.threshold = 0.25,
    min.pct = 0.1,
    test.use = "wilcox"
  )
  
  volcanoFUN(dataset = deg,
             title = ct,
             sample = "merge_", sampleoutpath = ".", 
             cut_off_geneLabel = 5,
             use_adj.pvalue = FALSE,
             cut_off_logFC = 0.25, cut_off_adj.pvalue = 0.5, cut_off_pvalue = 0.05,
             plotCol = c("#BC3C28", "#0072B5"),
             labelUp = "In-house Blastoid",
             labelDown = "Public Blastocyst",
             annoData.x=2,
             annoData.y=4,
             w=8, h=6)
  
  deg$cell_type <- ct
  deg$gene <- rownames(deg)
  
  deg_list[[ct]] <- deg
}

deg_all <- bind_rows(deg_list)
write.table(deg_all, "DEG_all_inte.txt", quote = F, row.names = T)
