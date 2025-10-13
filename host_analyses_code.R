#Code for host- and integrated host-microbial analyses for "Distinct Host and Microbial Biology Distinguishes Lower Respiratory Tract Infection from Incidental Pathogen Carriage"
#Includes code for generating the relevant panels of Figures 3 and 4

#loading necessary packages
library(dplyr)
library(limma)
library(ggplot2)
library(ggrepel)
library(vegan)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(tibble)
library(forcats)
library(ComplexHeatmap)
library(colorRamp2)
library(ggpubr)
library(ggsignif)
library(stringr)

#importing metadata and gene counts files
metadata <- read.csv("/Users/emilylydon/Library/CloudStorage/Box-Box/VAP_IvC_project/Code:Source Data File/metadata.csv", check.names = FALSE)
counts <- read.csv("/Users/emilylydon/Library/CloudStorage/Box-Box/VAP_IvC_project/Code:Source Data File/gene_counts.csv", check.names = FALSE, row.names=1)

#extracting key of ENSEMBL IDs/gene symbols
gene_symbols <- counts[, "gene_symbol", drop = FALSE] #extracting gene names for later
counts <- counts[, -1] #getting rid of gene symbol to have a clean gene counts table

metadata$group <- factor(
  metadata$group,
  levels = c("LRTI", "IPC", "CTRL")
)

#VST transformation of gene counts
dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts, 
                                      colData = metadata, 
                                      design = ~1)
vsd <- DESeq2::varianceStabilizingTransformation(dds)
counts_VST <- SummarizedExperiment::assay(vsd) %>%
  round(.,digits=2)
counts_VST <- as.data.frame(counts_VST)

####################################
#Figure 3
####################################

# PCA of all gene counts (figure 3A)
gene_matrix <- t(counts_VST)
pca <- prcomp(gene_matrix, center = TRUE, scale. = TRUE) #running PCA on transposed gene counts table
pca_df <- as.data.frame(pca$x[, 1:10])  # extracting data for the first 10 PCs
pca_df$Patient <- rownames(pca_df)
pca_df <- merge(pca_df, metadata, by = "Patient") #merging with metadata to assign groups for plotting

# Color palette for plotting
colors <- c(
  "IPC" = "#008080FF",
  "LRTI" = "#CA562CFF",
  "CTRL" = "gray"
)

#generating a PCA plot. can substitute different PCs (PC1 v PC3 in main text, PC1 v PC2 and PC2 vs PC3 in the supplement)
PCA_plot <- ggplot(pca_df, aes(x = PC1, y = PC3, color = group)) +
  geom_point(size = 1.5, alpha = 0.9) +
  stat_ellipse(level = 0.68, size = 1.2, alpha = 0.9) +
  scale_color_manual(values = colors) +
  theme_minimal(base_family = "Arial", base_size = 12) +
  labs(
    x = paste0("PC1"),
    y = paste0("PC3")
  ) +
  coord_cartesian(xlim = c(-200, 200), ylim = c(-65, 80))+
  theme(
    plot.title = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 12),
    legend.position = c(0.1, 0.8),
    legend.justification = c(0, 0),
    legend.key.size = unit(0.6, "lines"),  
    panel.grid = element_blank(),  
    axis.line.x = element_line(color = "black", linewidth = 0.4),
    axis.line.y = element_line(color = "black", linewidth = 0.4),
    axis.ticks = element_line(color = "black", linewidth = 0.4)
  )
PCA_plot

#statistical testing between each group using pairwise PERMANOVA with adonis and printing output
group_levels <- unique(metadata$group)
pairwise_comparisons <- combn(group_levels, 2, simplify = FALSE)

pairwise_results <- lapply(pairwise_comparisons, function(pair) {
  idx <- metadata$group %in% pair
  sub_gene_matrix <- gene_matrix[idx, ]
  sub_labels <- metadata[idx, ]
  dist_matrix <- dist(sub_gene_matrix, method = "euclidean")
  ad <- adonis2(dist_matrix ~ group, data = sub_labels)
  data.frame(
    Group1 = pair[1],
    Group2 = pair[2],
    F_value = ad$F[1],
    R2 = ad$R2[1],
    p_value = ad$`Pr(>F)`[1]
  )
})

pairwise_adonis_df <- do.call(rbind, pairwise_results)
pairwise_adonis_df$p_adj <- p.adjust(pairwise_adonis_df$p_value, method = "BH")
pairwise_adonis_df #prints output, can then be added to plot if desired

# Pairwise limma analyses between each of the three groups

#first comparison: LRTI versus controls
metadata_LRTI_CTRL <- metadata %>% filter(group == "LRTI" | group == "CTRL") #subsetting for the comparison
metadata_LRTI_CTRL$group <- factor(metadata_LRTI_CTRL$group, levels=c("CTRL", "LRTI")) #establishing factors/levels
counts_LRTI_CTRL <- counts[, metadata_LRTI_CTRL$Patient] #creating aligned counts table matching the metadata (non-VST transformed since limma will do transformation)

design <- model.matrix( ~group + X_AGE_YEARS + X_SEX,
                        data = metadata_LRTI_CTRL) #generating a model matrix, adding age and sex as covariates. For SDI adjustment (panel 4i, add +SDI as an additional covariate)

vwts <- voom(counts_LRTI_CTRL,
             design = design,
             normalize.method = "quantile",
             plot = TRUE) #running limma/voom
vfit <- lmFit(vwts)
vfit <- eBayes(vfit)
summary(decideTests(vfit)) #outputting DE gene table
DEresults_LRTI_CTRL <- topTable(
  vfit,
  coef = "groupLRTI",
  sort.by="none",
  number = Inf
)
DEresults_LRTI_CTRL$gene_symbol <- gene_symbols[rownames(DEresults_LRTI_CTRL), "gene_symbol"] #adding the gene symbol for reference

#second comparison: LRTI vs IPC
metadata_LRTI_IPC <- metadata %>% filter(group == "LRTI" | group == "IPC")
metadata_LRTI_IPC$group <- factor(metadata_LRTI_IPC$group, levels=c("IPC", "LRTI")) #establishing factors/levels
counts_LRTI_IPC <- counts[, metadata_LRTI_IPC$Patient] #creating aligned counts table matching the metadata (non-VST transformed since limma will do transformation)

design <- model.matrix( ~group + X_AGE_YEARS + X_SEX,
                        data = metadata_LRTI_IPC) #generating model matrix, adding age/sex as covariates. For SDI adjustment (panel 4i, add +SDI as an additional covariate

vwts <- voom(counts_LRTI_IPC,
             design = design,
             normalize.method = "quantile",
             plot = TRUE) #running limma-voom
vfit <- lmFit(vwts)
vfit <- eBayes(vfit)
summary(decideTests(vfit)) #outputting DE gene table
DEresults_LRTI_IPC <- topTable(
  vfit,
  coef = "groupLRTI",
  sort.by="none",
  number = Inf
)
DEresults_LRTI_IPC$gene_symbol <- gene_symbols[rownames(DEresults_LRTI_IPC), "gene_symbol"] #adding the gene symbol for reference

#third comparison: IPC vs CTRL
metadata_IPC_CTRL <- metadata %>% filter(group == "IPC" | group == "CTRL")
metadata_IPC_CTRL$group <- factor(metadata_IPC_CTRL$group, levels=c("CTRL", "IPC")) #establishing factors/levels for limma
counts_IPC_CTRL <- counts[, metadata_IPC_CTRL$Patient] #creating aligned counts table matching the metadata (non-VST transformed since limma will do transformation)
table(colnames(counts_IPC_CTRL) == metadata_IPC_CTRL$Patient) #sanity check to make sure things are lined up

design <- model.matrix( ~group + X_AGE_YEARS + X_SEX,
                        data = metadata_IPC_CTRL) #generating a model matrix, including age and sex. For SDI adjustment (panel 4i, add +SDI as an additional covariate

#running limma-voom
vwts <- voom(counts_IPC_CTRL,
             design = design,
             normalize.method = "quantile",
             plot = TRUE)
vfit <- lmFit(vwts)
vfit <- eBayes(vfit)
summary(decideTests(vfit)) #outputting DE gene table
DEresults_IPC_CTRL <- topTable(
  vfit,
  coef = "groupIPC",
  sort.by="none",
  number = Inf
)
DEresults_IPC_CTRL$gene_symbol <- gene_symbols[rownames(DEresults_IPC_CTRL), "gene_symbol"]

#generating the heat map (figure 3b)
top20_genes <- DEresults_LRTI_CTRL[order(DEresults_LRTI_CTRL$adj.P.Val), ][1:20, ] #selecting the top 20 DE genes based on adj.p.val between LRTI and CTRL
counts_heatmap <- counts_VST[rownames(top20_genes), ] #subsetting the VST gene counts matrix to just those 20 genes
gene_symbols_vector <- gene_symbols$gene_symbol #extracting the gene symbol for display purposes
names(gene_symbols_vector) <- rownames(gene_symbols)
rownames(counts_heatmap) <- gene_symbols_vector[rownames(counts_heatmap)]

class_colors <- c(
  "LRTI" = "#CA562CFF",
  "IPC" = "#008080FF",
  "CTRL" = "grey"
) #color mapping based on infection class

col_fun <- colorRamp2(
  c(-2, 0, 2),
  c("#6a51a3", "white", "#fff07a")
) #contrasting color scheme for the Z-scores

counts_heatmap_matrix <- as.matrix(counts_heatmap)
counts_heatmap_z <- t(scale(t(counts_heatmap_matrix))) #z-scoring for the heatmap

ha <- HeatmapAnnotation(
  df = data.frame(Group = metadata$group),
  col = list(Group = class_colors),
  show_legend = TRUE,
  annotation_name_gp = gpar(fontsize = 0),
  annotation_legend_param = list(
    Group = list(
      title = "", 
      direction = "horizontal",
      nrow = 1,
      labels_gp = gpar(fontsize = 10),
      at = c("LRTI", "IPC", "CTRL"), 
      labels = c("LRTI", "IPC", "CTRL")
    )
  )
) #annotations for the heat map

ht <- Heatmap(
  counts_heatmap_z,
  name = "Z score",
  col = col_fun,
  top_annotation = ha,
  show_row_names = TRUE,
  row_names_gp = gpar(fontsize = 9, fontface = "italic"),
  row_names_side = "left",
  cluster_rows = TRUE,
  show_row_dend = FALSE,
  cluster_columns = TRUE,
  show_column_dend = FALSE,
  clustering_distance_columns = function(x) as.dist(1 - cor(t(x), method = "pearson")),
  clustering_method_columns = "ward.D2",
  show_column_names = FALSE,
  heatmap_legend_param = list(
    title = "Z score",
    title_gp = gpar(fontsize = 10),
    labels_gp = gpar(fontsize = 10),
    direction = "horizontal",
    title_position = "topcenter",
    legend_height = unit(4, "cm"),
    grid_width = unit(0.6, "cm")
  )
)
draw(ht, heatmap_legend_side = "bottom", annotation_legend_side = "top") # plotting the heatmap

#volcano plots for DE analyses (figure 3c-3e)
#volcano plot #1: LRTI vs CTRL
DEresults_LRTI_CTRL <- DEresults_LRTI_CTRL %>% mutate(gene_type= case_when(logFC >0 & adj.P.Val <= 0.05 ~ "↑ in LRTI",
                                                                     logFC <=0 & adj.P.Val <= 0.05 ~ "↑ in CTRL",
                                                                     TRUE ~ "Not significant"))
n_DE <- sum(DEresults_LRTI_CTRL$gene_type != "Not significant")
cols <- c(
  "↑ in LRTI" = "#CA562CFF",
  "↑ in CTRL" = "grey",
  "Not significant" = "black"
)

DEresults_LRTI_CTRL$gene_type <- factor(
  DEresults_LRTI_CTRL$gene_type,
  levels = c("↑ in LRTI", "↑ in CTRL", "Not significant")
) # Set levels so they appear in desired order

volcano_LRTI_CTRL <- ggplot(DEresults_LRTI_CTRL, aes(x = logFC, y = -log10(adj.P.Val))) +
  geom_point(aes(color = gene_type), size = 1) +
  scale_color_manual(
    values = cols,
    breaks = c("↑ in LRTI", "↑ in CTRL"), # show only these in legend
    name = NULL
  ) +
  scale_x_continuous(limits = xlim_vals, expand = c(0, 0)) +
  scale_y_continuous(limits = ylim_vals, expand = c(0, 0)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +   # Significance line
  annotate("text",
           x = xlim_vals[2] - 2,                 
           y = ylim_vals[2] - 3,                
           label = paste0(n_DE, " DE genes"),
           size = 4,
           hjust = 0.5,             
           vjust = 1,               
           family = "Arial") +
  labs(
    title = "LRTI vs CTRL",
    x = "log2(fold change)",
    y = "-log10(adjusted P-value)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 12, face = "plain"),
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 12),
    text = element_text(size = 12, family = "Arial"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = unit(c(1, 1, 0.7, 0.5), "cm"),
    legend.position = "bottom",
    legend.direction = "vertical",
    legend.text = element_text(size = 12),
    legend.spacing.x = unit(0.4, 'cm'),
    axis.line.x = element_line(color = "black", linewidth = 0.4),
    axis.line.y = element_line(color = "black", linewidth = 0.4),
    axis.ticks = element_line(color = "black", linewidth = 0.4)
  ) +
  guides(color = guide_legend(override.aes = list(size = 2.5)))
volcano_LRTI_CTRL

#volcano plot #2: LRTI vs IPC
DEresults_LRTI_IPC <- DEresults_LRTI_IPC %>% mutate(gene_type= case_when(logFC >0 & adj.P.Val <= 0.05 ~ "↑ in LRTI",
                                                                 logFC <=0 & adj.P.Val <= 0.05 ~ "↑ in IPC",
                                                                 TRUE ~ "Not significant"))
n_DE <- sum(DEresults_LRTI_IPC$gene_type != "Not significant")
cols <- c(
  "↑ in LRTI" = "#CA562CFF",
  "↑ in IPC" = "#008080FF",
  "Not significant" = "black"
)
DEresults_LRTI_IPC$gene_type <- factor(
  DEresults_LRTI_IPC$gene_type,
  levels = c("↑ in LRTI", "↑ in IPC", "Not significant")
) # Set levels so they appear in desired order

xlim_vals <- c(-5, 5) #setting identical axes for panel 3c for easy visual comparison
ylim_vals <- c(0, 30) #setting identical axes for panel 3c for easy visual comparison

volcano_LRTI_IPC <- ggplot(DEresults_LRTI_IPC, aes(x = logFC, y = -log10(adj.P.Val))) +
  geom_point(aes(color = gene_type), size = 1) +
  scale_color_manual(
    values = cols,
    breaks = c("↑ in LRTI", "↑ in IPC"), 
    name = NULL
  ) +
  scale_x_continuous(limits = xlim_vals, expand = c(0, 0)) +
  scale_y_continuous(limits = ylim_vals, expand = c(0, 0)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +  # Significance line
  annotate("text",
           x = xlim_vals[2] - 2,                 
           y = ylim_vals[2] - 3,             
           label = paste0(n_DE, " DE genes"),
           size = 4,
           hjust = 0.5,          
           vjust = 1,        
           family = "Arial",
  ) +
  labs(
    title = "LRTI vs IPC",
    x = "log2(fold change)",
    y = "-log10(adjusted P-value)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 12, face = "plain"),
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 12),
    text = element_text(size = 12, family = "Arial"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = unit(c(1, 1, 0.7, 0.5), "cm"),
    legend.position = "bottom",
    legend.direction = "vertical",
    legend.text = element_text(size = 12),
    legend.spacing.x = unit(0.4, 'cm'),
    axis.line.x = element_line(color = "black", linewidth = 0.4),
    axis.line.y = element_line(color = "black", linewidth = 0.4),
    axis.ticks = element_line(color = "black", linewidth = 0.4),
  ) +
  guides(color = guide_legend(override.aes = list(size = 2.5)))
volcano_LRTI_IPC

#volcano plot #3: IPC vs CTRL
DEresults_IPC_CTRL <- DEresults_IPC_CTRL %>% mutate(gene_type= case_when(logFC >0 & adj.P.Val <= 0.05 ~ "↑ in IPC",
                                                                       logFC <=0 & adj.P.Val <= 0.05 ~ "↑ in CTRL",
                                                                       TRUE ~ "Not significant"))
n_DE <- sum(DEresults_IPC_CTRL$gene_type != "Not significant")
cols <- c(
  "↑ in IPC" = "#008080FF",
  "↑ in CTRL" = "grey",
  "Not significant" = "black"
)
DEresults_IPC_CTRL$gene_type <- factor(
  DEresults_IPC_CTRL$gene_type,
  levels = c("↑ in IPC", "↑ in CTRL", "Not significant")
) # set levels so they appear in desired order

volcano_IPC_CTRL <- ggplot(DEresults_IPC_CTRL, aes(x = logFC, y = -log10(adj.P.Val))) +
  geom_point(aes(color = gene_type), size = 1) +
  scale_color_manual(
    values = cols,
    breaks = c("↑ in IPC", "↑ in CTRL"), # show only these in legend
    name = NULL,
    drop = FALSE
  ) +
  scale_x_continuous(limits = xlim_vals, expand = c(0, 0)) + #identical to fig 3c for easy visual comparison
  scale_y_continuous(limits = ylim_vals, expand = c(0, 0)) + #identical to fig 3c for easy visual comparison
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  annotate("text",
           x = xlim_vals[2] - 2,                
           y = ylim_vals[2] - 3,        
           label = paste0(n_DE, " DE genes"),
           size = 4,
           hjust = 0.5,             
           vjust = 1,             
           family = "Arial") +
  labs(
    title = "IPC vs CTRL",
    x = "log2(fold change)",
    y = "-log10(adjusted P-value)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 12, face = "plain"),
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 12),
    text = element_text(size = 12, family = "Arial"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = unit(c(1, 1, 0.7, 0.5), "cm"),
    legend.position = "bottom",
    legend.direction = "vertical",
    legend.text = element_text(size = 12),
    legend.spacing.x = unit(0.4, 'cm'),
    axis.line.x = element_line(color = "black", linewidth = 0.4),
    axis.line.y = element_line(color = "black", linewidth = 0.4),
    axis.ticks = element_line(color = "black", linewidth = 0.4)
  ) +
  guides(color = guide_legend(override.aes = list(size = 2.5)))
volcano_IPC_CTRL

##GSEA for LRTI vs IPC (figure 3g)
set.seed(2)
gene_map <- bitr(rownames(DEresults_I_IC),
                 fromType = "ENSEMBL",
                 toType = "ENTREZID",
                 OrgDb = org.Hs.eg.db) #need entrez IDs for GSEA package (note: does lose some genes due to failure to map)
DEresults_LRTI_IPC_GSEA <- DEresults_LRTI_IPC %>%
  rownames_to_column("ENSEMBL")
DEresults_LRTI_IPC_GSEA <- inner_join(DEresults_LRTI_IPC_GSEA, gene_map, by = "ENSEMBL")

DEresults_LRTI_IPC_GSEA <- DEresults_LRTI_IPC_GSEA %>%
  filter(!is.na(t) & !is.na(ENTREZID)) %>%  
  group_by(ENTREZID) %>%
  arrange(adj.P.Val, .by_group = TRUE) %>%   
  dplyr::slice(1) %>%                                  
  ungroup() #getting rid of the very amount of duplicates (i.e. 1 entrez IDs mapping to 2 ENSEMBL IDs) by keeping the highest ranked one by adj p val

rankedgenes_LRTI_IPC <- DEresults_I_IC_GSEA$t #ranking all genes by T statistic for GSEA
names(rankedgenes_LRTI_IPC) <- DEresults_LRTI_IPC_GSEA$ENTREZID
rankedgenes_LRTI_IPC <- rankedgenes_LRTI_IPC + rnorm(length(rankedgenes_LRTI_IPC), mean = 0, sd = 1e-6) # adding a tiny jitter to break ties
rankedgenes_LRTI_IPC <- sort(rankedgenes_LRTI_IPC, decreasing = TRUE) #sorting in descending order for input into GSEA

gsea_LRTI_IPC <- gsePathway(
  geneList = rankedgenes_LRTI_IPC,
  organism = "human",
  pvalueCutoff = 0.05,
  minGSSize = 10,     
  maxGSSize = 1500, 
  eps = 0,
  verbose = TRUE
)
gsea_LRTI_IPC_df <- as.data.frame(gsea_LRTI_IPC)

##GSEA for IPC vs CTRL (figure 3g)
DEresults_IPC_CTRL_GSEA <- DEresults_IPC_CTRL %>%
  rownames_to_column("ENSEMBL")
DEresults_IPC_CTRL_GSEA <- inner_join(DEresults_IPC_CTRL_GSEA, gene_map, by = "ENSEMBL") #getting entrez IDs

DEresults_IPC_CTRL_GSEA <- DEresults_IPC_CTRL_GSEA %>%
  filter(!is.na(t) & !is.na(ENTREZID)) %>%    
  group_by(ENTREZID) %>%
  arrange(adj.P.Val, .by_group = TRUE) %>%      
  dplyr::slice(1) %>%                              
  ungroup() #getting rid of the small amount of duplicates (2 ENSEMBL IDs -> 1 entrez ID) by keeping the lowest p value

rankedgenes_IPC_CTRL <- DEresults_IPC_CTRL_GSEA$t #ranking all genes by T statistic for the GSEA
names(rankedgenes_IPC_CTRL) <- DEresults_IPC_CTRL_GSEA$ENTREZID
rankedgenes_IPC_CTRL <- rankedgenes_IPC_CTRL + rnorm(length(rankedgenes_IPC_CTRL), mean = 0, sd = 1e-6) #adding tiny jitter
rankedgenes_IPC_CTRL <- sort(rankedgenes_IPC_CTRL, decreasing = TRUE) #sorting in descending order

gsea_IC_CTRL <- gsePathway(
  geneList = rankedgenes_IPC_CTRL,
  organism = "human",
  pvalueCutoff = 0.05,
  minGSSize = 10,   
  maxGSSize = 1500,
  eps=0,
  verbose = TRUE
)
gsea_IPC_CTRL_df <- as.data.frame(gsea_IC_CTRL)

#Data visualization for overlapping GSEA - selecting top 20 pathways based on adj p value for the LRTI/IPC comparison, then overlaying IPC/CTRL comparison
top20_LRTI_IPC <- gsea_LRTI_IPC_df %>%
  arrange(p.adjust) %>%
  slice_head(n = 20) %>%
  dplyr::select(Description, NES, p.adjust) %>%
  mutate(Group = "LRTI vs IPC")

#matching pathways, if they exist, in the IC vs CTRL comparison
IPC_CTRL_matched_pathways <- gsea_IPC_CTRL_df %>%
  filter(Description %in% top20_LRTI_IPC$Description) %>%
  dplyr::select(Description, NES, p.adjust) %>%
  mutate(Group = "IPC vs CTRL")

missing_paths <- setdiff(top20_LRTI_IPC$Description, IPC_CTRL_matched_pathways$Description)
IPC_CTRL_missing_pathways <- data.frame(
  Description = missing_paths,
  NES = 0,
  p.adjust = NA,
  Group = "IPC vs CTRL"
) #if pathways don't exist, defaulting to NES = 0 and P adj to NA for IC vs CTRL (these won't be shown on the plot)
top20_IPC_CTRL <- bind_rows(IPC_CTRL_matched_pathways, IPC_CTRL_missing_pathways)

plot_df <- bind_rows(top20_LRTI_IPC, top20_IPC_CTRL) %>%
  mutate(
    Description = str_wrap(Description, width = 1000),
    Description = fct_reorder(Description, NES[top20_I$Group == "LRTI vs IPC"]),
    size = -log10(p.adjust),
    Group = factor(Group, levels = c("LRTI vs IPC", "IPC vs CTRL"))
  )

overlap_GSEA_plot <- ggplot(plot_df, aes(x = NES, y = Description, color = Group)) +
  geom_segment(aes(x = 0, xend = NES, y = Description, yend = Description), color = "grey50") +
  geom_vline(xintercept = 0, linetype = "solid", color = "grey50", linewidth = 0.5) +   # Vertical line at NES = 0
  geom_point(aes(size = size)) +
  scale_color_manual(
    values = c(
      "LRTI vs IPC" = "#6a51a3",
      "IPC vs CTRL" = "goldenrod2"
    )
  ) +
  scale_size_continuous(
    name = expression(P[adj]),
    range = c(1, 6),
    breaks = c(0, 10, 40),
    labels = c("1.0", "1e-10", "1e-40")
  ) +
  labs(
    x = "NES",
    y = NULL,
    color = NULL
  ) +
  xlim(-3.4, 3.4) +
  theme_bw(base_family = "Arial", base_size = 12) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 12),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.box.just = "left",  
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.spacing.y = unit(0.2, "cm"),
    plot.margin = margin(t = 10, r = 50, b = 20, l = 10),  
    plot.title = element_blank()
  ) +
  guides(
    color = guide_legend(
      title = NULL,
      direction = "horizontal",
      override.aes = list(size = 4),
      order = 1
    ),
    size = guide_legend(
      title.position = "left",
      direction = "horizontal",
      order = 2
    )
  )
overlap_GSEA_plot

####################################
#Figure 4
####################################

#pairwise limma comparisons for viral subgroup, starting limma #1: viral LRTI vs viral CTRL
metadata_viral_LRTI_CTRL <- metadata %>% filter(subgroup == "Viral LRTI" | subgroup == "Bacterial/viral LRTI" | subgroup == "CTRL")
metadata_viral_LRTI_CTRL$group <- factor(metadata_viral_LRTI_CTRL$group, levels=c("CTRL", "LRTI")) #establishing factors/levels for limma
counts_viral_LRTI_CTRL <- counts[, metadata_viral_LRTI_CTRL$Patient] #creating aligned counts table matching the metadata (non-VST transformed since limma will do transformation)

design <- model.matrix( ~group + X_AGE_YEARS + X_SEX,
                        data = metadata_viral_LRTI_CTRL) #generating a model matrix

vwts <- voom(counts_viral_LRTI_CTRL,
             design = design,
             normalize.method = "quantile",
             plot = TRUE) #running limma-voom
vfit <- lmFit(vwts)
vfit <- eBayes(vfit)
summary(decideTests(vfit))
DEresults_viral_LRTI_CTRL <- topTable(
  vfit,
  coef = "groupLRTI",
  sort.by="none",
  number = Inf
)
DEresults_viral_LRTI_CTRL$gene_symbol <- gene_symbols[rownames(DEresults_viral_LRTI_CTRL), "gene_symbol"]

#limma #2 - viral LRTI vs viral IPC
metadata_viral_LRTI_IPC <- metadata %>% filter(subgroup == "Viral LRTI" | subgroup == "Bacterial/viral LRTI" | subgroup == "Viral IPC" | subgroup == "Bacterial/viral IPC")
metadata_viral_LRTI_IPC$group <- factor(metadata_viral_LRTI_IPC$group, levels=c("IPC", "LRTI")) #establishing factors/levels for limma
counts_viral_LRTI_IPC <- counts[, metadata_viral_LRTI_IPC$Patient] #creating aligned counts table matching the metadata (non-VST transformed since limma will do transformation)

design <- model.matrix( ~group + X_AGE_YEARS + X_SEX,
                        data = metadata_viral_LRTI_IPC) #generating a model matrix

vwts <- voom(counts_viral_LRTI_IPC,
             design = design,
             normalize.method = "quantile",
             plot = TRUE) #running limma-voom
vfit <- lmFit(vwts)
vfit <- eBayes(vfit)
summary(decideTests(vfit))
DEresults_viral_LRTI_IPC <- topTable(
  vfit,
  coef = "groupLRTI",
  sort.by="none",
  number = Inf
)
DEresults_viral_LRTI_IPC$gene_symbol <- gene_symbols[rownames(DEresults_viral_LRTI_IPC), "gene_symbol"]

#limma #3 -  viral IPC vs CTRL
metadata_viral_IPC_CTRL <- metadata %>% filter(subgroup == "Viral IPC" | subgroup == "Bacterial/viral IPC" | subgroup == "CTRL")
metadata_viral_IPC_CTRL$group <- factor(metadata_viral_IPC_CTRL$group, levels=c("CTRL", "IPC")) #establishing factors/levels for limma
counts_viral_IPC_CTRL <- counts[, metadata_viral_IPC_CTRL$Patient] #creating aligned counts table matching the metadata (non-VST transformed since limma will do transformation)

design <- model.matrix( ~group + X_AGE_YEARS + X_SEX,
                        data = metadata_viral_IPC_CTRL) #generating a model matrix

vwts <- voom(counts_viral_IPC_CTRL,
             design = design,
             normalize.method = "quantile",
             plot = TRUE) #running limma-voom
vfit <- lmFit(vwts)
vfit <- eBayes(vfit)
summary(decideTests(vfit))
DEresults_viral_IPC_CTRL <- topTable(
  vfit,
  coef = "groupIPC",
  sort.by="none",
  number = Inf
)
DEresults_viral_IPC_CTRL$gene_symbol <- gene_symbols[rownames(DEresults_viral_IPC_CTRL), "gene_symbol"]

#now the same for bacterial subgroups with pairwise limma, starting with limma #1 - bacterial LRTI vs CTRL
metadata_bact_LRTI_CTRL <- metadata %>% filter(subgroup == "Bacterial LRTI" | subgroup == "Bacterial/viral LRTI" | subgroup == "CTRL")
metadata_bact_LRTI_CTRL$group <- factor(metadata_bact_LRTI_CTRL$group, levels=c("CTRL", "LRTI")) #establishing factors/levels for limma
counts_bact_LRTI_CTRL <- counts[, metadata_bact_LRTI_CTRL$Patient] #creating aligned counts table matching the metadata (non-VST transformed since limma will do transformation)

design <- model.matrix( ~group + X_AGE_YEARS + X_SEX,
                        data = metadata_bact_LRTI_CTRL) #generating a model matrix

vwts <- voom(counts_bact_LRTI_CTRL,
             design = design,
             normalize.method = "quantile",
             plot = TRUE) #running limma-voom
vfit <- lmFit(vwts)
vfit <- eBayes(vfit)
summary(decideTests(vfit))
DEresults_bact_LRTI_CTRL <- topTable(
  vfit,
  coef = "groupLRTI",
  sort.by="none",
  number = Inf
)
DEresults_bact_LRTI_CTRL$gene_symbol <- gene_symbols[rownames(DEresults_bact_LRTI_CTRL), "gene_symbol"]

#bacterial subgroup limma #2 - Bacterial LRTI vs bacterial IPC
metadata_bact_LRTI_IPC <- metadata %>% filter(subgroup == "Bacterial LRTI" | subgroup == "Bacterial/viral LRTI" | subgroup == "Bacterial IPC" | subgroup == "Bacterial/viral IPC")
metadata_bact_LRTI_IPC$group <- factor(metadata_bact_LRTI_IPC$group, levels=c("IPC", "LRTI")) #establishing factors/levels for limma
counts_bact_LRTI_IPC <- counts[, metadata_bact_LRTI_IPC$Patient] #creating aligned counts table matching the metadata (non-VST transformed since limma will do transformation)

design <- model.matrix( ~group + X_AGE_YEARS + X_SEX,
                        data = metadata_bact_LRTI_IPC) #generating a model matrix

vwts <- voom(counts_bact_LRTI_IPC,
             design = design,
             normalize.method = "quantile",
             plot = TRUE) #running limma-voom
vfit <- lmFit(vwts)
vfit <- eBayes(vfit)
summary(decideTests(vfit))
DEresults_bact_LRTI_IPC <- topTable(
  vfit,
  coef = "groupLRTI",
  sort.by="none",
  number = Inf
)
DEresults_bact_LRTI_IPC$gene_symbol <- gene_symbols[rownames(DEresults_bact_LRTI_IPC), "gene_symbol"]

#bacterial subgroup limma #3 - bacterial IPC vs CTRL
metadata_bact_IPC_CTRL <- metadata %>% filter(subgroup == "Bacterial IPC" | subgroup == "Bacterial/viral IPC" | subgroup == "CTRL")
metadata_bact_IPC_CTRL$group <- factor(metadata_bact_IPC_CTRL$group, levels=c("CTRL", "IPC")) #establishing factors/levels for limma
counts_bact_IPC_CTRL <- counts[, metadata_bact_IPC_CTRL$Patient] #creating aligned counts table matching the metadata (non-VST transformed since limma will do transformation)

design <- model.matrix( ~group + X_AGE_YEARS + X_SEX,
                        data = metadata_bact_IPC_CTRL) #generating a model matrix

vwts <- voom(counts_bact_IPC_CTRL,
             design = design,
             normalize.method = "quantile",
             plot = TRUE) #running limma-voom
vfit <- lmFit(vwts)
vfit <- eBayes(vfit)
summary(decideTests(vfit))
DEresults_bact_IPC_CTRL <- topTable(
  vfit,
  coef = "groupIPC",
  sort.by="none",
  number = Inf
)
DEresults_bact_IPC_CTRL$gene_symbol <- gene_symbols[rownames(DEresults_bact_IPC_CTRL), "gene_symbol"]

#volcano plot of viral IPC vs CTRL (figure 4b)
IFN_genes <- c(
  "SP140", "IFIH1", "OAS3", "RSAD2", "DDX60", "IFI44", "HERC6", "HERC5",
  "NLRC5", "GBP4", "MX1", "IFIT5", "DHX58", "ZBP1", "XAF1", "IFI44L",
  "DDX58", "IFIT2", "EPSTI1", "IFIT3", "PARP14", "IFIT1", "ISG15"
) #manually annotating the DE genes that are interferon-related genes per literature review

DEresults_viral_IPC_CTRL <- DEresults_viral_IPC_CTRL %>%
  mutate(
    gene_type = case_when(
      logFC > 0 & adj.P.Val <= 0.05 ~ "↑ in Viral IPC",
      logFC <= 0 & adj.P.Val <= 0.05 ~ "↑ in CTRL",
      TRUE ~ "Not significant"
    ),
    IFN_related = ifelse(gene_symbol %in% IFN_genes, "IFN-related", "Non-IFN")
  ) # Annotate DE results

DEresults_viral_IPC_CTRL$gene_type <- factor(
  DEresults_viral_IPC_CTRL$gene_type,
  levels = c("↑ in Viral IPC", "↑ in CTRL", "Not significant")
) # Set factor levels

DEresults_viral_IPC_CTRL$IFN_related <- factor(
  DEresults_viral_IPC_CTRL$IFN_related,
  levels = c("IFN-related", "Non-IFN")
) #factor for legend purposes

shape_map <- c("IFN-related" = 24, "Non-IFN" = 21)
cols <- c(
  "↑ in Viral IPC" = "#70A494FF",
  "↑ in CTRL" = "darkgrey",
  "Not significant" = "black"
) # aesthetic maps

xlim_vals <- c(-4.1, 4.1)
ylim_vals <- c(0, 2.9)

volcano_viral_IPC_CTRL <- ggplot(DEresults_viral_IPC_CTRL, aes(x = logFC, y = -log10(adj.P.Val))) + # code for volcano plot
  geom_point(
    aes(color = gene_type, fill = gene_type, shape = IFN_related),
    size = 2,
    stroke = 0.6
  ) +
  geom_point(  
    data = data.frame(x = as.numeric(NA), y = as.numeric(NA)),
    aes(x = x, y = y, shape = "IFN-related", color = "↑ in Viral IPC", fill = "↑ in Viral IPC"),
    size = 2,
    stroke = 0.6,
    inherit.aes = FALSE
  ) +   # Dummy point for IFN-related gene legend
  scale_color_manual(
    values = cols,
    breaks = c("↑ in Viral IPC", "↑ in CTRL"),
    name = NULL
  ) +
  scale_fill_manual(
    values = cols,
    guide = "none"
  ) +
  scale_shape_manual(
    values = shape_map,
    breaks = "IFN-related", 
    name = NULL
  ) +  # Only show IFN-related in legend
  guides(
    color = guide_legend(order = 1),
    shape = guide_legend(
      override.aes = list(
        shape = 24,
        fill = "#70A494FF",
        color = "#70A494FF",
        size = 2,
        stroke = 0.6
      ),
      order = 2
    )
  ) +
  scale_x_continuous(limits = xlim_vals, expand = c(0, 0)) +
  scale_y_continuous(limits = ylim_vals, expand = c(0, 0)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  labs(
    x = "log2(fold change)",
    y = "-log10(adjusted P-value)",
    title = "Viral IPC vs CTRL"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 12),
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 12),
    text = element_text(size = 12, family = "Arial"),
    panel.grid = element_blank(),
    axis.line.x = element_line(color = "black", linewidth = 0.4),
    axis.line.y = element_line(color = "black", linewidth = 0.4),
    legend.position = c(0.05, 0.95),
    legend.justification = c(0, 1),
    legend.direction = "vertical",
    legend.text = element_text(size = 12),
    legend.spacing.x = unit(0.4, 'cm'),
    legend.spacing.y = unit(0.01, "cm")
  )
volcano_viral_IPC_CTRL

#canonical viral gene analyses (Figure 3c, 3d, 3e, 3f, 3g, and s7)
gene_symbols_vector <- gene_symbols$gene_symbol
names(gene_symbols_vector) <- rownames(gene_symbols)
isg_ensembl <- names(gene_symbols_vector)[gene_symbols_vector == "IFIH1"] #can substitute other interferon gene (i.e. ISG15, IFI44 shown in supplement) 
isg_expr <- counts_VST[isg_ensembl, ] #obtaining gene counts for that specific ISG

isg_df <- data.frame(
  Patient = names(isg_expr),
  ISG = as.numeric(isg_expr)
) #compiling a dataframe with just the ISG values

isg_df <- isg_df %>%
  left_join(metadata, by = "Patient") %>%
  filter(subgroup %in% c("Viral LRTI", "Bacterial/viral LRTI", "Viral IPC", "Bacterial/viral IPC", "CTRL")) %>% #filtering just for viral cases + CTRLs
  mutate(
    group = factor(group,
                   levels = c("LRTI", "IPC", "CTRL"),
                   labels = c("Viral LRTI", "Viral IPC", "CTRL")) #relabeling for plotting purposes
  ) #merging ISG expression df and the clinical metadata

isg_df <- isg_df %>% filter(is.na(sum_nt_viral_rpm) == FALSE) #filtering the dataframe for the viral subjects who had a virus detected on mNGS (and thus have a viral load for regression analyses
isg_df$group <- factor(isg_df$group, levels = c("Viral IPC", "Viral LRTI"))

group_colors <- c(
  "Viral LRTI" = "#DE8A5AFF",       # Orange
  "Viral IPC" = "#70A494FF"  # Teal
)

#first, looking at regression of viral load and ISG expression (fig 4d)
#linear model of ISG15 by viral load, incorporating an interaction term. for sensitivity analyses for age, just add + X_AGE_YEARS to the model line
model_interaction <- lm(ISG ~ log10(sum_nt_viral_rpm) * group, data = isg_df)
summary(model_interaction)
adj_r2_int <- summary(model_interaction)$adj.r.squared #extracting the R squared value for plotting purposes

#anova testing whether the lines for viral LRTI and viral IPC are different (based on slope/intercept) using F-test
model_null <- lm(ISG ~ log10(sum_nt_viral_rpm), data = isg_df)
anova_result <- anova(model_null, model_interaction)
anova_result
f_test_p <- anova_result$`Pr(>F)`[2] #extracting the F-test p value for plotting

#plotting label
annot_label <- sprintf("F-test P = %.1e\nAdj R\u00B2 = %.2f",
                       f_test_p, adj_r2_int)

#obtaining model-based predictions for plotting purposes (plotting fit line)
xgrid <- tibble(sum_nt_viral_rpm = exp(seq(log(min(isg_df$sum_nt_viral_rpm)),
                                     log(max(isg_df$sum_nt_viral_rpm)),
                                     length.out = 200)))
newdata <- expand_grid(
  group = levels(isg_df$group),
  sum_nt_viral_rpm = xgrid$sum_nt_viral_rpm
)
pred <- cbind(
  newdata,
  as.data.frame(predict(model_interaction, newdata = newdata, interval = "confidence"))
)

#plotting the regression for viral RPM by ISG expression
ISG_viralrpm_plot <- ggplot(isg_df, aes(x = sum_nt_viral_rpm, y = ISG, color = group)) +
  geom_point(size = 2, alpha = 0.8) +
  # Model-implied confidence ribbons & lines
  geom_ribbon(data = pred,
              aes(x = sum_nt_viral_rpm, ymin = lwr, ymax = upr, fill = group),
              alpha = 0.15, color = NA, inherit.aes = FALSE) +
  geom_line(data = pred,
            aes(x = sum_nt_viral_rpm, y = fit, color = group),
            linewidth = 0.9, inherit.aes = FALSE) +
  scale_x_log10() +
  scale_color_manual(values = group_colors) +
  scale_fill_manual(values = group_colors) +
  scale_y_continuous(limits = c(1, 11), breaks = seq(2, 10, by = 2)) +
  labs(
    x = "∑ Viral pathogen RPM",
    y = expression(italic("ISGxx") * " expression"), #change labels depending on the gene you are using
    title = expression(italic("ISGxx") * " expression by viral abundance"), #change labels depending on the gene you are using
    color = NULL, fill = NULL
  ) +
  annotate("text", x = min(isg_df$sum_nt_viral_rpm)*1.05, y = 10,
           label = annot_label, hjust = 0, size = 4, family = "Arial") +
  theme_minimal(base_family = "Arial") +
  theme(
    axis.ticks = element_line(color = "black", linewidth = 0.4),
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 12),
    legend.position = c(0.98, 0.02),
    legend.justification = c(1, 0),
    legend.direction = "vertical",
    legend.box = "vertical",
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.background = element_blank(),
    plot.title = element_text(size = 12, hjust = 0.5),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line.x = element_line(color = "black", linewidth = 0.4),
    axis.line.y = element_line(color = "black", linewidth = 0.4)
  ) +
  guides(color = guide_legend(reverse = TRUE),
         fill  = guide_legend(reverse = TRUE))
ISG_viralrpm_plot

#next a regression of ISG and SDI (figure 4e)
#first an interaction model to determine if there are significant differences betweeen the viral LRTI vs viral IPC lines
model_interaction <- lm(ISG ~ SDI * group, data = isg_df)
summary(model_interaction)

#anova testing the interaction model - no significant difference between LRTI and IPC (p=0.09)
model_null <- lm(ISG ~ SDI, data = isg_df)
anova_result <- anova(model_null, model_interaction)
anova_result

#since no difference between groups, doing a simple model of ISG ~ SDI and showing that with the R2 value/P value in the plot
model_ISG_SDI <- lm(ISG ~ SDI, data = isg_df)
summary(model_ISG_SDI)
p_val <- summary(model_ISG_SDI)$coefficients["SDI", "Pr(>|t|)"] #extracting the p value for plotting
adj_r2_int <- summary(model_ISG_SDI)$adj.r.squared #extracting the R squared value for plotting

annot_label <- sprintf("P = %.1e\nAdj R\u00B2 = %.2f",
                       p_val, adj_r2_int) #plotting label

#building a prediction grid to plot model output
xgrid <- tibble(SDI = seq(min(isg_df$SDI, na.rm = TRUE),
                                  max(isg_df$SDI, na.rm = TRUE),
                                  length.out = 200))
newdata <- expand_grid(
  group = levels(isg_df$group),
  SDI = xgrid$SDI
)
pred <- cbind(
  newdata,
  as.data.frame(predict(model_ISG_SDI, newdata = newdata, interval = "confidence"))
)

ISG_SDI_plot <- ggplot(isg_df, aes(x = SDI, y = ISG)) +
  geom_point(size = 2, alpha = 0.8) +
  geom_ribbon(data = pred,
              aes(x = SDI, ymin = lwr, ymax = upr),
              alpha = 0.15, fill = "darkgrey", inherit.aes = FALSE) +
  geom_line(data = pred, aes(x = SDI, y = fit),
            color = "darkgrey", linewidth = 0.9, inherit.aes = FALSE) +
  scale_y_continuous(limits = c(3, 10), breaks = seq(2, 10, by = 2)) +
  labs(
    x = "Shannon Diversity Index (SDI)",
    y = expression(italic("ISGxx") * " expression"),  # annotate with whichever ISG you decided to use
    title = expression(italic("ISGxx") * " expression by SDI") # annotate with whichever ISG you decided to use
  ) +
  annotate("text",
           x = 0.65*max(isg_df$SDI, na.rm = TRUE) +
             0.05*diff(range(isg_df$SDI, na.rm = TRUE)),
           y = 9.0,
           label = annot_label, hjust = 0, size = 4, family = "Arial") +
  theme_minimal(base_family = "Arial") +
  theme(
    axis.ticks = element_line(color = "black", linewidth = 0.4),
    axis.text  = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 12),
    plot.title = element_text(size = 12, hjust = 0.5),
    legend.position   = "none",
    panel.grid = element_blank(),
    axis.line.x = element_line(color = "black", linewidth = 0.4),
    axis.line.y = element_line(color = "black", linewidth = 0.4)
  )
ISG_SDI_plot

#now plotting ISG expression by viral RPM (same as 4d) but adjusting for SDI (figure 4f)
isg_df$ISG_adj_div <- resid(
  lm(ISG ~ SDI, data = isg_df)
)

model_interaction <- lm(ISG_adj_div ~ log10(sum_nt_viral_rpm) * group, 
                        data = isg_df)
summary(model_interaction)
adj_r2_int <- summary(model_interaction)$adj.r.squared #extracting the R squared value for plotting

#anova testing the interaction model
model_null <- lm(ISG_adj_div ~ log10(sum_nt_viral_rpm), data = isg_df)
anova_result <- anova(model_null, model_interaction)
anova_result
f_test_p <- anova_result$`Pr(>F)`[2]

annot_label <- sprintf(
  "F-test P = %.2f\nAdj R\u00B2 = %.2f",
  f_test_p, adj_r2_int
)

#generating a prediction grid using the interaction model, to make fit lines
xmin <- min(isg_df$sum_nt_viral_rpm[isg_df$sum_nt_viral_rpm > 0], na.rm = TRUE)
xmax <- max(isg_df$sum_nt_viral_rpm, na.rm = TRUE)
xgrid <- tibble(sum_nt_viral_rpm = exp(seq(log(xmin), log(xmax), length.out = 200)))
newdata <- expand_grid(
  group = levels(isg_df$group),
  sum_nt_viral_rpm = xgrid$sum_nt_viral_rpm
)
pred <- cbind(
  newdata,
  as.data.frame(predict(model_interaction, newdata = newdata, interval = "confidence"))
)

#plotting
ISG_residuals_plot <- ggplot(isg_df, aes(x = sum_nt_viral_rpm, y = ISG_adj_div, color = group)) +
  geom_point(size = 2, alpha = 0.8) +
  geom_ribbon(data = pred,
              aes(x = sum_nt_viral_rpm, ymin = lwr, ymax = upr, fill = group),
              alpha = 0.15, color = NA, inherit.aes = FALSE) +
  geom_line(data = pred,
            aes(x = sum_nt_viral_rpm, y = fit, color = group),
            linewidth = 0.9, inherit.aes = FALSE) +
  scale_x_log10() +
  scale_color_manual(values = group_colors) +
  scale_fill_manual(values = group_colors) +
  scale_y_continuous(limits = c(-4, 3), breaks = seq(-4, 2, by = 2)) +
  labs(
    x = "∑ Viral pathogen RPM",
    y = expression(italic("ISGxx") * " residuals (adjusted for SDI)"), #insert whichever ISG you are looking at 
    title = expression(italic("ISGxx") * " residuals vs viral abundance (SDI-adjusted)"), #insert whichever ISG you are looking at 
    color = NULL, fill = NULL
  ) +
  annotate("text",
           x = xmin * 1.05,
           y = 2.5,
           label = annot_label,
           hjust = 0, size = 4, family = "Arial") +
  theme_minimal(base_family = "Arial") +
  theme(
    axis.ticks        = element_line(color = "black", linewidth = 0.4),
    axis.text         = element_text(size = 12, color = "black"),
    axis.title        = element_text(size = 12),
    legend.position   = c(0.98, 0.02),
    legend.justification = c(1, 0),
    legend.direction  = "vertical",
    legend.box        = "vertical",
    legend.text       = element_text(size = 12),
    legend.title      = element_text(size = 12),
    legend.background = element_blank(),
    plot.title        = element_text(size = 12, hjust = 0.5),
    panel.grid        = element_blank(),
    axis.line.x       = element_line(color = "black", linewidth = 0.4),
    axis.line.y       = element_line(color = "black", linewidth = 0.4)
  ) +
  guides(color = guide_legend(reverse = TRUE),
         fill  = guide_legend(reverse = TRUE))
ISG_residuals_plot

# Quasi-mediation analysis for viral LRTI vs IPC
set.seed(0)
model_mediator <- lm(SDI ~ group + log10(sum_nt_viral_rpm), data = isg_df)
summary(model_mediator)

# Outcome model: ISG15 ~ group + diversity + viral load
model_outcome <- lm(ISG ~ group + SDI + log10(sum_nt_viral_rpm), data = isg_df)
summary(model_outcome)

med_result <- mediate(
  model.m = model_mediator,
  model.y = model_outcome,
  treat = "group",
  mediator = "SDI",
  sims = 1000  # Number of bootstrap simulations
)
summary(med_result)

#analogous regression analysis for the bacterial subgroups (bacterial LRTI vs IPC)
gene_symbols_vector <- gene_symbols$gene_symbol
names(gene_symbols_vector) <- rownames(gene_symbols)
bact_ensembl <- names(gene_symbols_vector)[gene_symbols_vector == "GZMB"] #can adjust for different canonical anti-bacterial genes
bact_expr <- counts_VST[bact_ensembl, ]

bact_df <- data.frame(
  Patient = names(bact_expr),
  bact_gene = as.numeric(bact_expr)
)

bact_df <- bact_df %>%
  left_join(metadata, by = "Patient") %>%
  filter(subgroup %in% c("Bacterial LRTI", "Bacterial/viral LRTI", "Bacterial IPC", "Bacterial/viral IPC", "CTRL")) %>%
  mutate(
    group = factor(group,
                   levels = c("LRTI", "IPC", "CTRL"),
                   labels = c("Bacterial LRTI", "Bacterial IPC", "CTRL"))
  )

bact_df <- bact_df %>% filter(is.na(sum_nt_bacterial_rpm) == FALSE) #filtering the dataframe for the bacterial subjects who had a bacterial detected on mNGS RBM (and thus have a bacterial pathogen VL for regression analyses)
bact_df$group <- factor(bact_df$group, levels = c("Bacterial IPC", "Bacterial LRTI"))

group_colors <- c(
  "Bacterial LRTI" = "#DE8A5AFF",       # Orange
  "Bacterial IPC" = "#70A494FF"  # Teal
)

#generating an interaction model to plot bacterial RPM by bacterial gene expression (for sensitivity analyses, can add X_AGE_YEARS to end of model)
model_interaction <- lm(bact_gene ~ log10(sum_nt_bacterial_rpm) * group, data = bact_df)
summary(model_interaction)
adj_r2_int <- summary(model_interaction)$adj.r.squared #extracting the R squared value from the interaction model since that's what we're plotting

#anova testing the interaction model
model_null <- lm(bact_gene ~ log10(sum_nt_bacterial_rpm), data = bact_df)
anova_result <- anova(model_null, model_interaction)
anova_result
f_test_p <- anova_result$`Pr(>F)`[2] #storing F-test p value for plotting

#obtaining model-based predictions for plotting purposes
xgrid <- tibble(sum_nt_bacterial_rpm = exp(seq(log(min(bact_df$sum_nt_bacterial_rpm)),
                                     log(max(bact_df$sum_nt_bacterial_rpm)),
                                     length.out = 200)))
newdata <- expand_grid(
  group = levels(bact_df$group),
  sum_nt_bacterial_rpm = xgrid$sum_nt_bacterial_rpm
)
pred <- cbind(
  newdata,
  as.data.frame(predict(model_interaction, newdata = newdata, interval = "confidence"))
)

annot_label <- sprintf("F-test P = %.1e\nAdj R\u00B2 = %.2f",
                       f_test_p, adj_r2_int) #for labeling plot

bactgene_bactrpm_plot <- ggplot(bact_df, aes(x = sum_nt_bacterial_rpm, y = bact_gene, color = group)) +
  geom_point(size = 2, alpha = 0.8) +
  # Model-implied confidence ribbons & lines
  geom_ribbon(data = pred,
              aes(x = sum_nt_bacterial_rpm, ymin = lwr, ymax = upr, fill = group),
              alpha = 0.15, color = NA, inherit.aes = FALSE) +
  geom_line(data = pred,
            aes(x = sum_nt_bacterial_rpm, y = fit, color = group),
            linewidth = 0.9, inherit.aes = FALSE) +
  scale_x_log10(
    limits = c(1, max(bact_df$sum_nt_bacterial_rpm, na.rm = TRUE)),
    labels = scales::label_number(accuracy = 1, big.mark = "", decimal.mark = ".")
  ) +
  scale_color_manual(values = group_colors) +
  scale_fill_manual(values = group_colors) +
  scale_y_continuous(limits = c(0.5, 11), breaks = seq(2, 10, by = 2)) +
  labs(
    x = "∑ Bacterial pathogen RPM",
    y = expression(italic("Bacterial gene xx") * " expression"), #insert the bacterial gene you are using here
    title = expression(italic("Bacterial gene xx") * " expression by bacterial abundance"), #insert the bacterial gene you are using here
    color = NULL, fill = NULL
  ) +
  annotate("text", x = 1, y = 10,
           label = annot_label, hjust = 0, size = 4, family = "Arial") +
  theme_minimal(base_family = "Arial") +
  theme(
    axis.ticks = element_line(color = "black", linewidth = 0.4),
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 12),
    legend.position = c(0.98, 0.02),
    legend.justification = c(1, 0),
    legend.direction = "vertical",
    legend.box = "vertical",
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.background = element_blank(),
    plot.title = element_text(size = 12, hjust = 0.5),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line.x = element_line(color = "black", linewidth = 0.4),
    axis.line.y = element_line(color = "black", linewidth = 0.4)
  ) +
  guides(color = guide_legend(reverse = TRUE),
         fill  = guide_legend(reverse = TRUE))
bactgene_bactrpm_plot


#mediation analysis for bacterial gene - does diversity mediate the group difference seen in bacterial gene expression
set.seed(0)
model_mediator2 <- lm(SDI ~ group + log(sum_nt_bacterial_rpm), data = bact_df)
summary(model_mediator2)

# Outcome model: bact_gene ~ group + diversity + viral load
model_outcome2 <- lm(bact_gene ~ group + SDI + log(sum_nt_bacterial_rpm), data = bact_df)
summary(model_outcome2)

med_result2 <- mediate(
  model.m = model_mediator2,
  model.y = model_outcome2,
  treat = "group",
  mediator = "SDI",
  sims = 1000  # Number of bootstrap simulations
)
summary(med_result2)

