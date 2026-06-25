# library-----------------------------------------------------------------------
library(Seurat)
library(patchwork)
library(ggplot2)
library(readxl)
library(dplyr)

# Function ---------------------------------------------------------------------
LoadRData = function(file.name,
                     print.dim = T,
                     print.real.name = T) {
  file.real.name = load(file.name, verbose = print.real.name)
  if (print.dim) cat("dimensions: ", dim(get(file.real.name)))
  return(get(file.real.name))
}


# load data --------------------------------------------------------------------
gene_E_S <- read_xlsx(file.path("Supplementary_Table 1-4.xlsx"), sheet = "Table S3.2")
colnames(gene_E_S) <- gene_E_S[2, ]
gene_E_S <- gene_E_S[-c(1,2), ]  
gene_E <- gene_E_S$`Epithelioid markers`
gene_S <- gene_E_S$`Sarcomatoid markers`

MPM_merge <- readRDS(file.path("obj_snRNAseq_tumor.rds"))
MPM_merge <- MPM_merge %>% 
  subset(final_type == "tumor")

E_markers_list <- intersect(gene_E, rownames(MPM_merge@assays$RNA))
S_markers_list <- intersect(gene_S, rownames(MPM_merge@assays$RNA))

MPM_merge$epi_score <- apply(GetAssayData(object=MPM_merge, assay = "RNA", slot="data")[E_markers_list,], 2, mean)
MPM_merge$sarco_score <- apply(GetAssayData(object=MPM_merge, assay = "RNA", slot="data")[S_markers_list,], 2, mean)

cell_E_S_scores_df <- MPM_merge@meta.data %>%
  select(epi_mean = epi_score,
         sarco_mean = sarco_score,
         nCount_RNA,
         nFeature_RNA,
         percent.mt,
         orig.ident)

# Create table for normalize each Emean and Smean by the maxEmean of cells
maxEmean_cell <- cell_E_S_scores_df %>% pull(epi_mean) %>% max()
maxSmean_cell <- cell_E_S_scores_df %>% pull(sarco_mean) %>% max()

# For each cells
cell_E_S_scores_df <- cell_E_S_scores_df %>% 
  mutate(epi_mean_normalize = epi_mean/maxEmean_cell,
         sarco_mean_normalize = sarco_mean/maxSmean_cell)


color_tumors <- c("T267" = "#b15928", "T278_C" = "#cab2d6", "T343HP_C" = "#6a3d9a", "T043" = "#1f77b4",
                  "T093" = "#ff7f0e", "T094" = "#2ca02c", "T201_B" = "#d62728", "T038" = "#e377c2", "T161NE" = "#98df8a", 
                  "T277HP_A" = "#bcbd22", "T111NE" = "#ffbb78", "T325HP" = "#aec7e8", "T227LE-A" = "#17becf", "T255HP" = "#7f7f7f", 
                  "T265HP" = "#ff9896")

levels_tumors <- names(color_tumors)

E_Sscore_by_cells <- cell_E_S_scores_df %>% 
  .[sample(nrow(.)), ] %>%
  filter(orig.ident %in% c("T043", "T093", "T094", "T111NE", "T161NE", "T227LE-A", "T255HP")) %>% 
  ggplot(aes(x = epi_mean_normalize, y = sarco_mean_normalize, color = orig.ident)) +
  geom_point(size = 1) +
  scale_y_continuous(limits = c(0,1), breaks = seq(0, 2, 0.1)) +
  scale_x_continuous(limits = c(0,1), breaks = seq(0, 2, 0.1)) +
  scale_color_manual(limits = levels_tumors, values = color_tumors , "tumor") +
  theme_bw() +
  geom_vline(xintercept = 0.15, linetype="dashed") +
  geom_hline(yintercept = 0.25, linetype="dashed") +
  labs(x = "sc-Epithelioid score", y = "sc-Sarcomatoid score") +
  theme(axis.text.x = element_text(size = 12, color = "black"),
        axis.title.x = element_text(size = 20), 
        axis.text.y = element_text(size = 12, color = "black"),
        axis.title.y = element_text(size = 20),
        legend.text = element_blank(),
        legend.position = "none")
E_Sscore_by_cells

png(file.path(image_dir, "E_Sscore_by_cells.png"), width = 610, height = 600, res = 120)
print(E_Sscore_by_cells)
dev.off()


svg(file.path("E_Sscore_by_cells.svg"), width = 612/120, height = 600/120)
print(E_Sscore_by_cells)
dev.off()

E_Sscore_by_cells <- cell_E_S_scores_df %>% 
  .[sample(nrow(.)), ] %>%
  filter(orig.ident == "T038") %>% 
  ggplot(aes(x = epi_mean_normalize, y = sarco_mean_normalize, color = orig.ident)) +
  geom_point(size = 1) +
  scale_y_continuous(limits = c(0,1), breaks = seq(0, 2, 0.1)) +
  scale_x_continuous(limits = c(0,1), breaks = seq(0, 2, 0.1)) +
  scale_color_manual(limits = levels_tumors, values = color_tumors , "tumor") +
  theme_bw() +
  geom_vline(xintercept = 0.15, linetype="dashed") +
  geom_hline(yintercept = 0.25, linetype="dashed") +
  labs(x = "sc-Epithelioid score", y = "sc-Sarcomatoid score") +
  theme(axis.text.x = element_text(size = 12, color = "black"),
        axis.title.x = element_text(size = 20), 
        axis.text.y = element_text(size = 12, color = "black"),
        axis.title.y = element_text(size = 20),
        legend.text = element_blank(),
        legend.position = "none")
E_Sscore_by_cells

svg(file.path("E_Sscore_by_cells_T038.svg"), width = 612/120, height = 600/120)
print(E_Sscore_by_cells)
dev.off()
