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


# Load data --------------------------------------------------------------------
gene_E_S <- read_xlsx(file.path("Supplementary_Table 1-4.xlsx"), sheet = "Table S3.2")
  
colnames(gene_E_S) <- gene_E_S[2, ]
gene_E_S <- gene_E_S[-c(1,2), ]  
gene_E <- gene_E_S$`Epithelioid markers`
gene_S <- gene_E_S$`Sarcomatoid markers`

MPM_merge <- readRDS(file.path("obj_scRNAseq_cell_lines.rds"))

MPM_merge$epi_score <- apply(GetAssayData(object=MPM_merge, assay = "RNA", slot="data")[gene_E,], 2, mean)
MPM_merge$sarco_score <- apply(GetAssayData(object=MPM_merge, assay = "RNA", slot="data")[gene_S,], 2, mean)

cell_E_S_scores_df <- MPM_merge@meta.data %>%
  select(epi_mean = epi_score,
         sarco_mean = sarco_score,
         nCount_RNA,
         nFeature_RNA,
         percent.mt,
         orig.ident)
    
# Create table for normalize each Emean and Smean by the maxEmean and max Smean of cells
maxEmean_cell <- cell_E_S_scores_df %>% pull(epi_mean) %>% max()
maxSmean_cell <- cell_E_S_scores_df %>% pull(sarco_mean) %>% max()

# For each cells
cell_E_S_scores_df <- cell_E_S_scores_df %>% 
  mutate(epi_mean_normalize = epi_mean/maxEmean_cell,
         sarco_mean_normalize = sarco_mean/maxSmean_cell)

color_cell_lines = c("MPM05" = "#a6cee3", "MPM07" = "#1f78b4", "MPM12" = "#b2df8a", "MPM27" = "#33a02c", "MPM28" = "#FF5B6F",
                     "MPM37" = "#e31a1c", "MPM59" = "#5EBC88", "MPM60" = "#ff7f00", "MPM70" = "#cab2d6", "MPM86" = "#6a3d9a",
                     "MPM87" = "#FADF00", "MPM78" = "#b15928", "MPM16" = "#4230E2", "MPM82" = "#ff9896", "MPM83" = "#17becf")

level_color = names(color_cell_lines)

E_Sscore_by_cells <- cell_E_S_scores_df %>% 
  .[sample(nrow(.)), ] %>%
  filter(!orig.ident %in% c("MPM27", "MPM82", "MPM60", "MPM16")) %>% 
  ggplot(aes(x = epi_mean_normalize, y = sarco_mean_normalize, color = orig.ident)) +
  geom_point(size = 1) +
  scale_y_continuous(limits = c(0,1), breaks = seq(0, 2, 0.1)) +
  scale_x_continuous(limits = c(0,1), breaks = seq(0, 2, 0.1)) +
  scale_color_manual(limits = level_color, values = color_cell_lines , "cell_line") +
  theme_bw() +
  geom_vline(xintercept = 0.15, linetype="dashed") +
  geom_hline(yintercept = 0.25, linetype="dashed") +
  labs(x = "sc-Epithelioid score", y = "sc-Sarcomatoid score") +
  theme(axis.text.x = element_text(size = 12, color = "black"),
        axis.title.x = element_text(size = 20), 
        axis.text.y = element_text(size = 12, color = "black"),
        axis.title.y = element_text(size = 20),
        legend.position = "none")
E_Sscore_by_cells

png(file.path("E_Sscore_by_cells.png"), width = 612, height = 600, res = 120)
print(E_Sscore_by_cells)
dev.off()

svg(file.path("E_Sscore_by_cells.svg"), width = 612/120, height = 600/120)
print(E_Sscore_by_cells)
dev.off()

