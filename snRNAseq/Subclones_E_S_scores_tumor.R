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
         orig.ident,
         final_clone)

# Create table for normalize each Emean and Smean by the maxEmean of cells
maxEmean_cell <- cell_E_S_scores_df %>% pull(epi_mean) %>% max()
maxSmean_cell <- cell_E_S_scores_df %>% pull(sarco_mean) %>% max()

# For each cells
cell_E_S_scores_df <- cell_E_S_scores_df %>% 
  mutate(epi_mean_normalize = epi_mean/maxEmean_cell,
         sarco_mean_normalize = sarco_mean/maxSmean_cell)

# Subclones distribution -------------------------------------

name_tumor <- c("T267", "T278_C", "T343_C", "T043", "T093", "T094", "T201_B", "T038", 
                "T111NE", "T161NE", "T277HP-A", "T325HP", "T227LE-A", "T255HP", "T265HP")[1]
  
cell_E_S_scores_df <- cell_E_S_scores_df %>% 
  filter(orig.ident %in% name_tumor)

unique(cell_E_S_scores_df$final_clone)

# Couleur clones
"#ca0020"
"#f4a582"
"black"
"#33a02c"
"#92c5de"
"#0571b0"
"#6a3d9a"
"#fdbf6f"
"#2f5597"

# Depending of each tumors
color_clones <- c("clonal" = "#2f5597", "subclone_1" = "#ca0020", "subclonal_2" = "#fdbf6f")
levels_clones <- names(color_clones) 


E_Sscore_by_cells <- cell_E_S_scores_df %>%
  filter(!is.na(final_clone)) %>% 
  .[sample(nrow(.)), ] %>% 
  ggplot(aes(x = epi_mean_normalize, y = sarco_mean_normalize, color = final_clone)) +
  geom_point(size = 1) +
  scale_y_continuous(limits = c(0,1), breaks = seq(0, 2, 0.1)) +
  scale_x_continuous(limits = c(0,1), breaks = seq(0, 2, 0.1)) +
  scale_color_manual(limits = levels_clones, values = color_clones, "clones") +
  theme_bw() +
  labs(x = "sc-Epithelioid score", y = "sc-Sarcomatoid score") +
  theme(axis.text.x = element_text(size = 12, color = "black"),
        axis.title.x = element_text(size = 20), 
        axis.text.y = element_text(size = 12, color = "black"),
        axis.title.y = element_text(size = 20),
        legend.position = "none")
E_Sscore_by_cells

png(file.path(glue::glue("E_Sscore_cells_final_clones_{name_tumor}.png")), width = 610, height = 600, res = 120)
print(E_Sscore_by_cells)
dev.off()

svg(file.path(glue::glue("E_Sscore_cells_final_clones_{name_tumor}.svg")), width = 610/120, height = 600/120)
print(E_Sscore_by_cells)
dev.off()


# Depending of the tumor
color_clones <- c("clonal" = "#2f5597", "subclonal_2" = "#fdbf6f", "subclone_1" = "#ca0020")
levels_clones <- names(color_clones) 

# Rigde plot S
ridgeplot <- cell_E_S_scores_df %>%
  mutate(sarco_mean_normalize = -sarco_mean_normalize) %>%
  filter(!is.na(final_clone)) %>% 
  ggplot(aes(x = sarco_mean_normalize, y = final_clone,  fill = final_clone)) +
  ggridges::geom_density_ridges(scale = 5, size = 0.5, alpha = 0.8) +
  geom_segment(data = . %>% distinct(final_clone),
               aes(x = -0.87,
                   xend = -1, y = final_clone, yend = final_clone),
               inherit.aes = FALSE, color = "black", size = 0.5) +
  geom_segment(data = . %>% distinct(final_clone),
               aes(x = -0.054,
                   xend = 0, y = final_clone, yend = final_clone),
               inherit.aes = FALSE, color = "black", size = 0.5) +
  scale_x_continuous(limits = c(-1, 0)) +
  theme_void() +
  scale_y_discrete(limits = levels_clones) +
  scale_fill_manual(values = color_clones, guide = "none")
ridgeplot

png(file.path(glue::glue("ridge_plot_{name_tumor}_final.png")), width = 700, height = 200, res = 120)
print(ridgeplot)
dev.off()

svg(file.path(glue::glue("ridge_plot_{name_tumor}_final.svg")), width = 700/120, height = 200/120)
print(ridgeplot)
dev.off()

# For T265HP ----------------
# Rigde plot E
ridgeplot <- cell_E_S_scores_df %>%
  filter(!is.na(final_clone)) %>% 
  ggplot(aes(x = epi_mean_normalize, y = final_clone,  fill = final_clone)) +
  ggridges::geom_density_ridges(scale = 5, size = 0.5, alpha = 0.8) +
  geom_segment(data = . %>% distinct(final_clone),
               aes(x = 0.72,
                   xend = 1, y = final_clone, yend = final_clone),
               inherit.aes = FALSE, color = "black", size = 0.5) +
  scale_x_continuous(limits = c(0, 1)) +
  theme_void() +
  scale_y_discrete(limits = levels_clones) +
  scale_fill_manual(values = color_clones, guide = "none")
ridgeplot

png(file.path(glue::glue("ridge_plot_{name_tumor}_E.png")), width = 700, height = 200, res = 120)
print(ridgeplot)
dev.off()

svg(file.path(glue::glue("ridge_plot_{name_tumor}_E.svg")), width = 700/120, height = 200/120)
print(ridgeplot)
dev.off()
