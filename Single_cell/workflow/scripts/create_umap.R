library(Seurat)

# ============================================================
# 1. ARGUMENTS
# ============================================================

args <- commandArgs(trailingOnly = TRUE)

input_file <- args[1]
output_file <- args[2]

# ============================================================
# 2. CHARGER OBJET SEURAT
# ============================================================

obj <- readRDS(input_file)

cat(
  "\nObjet chargé :",
  nrow(obj), "gènes x",
  ncol(obj), "noyaux\n"
)

# ============================================================
# 3. NORMALISATION
# ============================================================

obj <- NormalizeData(
  obj,
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = FALSE
)

# ============================================================
# 4. GENES VARIABLES
# ============================================================

obj <- FindVariableFeatures(
  obj,
  selection.method = "vst",
  nfeatures = 2000,
  verbose = FALSE
)

# ============================================================
# 5. SCALE
# ============================================================

obj <- ScaleData(
  obj,
  features = VariableFeatures(obj),
  verbose = FALSE
)

# ============================================================
# 6. PCA
# ============================================================

obj <- RunPCA(
  obj,
  features = VariableFeatures(obj),
  npcs = 50,
  verbose = FALSE
)

# ============================================================
# 7. UMAP
# ============================================================

obj <- RunUMAP(
  obj,
  dims = 1:30,
  reduction = "pca",
  verbose = FALSE
)

# ============================================================
# 8. SAUVEGARDE
# ============================================================

saveRDS(
  obj,
  file = output_file,
  compress = FALSE
)

cat(
  "\nObjet Seurat avec UMAP sauvegardé dans :",
  output_file,
  "\n"
)
