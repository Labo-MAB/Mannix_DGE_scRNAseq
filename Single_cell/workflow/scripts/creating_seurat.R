library(Seurat)
library(Matrix)

# ============================================================
# 1. Chemins
# ============================================================

data_dir <- "data/single_cell"
output_dir <- "results/single_cell"

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

mtx_file <- file.path(data_dir, "GSE304446_matrix.mtx.gz")
features_file <- file.path(data_dir, "GSE304446_features.tsv.gz")
barcodes_file <- file.path(data_dir, "GSE304446_barcodes.tsv.gz")
metadata_file <- file.path(data_dir, "GSE304446_cell_metadata.tsv.gz")

output_file <- file.path(output_dir, "GSE304446_seurat.rds")

# ============================================================
# 2. Lire la matrice sparse
# ============================================================

counts <- readMM(mtx_file)

# ============================================================
# 3. Lire features et barcodes
# ============================================================

features <- read.delim(
  features_file,
  header = FALSE,
  stringsAsFactors = FALSE
)

barcodes <- read.delim(
  barcodes_file,
  header = FALSE,
  stringsAsFactors = FALSE
)

gene_names <- features[, 1]
cell_names <- barcodes[, 1]

# Vérifications
stopifnot(nrow(counts) == length(gene_names))
stopifnot(ncol(counts) == length(cell_names))

rownames(counts) <- make.unique(gene_names)
colnames(counts) <- cell_names

cat("Dimensions de la matrice :", dim(counts), "\n")

# ============================================================
# 4. Lire les metadata
# ============================================================

metadata <- read.delim(
  metadata_file,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

rownames(metadata) <- metadata$barcode

# Vérifier que tous les noyaux sont présents
matched <- sum(colnames(counts) %in% rownames(metadata))

cat(
  "Barcodes retrouvés dans metadata :",
  matched,
  "/",
  ncol(counts),
  "\n"
)

stopifnot(matched == ncol(counts))

# Remettre metadata dans exactement le même ordre que la matrice
metadata <- metadata[colnames(counts), , drop = FALSE]

# ============================================================
# 5. Créer objet Seurat
# ============================================================

seurat_obj <- CreateSeuratObject(
  counts = counts,
  meta.data = metadata,
  project = "GSE304446",
  min.cells = 0,
  min.features = 0
)

# ============================================================
# 6. Afficher quelques infos
# ============================================================

print(seurat_obj)

cat("\nCell types :\n")
print(table(seurat_obj$CellType))

cat("\nTreatments :\n")
print(table(seurat_obj$Treatment))

cat("\nBackgrounds :\n")
print(table(seurat_obj$Background))

cat("\nSamples :\n")
print(table(seurat_obj$SampleID))

# ============================================================
# 7. Sauvegarder
# ============================================================

saveRDS(
  seurat_obj,
  file = output_file,
  compress = FALSE
)

cat("\nObjet Seurat sauvegardé dans :", output_file, "\n")
