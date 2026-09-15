# ============================================================
# Wee1 snRNA-seq pseudobulk
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(DESeq2)
  library(ggplot2)
  library(dplyr)
  library(tibble)
})

# ============================================================
# 1. Arguments Snakemake
# ============================================================

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
  stop(
    "Usage: Rscript Wee1_pseudobulk.R <seurat.rds> <output_dir>"
  )
}

seurat_file <- args[1]
output_dir <- args[2]

dir.create(
  output_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

gene_of_interest <- "Wee1"


# ============================================================
# 2. Charger objet Seurat
# ============================================================

obj <- readRDS(
  seurat_file
)

cat(
  "\n============================================\n",
  "OBJET SEURAT CHARGÉ\n",
  "============================================\n",
  sep = ""
)

cat(
  "Nombre de noyaux : ",
  ncol(obj),
  "\n",
  sep = ""
)

cat(
  "Nombre de gènes : ",
  nrow(obj),
  "\n",
  sep = ""
)


# ============================================================
# 3. Vérifications
# ============================================================

if (!gene_of_interest %in% rownames(obj)) {

  stop(
    paste0(
      gene_of_interest,
      " absent de l'objet Seurat."
    )
  )
}


required_metadata <- c(
  "SampleID",
  "CellType",
  "Treatment",
  "Background"
)

missing_metadata <- setdiff(
  required_metadata,
  colnames(obj@meta.data)
)

if (length(missing_metadata) > 0) {

  stop(
    paste(
      "Metadata manquante :",
      paste(
        missing_metadata,
        collapse = ", "
      )
    )
  )
}


# ============================================================
# 4. Ordre des groupes
# ============================================================

obj$Treatment <- factor(
  obj$Treatment,
  levels = c(
    "Ctrl",
    "HFpEF"
  )
)

obj$Background <- factor(
  obj$Background,
  levels = c(
    "J",
    "N"
  )
)


cat(
  "\n============================================\n",
  "TRAITEMENTS\n",
  "============================================\n",
  sep = ""
)

print(
  table(
    obj$Treatment,
    obj$Background
  )
)


cat(
  "\n============================================\n",
  "CELL TYPES\n",
  "============================================\n",
  sep = ""
)

print(
  table(
    obj$CellType
  )
)


# ============================================================
# 5. Couleurs
# ============================================================

celltype_colors <- c(
  "Cardiomyocyte" = "#E53935",
  "EC"             = "#1565C0",
  "Fibroblast"     = "#F57C00",
  "Mural_Cell"     = "#8E44AD",
  "Mast_Cell"      = "#43A047",
  "Myeloid"        = "#00ACC1",
  "Lymphoid"       = "#1E88E5",
  "Neural_Cell"    = "#C99000"
)

treatment_colors <- c(
  "Ctrl"  = "#F8766D",
  "HFpEF" = "#00BFC4"
)


# ============================================================
# 6. Extraire raw counts Wee1
# ============================================================

raw_wee1 <- FetchData(
  obj,
  vars = gene_of_interest,
  layer = "counts"
)

raw_wee1$barcode <- rownames(
  raw_wee1
)

colnames(raw_wee1)[1] <- "raw_count"


gene_data <- obj@meta.data %>%

  mutate(
    barcode = rownames(.)
  ) %>%

  select(
    barcode,
    SampleID,
    CellType,
    Treatment,
    Background
  ) %>%

  left_join(
    raw_wee1,
    by = "barcode"
  )


# ============================================================
# 7. FIGURE :
#    abondance vs signal Wee1
#    C57BL/6J seulement
# ============================================================

gene_data_J <- gene_data %>%

  filter(
    Background == "J"
  )


summary_signal_J <- gene_data_J %>%

  group_by(
    Treatment,
    CellType
  ) %>%

  summarise(

    n_nuclei = n(),

    total_Wee1_counts =
      sum(
        raw_count,
        na.rm = TRUE
      ),

    .groups = "drop"
  ) %>%

  group_by(
    Treatment
  ) %>%

  mutate(

    percent_of_nuclei =
      100 *
      n_nuclei /
      sum(n_nuclei),

    percent_of_Wee1_signal =
      100 *
      total_Wee1_counts /
      sum(total_Wee1_counts)
  ) %>%

  ungroup()


# ordre selon contribution Wee1

cell_order_J <- summary_signal_J %>%

  group_by(
    CellType
  ) %>%

  summarise(

    average_signal =
      mean(
        percent_of_Wee1_signal
      ),

    .groups = "drop"
  ) %>%

  arrange(
    desc(average_signal)
  ) %>%

  pull(
    CellType
  )


summary_signal_J$CellType <- factor(
  summary_signal_J$CellType,
  levels = rev(cell_order_J)
)


plot_signal_J <- ggplot(

  summary_signal_J,

  aes(
    y = CellType
  )
) +

  geom_segment(

    aes(
      x = percent_of_nuclei,
      xend = percent_of_Wee1_signal,
      yend = CellType
    ),

    color = "grey72",
    linewidth = 1.2
  ) +

  geom_point(

    aes(
      x = percent_of_nuclei
    ),

    color = "grey40",
    size = 4
  ) +

  geom_point(

    aes(
      x = percent_of_Wee1_signal,
      color = CellType
    ),

    size = 5
  ) +

  geom_text(

    aes(
      x = percent_of_nuclei,

      label = paste0(
        round(
          percent_of_nuclei,
          1
        ),
        "%"
      )
    ),

    nudge_y = -0.20,
    color = "grey40",
    size = 4
  ) +

  geom_text(

    aes(
      x = percent_of_Wee1_signal,

      label = paste0(
        round(
          percent_of_Wee1_signal,
          1
        ),
        "%"
      ),

      color = CellType
    ),

    nudge_y = 0.20,
    size = 4.2,
    fontface = "bold"
  ) +

  facet_wrap(
    ~ Treatment,
    ncol = 2
  ) +

  scale_color_manual(
    values = celltype_colors
  ) +

  labs(

    title =
      "Abondance cellulaire vs signal Wee1: C57BL/6J",

    subtitle =
      "Gris = % des noyaux | Couleur = % des raw counts Wee1",

    x =
      "Pourcentage (%)",

    y =
      NULL
  ) +

  theme_classic(
    base_size = 15
  ) +

  theme(

    plot.title =
      element_text(
        face = "bold",
        size = 18
      ),

    axis.text.y =
      element_text(
        face = "bold",
        size = 12
      ),

    strip.text =
      element_text(
        face = "bold",
        size = 14
      ),

    legend.position =
      "none"
  )


ggsave(

  filename = file.path(
    output_dir,
    "01_Wee1_J_abundance_vs_signal.png"
  ),

  plot = plot_signal_J,

  width = 14,
  height = 7,
  dpi = 300
)


# ============================================================
# 8. FIGURE :
#    abondance vs signal Wee1
#    C57BL/6J + C57BL/6N COMBINÉS
# ============================================================

summary_signal_NJ <- gene_data %>%

  group_by(
    Treatment,
    CellType
  ) %>%

  summarise(
    n_nuclei = n(),

    total_Wee1_counts = sum(
      raw_count,
      na.rm = TRUE
    ),

    .groups = "drop"
  ) %>%

  group_by(
    Treatment
  ) %>%

  mutate(
    percent_of_nuclei =
      100 * n_nuclei / sum(n_nuclei),

    percent_of_Wee1_signal =
      100 * total_Wee1_counts / sum(total_Wee1_counts)
  ) %>%

  ungroup()


# Ordre des cell types selon contribution moyenne au signal Wee1
cell_order_NJ <- summary_signal_NJ %>%

  group_by(
    CellType
  ) %>%

  summarise(
    average_signal = mean(
      percent_of_Wee1_signal,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%

  arrange(
    desc(average_signal)
  ) %>%

  pull(
    CellType
  )


summary_signal_NJ$CellType <- factor(
  summary_signal_NJ$CellType,
  levels = rev(cell_order_NJ)
)


plot_signal_NJ <- ggplot(
  summary_signal_NJ,
  aes(
    y = CellType
  )
) +

  geom_segment(
    aes(
      x = percent_of_nuclei,
      xend = percent_of_Wee1_signal,
      yend = CellType
    ),
    color = "grey72",
    linewidth = 1.2
  ) +

  geom_point(
    aes(
      x = percent_of_nuclei
    ),
    color = "grey40",
    size = 4
  ) +

  geom_point(
    aes(
      x = percent_of_Wee1_signal,
      color = CellType
    ),
    size = 5
  ) +

  geom_text(
    aes(
      x = percent_of_nuclei,
      label = paste0(
        round(percent_of_nuclei, 1),
        "%"
      )
    ),
    nudge_y = -0.20,
    color = "grey40",
    size = 4
  ) +

  geom_text(
    aes(
      x = percent_of_Wee1_signal,
      label = paste0(
        round(percent_of_Wee1_signal, 1),
        "%"
      ),
      color = CellType
    ),
    nudge_y = 0.20,
    size = 4.2,
    fontface = "bold"
  ) +

  facet_wrap(
    ~ Treatment,
    ncol = 2
  ) +

  scale_color_manual(
    values = celltype_colors
  ) +

  labs(
    title =
      "Abondance cellulaire vs signal Wee1: C57BL/6J +/6N",

    subtitle =
      " Gris = % des noyaux | Couleur = % des raw counts Wee1",

    x =
      "Pourcentage (%)",

    y =
      NULL
  ) +

  theme_classic(
    base_size = 15
  ) +

  theme(
    plot.title =
      element_text(
        face = "bold",
        size = 18
      ),

    axis.text.y =
      element_text(
        face = "bold",
        size = 12
      ),

    strip.text =
      element_text(
        face = "bold",
        size = 14
      ),

    legend.position =
      "none"
  )


ggsave(
  filename = file.path(
    output_dir,
    "02_Wee1_NJ_abundance_vs_signal.png"
  ),

  plot = plot_signal_NJ,

  width = 14,
  height = 7,
  dpi = 300
)

# ============================================================
# 9. DOTPLOT Wee1
# ============================================================

p_dotplot <- DotPlot(

  obj,

  features = gene_of_interest,

  group.by = "CellType"
) +

  RotatedAxis() +

  labs(

    title =
      "Wee1 expression by cell type",

    x =
      NULL,

    y =
      NULL
  ) +

  theme_classic(
    base_size = 14
  ) +

  theme(

    plot.title =
      element_text(
        face = "bold",
        hjust = 0.5
      ),

    axis.text.y =
      element_text(
        face = "bold"
      )
  )


ggsave(

  filename = file.path(
    output_dir,
    "03_Wee1_DotPlot_CellType.png"
  ),

  plot = p_dotplot,

  width = 10,
  height = 6,
  dpi = 300
)


# ============================================================
# 10. Créer ID pseudobulk
# ============================================================

obj$PB_ID <- paste0(

  "PB",
  obj$SampleID,
  ".",

  gsub(
    "_",
    "-",
    obj$CellType
  )
)


# ============================================================
# 11. Metadata pseudobulk
# ============================================================

pb_meta <- obj@meta.data %>%

  select(
    PB_ID,
    SampleID,
    CellType,
    Treatment,
    Background
  ) %>%

  distinct()


check_pb <- pb_meta %>%

  count(
    PB_ID
  ) %>%

  filter(
    n != 1
  )


if (nrow(check_pb) > 0) {

  print(
    check_pb
  )

  stop(
    "Certains PB_ID possèdent plusieurs annotations."
  )
}


# ============================================================
# 12. Construire pseudobulk
# ============================================================

cat(
  "\n============================================\n",
  "CRÉATION PSEUDOBULK\n",
  "============================================\n",
  sep = ""
)


pb <- AggregateExpression(

  obj,

  assays = "RNA",

  group.by = "PB_ID",

  return.seurat = FALSE
)


pb_counts <- pb$RNA


cat(
  "Dimensions pseudobulk :\n"
)

print(
  dim(pb_counts)
)


cat(
  "\nNombre de pseudobulks : ",
  ncol(pb_counts),
  "\n",
  sep = ""
)


# ============================================================
# 13. Aligner metadata / counts
# ============================================================

rownames(pb_meta) <- pb_meta$PB_ID


missing_pb <- setdiff(
  colnames(pb_counts),
  rownames(pb_meta)
)


if (length(missing_pb) > 0) {

  cat(
    "\nColonnes pseudobulk non retrouvées :\n"
  )

  print(
    missing_pb
  )

  stop(
    "Les noms de colonnes pseudobulk ne correspondent pas."
  )
}


pb_meta <- pb_meta[
  colnames(pb_counts),
  ,
  drop = FALSE
]


stopifnot(
  all(
    rownames(pb_meta) ==
      colnames(pb_counts)
  )
)


# ============================================================
# 14. Nombre de noyaux par pseudobulk
# ============================================================

pb_nuclei <- obj@meta.data %>%

  count(
    PB_ID,
    name = "n_nuclei"
  )


pb_meta <- pb_meta %>%

  rownames_to_column(
    "PB_ID_row"
  ) %>%

  left_join(

    pb_nuclei,

    by = c(
      "PB_ID_row" = "PB_ID"
    )
  )


rownames(pb_meta) <- pb_meta$PB_ID_row

pb_meta$PB_ID_row <- NULL


# ============================================================
# 15. Fonction DESeq2 pseudobulk
# ============================================================

run_pseudobulk <- function(

  celltype,

  counts_matrix,

  metadata,

  target_gene,

  mode = "NJ"
) {


  cat(
    "\n--------------------------------------------\n"
  )

  cat(
    "Cell type : ",
    celltype,
    " | Mode : ",
    mode,
    "\n",
    sep = ""
  )


  # ----------------------------------------------------------
  # Sélection CellType
  # ----------------------------------------------------------

  keep_samples <-
    metadata$CellType == celltype


  if (mode == "J") {

    keep_samples <-
      keep_samples &
      metadata$Background == "J"
  }


  counts_ct <- counts_matrix[
    ,
    keep_samples,
    drop = FALSE
  ]


  meta_ct <- metadata[
    keep_samples,
    ,
    drop = FALSE
  ]


  cat(
    "Nombre de pseudobulks : ",
    ncol(counts_ct),
    "\n",
    sep = ""
  )

  print(
    table(
      meta_ct$Treatment
    )
  )


  # ----------------------------------------------------------
  # Vérifier Ctrl/HFpEF
  # ----------------------------------------------------------

  if (
    !all(
      c(
        "Ctrl",
        "HFpEF"
      ) %in%
        meta_ct$Treatment
    )
  ) {

    warning(
      paste(
        "Ctrl ou HFpEF absent pour",
        celltype,
        mode
      )
    )

    return(
      NULL
    )
  }


  # ----------------------------------------------------------
  # Retirer gènes totalement absents
  # ----------------------------------------------------------

  counts_ct <- counts_ct[
    rowSums(counts_ct) > 0,
    ,
    drop = FALSE
  ]


  if (!target_gene %in% rownames(counts_ct)) {

    warning(
      paste(
        target_gene,
        "sans aucun count dans",
        celltype,
        mode
      )
    )

    return(
      NULL
    )
  }


  # ----------------------------------------------------------
  # Facteurs
  # ----------------------------------------------------------

  meta_ct$Treatment <- factor(

    meta_ct$Treatment,

    levels = c(
      "Ctrl",
      "HFpEF"
    )
  )


  meta_ct$Background <- factor(

    meta_ct$Background,

    levels = c(
      "J",
      "N"
    )
  )


  # ----------------------------------------------------------
  # Créer DESeq2
  # ----------------------------------------------------------

  if (mode == "J") {

    dds <- DESeqDataSetFromMatrix(

      countData =
        round(
          counts_ct
        ),

      colData =
        meta_ct,

      design =
        ~ Treatment
    )

  } else {

    dds <- DESeqDataSetFromMatrix(

      countData =
        round(
          counts_ct
        ),

      colData =
        meta_ct,

      design =
        ~ Background + Treatment
    )
  }


  # ----------------------------------------------------------
  # Filtre léger
  #
  # >= 10 counts dans >= 2 souris
  # Wee1 toujours conservé
  # ----------------------------------------------------------

  keep_genes <-

    rowSums(
      counts(dds) >= 10
    ) >= 2 |

    rownames(dds) == target_gene


  dds <- dds[
    keep_genes,
  ]


  cat(
    "Gènes gardés : ",
    nrow(dds),
    "\n",
    sep = ""
  )


  # ----------------------------------------------------------
  # DESeq2
  # ----------------------------------------------------------

  dds <- DESeq(
    dds,
    quiet = TRUE
  )


  # ----------------------------------------------------------
  # HFpEF vs Ctrl
  # ----------------------------------------------------------

  res <- results(

    dds,

    contrast = c(
      "Treatment",
      "HFpEF",
      "Ctrl"
    )
  )


  res_df <- as.data.frame(
    res
  )

  res_df$gene <- rownames(
    res_df
  )


  # ==========================================================
  # IMPORTANT :
  # garder SEULEMENT Wee1
  #
  # Correction du bug gene == gene
  # ==========================================================

  wee1_result <- res_df %>%

    filter(
      .data$gene == .env$target_gene
    ) %>%

    mutate(

      CellType =
        celltype,

      Analysis =
        mode
    )


  cat(
    "\nRésultat ",
    target_gene,
    " :\n",
    sep = ""
  )

  print(
    wee1_result
  )


  # ----------------------------------------------------------
  # Normalized counts DESeq2
  # ----------------------------------------------------------

  normalized_counts <- counts(
    dds,
    normalized = TRUE
  )


  wee1_norm <- normalized_counts[
    target_gene,
    ,
    drop = TRUE
  ]


  plot_data <- data.frame(

    PB_ID =
      rownames(meta_ct),

    SampleID =
      meta_ct$SampleID,

    CellType =
      meta_ct$CellType,

    Treatment =
      meta_ct$Treatment,

    Background =
      meta_ct$Background,

    n_nuclei =
      meta_ct$n_nuclei,

    Wee1_normalized_count =
      as.numeric(
        wee1_norm
      )
  )


  plot_data <- plot_data %>%

    mutate(

      log2_Wee1 =
        log2(
          Wee1_normalized_count + 1
        )
    )


  return(

    list(

      dds =
        dds,

      result =
        wee1_result,

      plot_data =
        plot_data,

      all_results =
        res_df
    )
  )
}


# ============================================================
# 16. Cell types
# ============================================================

celltypes <- sort(
  unique(
    pb_meta$CellType
  )
)


# ============================================================
# 17. DESeq2 J seulement
# ============================================================

cat(
  "\n============================================\n",
  "DESEQ2 PSEUDOBULK : J SEULEMENT\n",
  "design = ~ Treatment\n",
  "============================================\n",
  sep = ""
)


results_J_list <- lapply(

  celltypes,

  function(ct) {

    run_pseudobulk(

      celltype =
        ct,

      counts_matrix =
        pb_counts,

      metadata =
        pb_meta,

      target_gene =
        gene_of_interest,

      mode =
        "J"
    )
  }
)


names(
  results_J_list
) <- celltypes


results_J_valid <- results_J_list[

  !vapply(

    results_J_list,

    is.null,

    logical(1)
  )
]


# ============================================================
# 18. Stats Wee1 J
# ============================================================

wee1_stats_J <- bind_rows(

  lapply(

    results_J_valid,

    function(x) {
      x$result
    }
  )
)


wee1_stats_J <- wee1_stats_J %>%

  select(
    CellType,
    Analysis,
    gene,
    baseMean,
    log2FoldChange,
    lfcSE,
    stat,
    pvalue,
    padj
  ) %>%

  arrange(
    pvalue
  )


write.csv(

  wee1_stats_J,

  file = file.path(
    output_dir,
    "04_Wee1_pseudobulk_DESeq2_J_statistics.csv"
  ),

  row.names = FALSE
)


# ============================================================
# 19. Plot data J
# ============================================================

plot_data_J <- bind_rows(

  lapply(

    results_J_valid,

    function(x) {
      x$plot_data
    }
  )
)


plot_data_J$Treatment <- factor(

  plot_data_J$Treatment,

  levels = c(
    "Ctrl",
    "HFpEF"
  )
)


# ============================================================
# 20. Labels J
# ============================================================

labels_J <- wee1_stats_J %>%

  mutate(

    label =
      case_when(

        is.na(padj) ~ "NA",

        padj < 0.001 ~ "***",

        padj < 0.01 ~ "**",

        padj < 0.05 ~ "*",

        TRUE ~ "ns"
      )
  ) %>%

  left_join(

    plot_data_J %>%

      group_by(
        CellType
      ) %>%

      summarise(

        ymax =
          max(
            log2_Wee1,
            na.rm = TRUE
          ),

        ymin =
          min(
            log2_Wee1,
            na.rm = TRUE
          ),

        .groups = "drop"
      ),

    by = "CellType"
  ) %>%

  mutate(

    range_y =
      pmax(
        ymax - ymin,
        0.2
      ),

    y_position =
      ymax +
      0.15 *
      range_y
  )


# ============================================================
# 21. Boxplot pseudobulk J
# ============================================================

plot_J <- ggplot(

  plot_data_J,

  aes(

    x =
      Treatment,

    y =
      log2_Wee1,

    color =
      Treatment
  )
) +

  geom_boxplot(

    outlier.shape = NA,

    width = 0.50,

    color = "grey40"
  ) +

  geom_jitter(

    width = 0.08,

    size = 3.2,

    alpha = 0.95
  ) +

  facet_wrap(

    ~ CellType,

    scales = "free_y"
  ) +

  geom_text(

    data =
      labels_J,

    aes(

      x =
        1.5,

      y =
        y_position,

      label =
        label
    ),

    inherit.aes = FALSE,

    size = 5,

    fontface = "bold"
  ) +

  scale_color_manual(
    values = treatment_colors
  ) +

  labs(

    title =
      "Wee1 expression: C57BL/6J pseudobulk",

    x =
      NULL,

    y =
      "log2(DESeq2 normalized Wee1 count + 1)"
  ) +

  theme_classic(
    base_size = 14
  ) +

  theme(

    legend.position =
      "none",

    strip.text =
      element_text(
        face = "bold"
      ),

    axis.text.x =
      element_text(
        face = "bold"
      )
  )


ggsave(

  filename = file.path(
    output_dir,
    "05_Wee1_pseudobulk_J_boxplot.png"
  ),

  plot = plot_J,

  width = 12,
  height = 9,
  dpi = 300
)


# ============================================================
# 22. DESeq2 N + J
# ============================================================

cat(
  "\n============================================\n",
  "DESEQ2 PSEUDOBULK : N + J\n",
  "design = ~ Background + Treatment\n",
  "============================================\n",
  sep = ""
)


results_NJ_list <- lapply(

  celltypes,

  function(ct) {

    run_pseudobulk(

      celltype =
        ct,

      counts_matrix =
        pb_counts,

      metadata =
        pb_meta,

      target_gene =
        gene_of_interest,

      mode =
        "NJ"
    )
  }
)


names(
  results_NJ_list
) <- celltypes


results_NJ_valid <- results_NJ_list[

  !vapply(

    results_NJ_list,

    is.null,

    logical(1)
  )
]


# ============================================================
# 23. Stats Wee1 N + J
# ============================================================

wee1_stats_NJ <- bind_rows(

  lapply(

    results_NJ_valid,

    function(x) {
      x$result
    }
  )
)


wee1_stats_NJ <- wee1_stats_NJ %>%

  select(
    CellType,
    Analysis,
    gene,
    baseMean,
    log2FoldChange,
    lfcSE,
    stat,
    pvalue,
    padj
  ) %>%

  arrange(
    pvalue
  ) %>%

  mutate(

    Wee1_padj_across_CellTypes =
      p.adjust(
        pvalue,
        method = "BH"
      )
  )


write.csv(

  wee1_stats_NJ,

  file = file.path(
    output_dir,
    "06_Wee1_pseudobulk_DESeq2_NJ_statistics.csv"
  ),

  row.names = FALSE
)


# ============================================================
# 24. Plot data N + J
# ============================================================

plot_data_NJ <- bind_rows(

  lapply(

    results_NJ_valid,

    function(x) {
      x$plot_data
    }
  )
)


plot_data_NJ$Treatment <- factor(

  plot_data_NJ$Treatment,

  levels = c(
    "Ctrl",
    "HFpEF"
  )
)


plot_data_NJ$Background <- factor(

  plot_data_NJ$Background,

  levels = c(
    "J",
    "N"
  )
)


# ============================================================
# 25. Labels N + J
# ============================================================

labels_NJ <- wee1_stats_NJ %>%

  mutate(

    label =
      case_when(

        is.na(padj) ~ "NA",

        padj < 0.001 ~ "***",

        padj < 0.01 ~ "**",

        padj < 0.05 ~ "*",

        TRUE ~ "ns"
      )
  ) %>%

  left_join(

    plot_data_NJ %>%

      group_by(
        CellType
      ) %>%

      summarise(

        ymax =
          max(
            log2_Wee1,
            na.rm = TRUE
          ),

        ymin =
          min(
            log2_Wee1,
            na.rm = TRUE
          ),

        .groups =
          "drop"
      ),

    by =
      "CellType"
  ) %>%

  mutate(

    range_y =
      pmax(
        ymax - ymin,
        0.2
      ),

    y_position =
      ymax +
      0.15 *
      range_y
  )


# ============================================================
# 26. Boxplot pseudobulk N + J
# ============================================================

plot_NJ <- ggplot(

  plot_data_NJ,

  aes(

    x =
      Treatment,

    y =
      log2_Wee1,

    color =
      Treatment
  )
) +

  geom_boxplot(

    outlier.shape =
      NA,

    width =
      0.50,

    color =
      "grey40"
  ) +

  geom_jitter(

    aes(
      shape =
        Background
    ),

    width =
      0.08,

    size =
      3.2,

    alpha =
      0.95
  ) +

  facet_wrap(

    ~ CellType,

    scales =
      "free_y"
  ) +

  geom_text(

    data =
      labels_NJ,

    aes(

      x =
        1.5,

      y =
        y_position,

      label =
        label
    ),

    inherit.aes =
      FALSE,

    size =
      5,

    fontface =
      "bold"
  ) +

  scale_color_manual(
    values =
      treatment_colors
  ) +

  scale_shape_manual(

    values = c(
      "J" = 16,
      "N" = 17
    )
  ) +

  labs(

    title =
      "Wee1 expression: C57BL/6J + C57BL/6N pseudobulk",

    x =
      NULL,

    y =
      "log2(DESeq2 normalized Wee1 count + 1)",

    shape =
      "Background",

    color =
      "Treatment"
  ) +

  theme_classic(
    base_size = 14
  ) +

  theme(

    strip.text =
      element_text(
        face = "bold"
      ),

    axis.text.x =
      element_text(
        face = "bold"
      ),

    legend.position =
      "right"
  )


ggsave(

  filename = file.path(
    output_dir,
    "07_Wee1_pseudobulk_NJ_boxplot.png"
  ),

  plot = plot_NJ,

  width = 13,
  height = 9,
  dpi = 300
)


# ============================================================
# 27. Cardiomyocytes
# ============================================================

cat(
  "\n============================================\n",
  "CARDIOMYOCYTES : J SEULEMENT\n",
  "============================================\n",
  sep = ""
)

print(

  wee1_stats_J %>%

    filter(
      CellType == "Cardiomyocyte"
    )
)


cat(
  "\n============================================\n",
  "CARDIOMYOCYTES : N + J\n",
  "============================================\n",
  sep = ""
)

print(

  wee1_stats_NJ %>%

    filter(
      CellType == "Cardiomyocyte"
    )
)


# ============================================================
# 28. Résumé direction N + J
# ============================================================

direction_NJ <- wee1_stats_NJ %>%

  mutate(

    Direction =
      case_when(

        is.na(log2FoldChange) ~
          "NA",

        log2FoldChange > 0 ~
          "HFpEF > Ctrl",

        log2FoldChange < 0 ~
          "HFpEF < Ctrl",

        TRUE ~
          "No change"
      )
  ) %>%

  select(

    CellType,

    baseMean,

    log2FoldChange,

    pvalue,

    padj,

    Wee1_padj_across_CellTypes,

    Direction
  )


write.csv(

  direction_NJ,

  file = file.path(
    output_dir,
    "08_Wee1_pseudobulk_NJ_summary.csv"
  ),

  row.names = FALSE
)


# ============================================================
# 29. Fin
# ============================================================

cat(
  "\n============================================\n",
  "ANALYSE TERMINÉE\n",
  "============================================\n",
  sep = ""
)

cat(
  "\nFichiers générés :\n",
  "01_Wee1_J_abundance_vs_signal.png\n",
  "02_Wee1_NJ_abundance_vs_signal.png\n",
  "03_Wee1_DotPlot_CellType.png\n",
  "04_Wee1_pseudobulk_DESeq2_J_statistics.csv\n",
  "05_Wee1_pseudobulk_J_boxplot.png\n",
  "06_Wee1_pseudobulk_DESeq2_NJ_statistics.csv\n",
  "07_Wee1_pseudobulk_NJ_boxplot.png\n",
  "08_Wee1_pseudobulk_NJ_summary.csv\n",
  sep = ""
)
