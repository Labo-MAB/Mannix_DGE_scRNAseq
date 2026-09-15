library(Seurat)
library(ggplot2)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)

umap_file <- args[1]
output_dir <- args[2]

dir.create(
  output_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

gene <- "Aplnr"

obj <- readRDS(umap_file)

# ============================================================
# 2. CHARGER L'OBJET
# ============================================================

if (file.exists(umap_file)) {

  obj <- readRDS(umap_file)

} else {

  obj <- readRDS(base_file)
}

cat(
  "\nObjet chargé :",
  ncol(obj),
  "noyaux\n"
)

if (!gene %in% rownames(obj)) {
  stop("Aplnr absent de l'objet Seurat.")
}

if (!"umap" %in% names(obj@reductions)) {
  stop(
    "Le UMAP n'est pas présent dans l'objet Seurat."
  )
}


# ============================================================
# 3. EXTRAIRE Aplnr
# ============================================================

aplnr_count <- LayerData(
  obj,
  assay = "RNA",
  layer = "counts"
)[gene, ]

aplnr_count <- as.numeric(aplnr_count)

obj$Aplnr_count <- aplnr_count


# ============================================================
# 4. NORMALISATION DE Aplnr
#
# log(1 + Aplnr counts / total UMI * 10000)
# ============================================================

obj$Aplnr_norm <- log1p(
  obj$Aplnr_count /
    obj$nCount_RNA *
    10000
)


# ============================================================
# 5. DATAFRAME
# ============================================================

df <- obj@meta.data %>%
  mutate(
    barcode = rownames(.)
  ) %>%
  select(
    barcode,
    SampleID,
    CellType,
    Treatment,
    Background,
    nCount_RNA,
    Aplnr_count,
    Aplnr_norm
  )


# ============================================================
# 6. ORDRE DES GROUPES
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

df$Treatment <- factor(
  df$Treatment,
  levels = c(
    "Ctrl",
    "HFpEF"
  )
)

df$Background <- factor(
  df$Background,
  levels = c(
    "J",
    "N"
  )
)


# ============================================================
# 7. COULEURS
# ============================================================

celltype_colors <- c(
  "Cardiomyocyte" = "#E53935",
  "EC" = "#1565C0",
  "Fibroblast" = "#F57C00",
  "Mural_Cell" = "#8E44AD",
  "Mast_Cell" = "#43A047",
  "Myeloid" = "#00ACC1",
  "Lymphoid" = "#6A5ACD",
  "Neural_Cell" = "#C99000"
)

treatment_colors <- c(
  "Ctrl" = "#4C78A8",
  "HFpEF" = "#E45756"
)


# ============================================================
# 8. GARDER LES C57BL/6J
# ============================================================

obj_J <- subset(
  obj,
  subset = Background == "J"
)

df_J <- df %>%
  filter(
    Background == "J"
  )

cat(
  "\n6J seulement :\n"
)

print(
  table(
    df_J$Treatment
  )
)


# ============================================================
# 9. UMAP CELL TYPES
#    Ctrl J vs HFpEF J
# ============================================================

p_umap_celltype <- DimPlot(
  obj_J,
  reduction = "umap",
  group.by = "CellType",
  split.by = "Treatment",
  cols = celltype_colors,
  pt.size = 0.15,
  label = FALSE,
  raster = TRUE
) +
  labs(
    title = "Cell types — C57BL/6J"
  ) +
  theme_classic()

ggsave(
  file.path(
    output_dir,
    "01_UMAP_CellType_J_Ctrl_vs_HFpEF.png"
  ),
  p_umap_celltype,
  width = 14,
  height = 7,
  dpi = 300
)


# ============================================================
# 10. UMAP Aplnr
#     Ctrl J vs HFpEF J
#
# même échelle entre Ctrl et HFpEF
# ============================================================

p_umap_aplnr <- FeaturePlot(
  obj_J,
  features = "Aplnr_norm",
  reduction = "umap",
  split.by = "Treatment",
  pt.size = 0.20,
  order = TRUE,
  raster = TRUE,
  min.cutoff = 0,
  max.cutoff = "q95",
  cols = c(
    "grey90",
    "red"
  ),
  keep.scale = "feature"
)

ggsave(
  file.path(
    output_dir,
    "02_UMAP_Aplnr_J_Ctrl_vs_HFpEF.png"
  ),
  p_umap_aplnr,
  width = 14,
  height = 7,
  dpi = 300
)


# ============================================================
# 11. ABONDANCE CELLULAIRE VS SIGNAL Aplnr
#     J SEULEMENT
# ============================================================

signal_J <- df_J %>%
  group_by(
    Treatment,
    CellType
  ) %>%
  summarise(
    n_nuclei = n(),

    total_Aplnr =
      sum(Aplnr_count),

    .groups = "drop"
  ) %>%
  group_by(
    Treatment
  ) %>%
  mutate(
    percent_nuclei =
      100 *
      n_nuclei /
      sum(n_nuclei),

    percent_Aplnr =
      100 *
      total_Aplnr /
      sum(total_Aplnr)
  ) %>%
  ungroup()


# ordre identique dans les deux panels

cell_order <- signal_J %>%
  group_by(CellType) %>%
  summarise(
    signal =
      mean(percent_Aplnr),

    .groups = "drop"
  ) %>%
  arrange(
    desc(signal)
  ) %>%
  pull(CellType)

signal_J$CellType <- factor(
  signal_J$CellType,
  levels = rev(cell_order)
)


p_signal <- ggplot(
  signal_J,
  aes(
    y = CellType
  )
) +

  geom_segment(
    aes(
      x = percent_nuclei,
      xend = percent_Aplnr,
      yend = CellType
    ),
    color = "grey70",
    linewidth = 1.1
  ) +

  geom_point(
    aes(
      x = percent_nuclei
    ),
    color = "grey40",
    size = 4
  ) +

  geom_point(
    aes(
      x = percent_Aplnr,
      color = CellType
    ),
    size = 5
  ) +

  geom_text(
    aes(
      x = percent_nuclei,
      label = paste0(
        round(percent_nuclei, 1),
        "%"
      )
    ),
    nudge_y = -0.2,
    color = "grey40",
    size = 3.5
  ) +

  geom_text(
    aes(
      x = percent_Aplnr,
      label = paste0(
        round(percent_Aplnr, 1),
        "%"
      ),
      color = CellType
    ),
    nudge_y = 0.2,
    fontface = "bold",
    size = 3.8
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
      "Abondance cellulaire vs signal Aplnr — C57BL/6J",

    subtitle =
      "Gris = % des noyaux | Couleur = % des counts Aplnr",

    x =
      "Pourcentage (%)",

    y =
      NULL
  ) +

  theme_classic(
    base_size = 14
  ) +

  theme(
    legend.position = "none",
    strip.text = element_text(
      face = "bold"
    ),
    axis.text.y = element_text(
      face = "bold"
    )
  )


ggsave(
  file.path(
    output_dir,
    "03_Aplnr_J_abundance_vs_signal.png"
  ),
  p_signal,
  width = 13,
  height = 7,
  dpi = 300
)


# ============================================================
# 12. MOYENNE NORMALISEE PAR SOURIS
# ============================================================

make_summary <- function(data) {

  data %>%
    group_by(
      SampleID,
      CellType,
      Treatment,
      Background
    ) %>%
    summarise(
      n_nuclei = n(),

      mean_Aplnr =
        mean(
          Aplnr_norm,
          na.rm = TRUE
        ),

      percent_positive =
        100 *
        mean(
          Aplnr_count > 0
        ),

      .groups = "drop"
    )
}


summary_J <- make_summary(
  df_J
)

summary_NJ <- make_summary(
  df
)


# ============================================================
# 13. STATISTIQUE
#     J SEULEMENT
#
# un point = une souris
# ============================================================

stats_J <- summary_J %>%
  group_by(
    CellType
  ) %>%
  group_modify(
    ~ {

      x <- .x

      test <- wilcox.test(
        mean_Aplnr ~ Treatment,
        data = x,
        exact = FALSE
      )

      data.frame(
        n_Ctrl =
          sum(
            x$Treatment == "Ctrl"
          ),

        n_HFpEF =
          sum(
            x$Treatment == "HFpEF"
          ),

        mean_Ctrl =
          mean(
            x$mean_Aplnr[
              x$Treatment == "Ctrl"
            ],
            na.rm = TRUE
          ),

        mean_HFpEF =
          mean(
            x$mean_Aplnr[
              x$Treatment == "HFpEF"
            ],
            na.rm = TRUE
          ),

        difference =
          mean(
            x$mean_Aplnr[
              x$Treatment == "HFpEF"
            ],
            na.rm = TRUE
          ) -
          mean(
            x$mean_Aplnr[
              x$Treatment == "Ctrl"
            ],
            na.rm = TRUE
          ),

        pvalue =
          test$p.value
      )
    }
  ) %>%
  ungroup() %>%
  mutate(
    padj =
      p.adjust(
        pvalue,
        method = "BH"
      ),

    analysis =
      "J_only"
  )


# ============================================================
# 14. STATISTIQUE N + J
#     avec ajustement du Background
# ============================================================

stats_NJ <- summary_NJ %>%
  group_by(
    CellType
  ) %>%
  group_modify(
    ~ {

      x <- .x

      x$Treatment <- factor(
        x$Treatment,
        levels = c(
          "Ctrl",
          "HFpEF"
        )
      )

      x$Background <- factor(
        x$Background,
        levels = c(
          "J",
          "N"
        )
      )

      fit <- lm(
        mean_Aplnr ~
          Treatment +
          Background,
        data = x
      )

      ct <- summary(fit)$coefficients

      pvalue <-
        ct[
          "TreatmentHFpEF",
          "Pr(>|t|)"
        ]

      data.frame(
        n_Ctrl =
          sum(
            x$Treatment == "Ctrl"
          ),

        n_HFpEF =
          sum(
            x$Treatment == "HFpEF"
          ),

        mean_Ctrl =
          mean(
            x$mean_Aplnr[
              x$Treatment == "Ctrl"
            ],
            na.rm = TRUE
          ),

        mean_HFpEF =
          mean(
            x$mean_Aplnr[
              x$Treatment == "HFpEF"
            ],
            na.rm = TRUE
          ),

        difference =
          mean(
            x$mean_Aplnr[
              x$Treatment == "HFpEF"
            ],
            na.rm = TRUE
          ) -
          mean(
            x$mean_Aplnr[
              x$Treatment == "Ctrl"
            ],
            na.rm = TRUE
          ),

        pvalue =
          pvalue
      )
    }
  ) %>%
  ungroup() %>%
  mutate(
    padj =
      p.adjust(
        pvalue,
        method = "BH"
      ),

    analysis =
      "N_plus_J_adjusted"
  )


# ============================================================
# 15. UN SEUL CSV
# ============================================================

stats_all <- bind_rows(
  stats_J,
  stats_NJ
)

write.csv(
  stats_all,
  file.path(
    output_dir,
    "06_Aplnr_statistics.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 16. LABELS STATISTIQUES J
# ============================================================

labels_J <- summary_J %>%
  group_by(CellType) %>%
  summarise(
    ymax =
      max(
        mean_Aplnr,
        na.rm = TRUE
      ),

    .groups = "drop"
  ) %>%
  left_join(
    stats_J %>%
      select(
        CellType,
        padj
      ),
    by = "CellType"
  ) %>%
  mutate(
    label =
      case_when(
        padj < 0.001 ~ "***",
        padj < 0.01  ~ "**",
        padj < 0.05  ~ "*",
        TRUE         ~ "ns"
      ),

    y_position =
      ifelse(
        ymax == 0,
        0.05,
        ymax * 1.15
      )
  )


# ============================================================
# 17. BOXPLOT J
# ============================================================

p_J <- ggplot(
  summary_J,
  aes(
    x = Treatment,
    y = mean_Aplnr,
    color = Treatment
  )
) +

  geom_boxplot(
    outlier.shape = NA,
    width = 0.5,
    color = "grey40"
  ) +

  geom_jitter(
    width = 0.08,
    size = 3.2
  ) +

  geom_text(
    data = labels_J,
    aes(
      x = 1.5,
      y = y_position,
      label = label
    ),
    inherit.aes = FALSE,
    size = 5,
    fontface = "bold"
  ) +

  facet_wrap(
    ~ CellType,
    scales = "free_y"
  ) +

  scale_color_manual(
    values = treatment_colors
  ) +

  labs(
    title =
      "Aplnr expression — C57BL/6J",

    subtitle =
      "Ctrl vs HFpEF | chaque point = une souris",

    x =
      NULL,

    y =
      "Mean normalized Aplnr expression"
  ) +

  theme_classic(
    base_size = 14
  ) +

  theme(
    legend.position = "none",
    strip.text = element_text(
      face = "bold"
    )
  )


ggsave(
  file.path(
    output_dir,
    "04_Aplnr_J_boxplot.png"
  ),
  p_J,
  width = 12,
  height = 8,
  dpi = 300
)


# ============================================================
# 18. LABELS N + J
# ============================================================

labels_NJ <- summary_NJ %>%
  group_by(CellType) %>%
  summarise(
    ymax =
      max(
        mean_Aplnr,
        na.rm = TRUE
      ),

    .groups = "drop"
  ) %>%
  left_join(
    stats_NJ %>%
      select(
        CellType,
        padj
      ),
    by = "CellType"
  ) %>%
  mutate(
    label =
      case_when(
        padj < 0.001 ~ "***",
        padj < 0.01  ~ "**",
        padj < 0.05  ~ "*",
        TRUE         ~ "ns"
      ),

    y_position =
      ifelse(
        ymax == 0,
        0.05,
        ymax * 1.15
      )
  )


# ============================================================
# 19. BOXPLOT N + J
# ============================================================

p_NJ <- ggplot(
  summary_NJ,
  aes(
    x = Treatment,
    y = mean_Aplnr,
    color = Treatment
  )
) +

  geom_boxplot(
    aes(
      group = Treatment
    ),
    outlier.shape = NA,
    width = 0.5,
    color = "grey40"
  ) +

  geom_jitter(
    aes(
      shape = Background
    ),
    width = 0.08,
    size = 3.2
  ) +

  geom_text(
    data = labels_NJ,
    aes(
      x = 1.5,
      y = y_position,
      label = label
    ),
    inherit.aes = FALSE,
    size = 5,
    fontface = "bold"
  ) +

  facet_wrap(
    ~ CellType,
    scales = "free_y"
  ) +

  scale_color_manual(
    values = treatment_colors
  ) +

  labs(
    title =
      "Aplnr expression — C57BL/6J + C57BL/6N",

    subtitle =
      "Chaque point = une souris | forme = background",

    x =
      NULL,

    y =
      "Mean normalized Aplnr expression",

    shape =
      "Background"
  ) +

  theme_classic(
    base_size = 14
  ) +

  theme(
    strip.text = element_text(
      face = "bold"
    )
  )


ggsave(
  file.path(
    output_dir,
    "05_Aplnr_NJ_boxplot.png"
  ),
  p_NJ,
  width = 13,
  height = 8,
  dpi = 300
)


# ============================================================
# 20. CARDIOMYOCYTES
# ============================================================

cat(
  "\n========================================\n",
  "Aplnr — CARDIOMYOCYTES — J SEULEMENT\n",
  "========================================\n"
)

print(
  stats_J %>%
    filter(
      CellType == "Cardiomyocyte"
    )
)


cat(
  "\n========================================\n",
  "Aplnr — CARDIOMYOCYTES — N + J\n",
  "========================================\n"
)

print(
  stats_NJ %>%
    filter(
      CellType == "Cardiomyocyte"
    )
)


# ============================================================
# 21. FIN
# ============================================================

cat(
  "\n6 fichiers créés :\n",
  "01_UMAP_CellType_J_Ctrl_vs_HFpEF.png\n",
  "02_UMAP_Aplnr_J_Ctrl_vs_HFpEF.png\n",
  "03_Aplnr_J_abundance_vs_signal.png\n",
  "04_Aplnr_J_boxplot.png\n",
  "05_Aplnr_NJ_boxplot.png\n",
  "06_Aplnr_statistics.csv\n"
)
