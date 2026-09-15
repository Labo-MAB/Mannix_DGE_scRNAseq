rule create_seurat:
    input:
        mtx="data/single_cell/GSE304446_matrix.mtx.gz",
        features="data/single_cell/GSE304446_features.tsv.gz",
        barcodes="data/single_cell/GSE304446_barcodes.tsv.gz",
        metadata="data/single_cell/GSE304446_cell_metadata.tsv.gz"
    output:
        rds="results/single_cell/GSE304446_seurat.rds"
    params:
        script="scripts/create_seurat.R"
    shell:
        """
        mkdir -p results/single_cell
        Rscript {params.script}
        """

rule create_umap:
    input:
        rds=rules.create_seurat.output.rds
    output:
        rds="results/single_cell/umap/GSE304446_seurat_UMAP.rds"
    params:
        script="scripts/create_umap.R"
    shell:
        """
        mkdir -p results/single_cell/umap

        Rscript {params.script} \
            {input.rds} \
            {output.rds}
        """

rule wee1_pseudobulk:
    input:
        rds=rules.create_seurat.output.rds
    output:
        abundance_J="results/single_cell/gene_expression/01_Wee1_J_abundance_vs_signal.png",
        abundance_NJ="results/single_cell/gene_expression/02_Wee1_NJ_abundance_vs_signal.png",
        dotplot="results/single_cell/gene_expression/03_Wee1_DotPlot_CellType.png",
        stats_J="results/single_cell/gene_expression/04_Wee1_pseudobulk_DESeq2_J_statistics.csv",
        boxplot_J="results/single_cell/gene_expression/05_Wee1_pseudobulk_J_boxplot.png",
        stats_NJ="results/single_cell/gene_expression/06_Wee1_pseudobulk_DESeq2_NJ_statistics.csv",
        boxplot_NJ="results/single_cell/gene_expression/07_Wee1_pseudobulk_NJ_boxplot.png",
        summary_NJ="results/single_cell/gene_expression/08_Wee1_pseudobulk_NJ_summary.csv"
    params:
        script="scripts/Wee1_pseudobulk.R",
        output_dir="results/single_cell/gene_expression"
    shell:
        """
        Rscript {params.script} \
            {input.rds} \
            {params.output_dir}
        """
