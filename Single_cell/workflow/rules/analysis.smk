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
