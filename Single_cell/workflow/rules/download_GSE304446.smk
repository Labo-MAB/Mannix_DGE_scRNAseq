rule download_GSE304446:
    output:
        mtx="data/single_cell/GSE304446_matrix.mtx.gz",
        features="data/single_cell/GSE304446_features.tsv.gz",
        barcodes="data/single_cell/GSE304446_barcodes.tsv.gz",
        metadata="data/single_cell/GSE304446_cell_metadata.tsv.gz"
    params:
        base_url="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE304nnn/GSE304446/suppl"
    shell:
        r"""
        mkdir -p data/single_cell

        wget -O {output.mtx} \
            {params.base_url}/GSE304446_matrix.mtx.gz

        wget -O {output.features} \
            {params.base_url}/GSE304446_features.tsv.gz

        wget -O {output.barcodes} \
            {params.base_url}/GSE304446_barcodes.tsv.gz

        wget -O {output.metadata} \
            {params.base_url}/GSE304446_cell_metadata.tsv.gz
        """
