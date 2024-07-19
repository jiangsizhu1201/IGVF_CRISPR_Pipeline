
nextflow.enable.dsl=2

include { PreprocessAnnData } from './processes/PreprocessAnnData.nf'
include { CreateMuData } from './processes/CreateMuData.nf'

workflow preprocess_adata_pipeline {

    take:
    trans_adata
    trans_gnames
    guide_adata

    main:
    Preprocessed_AnnData = PreprocessAnnData(
        trans_adata,
        trans_gnames,
        params.min_genes,
        params.min_cells,
        params.TRANSCRIPTOME_REFERENCE)

    MuData = CreateMuData(
        Preprocessed_AnnData.filtered_anndata_rna,
        guide_adata,
        file(params.guide_features))
}
