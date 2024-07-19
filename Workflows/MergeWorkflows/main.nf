
nextflow.enable.dsl=2

include { mapping_rna_pipeline } from './mapping_rna_pipeline.nf'
include { mapping_guide_pipeline } from './mapping_guide_pipeline.nf'
include { preprocess_adata_pipeline } from './preprocess_adata_pipeline.nf'

workflow {
  mapping_rna_pipeline()
  mapping_guide_pipeline()

  preprocess_adata_pipeline(
    mapping_rna_pipeline.out.trans_adata, 
    mapping_rna_pipeline.out.trans_gnames, 
    mapping_guide_pipeline.out.guide_adata
    )
}

