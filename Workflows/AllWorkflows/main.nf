nextflow.enable.dsl=2

include { mapping_rna_pipeline } from './mapping_rna_pipeline.nf'
include { mapping_guide_pipeline } from './mapping_guide_pipeline.nf'
include { process_mudata_pipeline } from './process_mudata_pipeline.nf'
include { evaluation } from './evaluation.nf'

workflow {
  mapping_rna_pipeline()
  mapping_guide_pipeline()

  process_mudata_pipeline(
    mapping_rna_pipeline.out.concat_anndata_rna,
    mapping_rna_pipeline.out.trans_out_dir,
    mapping_guide_pipeline.out.concat_anndata_guide,
    mapping_guide_pipeline.out.guide_out_dir,
    mapping_rna_pipeline.out.covariate_list
  )

  evaluation(process_mudata_pipeline.out.sceptre_inference_mudata)
}
