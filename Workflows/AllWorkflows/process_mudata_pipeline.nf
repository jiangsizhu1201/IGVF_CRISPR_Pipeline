
nextflow.enable.dsl=2

include { PreprocessAnnData } from './processes/PreprocessAnnData.nf'
include { CreateMuData } from './processes/CreateMuData.nf'
include { doublets_scrub } from './processes/doublets_scrub.nf'
include { guide_assignment_cleanser } from './processes/guide_assignment_cleanser.nf'
include { guide_assignment_sceptre } from './processes/guide_assignment_sceptre.nf'
include { downloadGTF } from './processes/downloadGTF.nf'
include { prepare_guide_inference } from './processes/prepare_guide_inference.nf'
include { prepare_inference } from './processes/prepare_inference.nf'
include { prepare_user_guide_inference } from './processes/prepare_user_guide_inference.nf'
include { inference_sceptre } from './processes/inference_sceptre_v1.nf'

workflow process_mudata_pipeline {

    take:
    concat_anndata_rna
    trans_out_dir
    concat_anndata_guide
    guide_out_dir
    covariate_list

    main:

    Preprocessed_AnnData = PreprocessAnnData(
        concat_anndata_rna,
        trans_out_dir.flatten().first(),
        params.min_genes,
        params.min_cells,
        params.TRANSCRIPTOME_REFERENCE
        )
    
    GTF_Reference = downloadGTF(params.gtf_path)

    MuData = CreateMuData(
        Preprocessed_AnnData.filtered_anndata_rna,
        concat_anndata_guide, 
        file(params.guide_metadata),
        GTF_Reference.gtf
        )

    MuData_Doublets = doublets_scrub(MuData.mudata)
    
    //adding cov
    
    if (params.assignment_method == "cleanser") {
    Guide_Assignment = guide_assignment_cleanser(MuData_Doublets.mudata_doublet, params.THRESHOLD)}
    else if (params.assignment_method == "sceptre") {
    Guide_Assignment = guide_assignment_sceptre(MuData_Doublets.mudata_doublet)}

    if (file(params.user_inference).exists()) {
        PrepareInference = prepare_user_guide_inference(
            Guide_Assignment.guide_assignment_mudata_output,
            file(params.user_inference)
        )}
    else {
        PrepareInference = prepare_guide_inference(
            Guide_Assignment.guide_assignment_mudata_output,
            GTF_Reference.gtf,
            params.limit
            )}

    GuideInference = inference_sceptre(PrepareInference.mudata_inference_input, covariate_list)

    emit:
    sceptre_inference_mudata = GuideInference.sceptre_inference_mudata

}
