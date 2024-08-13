
nextflow.enable.dsl=2

include { PreprocessAnnData } from './processes/PreprocessAnnData.nf'
include { CreateMuData } from './processes/CreateMuData.nf'
include { doublets_scrub } from './processes/doublets_scrub.nf'
include { guide_assignment_cleanser } from './processes/guide_assignment_cleanser.nf'
include { guide_assignment_sceptre } from './processes/guide_assignment_sceptre.nf'
include { downloadGTF } from './processes/downloadGTF.nf'

workflow {

    Preprocessed_AnnData = PreprocessAnnData(
        file(params.adata_rna),
        file(params.gname_rna),
        params.min_genes,
        params.min_cells,
        params.TRANSCRIPTOME_REFERENCE)
    
    GTF_Reference = downloadGTF(params.gtf_path)

    MuData = CreateMuData(
        Preprocessed_AnnData.filtered_anndata_rna,
        file(params.adata_guide), 
        file(params.guide_features),
        GTF_Reference.gtf)
// emit:adata_guide

    MuData_Doublets = doublets_scrub(MuData.mudata)
// multi-seq doublets 
    if (params.assignment_method == "cleanser") {
    Guide_Assignment = guide_assignment_cleanser(MuData_Doublets.mudata_doublet, params.THRESHOLD)}
    else if (params.assignment_method == "sceptre") {
    Guide_Assignment = guide_assignment_sceptre(MuData_Doublets.mudata_doublet)}
// define default 
}
