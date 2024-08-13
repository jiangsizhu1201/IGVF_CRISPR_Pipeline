
nextflow.enable.dsl=2

include { seqSpecParser } from './processes/seqSpecParser.nf'
include { downloadReference } from './processes/downloadReference.nf'
include { prepare_inference } from "./processes/prepare_inference"
include { mappingscRNA } from './processes/mappingscRNA.nf'
include { anndata_concat } from './processes/anndata_concat.nf'

workflow mapping_rna_pipeline {
    main:
    SeqSpecResult = seqSpecParser(
        file("${params.seqSpecDirectory}/${params.seqSpec_yaml_rna}"),
        file(params.seqSpecDirectory),
        'rna'
        )

    DownloadRefResult = downloadReference(params.TRANSCRIPTOME_REFERENCE)
    
    Prepare_covariate = prepare_inference(params.covariate_list_rna)

    MappingOut = mappingscRNA(
        file(params.fastq_files_rna_dir),
        params.fastq_files_rna,
        Prepare_covariate.parse_covariate_file,
        DownloadRefResult.transcriptome_idx,
        DownloadRefResult.t2g_transcriptome_index,
        SeqSpecResult.parsed_seqspec,
        SeqSpecResult.whitelist
        )

    AnndataConcatenate = anndata_concat(
        Prepare_covariate.parse_covariate_file,
        MappingOut.ks_transcripts_out_dir
    )

    emit:
    trans_out_dir = MappingOut.ks_transcripts_out_dir
    concat_anndata_rna = AnndataConcatenate.concat_anndata
    covariate_list =  Prepare_covariate.cov_string_file
}
