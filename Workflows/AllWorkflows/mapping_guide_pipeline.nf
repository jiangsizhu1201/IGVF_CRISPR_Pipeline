
nextflow.enable.dsl=2

include { seqSpecParser } from './processes/seqSpecParser.nf'
include { downloadGenome } from './processes/downloadGenome.nf'
include { creatingGuideRef } from './processes/createGuideIndex.nf'
include { prepare_inference } from "./processes/prepare_inference"
include { mappingGuide } from './processes/mappingGuide.nf'
include { anndata_concat } from './processes/anndata_concat.nf'

workflow mapping_guide_pipeline {
    main:
    SeqSpecResult = seqSpecParser(
        file("${params.seqSpecDirectory}/${params.seqSpec_yaml_guide}"),
        file(params.seqSpecDirectory),
        'guide')

    Genome = downloadGenome(params.genome_path)
    GuideRef = creatingGuideRef(Genome.genome, file(params.guide_metadata))

    Prepare_covariate = prepare_inference(params.covariate_list_guide)

    MappingOut = mappingGuide(
        file(params.fastq_files_guide_dir),
        params.fastq_files_guide,
        Prepare_covariate.parse_covariate_file,
        GuideRef.guide_index,
        GuideRef.t2guide,
        SeqSpecResult.parsed_seqspec,
        SeqSpecResult.whitelist
        )
    
    AnndataConcatenate = anndata_concat(
        Prepare_covariate.parse_covariate_file,
        MappingOut.ks_guide_out_dir
    )

    emit:
    guide_out_dir = MappingOut.ks_guide_out_dir
    concat_anndata_guide = AnndataConcatenate.concat_anndata
}
