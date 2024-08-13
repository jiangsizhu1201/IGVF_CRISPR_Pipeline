nextflow.enable.dsl=2

include { prepare_inference } from './processes/prepare_inference.nf'
include { inference_sceptre } from './processes/inference_sceptre_v1.nf'

workflow {
    covariate_list = prepare_inference(params.covariate_list)
    sceptre_inference_mudata = inference_sceptre(
        file(params.mudata_input),
        covariate_list.cov_string_file
        )
}
