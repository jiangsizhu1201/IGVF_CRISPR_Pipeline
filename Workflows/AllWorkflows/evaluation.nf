
nextflow.enable.dsl=2
include { evaluation_plot } from './processes/evaluation_plot.nf'

workflow evaluation {
    take:
    sceptre_inference_mudata

    main:
    Evaluation = evaluation_plot(sceptre_inference_mudata, params.central_node, params.min_weight)
}
