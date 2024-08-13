
process inference_sceptre {
    container 'sceptre_renv'

    input:
    path mudata_fp

    output:
    path "sceptre_inference_mudata.h5mu", emit: sceptre_inference_mudata

    script:
    """
    inference_sceptre_v1.R \
        $mudata_fp \
        $params.side \
        $params.grna_integration_strategy \
        $params.resampling_approximation \
        $params.control_group \
        $params.resampling_mechanism \
        '${params.formula_object}'
    """
}
