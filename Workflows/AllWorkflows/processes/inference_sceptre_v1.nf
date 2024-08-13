
process inference_sceptre {
container 'igvf/sceptre-igvf:v0.2'

    input:
    path mudata_fp
    path cov_string_file

    output:
    path "sceptre_inference_mudata.h5mu", emit: sceptre_inference_mudata

    script:
    """
    cov_string=\$(cat $cov_string_file)

    cat <<EOF > args.txt
$mudata_fp
$params.side
$params.grna_integration_strategy
$params.resampling_approximation
$params.control_group
$params.resampling_mechanism
${params.formula_object}
\$cov_string
EOF

    inference_sceptre_v1.R args.txt
    """
}
