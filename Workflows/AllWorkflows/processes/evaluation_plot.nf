
process evaluation_plot {
    conda "./conda_envs/evaluation_plot.yaml"

    input:

    path mdata
    val central_node
    val min_weight

    output:
    path "evaluation_output" , emit: evaluation_output

    script:
            """
            network_plot.py ${mdata} --central_node ${central_node}  --min_weight=${min_weight}
            volcano_plot.py ${mdata}
            igv.py ${mdata}
            """
}
