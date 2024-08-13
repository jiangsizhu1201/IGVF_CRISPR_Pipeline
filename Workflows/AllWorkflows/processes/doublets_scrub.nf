
process doublets_scrub {
    conda './conda_envs/doublets_scrub.yaml'
    input:
        path mudata

    output:
        path "mdata_doublets.h5mu", emit: mudata_doublet

    script:
        """
        doublets.py --input ${mudata}
        """
}
