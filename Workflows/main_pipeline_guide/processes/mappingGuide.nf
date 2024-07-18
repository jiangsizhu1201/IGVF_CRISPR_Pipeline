process mappingGuide {
    cache 'lenient'
    cpus 6
    debug true
    conda "./conda_envs/kallisto_ref.yaml"

    input:
    tuple path(fastq_file1), path(fastq_file2)
    path guide_index
    path t2guide
    path parsed_seqSpec_file
    path whitelist

    output:
    path("ks_guide_out"), emit: ks_guide_out
    path("ks_guide_out/counts_unfiltered/adata.h5ad"), emit: ks_guide_out_adata_h5ad

    script:
        """
        k_bin=\$(which kallisto)
        chemistry=\$(extract_parsed_seqspec.py --file ${parsed_seqSpec_file})
        kb count -i ${guide_index} -g ${t2guide} --verbose -w ${whitelist} --h5ad --kallisto \$k_bin -x \$chemistry -o ks_guide_out -t ${task.cpus} ${fastq_file1} ${fastq_file2} --overwrite -m 48G
        """
}
