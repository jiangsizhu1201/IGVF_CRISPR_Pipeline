
process creatingGuideRef {
    cache 'lenient'
    conda "./conda_envs/kallisto_ref.yaml"

    input:
    path genome
    path guide_metadata

    output:
    path "guide_index.idx" ,  emit: guide_index
    path "t2guide.txt" , emit: t2guide

    script:

    """
        k_bin=\$(which kallisto)
        guide_features_table=\$(guide_features.py --guide_table ${guide_metadata})
        kb ref -i guide_index.idx -f1 $genome -g t2guide.txt --kallisto \$k_bin  --workflow kite guide_features.txt
    """

}
