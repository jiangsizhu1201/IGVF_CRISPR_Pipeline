
process downloadGTF {
    input:
    val gtf_path

    output:
    path "hg38.ensGene.gtf" , emit: gtf

    script:
    """
        wget $gtf_path
        gunzip hg38.ensGene.gtf.gz
    """
}
