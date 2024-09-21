
process mappingGuide {
   
    cache 'lenient'
    cpus 2
    debug true

    input:
    each batch_fastq
    path guide_index
    path t2g_guide
    path parsed_seqSpec_file
    path barcode_file

    output:
    path "*_ks_guide_out", emit: ks_guide_out_dir
    path "*_ks_guide_out/counts_unfiltered/adata.h5ad", emit: ks_guide_out_adata

    script:
        batch = batch_fastq[0]
        fastq_files = batch_fastq[1]
        """
        echo "Processing $batch with $fastq_files"
        
        k_bin=\$(which kallisto)
        bustools_bin=\$(which bustools)
        chemistry=\$(extract_parsed_seqspec.py --file ${parsed_seqSpec_file})

        kb count -i ${guide_index} -g ${t2g_guide} --verbose -w ${barcode_file} \\
                --h5ad --kallisto \$k_bin --bustools \$bustools_bin -x \$chemistry -o ${batch}_ks_guide_out -t ${task.cpus} \\
                ${fastq_files} --overwrite -m 30G

        echo "KB mapping Complete"
        """
}
