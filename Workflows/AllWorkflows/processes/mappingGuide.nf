
process mappingGuide {
    cache 'lenient'
    cpus 2
    debug true
    conda "./conda_envs/kallisto_ref.yaml"

    input:
    path fastq_files_guide_dir
    val  fastq_files_guide
    path parsed_covariate_df
    path guide_index
    path t2guide
    path parsed_seqSpec_file
    path whitelist

    output:
    path "*_ks_guide_out", emit: ks_guide_out_dir

    script:
        """
        k_bin=\$(which kallisto)
        chemistry=\$(extract_parsed_seqspec.py --file ${parsed_seqSpec_file})
        processed_batches=\$(process_batches.py --dir ${fastq_files_guide_dir} --fastq '${fastq_files_guide}')
        echo "\${processed_batches}"
        
        awk -F',' 'NR>1 {print \$1}' ${parsed_covariate_df} > batch_array.txt

        while IFS= read -r line; do
            index=\$(echo \$line | cut -d' ' -f1)
            num=\$(sed -n "\$((index+1))p" batch_array.txt)
            batch_files=\$(echo \$line | cut -d' ' -f2-)
            output_dir="\${num}_ks_guide_out"
            echo "\${output_dir}"
            echo "\${batch_files}"
            kb count -i ${guide_index} -g ${t2guide} --verbose -w ${whitelist} \\
                --h5ad --kallisto \$k_bin -x \$chemistry -o \${output_dir} -t ${task.cpus} \\
                \$batch_files --overwrite -m 50G
            echo "KB mapping Complete"
        done <<< "\${processed_batches}"

        # Clean up
        rm batch_array.txt
        """
}
