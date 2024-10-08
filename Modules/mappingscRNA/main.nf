nextflow.enable.dsl=2

include { mappingscRNA } from './processes/mappingscRNA.nf'

workflow {
  map_rna = mappingscRNA(
        tuple(file(params.fastq_file1), file(params.fastq_file2)),
        file(params.transcriptome_idx),
        file(params.transcriptome_t2g),
        file(params.parsed_seqSpec_file),
        file(params.whitelist)
     )
}