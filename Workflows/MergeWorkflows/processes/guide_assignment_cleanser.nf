process guide_assignment_cleanser {
    container 'cleanser_env'

    input:
        path mudata_input
        val threshold

    output:
        path "cleanser_mudata_output.h5mu", emit: cleanser_mudata_output

    script:
        """
            igvf_guide_assignment.py  -i ${mudata_input} -o cleanser_mudata_output.h5mu -t ${threshold} --cleanser
        """
}