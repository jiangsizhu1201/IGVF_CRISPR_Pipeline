process guide_assignment_sceptre {
  
    input:
    path mudata_input

    output:
    path "sceptre_assignment_mudata.h5mu", emit: guide_assignment_mudata_output

    script:
    """
      assign_grnas_sceptre_v1.R ${mudata_input}
    """
}
