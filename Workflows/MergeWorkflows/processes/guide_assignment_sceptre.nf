process guide_assignment_sceptre {
    container 'sceptre_renv'

    input:
    path mudata_input

    output:
    path "sceptre_mudata_output.h5mu", emit: sceptre_mudata_output

    script:
    """
    Rscript -e \"
      mudata_in <- MuData::readH5MU('${mudata_input}');
      mudata_out <- sceptreIGVF::assign_grnas_sceptre(mudata = mudata_in);
      MuData::writeH5MU(mudata_out, 'sceptre_mudata_output.h5mu');
    \"
    """
}
