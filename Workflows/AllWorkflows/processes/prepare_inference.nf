
process prepare_inference {
    conda './conda_envs/pyenv.yaml'

    input:
    val covariate_list

    output:
    path "cov_string.txt", emit: cov_string_file
    path "parse_covariate.csv", emit: parse_covariate_file
    
    script:
    def jsonString = groovy.json.JsonOutput.toJson(covariate_list).replaceAll('"', '\\\\"')
    """
    parse_covariate.py \"$jsonString\"
    prepare_formula.py parse_covariate.csv
    """

}
