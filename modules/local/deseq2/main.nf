process DESEQ2 {
    tag "$meta.id"
    label 'process_medium'

    container 'quay.io/biocontainers/bioconductor-deseq2:1.50.2--r45ha27e39d_0'

    input:
    tuple val(meta), path(counts)
    path samplesheet

    output:
    tuple val(meta), path("*.de_results.tsv"), emit: results
    tuple val(meta), path("*.pca.png")       , emit: pca
    tuple val(meta), path("*.volcano.png")    , emit: volcano
    tuple val(meta), path("*.de_report.html") , emit: report
    path "versions.yml"                        , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'deseq2.R'
}
