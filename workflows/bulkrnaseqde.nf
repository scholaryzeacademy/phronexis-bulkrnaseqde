/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FASTQC                 } from '../modules/nf-core/fastqc/main'
include { FASTP                  } from '../modules/nf-core/fastp/main'
include { SALMON_INDEX           } from '../modules/nf-core/salmon/index/main'
include { SALMON_QUANT           } from '../modules/nf-core/salmon/quant/main'
include { CUSTOM_TX2GENE         } from '../modules/nf-core/custom/tx2gene/main'
include { TXIMETA_TXIMPORT       } from '../modules/nf-core/tximeta/tximport/main'
include { GUNZIP                 } from '../modules/nf-core/gunzip/main'
include { DESEQ2                 } from '../modules/local/deseq2/main'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_bulkrnaseqde_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow BULKRNASEQDE {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    multiqc_config
    multiqc_logo
    multiqc_methods_description
    outdir

    main:

    def ch_versions = channel.empty()
    def ch_multiqc_files = channel.empty()

    //
    // Reference channels (value channels so they are reused across all samples)
    //
    def ch_transcript_fasta = channel.value(file(params.transcript_fasta))
    def ch_gtf
    if (params.gtf.toString().endsWith('.gz')) {
        GUNZIP(channel.value([ [id: 'gtf'], file(params.gtf) ]))
        ch_gtf = GUNZIP.out.gunzip.first()
    } else {
        ch_gtf = channel.value([ [id: 'gtf'], file(params.gtf) ])
    }

    //
    // MODULE: FastQC on raw reads
    //
    FASTQC(ch_samplesheet)
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.map{ _meta, file -> file })

    //
    // MODULE: fastp adapter/quality trimming
    //
    FASTP(
        ch_samplesheet.map { meta, reads -> [ meta, reads, [] ] }, // no adapter fasta
        false,  // discard_trimmed_pass
        false,  // save_trimmed_fail
        false   // save_merged
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTP.out.json.map{ _meta, file -> file })

    //
    // MODULE: Build Salmon index (once, reused for all samples)
    //
    SALMON_INDEX(
        [],                       // genome_fasta (empty = lightweight, non-decoy index)
        params.transcript_fasta
    )
    def ch_index = SALMON_INDEX.out.index.collect()

    //
    // MODULE: Salmon quantification per sample
    //
    SALMON_QUANT(
        FASTP.out.reads,
        ch_index,
        params.gtf,
        params.transcript_fasta,
        false,   // alignment_mode (false = mapping-based)
        ''       // lib_type ('' lets Salmon auto-detect with -l A)
    )
    ch_multiqc_files = ch_multiqc_files.mix(SALMON_QUANT.out.json_info.map{ _meta, file -> file })

    //
    // MODULE: Build transcript-to-gene map from GTF + collected quant dirs
    //
    def ch_quants = SALMON_QUANT.out.results.collect{ _meta, dir -> dir }.map { dirs -> [ [id: 'bulkrnaseqde'], dirs ] }

    CUSTOM_TX2GENE(
        ch_gtf,
        ch_quants,
        'salmon',
        'gene_id',
        'gene_name'
    )

    //
    // MODULE: tximport - merge per-sample quants into gene-level matrices
    //
    TXIMETA_TXIMPORT(
        ch_quants,
        CUSTOM_TX2GENE.out.tx2gene,
        'salmon'
    )

    //
    // MODULE: DESeq2 differential expression + HTML report (local module)
    //
    DESEQ2(
        TXIMETA_TXIMPORT.out.counts_gene,
        file(params.input)
    )

    //
    // Collate and save software versions
    //
    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    def ch_collated_versions = softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${outdir}/pipeline_info",
            name:  'bulkrnaseqde_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        )

    //
    // MODULE: MultiQC
    //
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    def ch_summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def ch_workflow_summary = channel.value(paramsSummaryMultiqc(ch_summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    def ch_multiqc_custom_methods_description = multiqc_methods_description
        ? file(multiqc_methods_description, checkIfExists: true)
        : file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true)
    def ch_methods_description = channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: true))
    MULTIQC(
        ch_multiqc_files.flatten().collect().map { files ->
            [
                [id: 'bulkrnaseqde'],
                files,
                multiqc_config
                    ? file(multiqc_config, checkIfExists: true)
                    : file("${projectDir}/assets/multiqc_config.yml", checkIfExists: true),
                multiqc_logo ? file(multiqc_logo, checkIfExists: true) : [],
                [],
                [],
            ]
        }
    )
    emit:multiqc_report = MULTIQC.out.report.map { _meta, report -> [report] }.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
