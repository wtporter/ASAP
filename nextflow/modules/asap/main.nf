#! /usr/bin/env nextflow

process PREPARE_ASAP_JSON {
    tag "Preparing ASAP JSON"

    input:
    path input_files

    output:
    path "assay_input.json", emit: json

    script:
    def args = ""
    def file_list = input_files instanceof List ? input_files : [input_files]
    def num_files = file_list.size()
    def first_file = file_list[0].name.toLowerCase()

    // 1. Identify Format
    def is_gb = first_file.endsWith('.gb') || first_file.endsWith('.gbk') || first_file.endsWith('.gbb') || first_file.endsWith('.gbf') || first_file.endsWith('.gbff') || first_file.endsWith('.genbank')
    def is_fasta = first_file.endsWith('.fasta') || first_file.endsWith('.fa')
    def is_excel = first_file.endsWith('.xlsx') || first_file.endsWith('.xls')

    // 2. Validate Multi-file Rule
    if (num_files > 1 && !is_gb) {
        error """
        ERROR: Multiple files detected for non-GenBank input.
        Format detected: ${is_fasta ? 'FASTA' : is_excel ? 'Excel' : 'Unknown'}
        Number of files: ${num_files}
        
        ASAP only supports multiple reference files when using GenBank (.gb, .gbf, .gbb, .gbk, .gbff) format.
        Please provide only one file for FASTA or Excel inputs.
        """.stripIndent()
    }

    // 3. Construct Arguments & Log Status
    if (is_fasta) {
        log.info "[ASAP] Preparing JSON from single FASTA: ${file_list[0].name}"
        args = "-f ${file_list[0]}"
    } else if (is_gb) {
        log.info "[ASAP] Preparing JSON from ${num_files} GenBank file(s): ${file_list*.name.join(', ')}"
        args = "-g ${file_list.join(' ')}"
    } else if (is_excel) {
        log.info "[ASAP] Preparing JSON from single Excel sheet: ${file_list[0].name}"
        args = "-x ${file_list[0]}"
    } else {
        error "[ASAP] Unsupported reference format: ${first_file}. Expected FASTA, GB, or Excel."
    }

    """
    prepareJSONInput_nextflow.py \\
        ${args} \\
        -o assay_input.json
    """
}

process GENERATE_REFERENCE_FASTA {
    tag "generate_reference"
    publishDir "${params.outdir}/reference", mode: 'copy'

    input:
    path assay_json

    output:
    path "reference.fasta", emit: ref_fasta

    script:
    """
    assayInfo.py ${assay_json} > reference.fasta
    """
}

process MASK_PRIMERS {
    tag "$sample_id"
    publishDir "${params.outdir}/sample_info/${sample_id}/mask_primers", mode: 'copy'

    input:
    tuple val(sample_id), path(bamfile), path(bamindex), path(primer_file)

    output:
    tuple val(sample_id), path("${bamfile.getBaseName()}_primerMasked.bam"), path("${bamfile.getBaseName()}_primerMasked.bam.bai"), emit: mask_primers_output
    tuple val(sample_id), path("primer_masking.tsv"), path("primer_masking.log"), emit: mask_primers_logging
    tuple val(sample_id), path("primer_masking_stats.tsv"), emit: mask_primers_stats
    tuple val(sample_id), path("${sample_id}_masked_reads_per_primer.tsv"), emit: mask_primers_primer_stats

    script:
    def mask_bam_string = params.mask_bam ? "--mask-bam" : "--no-mask-bam"
    def ponly_string = params.primer_only ? "--primer-only" : "--no-primer-only"
    """
    maskPrimers.py -b ${bamfile} -p ${primer_file} --wiggle ${params.wiggle} ${mask_bam_string} ${ponly_string}

    # Per-sample masked-reads-per-primer report (prepend sample_id so the same
    # file also feeds the combined cross-sample report via collectFile).
    awk -v s="${sample_id}" 'BEGIN{OFS="\\t"} NR==1{print "sample_id", \$0; next} {print s, \$0}' \\
        primer_masking_primer_stats.tsv > ${sample_id}_masked_reads_per_primer.tsv
    """
}

process IDENTITY_FILTER {
    tag "$sample_id"
    publishDir "${params.outdir}/sample_info/${sample_id}/identity_filter", mode: 'copy'

    input:
    tuple val(sample_id), path(bamfile), path(bamindex)
    
    output:
    tuple val(sample_id), path("${bamfile.getBaseName()}_identityFiltered.bam"), path("${bamfile.getBaseName()}_identityFiltered.bam.bai"), emit: identity_filter_output
    tuple val(sample_id), path("identity_filtering.log"), emit: identity_filter_logging
    tuple val(sample_id), path("identity_filter_stats.tsv"), emit: identity_filter_stats

    script:
    def filter_pairs_flag = params.filter_pairs ? "" : "--no-filter-pairs"
    """
    identityFilter.py -b ${bamfile} -i ${params.identity} ${filter_pairs_flag}
    """
}

process SMOR {
    tag "$sample_id"
    publishDir "${params.outdir}/sample_info/${sample_id}/smor", mode: 'copy'

    input:
    tuple val(sample_id), path(bamfile), path(bamindex)

    output:
    tuple val(sample_id), path("${bamfile.getBaseName()}_SMOR.bam"), path("${bamfile.getBaseName()}_SMOR.bam.bai"), emit: smor_output
    tuple val(sample_id), path("smor_processing.log"), emit: smor_logging
    tuple val(sample_id), path("smor_stats.tsv"), emit: smor_stats

    script:
    """
    generateSMORbam.py -b ${bamfile} -c ${params.fill_character} 
    """
}

process SMOR_CORRECTION {
    tag "$sample_id"
    publishDir "${params.outdir}/sample_info/${sample_id}/smor_correction/", mode: 'copy'

    input:
    tuple val(sample_id), path(bamfile), path(bamindex)

    output:
    tuple val(sample_id), path("${bamfile.getBaseName()}_SMOR.bam"), path("${bamfile.getBaseName()}_SMOR.bam.bai"), emit: smor_output
    tuple val(sample_id), path("smor_processing.log"), emit: smor_logging
    tuple val(sample_id), path("smor_stats.tsv"), emit: smor_stats

    script:
    """
    generateSMORbam_correction.py -b ${bamfile} -c ${params.fill_character} -q ${params.smor_correction_qual_diff_threshold} -a ${params.smor_correction_agreement_method}
    """
}

process PROCESS_BAM {
    tag "$sample_id"
    publishDir "${params.outdir}/sample_info/${sample_id}/xml", mode: 'copy'

    input:
    tuple val(sample_id), path(bamfile), path(bamindex),
          path(original_bam,   stageAs: 'pre_filter.bam'),
          path(original_bai,   stageAs: 'pre_filter.bam.bai'),
          path(fastp_json),
          path(primer_stats,   stageAs: 'primer_stats'),
          path(identity_stats, stageAs: 'identity_stats'),
          path(smor_stats,     stageAs: 'smor_stats'),
          path(assay_json),
          path("genbank_input/*")

    output:
    tuple val(sample_id), path("${sample_id}.xml"), emit: xml_output

    script:
    def wg_flag = params.suppress_per_base ? "--suppress-per-base" : ""
    def prune_flag = params.prune_per_base ? "--prune-per-base" : ""
    def primer_flag   = primer_stats.size()   > 0 ? "--primer-stats ${primer_stats}"     : ""
    def identity_flag = identity_stats.size() > 0 ? "--identity-stats ${identity_stats}" : ""
    def smor_flag     = smor_stats.size()     > 0 ? "--smor-stats ${smor_stats}"         : ""
    def codon_flag        = params.codon_correction ? "--codon-correction" : ""
    def codon_gb_flag     = params.codon_correction ? "--codon-correction-genbank genbank_input/*" : ""
    def codon_err_flag    = params.codon_correction
                              ? "--codon-correction-error ${params.codon_correction_error}" : ""
    def codon_min_flag    = params.codon_correction
                              ? "--codon-correction-min-reads ${params.codon_correction_min_reads}" : ""
    def droi_flag         = params.discover_roi ? "--discover-roi" : ""
    def droi_perc_flag    = params.discover_roi
                              ? "--discover-roi-min-perc ${params.discover_roi_min_perc}" : ""
    def droi_min_flag     = params.discover_roi
                              ? "--discover-roi-min-reads ${params.discover_roi_min_reads}" : ""
    def droi_min_snp_flag = params.discover_roi
                              ? "--discover-roi-min-snp-perc ${params.discover_roi_min_snp_perc}" : ""

    """
    newBamProcessor.py \\
        -j ${assay_json} \\
        -b ${bamfile} \\
        -d ${params.depth} \\
        --breadth ${params.breadth} \\
        -p ${params.proportion} \\
        -m ${params.mutation_depth} \\
        --min-base-qual ${params.min_base_qual} \\
        --consensus-proportion ${params.consensus_proportion} \\
        --fill-gaps ${params.fill_gaps} \\
        --mark-deletions ${params.mark_deletions} \\
        --original-bam ${original_bam} \\
        --fastp-json ${fastp_json} \\
        ${primer_flag} \\
        ${identity_flag} \\
        ${smor_flag} \\
        ${wg_flag} \\
        ${prune_flag} \\
        ${codon_flag} \\
        ${codon_gb_flag} \\
        ${codon_err_flag} \\
        ${codon_min_flag} \\
        ${droi_flag} \\
        ${droi_perc_flag} \\
        ${droi_min_flag} \\
        ${droi_min_snp_flag} \\
        -o ${sample_id}.xml
    """
}

process OUTPUT_COMBINER {
    tag "output_combiner"
    publishDir "${params.outdir}/sample_reports/general_reports", mode: 'copy'

    input:
    path xml_files
    
    output:
    path("${params.file_name}_analysis.xml"), emit: final_xml

    script:
    """
    outputCombiner.py -x . -n ${params.file_name}
    """
}

process FORMAT_OUTPUT {
    tag "format_output"
    publishDir "${params.outdir}/sample_reports/general_reports", mode: 'copy'
    stageInMode = 'copy'

    def out_file = params.out_file ? params.out_file : "${params.file_name}_report.html"

    input:
    path final_xml
    path stylesheet
    
    output:
    path("*.html"), emit: asap_output
    path("${params.file_name}/"), emit : extra_output, optional: true

    script:
    """
    formatOutput.py -x ${final_xml} -s ${stylesheet} -o ${out_file}
    """
}
