#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
============================================================================
  BatchITSx: Nextflow-based ITSx pipeline
============================================================================
  Version: v0.1
  License: Apache-2.0
  Github : https://github.com/kessya/BatchITSx
  Author : Kessy Abarenkov
----------------------------------------------------------------------------
*/

// =============================
// Parameters
// =============================
version = '0.1'

params.input      = params.input ?: "../test.fasta"
params.outdir     = params.outdir ?: "${launchDir}/results"
params.method     = "itsx"
params.itsx_bin   = params.itsx_bin ?: "ITSx"
params.vsearch_bin = params.vsearch_bin ?: "vsearch"
params.region     = params.region ?: "itsfull"

params.itsx_target     = params.itsx_target ?: "F,All"
params.itsx_chunksize  = params.itsx_chunksize ?: 100
params.itsx_regions    = "ITS1,5.8S,ITS2"
params.itsx_partial    = "50"
params.itsx_detailed   = "T"
params.itsx_concat     = "T"
params.itsx_preserve   = "T"
params.itsx_eval       = 0.1
params.itsx_e          = 0.1
params.itsx_graphical  = "F"
params.itsx_complement = "T"

// Check if input path was provided
if( !params.input ) {
    error "Please provide --input file (FASTA)"
}

// Convert ITSx targets string to list if necessary
def itsx_targets_list = params.itsx_target instanceof String ?
    params.itsx_target.split(',').collect { it.trim() } :
    params.itsx_target

// Help message
def helpMsg() {
    log.info """
    =====================================================================
    BatchITSx ${version}
    =====================================================================

    Usage:
      nextflow run kessya/batchitsx --input ... --outdir ...

    Options:
      --input         Input file with sequences (FASTA)
      --outdir        Output directory
      --itsx_target   Organism groups (default: F,All)
      --itsx_regions  ITSx extracted regions (default: ITS1,5.8S,ITS2)
      --itsx_partial  Partial cutoff (default: 50)

    Nextflow:
      -resume         Resume pipeline
    """.stripIndent()
}
if( params.help ) {
    helpMsg()
    System.exit(0)
}

// Log parameters
log.info """
=======================================================================
BatchITSx ${version}
=======================================================================
Input data path : ${params.input}
Output path     : ${params.outdir}
Method          : ${params.method}
Target group(s) : ${params.itsx_target}
Target region   : ${params.region}
""".stripIndent()

// =============================
// Channels
// =============================
ch_fasta   = Channel.fromPath(params.input)
ch_targets = Channel.fromList(itsx_targets_list)

// =============================
// Processes
// =============================

process dereplicate_vsearch {
    input:
    path fasta_file

    output:
    tuple path("${fasta_file}.unique"), path("${fasta_file}.uc")

    script:
    """
    ${params.vsearch_bin} \
      --fastx_uniques "${fasta_file}" \
      --fastaout "${fasta_file}.unique" \
      --uc "${fasta_file}.uc"
    """
}

process split_fasta {
    tag { "${derep_fasta_file.simpleName}" }

    input:
    path derep_fasta_file

    output:
    path "work_chunks/*.fasta", emit: split_chunks

    script:
    """
    mkdir -p work_chunks

    awk -v chunk_size=${params.itsx_chunksize} -v prefix="work_chunks/${derep_fasta_file.simpleName}" '
    BEGIN {
        file_index = 1
        seq_in_chunk = 0
        filename = prefix "." file_index ".fasta"
    }
    /^>/ {
        if (seq_in_chunk >= chunk_size) {
            close(filename)
            file_index++
            filename = prefix "." file_index ".fasta"
            seq_in_chunk = 0
        }
        seq_in_chunk++
    }
    { print >> filename }
    ' ${derep_fasta_file}

    echo "Created files:"
    ls -lh work_chunks
    echo "Sequence counts per file:"
    grep -c "^>" work_chunks/*.fasta
    """
}

process itsx {
    tag { "${fasta_chunk.simpleName}_${target}" }

    input:
    tuple path(fasta_chunk), val(target)

    output:
    tuple val(target),
          path("${fasta_chunk.baseName}.${target}.ITS2.full_and_partial.fasta"),
          path("${fasta_chunk.baseName}.${target}.positions.txt"),
          path("${fasta_chunk.baseName}.${target}.full_and_partial.fasta"),
          path("${fasta_chunk.baseName}.${target}_no_detections.txt", optional: true)
    
    script:
    """
    echo "Running ITSx on: ${fasta_chunk} (target: ${target})"

    ${params.itsx_bin} \
      -i "${fasta_chunk}" \
      -o "${fasta_chunk.baseName}.${target}" \
      -t "${target}" \
      --cpu ${task.cpus} \
      --save_regions "${params.itsx_regions}" \
      --partial "${params.itsx_partial}" \
      --detailed_results "${params.itsx_detailed}" \
      --concat "${params.itsx_concat}" \
      --preserve "${params.itsx_preserve}" \
      --graphical "${params.itsx_graphical}" \
      --complement "${params.itsx_complement}"

    # ensure file exists, even if ITSx didn’t create it
    touch "${fasta_chunk.baseName}.${target}_no_detections.txt"
    """
}

process itsx_merge {
    tag { target }

    input:
    tuple val(target),
          path(its2_fastas),
          path(position_files),
          path(full_fastas),
          path(no_detect_files)

    output:
    tuple(
        val(target),
        path("${target}.ITS2.merged.fasta"),
        path("${target}.positions.txt"),
        path("${target}.merged.fasta"),
        path("${target}.no_detections.txt")
    )

    script:
    """
    echo "Merging ITSx outputs for target: ${target}"

    cat ${its2_fastas}   > ${target}.ITS2.merged.fasta
    cat ${full_fastas}   > ${target}.merged.fasta
    cat ${position_files} > ${target}.positions.txt
    cat ${no_detect_files} > ${target}.no_detections.txt
    """
}

process parse_itsx_results_F {
    tag { "parse_F" }
    publishDir "${params.outdir}", mode: 'copy'

    input:
    tuple(
        val(target),
        path(its2_fasta),
        path(positions),
        path(full_fasta),
        path(no_detections)
    )

    output:
    path "parsed_results_f_${params.region}.fasta", emit: itsx_f_fasta
    path "parsed_results_f_${params.region}.log",   emit: itsx_f_log

    script:
    """
    parse_itsx_f.py \
      --its2 ${its2_fasta} \
      --positions ${positions} \
      --full ${full_fasta} \
      --output parsed_results_f_${params.region}.fasta \
      --region ${params.region}
    """
}

process parse_itsx_results_All {
    tag { "parse_All" }
    publishDir "${params.outdir}", mode: 'copy'

    input:
    tuple(
        path(parsed_fungi),
        path(its2_other),
        path(positions_fungi),
        path(positions_other),
        path(full_fasta_other),
        path(no_detect_fungi),
        path(no_detect_other)
    )

    output:
    path "parsed_results_all_${params.region}.fasta", emit: itsx_all_fasta
    path "parsed_results_all_${params.region}.log",   emit: itsx_all_log

    script:
    """
    parse_itsx_all.py \
      --its2_fungi ${parsed_fungi} \
      --its2_other ${its2_other} \
      --positions_fungi ${positions_fungi} \
      --positions_other ${positions_other} \
      --full ${full_fasta_other} \
      --no_detect_fungi ${no_detect_fungi} \
      --no_detect_other ${no_detect_other} \
      --output parsed_results_all_${params.region}.fasta \
      --region ${params.region}
    """
}

process concatenate_fastas {
    tag { "concat_${params.region}" }
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path f_fasta
    path all_fasta

    output:
    path "parsed_results_combined_${params.region}.fasta"

    script:
    """
    cat ${f_fasta} ${all_fasta} > parsed_results_combined_${params.region}.fasta
    """
}

process restore_duplicates {
    tag { "restore_${params.region}:${fasta_file.simpleName}" }
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path fasta_file
    path fasta_file_uc

    output:
    path "parsed_results_combined_${params.region}.full.fasta"

    script:
    """
    restore_duplicates.py \
        --infile "${fasta_file}" \
        --outfile "parsed_results_combined_${params.region}.full.fasta" \
        --uc "${fasta_file_uc}"
    """
}

process iupac_filter {
    tag { "filter_${params.region}" }
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path restored_fasta

    output:
    path "final_${params.region}.fasta"
    path "excluded_${params.region}.txt"

    script:
    """
    filter_non_iupac.py \
        --infile "${restored_fasta}" \
        --outfile "final_${params.region}.fasta" \
        --excluded "excluded_${params.region}.txt" \
        --maxN 6 \
        --maxIUPAC 16
    """
}

// =============================
// Workflow
// =============================
workflow {
    // create logs folder
    new File("${params.outdir}/logs").mkdirs()

    //
    // 1. Dereplicate and split input
    //
    // derep_fasta_ch = dereplicate_vsearch(ch_fasta)
    //                 .map { unique, uc -> unique }
    derep_out_ch = dereplicate_vsearch(ch_fasta)
    derep_fasta_ch = derep_out_ch.map { unique, uc -> unique }
    derep_uc_ch    = derep_out_ch.map { unique, uc -> uc }

    chunked_fasta_ch = split_fasta(derep_fasta_ch)
                       .split_chunks
                       .flatten()

    //
    // 2. Pair each chunk with ITSx targets
    //
    fasta_target_pairs = chunked_fasta_ch
                         .combine(ch_targets)
                         .map { f, t -> tuple(f, t) }

    //
    // 3. Run ITSx on each chunk/target
    //
    itsx_out = itsx(fasta_target_pairs)

    //
    // 4. Merge results per target
    //
    merged_ch = itsx_out
        .groupTuple(by: 0)   // group by target (F, All, etc.)
        .map { target, its2s, poss, fulls, nodets ->
            tuple(target, its2s, poss, fulls, nodets)
        }
        | itsx_merge

    //
    // 5. Split into fungi vs other
    //
    fungi_ch = merged_ch.filter { it[0] == "F" }
    other_ch = merged_ch.filter { it[0] == "All" }

    //
    // 6. Parse fungal separately
    //
    fungi_out = parse_itsx_results_F(fungi_ch)

    //
    // 7. Join fungi + other for All parser
    //
    paired_ch = fungi_out.itsx_f_fasta
        .combine(fungi_ch.toList())
        .combine(other_ch.toList())
        .map { f_parsed, f_raw, o ->
            def (targetF, its2F, posF, mergedF, noDetF) = f_raw
            def (targetO, its2O, posO, mergedO, noDetO) = o

            tuple(
                f_parsed,   // parsed_results_f_${params.region}.fasta
                its2O,
                posF,
                posO,
                mergedO,
                noDetF,
                noDetO
            )
        }

    all_out = parse_itsx_results_All(paired_ch)

    //
    // 8. Concatenate fungi + other for the final FASTA
    //
    concatenated = concatenate_fastas(
        fungi_out.itsx_f_fasta,
        all_out.itsx_all_fasta
    )

    //
    // 9. Restore duplicates using concatenated fasta + UC file
    //
    restore_out = restore_duplicates(
        concatenated,
        derep_uc_ch
    )

    //
    // 10. Filter sequences with too many ambiguous bases (N/IUPAC)
    //
    iupac_out = iupac_filter(
        restore_out
    )
}

// =============================
// Completion hooks
// =============================
workflow.onComplete {
    println "Pipeline completed at : $workflow.complete"
    println "Duration              : ${workflow.duration}"
    println "Execution status      : ${workflow.success ? 'All done!' : 'Failed'}"
}

workflow.onError {
    println "Pipeline execution stopped with the following message: ${workflow.errorMessage}"
}
