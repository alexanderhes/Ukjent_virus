/*
 * SPADES_ASSEMBLY
 *
 * De novo assembly of per-species read sets using SPAdes.
 *
 * Two paths:
 *   ≥ params.validate_min_reads → SPAdes assembly (mode driven by spades_mode input);
 *                                  contigs shorter than params.validate_min_contig_len
 *                                  are filtered out. If SPAdes produces no usable
 *                                  contigs the output is an empty FASTA.
 *   < params.validate_min_reads → empty query FASTA; BLASTN_VALIDATE will skip
 *                                  BLAST and write a header-only result.
 *
 * In all cases the output file is named {id}_query.fasta so
 * the downstream BLASTN_VALIDATE module receives a consistent input.
 */

process SPADES_ASSEMBLY {
    tag "${meta.id}"
    label 'process_high'

    publishDir "${params.outdir}/validation/${meta.id}/assembly", mode: 'copy'

    input:
    tuple val(meta), path(r1), path(r2)
    val spades_mode   // assembly mode: rnaviral | meta | isolate (passed from main.nf)

    output:
    tuple val(meta), path("${meta.id}_query.fasta"), emit: query

    script:
    def min_reads  = params.validate_min_reads
    def min_contig = params.validate_min_contig_len
    def mem_gb     = task.memory.toGiga()
    """
    # Count reads in R1 using seqkit (consistent with HOST_FILTER)
    count=\$(seqkit stats -T -j ${task.cpus} ${r1} | awk 'NR==2{print \$4}')

    if [ "\${count}" -ge "${min_reads}" ]; then
        spades.py \\
            --${spades_mode} \\
            -1 ${r1} \\
            -2 ${r2} \\
            -o spades_out \\
            -t ${task.cpus} \\
            -m ${mem_gb}

        # --meta normally produces contigs.fasta; guard in case assembly fails
        if [ -f spades_out/contigs.fasta ]; then
            seqkit seq --min-len ${min_contig} spades_out/contigs.fasta \
                > ${meta.id}_query.fasta
        else
            touch ${meta.id}_query.fasta
        fi

        # If SPAdes produced no contigs passing the length filter, leave the
        # output empty — downstream BLASTN_VALIDATE will handle it gracefully.
        if [ ! -s "${meta.id}_query.fasta" ]; then
            echo "WARNING: SPAdes produced no contigs >= ${min_contig} bp for ${meta.id} -- outputting empty file"
        fi
    else
        # Too few reads for assembly — output empty file, skip BLAST
        echo "WARNING: only \${count} reads for ${meta.id} (< ${min_reads}) -- skipping assembly and BLAST"
        touch ${meta.id}_query.fasta
    fi
    """
}
