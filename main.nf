nextflow.enable.dsl=2

workflow {
// Load your samplesheet and emit tuples: (sample_id, R1, R2)
    reads_ch = Channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true)
        .map { row -> tuple(row.sample, file(row.R1), file(row.R2)) }
        .ifEmpty { error "No reads found in samplesheet: ${params.samplesheet}" }

    // QC & merge
    merged_reads_ch = qc_reads(reads_ch)

    // Run HUMAnN3 with database paths
    run_humann(merged_reads_ch)
}


// --- PROCESS: Run fastp for QC and merge paired reads ---

process qc_reads {
  label 'qc'
  publishDir "${params.outdir}/qc", mode: 'copy'
  container 'quay.io/biocontainers/fastp:0.24.0--heae3180_1'

  input:
  tuple val(sample_id), path(read1), path(read2)

  output:
  tuple val(sample_id), path("${sample_id}_merged.fastq.gz")

  script:
  """
  fastp \\
    --in1 ${read1} \\
    --in2 ${read2} \\
    --merge \\
    --merged_out ${sample_id}_merged.fastq.gz \\
    --thread ${params.fastp_threads} \\
    --detect_adapter_for_pe \\
    --html ${sample_id}_fastp.html \\
    --json ${sample_id}_fastp.json
  """
}

// --- PROCESS: Run HUMAnN3 (with gcsfuse) ---
process run_humann {
  label 'humann'
  publishDir "${params.outdir}/humann", mode: 'copy'
  container 'robmbio/humann3-gcsfuse:latest'

  input:
  tuple val(sample_id), path(merged_read)

  output:
  path "*", emit: results

  script:
  """
  set -euo pipefail

  mkdir -p /mnt/gcsdb/chocophlan /mnt/gcsdb/uniref /mnt/gcsdb/metaphlan

  gcsfuse --implicit-dirs --only-dir humann/chocophlan/chocophlan nextflow-bucket-batch /mnt/gcsdb/chocophlan
  gcsfuse --implicit-dirs --only-dir humann/uniref nextflow-bucket-batch /mnt/gcsdb/uniref
  gcsfuse --implicit-dirs --only-dir metaphlan-db nextflow-bucket-batch /mnt/gcsdb/metaphlan

  # Extract MetaPhlAn index base name
  BT2_FILE=\$(find /mnt/gcsdb/metaphlan -name "*.1.bt2l" | head -n 1)
  BT2_DB=\$(dirname "\$BT2_FILE")
  BT2_INDEX=\$(basename "\$BT2_FILE" | sed 's/\\.1\\.bt2l\$//')

  humann \\
    --input ${merged_read} \\
    --output ${sample_id} \\
    --threads ${params.humann_threads} \\
    --nucleotide-database /mnt/gcsdb/chocophlan \\
    --protein-database /mnt/gcsdb/uniref \\
    --metaphlan-options "--offline --bowtie2db /mnt/gcsdb/metaphlan --index mpa_vJun23_CHOCOPhlAnSGB_202403" \\
    --output-basename ${sample_id} \\
    --verbose
   """
}
