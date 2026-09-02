version 1.1

task PrepareExpression {
  input {
    File expression
    File gtf
    String docker_image
    Int cpu = 4
    String memory = "64 GB"
    Int disk_gb = 400
    Int preemptible_attempts = 2
    Int max_retries = 2
  }

  command <<<
    set -euo pipefail
    stage="prepare_expression"
    log="$stage.log"
    status=0
    printf 'stage=%s start_time=%s\n' "$stage" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" | tee -a "$log"
    trap 'status=$?; printf "stage=%s error_status=%s time=%s\\n" "$stage" "$status" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" | tee -a "$log"; exit "$status"' ERR
    Rscript /opt/celltype/scripts/prepare_expression.R \
      --expression '~{expression}' \
      --gtf '~{gtf}' \
      --output-dir outputs 2>&1 | tee -a "$log"
    printf 'stage=%s dimensions=%s outputs=%s completion_time=%s\n' \
      "$stage" "$(wc -l < outputs/gene_mapping_report.tsv)" \
      "prepared_cpm,prepared_log2_cpm,mapping_report,excluded_genes" \
      "$(date -u +%Y-%m-%dT%H:%M:%SZ)" | tee -a "$log"
  >>>

  output {
    File prepared_cpm = "outputs/prepared_cpm.tsv.gz"
    File prepared_log2_cpm = "outputs/prepared_log2_cpm.tsv.gz"
    File mapping_report = "outputs/gene_mapping_report.tsv"
    File excluded_genes = "outputs/excluded_genes.tsv"
    File log = "prepare_expression.log"
  }

  runtime {
    docker: docker_image
    cpu: cpu
    memory: memory
    disks: "local-disk ~{disk_gb} HDD"
    preemptible: preemptible_attempts
    maxRetries: max_retries
  }
}
