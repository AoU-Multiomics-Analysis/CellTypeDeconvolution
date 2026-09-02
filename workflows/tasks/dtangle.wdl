version 1.1

task RunDtangle {
  input {
    File prepared_log2_cpm
    File lm22
    Float min_overlap = 0.80
    Float marker_fraction = 0.10
    String marker_method = "ratio"
    Boolean quantile_normalize = false
    String docker_image
    Int cpu = 4
    String memory = "32 GB"
    Int disk_gb = 100
    Int preemptible_attempts = 2
    Int max_retries = 2
  }

  String quantile_normalize_argument = if quantile_normalize then "--quantile-normalize" else ""

  command <<<
    set -euo pipefail
    stage="run_dtangle"
    log="$stage.log"
    status=0
    printf 'stage=%s start_time=%s\n' "$stage" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" | tee -a "$log"
    trap 'status=$?; printf "stage=%s error_status=%s time=%s\\n" "$stage" "$status" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" | tee -a "$log"; exit "$status"' ERR
    Rscript /opt/celltype/scripts/run_dtangle.R \
      --bulk-log '~{prepared_log2_cpm}' \
      --lm22 '~{lm22}' \
      --min-overlap '~{min_overlap}' \
      --marker-fraction '~{marker_fraction}' \
      --marker-method '~{marker_method}' \
      ~{quantile_normalize_argument} \
      --output-dir outputs 2>&1 | tee -a "$log"
    printf 'stage=%s dimensions=%s outputs=%s completion_time=%s\n' \
      "$stage" "$(wc -l < outputs/dtangle_proportions.tsv)" \
      "proportions,markers,metadata,overlap_report,transformed_lm22,shared_bulk" \
      "$(date -u +%Y-%m-%dT%H:%M:%SZ)" | tee -a "$log"
  >>>

  output {
    File proportions = "outputs/dtangle_proportions.tsv"
    File markers = "outputs/dtangle_markers.tsv"
    File metadata = "outputs/dtangle_metadata.json"
    File overlap_report = "outputs/dtangle_overlap.tsv"
    File transformed_lm22 = "outputs/dtangle_lm22_log.tsv.gz"
    File shared_bulk = "outputs/dtangle_shared_bulk.tsv.gz"
    File log = "run_dtangle.log"
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
