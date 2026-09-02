version 1.1

task ProcessProportions {
  input {
    File proportions
    Float mean_threshold = 0.0001
    Float zero_floor = 0.000001
    String docker_image
    Int cpu = 2
    String memory = "16 GB"
    Int disk_gb = 50
    Int preemptible_attempts = 2
    Int max_retries = 2
  }

  command <<<
    set -euo pipefail
    stage="process_proportions"
    log="$stage.log"
    status=0
    printf 'stage=%s start_time=%s\n' "$stage" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" | tee -a "$log"
    trap 'status=$?; printf "stage=%s error_status=%s time=%s\\n" "$stage" "$status" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" | tee -a "$log"; exit "$status"' ERR
    Rscript /opt/celltype/scripts/process_proportions.R \
      --proportions '~{proportions}' \
      --mean-threshold '~{mean_threshold}' \
      --zero-floor '~{zero_floor}' \
      --output-dir outputs 2>&1 | tee -a "$log"
    printf 'stage=%s dimensions=%s outputs=%s completion_time=%s\n' \
      "$stage" "$(wc -l < outputs/proportions_tca_weights.tsv)" \
      "original,combined,tca_weights,filter_report" \
      "$(date -u +%Y-%m-%dT%H:%M:%SZ)" | tee -a "$log"
  >>>

  output {
    File original = "outputs/proportions_lm22.tsv"
    File combined = "outputs/proportions_combined.tsv"
    File tca_weights = "outputs/proportions_tca_weights.tsv"
    File filter_report = "outputs/cell_group_filter_report.tsv"
    File log = "process_proportions.log"
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
