version 1.1

task FitTca {
  input {
    File prepared_log2_cpm
    File tca_weights
    File? covariates
    Int max_iters = 10
    Int random_seed = 20260901
    String docker_image
    Int cpu = 16
    String memory = "192 GB"
    Int disk_gb = 750
    Int preemptible_attempts = 0
    Int max_retries = 1
  }

  String covariates_argument = if defined(covariates) then "--covariates " + select_first([covariates]) else ""

  command <<<
    set -euo pipefail
    stage="fit_tca"
    log="$stage.log"
    status=0
    printf 'stage=%s start_time=%s\n' "$stage" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" | tee -a "$log"
    trap 'status=$?; printf "stage=%s error_status=%s time=%s\\n" "$stage" "$status" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" | tee -a "$log"; exit "$status"' ERR
    Rscript /opt/celltype/scripts/fit_tca.R \
      --expression-log '~{prepared_log2_cpm}' \
      --weights '~{tca_weights}' \
      ~{covariates_argument} \
      --num-cores '~{cpu}' \
      --max-iters '~{max_iters}' \
      --random-seed '~{random_seed}' \
      --output-dir outputs 2>&1 | tee -a "$log"
    printf 'stage=%s dimensions=%s outputs=%s completion_time=%s\n' \
      "$stage" "$(wc -l < outputs/tca_expression.tsv.gz)" \
      "model,model_log,tca_expression,excluded_genes" \
      "$(date -u +%Y-%m-%dT%H:%M:%SZ)" | tee -a "$log"
  >>>

  output {
    File model = "outputs/tca_model.rds"
    File model_log = "outputs/tca_model.log"
    File tca_expression = "outputs/tca_expression.tsv.gz"
    File excluded_genes = "outputs/tca_excluded_genes.tsv"
    File log = "fit_tca.log"
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

task ExportTcaBeds {
  input {
    File tca_expression
    File coordinates
    File model
    File tca_weights
    File? covariates
    String docker_image
    Int cpu = 8
    String memory = "128 GB"
    Int disk_gb = 500
    Int preemptible_attempts = 0
    Int max_retries = 1
  }

  String covariates_argument = if defined(covariates) then "--covariates " + select_first([covariates]) else ""

  command <<<
    set -euo pipefail
    stage="export_tca_beds"
    log="$stage.log"
    status=0
    printf 'stage=%s start_time=%s\n' "$stage" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" | tee -a "$log"
    trap 'status=$?; printf "stage=%s error_status=%s time=%s\\n" "$stage" "$status" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" | tee -a "$log"; exit "$status"' ERR
    Rscript /opt/celltype/scripts/export_tca_beds.R \
      --expression-log '~{tca_expression}' \
      --coordinates '~{coordinates}' \
      --model '~{model}' \
      --weights '~{tca_weights}' \
      ~{covariates_argument} \
      --num-cores '~{cpu}' \
      --output-dir outputs \
      --log-file outputs/export_tca_beds.log 2>&1 | tee -a "$log"
    printf 'stage=%s dimensions=%s outputs=%s completion_time=%s\n' \
      "$stage" "$(wc -l < outputs/cell_type_bed_inventory.tsv)" \
      "cell_type_beds,inventory,reconstruction,qc_summary,qc_plots" \
      "$(date -u +%Y-%m-%dT%H:%M:%SZ)" | tee -a "$log"
  >>>

  output {
    Array[File] cell_type_beds = glob("outputs/*.bed.gz")
    File cell_type_bed_inventory = "outputs/cell_type_bed_inventory.tsv"
    File reconstruction_by_sample = "outputs/reconstruction_by_sample.tsv"
    File qc_summary = "outputs/qc_summary.tsv"
    File qc_plots = "outputs/qc_plots.pdf"
    File export_detail_log = "outputs/export_tca_beds.log"
    File log = "export_tca_beds.log"
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
