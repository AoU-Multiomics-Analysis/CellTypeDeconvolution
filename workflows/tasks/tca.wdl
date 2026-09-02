version 1.1

task FitTca {
  input {
    File prepared_log2_cpm
    File tca_weights
    File? covariates
    Int shard_size = 500
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
      --shard-size '~{shard_size}' \
      --max-iters '~{max_iters}' \
      --random-seed '~{random_seed}' \
      --output-dir outputs 2>&1 | tee -a "$log"
    printf 'stage=%s dimensions=%s outputs=%s completion_time=%s\n' \
      "$stage" "$(wc -l < outputs/gene_shard_manifest.tsv)" \
      "model,model_log,tca_expression,excluded_genes,shard_manifest,shards" \
      "$(date -u +%Y-%m-%dT%H:%M:%SZ)" | tee -a "$log"
  >>>

  output {
    File model = "outputs/tca_model.rds"
    File model_log = "outputs/tca_model.log"
    File tca_expression = "outputs/tca_expression.tsv.gz"
    File excluded_genes = "outputs/tca_excluded_genes.tsv"
    File shard_manifest = "outputs/gene_shard_manifest.tsv"
    Array[File] shards = glob("outputs/shard_*.txt")
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

task ExtractTcaShard {
  input {
    File tca_expression
    File model
    File shard
    String docker_image
    Int cpu = 8
    String memory = "64 GB"
    Int disk_gb = 200
    Int preemptible_attempts = 2
    Int max_retries = 2
  }

  command <<<
    set -euo pipefail
    stage="extract_tca_shard"
    log="$stage.log"
    status=0
    shard_name="$(basename '~{shard}')"
    shard_id="${shard_name#shard_}"
    shard_id="${shard_id%.txt}"
    shard_id="$((10#$shard_id))"
    printf 'stage=%s start_time=%s\n' "$stage" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" | tee -a "$log"
    trap 'status=$?; printf "stage=%s error_status=%s time=%s\\n" "$stage" "$status" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" | tee -a "$log"; exit "$status"' ERR
    Rscript /opt/celltype/scripts/extract_tca_shard.R \
      --expression-log '~{tca_expression}' \
      --model '~{model}' \
      --genes '~{shard}' \
      --shard-id "$shard_id" \
      --num-cores '~{cpu}' \
      --output outputs/tca_shard.h5 \
      --log-file outputs/extract_tca_shard.log 2>&1 | tee -a "$log"
    printf 'stage=%s dimensions=%s outputs=%s completion_time=%s\n' \
      "$stage" "$(wc -l < '~{shard}')" "shard_hdf5" \
      "$(date -u +%Y-%m-%dT%H:%M:%SZ)" | tee -a "$log"
  >>>

  output {
    File shard_hdf5 = "outputs/tca_shard.h5"
    File log = "extract_tca_shard.log"
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
