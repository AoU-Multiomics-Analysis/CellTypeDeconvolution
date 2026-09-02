version 1.1

task AssembleTca {
  input {
    Array[File] shard_hdf5
    File shard_manifest
    File tca_expression
    File model
    File tca_weights
    File? covariates
    Boolean write_tsv = false
    String pipeline_version
    String docker_image
    Int cpu = 8
    String memory = "128 GB"
    Int disk_gb = 500
    Int preemptible_attempts = 0
    Int max_retries = 1
  }

  String covariates_argument = if defined(covariates) then "--covariates " + select_first([covariates]) else ""
  String write_tsv_argument = if write_tsv then "--write-tsv" else ""

  command <<<
    set -euo pipefail
    stage="assemble_tca"
    log="$stage.log"
    status=0
    printf 'stage=%s start_time=%s\n' "$stage" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" | tee -a "$log"
    trap 'status=$?; printf "stage=%s error_status=%s time=%s\\n" "$stage" "$status" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" | tee -a "$log"; exit "$status"' ERR
    printf '%s\n' ~{sep(' ', shard_hdf5)} > shards.txt
    Rscript /opt/celltype/scripts/assemble_tca_outputs.R \
      --shard-list shards.txt \
      --manifest '~{shard_manifest}' \
      --expression-log '~{tca_expression}' \
      --model '~{model}' \
      --weights '~{tca_weights}' \
      ~{covariates_argument} \
      --pipeline-version '~{pipeline_version}' \
      ~{write_tsv_argument} \
      --output-dir outputs \
      --log-file outputs/assemble_tca_outputs.log 2>&1 | tee -a "$log"
    printf 'stage=%s dimensions=%s outputs=%s completion_time=%s\n' \
      "$stage" "$(wc -l < outputs/assembled_outputs.tsv)" \
      "group_hdf5,group_tsv,reconstruction_by_sample,assembly_qc" \
      "$(date -u +%Y-%m-%dT%H:%M:%SZ)" | tee -a "$log"
  >>>

  output {
    Array[File] group_hdf5 = glob("outputs/*.h5")
    Array[File] group_tsv = glob("outputs/*.tsv.gz")
    File reconstruction_by_sample = "outputs/reconstruction_by_sample.tsv"
    File assembly_qc = "outputs/qc_summary.tsv"
    File log = "assemble_tca.log"
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

task BuildManifest {
  input {
    File output_inventory
    File reconstruction_by_sample
    File assembly_qc_summary
    File assembly_qc_plots
    Array[File] group_hdf5
    Array[File] group_tsv
    File model
    File model_log
    File mapping_report
    File excluded_genes
    File original_proportions
    File combined_proportions
    File tca_weights
    File filter_report
    File? parameters_json
    String pipeline_version
    String tca_version = "1.2.1"
    String container_image
    String docker_image
    Int cpu = 4
    String memory = "32 GB"
    Int disk_gb = 100
    Int preemptible_attempts = 1
    Int max_retries = 2
  }

  String parameters_argument = if defined(parameters_json) then "--parameters-json " + select_first([parameters_json]) else ""

  command <<<
    set -euo pipefail
    stage="build_manifest"
    log="$stage.log"
    status=0
    printf 'stage=%s start_time=%s\n' "$stage" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" | tee -a "$log"
    trap 'status=$?; printf "stage=%s error_status=%s time=%s\\n" "$stage" "$status" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" | tee -a "$log"; exit "$status"' ERR
    printf '%s\n' \
      '~{reconstruction_by_sample}' \
      ~{sep(' ', group_hdf5)} \
      ~{sep(' ', group_tsv)} \
      '~{model}' \
      '~{model_log}' \
      '~{mapping_report}' \
      '~{excluded_genes}' \
      '~{original_proportions}' \
      '~{combined_proportions}' \
      '~{tca_weights}' \
      '~{filter_report}' > supporting_inputs.txt
    Rscript /opt/celltype/scripts/build_manifest.R \
      --outputs '~{output_inventory}' \
      --pipeline-version '~{pipeline_version}' \
      --tca-version '~{tca_version}' \
      ~{parameters_argument} \
      --container-image '~{container_image}' \
      --output outputs/output_manifest.json \
      --log-file outputs/output_manifest.log 2>&1 | tee -a "$log"
    printf 'stage=%s dimensions=%s outputs=%s completion_time=%s\n' \
      "$stage" "$(wc -l < supporting_inputs.txt)" \
      "output_manifest,qc_summary,qc_plots,provenance" \
      "$(date -u +%Y-%m-%dT%H:%M:%SZ)" | tee -a "$log"
  >>>

  output {
    File output_manifest = "outputs/output_manifest.json"
    File qc_summary = assembly_qc_summary
    File qc_plots = assembly_qc_plots
    File provenance = output_inventory
    File log = "build_manifest.log"
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
