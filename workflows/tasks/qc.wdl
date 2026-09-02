version 1.1

task BuildManifest {
  input {
    File cell_type_bed_inventory
    Array[File] cell_type_beds
    File export_qc_summary
    File export_qc_plots
    File model
    File model_log
    File mapping_report
    File excluded_genes
    File original_proportions
    File combined_proportions
    File tca_weights
    File filter_report
    File? dtangle_metadata
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

  String parameters_path = if defined(parameters_json) then select_first([parameters_json]) else ""
  String dtangle_metadata_path = if defined(dtangle_metadata) then select_first([dtangle_metadata]) else ""

  command <<<
    set -euo pipefail
    stage="build_manifest"
    log="$stage.log"
    status=0
    printf 'stage=%s start_time=%s\n' "$stage" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" | tee -a "$log"
    trap 'status=$?; printf "stage=%s error_status=%s time=%s\\n" "$stage" "$status" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" | tee -a "$log"; exit "$status"' ERR
    localized_bed_files="$PWD/localized_bed_files"
    checksum_inventory=checksum_inventory.localized.tsv
    public_inventory=outputs/output_inventory.tsv
    parameters_path="~{parameters_path}"
    parameters_arguments=()
    if [[ -n "$parameters_path" ]]; then
      parameters_arguments=(--parameters-json "$parameters_path")
    fi
    dtangle_metadata_path="~{dtangle_metadata_path}"
    dtangle_metadata_arguments=()
    if [[ -n "$dtangle_metadata_path" ]]; then
      dtangle_metadata_arguments=(--dtangle-metadata "$dtangle_metadata_path")
    fi
    mkdir -p "$localized_bed_files" outputs
    stage_localized_file() {
      source_path="$1"
      source_name="$(basename -- "$source_path")"
      if [[ -z "$source_name" || "$source_name" == "." || "$source_name" == ".." ]]; then
        printf 'stage=%s error_status=1 time=%s message=unsafe_output_basename\n' \
          "$stage" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" | tee -a "$log"
        exit 1
      fi
      source_dir="$(dirname -- "$source_path")"
      source_absolute="$(cd -P -- "$source_dir" && pwd -P)/$source_name"
      destination="$localized_bed_files/$source_name"
      if [[ ! -e "$source_absolute" ]]; then
        printf 'stage=%s error_status=1 time=%s message=missing_source_output\n' \
          "$stage" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" | tee -a "$log"
        exit 1
      fi
      if [[ -e "$destination" || -L "$destination" ]]; then
        printf 'stage=%s error_status=1 time=%s message=duplicate_output_basename\n' \
          "$stage" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" | tee -a "$log"
        exit 1
      fi
      ln -s "$source_absolute" "$destination"
      if [[ ! -e "$destination" ]]; then
        printf 'stage=%s error_status=1 time=%s message=failed_to_stage_output\n' \
          "$stage" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" | tee -a "$log"
        exit 1
      fi
    }
    cat > bed_output_sources.txt <<'BED_OUTPUT_SOURCES'
~{sep('\n', cell_type_beds)}
BED_OUTPUT_SOURCES
    : > array_basenames.txt
    while IFS= read -r source_path; do
      [[ -z "$source_path" ]] && continue
      stage_localized_file "$source_path"
      basename -- "$source_path" >> array_basenames.txt
    done < bed_output_sources.txt
    cp '~{cell_type_bed_inventory}' "$public_inventory"
    awk -F '\t' -v OFS='\t' -v localized_dir="$localized_bed_files" '
      NR == 1 {
        path_column = 0
        for (column = 1; column <= NF; column++) {
          if ($column == "path") {
            path_column = column
            break
          }
        }
        if (path_column == 0) {
          print "Missing path column in cell-type BED inventory" > "/dev/stderr"
          exit 1
        }
        print
        next
      }
      {
        output_basename = $path_column
        sub(/^.*\//, "", output_basename)
        if (output_basename == "") {
          print "Empty path in cell-type BED inventory" > "/dev/stderr"
          exit 1
        }
        if ($path_column != output_basename) {
          print "Cell-type BED inventory paths must be basenames" > "/dev/stderr"
          exit 1
        }
        $path_column = localized_dir "/" output_basename
        print
      }
    ' "$public_inventory" > "$checksum_inventory"
    awk -F '\t' '
      NR == 1 {
        path_column = 0
        for (column = 1; column <= NF; column++) {
          if ($column == "path") {
            path_column = column
            break
          }
        }
        if (path_column == 0) {
          print "Missing path column in localized BED inventory" > "/dev/stderr"
          exit 1
        }
        next
      }
      {
        print $path_column
        output_basename = $path_column
        sub(/^.*\//, "", output_basename)
        print output_basename > "inventory_basenames.txt"
      }
    ' "$checksum_inventory" > rewritten_paths.txt
    while IFS= read -r rewritten_path; do
      if [[ ! -e "$rewritten_path" ]]; then
        printf 'stage=%s error_status=1 time=%s message=missing_localized_output\n' \
          "$stage" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" | tee -a "$log"
        exit 1
      fi
    done < rewritten_paths.txt
    LC_ALL=C sort array_basenames.txt > array_basenames.sorted.txt
    LC_ALL=C sort inventory_basenames.txt > inventory_basenames.sorted.txt
    if ! cmp -s array_basenames.sorted.txt inventory_basenames.sorted.txt; then
      printf 'stage=%s error_status=1 time=%s message=inventory_array_mismatch\n' \
        "$stage" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" | tee -a "$log"
      exit 1
    fi
    printf '%s\n' \
      '~{export_qc_summary}' \
      '~{export_qc_plots}' \
      '~{model}' \
      '~{model_log}' \
      '~{mapping_report}' \
      '~{excluded_genes}' \
      '~{original_proportions}' \
      '~{combined_proportions}' \
      '~{tca_weights}' \
      '~{filter_report}' > supporting_inputs.txt
    Rscript /opt/celltype/scripts/build_manifest.R \
      --outputs "$checksum_inventory" \
      --export-qc-summary '~{export_qc_summary}' \
      --mapping-report '~{mapping_report}' \
      --original-proportions '~{original_proportions}' \
      --combined-proportions '~{combined_proportions}' \
      --tca-weights '~{tca_weights}' \
      --filter-report '~{filter_report}' \
      --model '~{model}' \
      --model-log '~{model_log}' \
      "${dtangle_metadata_arguments[@]}" \
      --pipeline-version '~{pipeline_version}' \
      --tca-version '~{tca_version}' \
      "${parameters_arguments[@]}" \
      --container-image '~{container_image}' \
      --output outputs/output_manifest.json \
      --qc-output outputs/qc_summary.tsv \
      --log-file outputs/output_manifest.log 2>&1 | tee -a "$log"
    bed_count="$(awk 'END { print NR - 1 }' "$public_inventory")"
    printf 'stage=%s dimensions=beds:%s outputs=%s completion_time=%s\n' \
      "$stage" "$bed_count" \
      "output_manifest,qc_summary,qc_plots,provenance" \
      "$(date -u +%Y-%m-%dT%H:%M:%SZ)" | tee -a "$log"
  >>>

  output {
    File output_manifest = "outputs/output_manifest.json"
    File qc_summary = "outputs/qc_summary.tsv"
    File qc_plots = export_qc_plots
    File provenance = "outputs/output_inventory.tsv"
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
