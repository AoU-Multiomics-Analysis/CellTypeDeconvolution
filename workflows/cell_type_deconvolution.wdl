version 1.0

import "tasks/dtangle.wdl" as dtangle_tasks
import "tasks/proportions.wdl" as proportion_tasks
import "tasks/tca.wdl" as tca_tasks
import "tasks/qc.wdl" as qc_tasks

struct EffectiveParameters {
  String proportion_mode
  Float min_lm22_overlap
  Float dtangle_marker_fraction
  String dtangle_marker_method
  Boolean dtangle_quantile_normalize
  Float group_mean_threshold
  Float zero_floor
  Int tca_max_iters
  Int random_seed
  String scale
}

workflow cell_type_deconvolution {
  input {
    File expression
    File gtf
    File? lm22
    File? precomputed_proportions
    File? covariates
    String docker_image = "ghcr.io/aou-multiomics-analysis/celltypedeconvolution:latest"
    Int preemptible_attempts = 2
    Int max_retries = 2
    Float min_lm22_overlap = 0.80
    Float dtangle_marker_fraction = 0.10
    Boolean dtangle_quantile_normalize = false
    Float group_mean_threshold = 0.0001
    Float zero_floor = 0.000001
    Int tca_max_iters = 10
    Int random_seed = 20260901

    Int dtangle_cpu = 4
    String dtangle_memory = "32 GB"
    Int dtangle_disk_gb = 100
    Int proportions_cpu = 2
    String proportions_memory = "16 GB"
    Int proportions_disk_gb = 50
    Int fit_cpu = 16
    String fit_memory = "192 GB"
    Int fit_disk_gb = 750
    Int export_cpu = 8
    String export_memory = "128 GB"
    Int export_disk_gb = 500
    Int manifest_cpu = 4
    String manifest_memory = "32 GB"
    Int manifest_disk_gb = 100
  }

  String tca_version = "1.2.1"
  call proportion_tasks.ValidateProportionMode {
    input:
      lm22 = lm22,
      precomputed_proportions = precomputed_proportions,
      docker_image = docker_image,
      preemptible_attempts = preemptible_attempts,
      max_retries = max_retries
  }

  String proportion_mode = ValidateProportionMode.selected_mode
  String dtangle_marker_method = "ratio"
  EffectiveParameters effective_parameters = object {
    proportion_mode: proportion_mode,
    min_lm22_overlap: min_lm22_overlap,
    dtangle_marker_fraction: dtangle_marker_fraction,
    dtangle_marker_method: dtangle_marker_method,
    dtangle_quantile_normalize: dtangle_quantile_normalize,
    group_mean_threshold: group_mean_threshold,
    zero_floor: zero_floor,
    tca_max_iters: tca_max_iters,
    random_seed: random_seed,
    scale: "log2_cpm"
  }
  File effective_parameters_json = write_json(effective_parameters)

  if (ValidateProportionMode.estimate_proportions) {
    call dtangle_tasks.RunDtangle {
      input:
        expression = expression,
        gtf = gtf,
        lm22 = select_first([lm22]),
        min_overlap = min_lm22_overlap,
        marker_fraction = dtangle_marker_fraction,
        marker_method = dtangle_marker_method,
        quantile_normalize = dtangle_quantile_normalize,
        docker_image = docker_image,
        cpu = dtangle_cpu,
        memory = dtangle_memory,
        disk_gb = dtangle_disk_gb,
        preemptible_attempts = preemptible_attempts,
        max_retries = max_retries
    }
  }

  File proportions_for_processing = if ValidateProportionMode.estimate_proportions
    then select_first([RunDtangle.proportions])
    else select_first([precomputed_proportions])

  call proportion_tasks.ProcessProportions {
    input:
      proportions = proportions_for_processing,
      mean_threshold = group_mean_threshold,
      zero_floor = zero_floor,
      docker_image = docker_image,
      cpu = proportions_cpu,
      memory = proportions_memory,
      disk_gb = proportions_disk_gb,
      preemptible_attempts = preemptible_attempts,
      max_retries = max_retries
  }

  call tca_tasks.FitTca {
    input:
      expression = expression,
      tca_weights = ProcessProportions.tca_weights,
      covariates = covariates,
      max_iters = tca_max_iters,
      random_seed = random_seed,
      docker_image = docker_image,
      cpu = fit_cpu,
      memory = fit_memory,
      disk_gb = fit_disk_gb,
      preemptible_attempts = preemptible_attempts,
      max_retries = max_retries
  }

  call tca_tasks.ExportTcaBeds {
    input:
      tca_expression = FitTca.tca_expression,
      expression = expression,
      model = FitTca.model,
      tca_weights = ProcessProportions.tca_weights,
      covariates = covariates,
      docker_image = docker_image,
      cpu = export_cpu,
      memory = export_memory,
      disk_gb = export_disk_gb,
      preemptible_attempts = preemptible_attempts,
      max_retries = max_retries
  }

  call qc_tasks.BuildManifest {
    input:
      cell_type_bed_inventory = ExportTcaBeds.cell_type_bed_inventory,
      cell_type_beds = ExportTcaBeds.cell_type_beds,
      export_qc_summary = ExportTcaBeds.qc_summary,
      export_qc_plots = ExportTcaBeds.qc_plots,
      model = FitTca.model,
      model_log = FitTca.model_log,
      original_proportions = ProcessProportions.original,
      combined_proportions = ProcessProportions.combined,
      tca_weights = ProcessProportions.tca_weights,
      filter_report = ProcessProportions.filter_report,
      dtangle_metadata = RunDtangle.metadata,
      parameters_json = effective_parameters_json,
      tca_version = tca_version,
      container_image = docker_image,
      docker_image = docker_image,
      cpu = manifest_cpu,
      memory = manifest_memory,
      disk_gb = manifest_disk_gb,
      preemptible_attempts = preemptible_attempts,
      max_retries = max_retries
  }

  output {
    File proportion_mode_validation_log = ValidateProportionMode.log

    File? estimated_proportions = RunDtangle.proportions
    File? dtangle_markers = RunDtangle.markers
    File? dtangle_metadata = RunDtangle.metadata
    File? dtangle_overlap_report = RunDtangle.overlap_report
    File? transformed_lm22 = RunDtangle.transformed_lm22
    File? dtangle_shared_bulk = RunDtangle.shared_bulk
    File? dtangle_log = RunDtangle.log

    File proportions_lm22 = ProcessProportions.original
    File proportions_combined = ProcessProportions.combined
    File tca_weights = ProcessProportions.tca_weights
    File cell_group_filter_report = ProcessProportions.filter_report
    File proportions_log = ProcessProportions.log

    File tca_model = FitTca.model
    File tca_model_log = FitTca.model_log
    File tca_expression = FitTca.tca_expression
    File tca_excluded_genes = FitTca.excluded_genes
    File fit_tca_log = FitTca.log

    Array[File] cell_type_beds = ExportTcaBeds.cell_type_beds
    File cell_type_bed_inventory = ExportTcaBeds.cell_type_bed_inventory
    File reconstruction_by_sample = ExportTcaBeds.reconstruction_by_sample
    File qc_summary = BuildManifest.qc_summary
    File qc_plots = ExportTcaBeds.qc_plots
    File export_log = ExportTcaBeds.log
    File export_detail_log = ExportTcaBeds.export_detail_log
    File output_manifest = BuildManifest.output_manifest
    File output_inventory = BuildManifest.provenance
    File manifest_log = BuildManifest.log
    File effective_parameters_file = effective_parameters_json
  }
}
