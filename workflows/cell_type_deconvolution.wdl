version 1.1

import "tasks/expression.wdl" as expression_tasks
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
    String docker_image
    String pipeline_version = "0.1.0"
    Float min_lm22_overlap = 0.80
    Float dtangle_marker_fraction = 0.10
    Boolean dtangle_quantile_normalize = false
    Float group_mean_threshold = 0.0001
    Float zero_floor = 0.000001
    Int tca_max_iters = 10
    Int random_seed = 20260901

    Int prepare_cpu = 4
    String prepare_memory = "64 GB"
    Int prepare_disk_gb = 400
    Int prepare_preemptible_attempts = 2
    Int prepare_max_retries = 2
    Int dtangle_cpu = 4
    String dtangle_memory = "32 GB"
    Int dtangle_disk_gb = 100
    Int dtangle_preemptible_attempts = 2
    Int dtangle_max_retries = 2
    Int proportions_cpu = 2
    String proportions_memory = "16 GB"
    Int proportions_disk_gb = 50
    Int proportions_preemptible_attempts = 2
    Int proportions_max_retries = 2
    Int fit_cpu = 16
    String fit_memory = "192 GB"
    Int fit_disk_gb = 750
    Int fit_preemptible_attempts = 0
    Int fit_max_retries = 1
    Int export_cpu = 8
    String export_memory = "128 GB"
    Int export_disk_gb = 500
    Int export_preemptible_attempts = 0
    Int export_max_retries = 1
    Int manifest_cpu = 4
    String manifest_memory = "32 GB"
    Int manifest_disk_gb = 100
    Int manifest_preemptible_attempts = 1
    Int manifest_max_retries = 2
  }

  String tca_version = "1.2.1"
  call proportion_tasks.ValidateProportionMode {
    input:
      lm22 = lm22,
      precomputed_proportions = precomputed_proportions,
      docker_image = docker_image
  }

  String proportion_mode = ValidateProportionMode.selected_mode
  String dtangle_marker_method = "ratio"
  EffectiveParameters effective_parameters = EffectiveParameters {
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

  call expression_tasks.PrepareExpression {
    input:
      expression = expression,
      gtf = gtf,
      docker_image = docker_image,
      cpu = prepare_cpu,
      memory = prepare_memory,
      disk_gb = prepare_disk_gb,
      preemptible_attempts = prepare_preemptible_attempts,
      max_retries = prepare_max_retries
  }

  if (ValidateProportionMode.estimate_proportions) {
    call dtangle_tasks.RunDtangle {
      input:
        prepared_log2_cpm = PrepareExpression.prepared_dtangle_log2_cpm,
        lm22 = select_first([lm22]),
        min_overlap = min_lm22_overlap,
        marker_fraction = dtangle_marker_fraction,
        marker_method = dtangle_marker_method,
        quantile_normalize = dtangle_quantile_normalize,
        docker_image = docker_image,
        cpu = dtangle_cpu,
        memory = dtangle_memory,
        disk_gb = dtangle_disk_gb,
        preemptible_attempts = dtangle_preemptible_attempts,
        max_retries = dtangle_max_retries
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
      preemptible_attempts = proportions_preemptible_attempts,
      max_retries = proportions_max_retries
  }

  call tca_tasks.FitTca {
    input:
      prepared_log2_cpm = PrepareExpression.prepared_tca_log2_cpm,
      tca_weights = ProcessProportions.tca_weights,
      covariates = covariates,
      max_iters = tca_max_iters,
      random_seed = random_seed,
      docker_image = docker_image,
      cpu = fit_cpu,
      memory = fit_memory,
      disk_gb = fit_disk_gb,
      preemptible_attempts = fit_preemptible_attempts,
      max_retries = fit_max_retries
  }

  call tca_tasks.ExportTcaBeds {
    input:
      tca_expression = FitTca.tca_expression,
      coordinates = PrepareExpression.prepared_coordinates,
      model = FitTca.model,
      tca_weights = ProcessProportions.tca_weights,
      covariates = covariates,
      docker_image = docker_image,
      cpu = export_cpu,
      memory = export_memory,
      disk_gb = export_disk_gb,
      preemptible_attempts = export_preemptible_attempts,
      max_retries = export_max_retries
  }

  call qc_tasks.BuildManifest {
    input:
      cell_type_bed_inventory = ExportTcaBeds.cell_type_bed_inventory,
      cell_type_beds = ExportTcaBeds.cell_type_beds,
      export_qc_summary = ExportTcaBeds.qc_summary,
      export_qc_plots = ExportTcaBeds.qc_plots,
      model = FitTca.model,
      model_log = FitTca.model_log,
      mapping_report = PrepareExpression.mapping_report,
      excluded_genes = PrepareExpression.excluded_genes,
      original_proportions = ProcessProportions.original,
      combined_proportions = ProcessProportions.combined,
      tca_weights = ProcessProportions.tca_weights,
      filter_report = ProcessProportions.filter_report,
      dtangle_metadata = RunDtangle.metadata,
      parameters_json = effective_parameters_json,
      pipeline_version = pipeline_version,
      tca_version = tca_version,
      container_image = docker_image,
      docker_image = docker_image,
      cpu = manifest_cpu,
      memory = manifest_memory,
      disk_gb = manifest_disk_gb,
      preemptible_attempts = manifest_preemptible_attempts,
      max_retries = manifest_max_retries
  }

  output {
    File prepared_tca_cpm = PrepareExpression.prepared_tca_cpm
    File prepared_tca_log2_cpm = PrepareExpression.prepared_tca_log2_cpm
    File prepared_dtangle_cpm = PrepareExpression.prepared_dtangle_cpm
    File prepared_dtangle_log2_cpm = PrepareExpression.prepared_dtangle_log2_cpm
    File prepared_coordinates = PrepareExpression.prepared_coordinates
    File mapping_report = PrepareExpression.mapping_report
    File prepare_excluded_genes = PrepareExpression.excluded_genes
    File prepare_log = PrepareExpression.log
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
