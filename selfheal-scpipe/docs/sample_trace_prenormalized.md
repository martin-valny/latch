# Run trace: `prenormalized_input` (planner: rules)

_Input is already log-normalized, but the run config says raw counts._

**[  0.04s] pipeline** from `load` -> **FAILED**
  - symptom: `load.value_scale`: expected raw integer counts but 100.0% of non-zero values are non-integer

**[  0.04s] agent** -> `get_run_status()`
  - hypothesis: start from the failed check

  - result: {"config": {"species": "human", "input_scale": "counts", "orientation": "cells_x_genes", "min_genes": 200, "max_pct_mito": 5.0, "min_cells_per_gene": 3, "target_sum": 10000.0}, "agent_editable": {"species": ["human", "mouse"], "input_scale": ["counts", "log1p"], "orientation": ["cells_x_genes", "genes_x_cells"]}, "last_run": {"status": "failed", "failed_step": "load", "steps_run": ["load"], "faile

**[  0.04s] agent** -> `inspect_matrix()`
  - hypothesis: non-integer values: normalized input or corruption?

  - result: {"shape_on_disk": [2638, 13714], "obs_name_examples": ["AAACATACAACCAC-1", "AAACATTGAGCTAC-1", "AAACATTGATCAGC-1", "AAACCGTGCTTCCG-1"], "var_name_examples": ["AL627309.1", "AP006222.2", "RP11-206L10.2", "RP11-206L10.9"], "obs_barcode_like_fraction": 1.0, "var_barcode_like_fraction": 0.0, "density": 0.0619, "n_nan": 0, "n_negative": 0, "min_nonzero": 0.7545999884605408, "max": 7.46959924697876, "fr

**[  0.08s] agent** -> `set_config("key": "input_scale", "value": "log1p")`
  - hypothesis: root cause of load.value_scale

  - result: {"changed": {"input_scale": ["counts", "log1p"]}, "must_relaunch_from_or_before": "load"}

**[  0.08s] agent** -> `relaunch("from_step": "load")`
  - hypothesis: verify the input_scale fix

**[  0.39s] pipeline** from `load` -> **PASSED**

  - result: {"status": "passed", "failed_step": null, "steps_run": ["load", "qc", "normalize", "annotate"], "failed_checks": [], "passed_checks": ["load.finite_nonnegative", "load.orientation", "load.value_scale", "qc.mito_genes_detected", "qc.cells_retained", "normalize.log_scale_range", "annotate.marker_overlap", "annotate.confident_fraction", "annotate.not_collapsed"], "relaunches_left": 2}


**Outcome: FIXED**: config changes: [{'key': 'input_scale', 'old': 'counts', 'new': 'log1p'}]

_agreement with published labels: 79.3%_
