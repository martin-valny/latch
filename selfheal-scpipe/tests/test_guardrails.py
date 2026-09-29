from healpipe.faults import SCENARIOS, materialize
from healpipe.pipeline import Pipeline
from healpipe.tools import Toolbox
from healpipe.trace import Trace


def _box(name, clean_h5ad, tmp_path):
    s = SCENARIOS[name]
    pipe = Pipeline(materialize(s, clean_h5ad, tmp_path), s.config)
    pipe.run()
    return Toolbox(pipe, Trace(name, "test"))


def test_qc_thresholds_are_not_editable(clean_h5ad, tmp_path):
    box = _box("shallow_sequencing", clean_h5ad, tmp_path)
    out = box.call("set_config", {"rationale": "x", "key": "min_genes", "value": "10"})
    assert "not agent-editable" in out["error"]
    assert box.pipeline.cfg.min_genes == 200


def test_relaunch_without_change_is_refused(clean_h5ad, tmp_path):
    box = _box("species_mislabel", clean_h5ad, tmp_path)
    out = box.call("relaunch", {"rationale": "x", "from_step": "load"})
    assert "no config change" in out["error"]
    assert box.relaunches == 0


def test_relaunch_must_cover_invalidated_steps(clean_h5ad, tmp_path):
    box = _box("transposed_matrix", clean_h5ad, tmp_path)
    box.call("set_config", {"rationale": "x", "key": "orientation", "value": "genes_x_cells"})
    out = box.call("relaunch", {"rationale": "x", "from_step": "qc"})
    assert "Relaunch from 'load'" in out["error"]
    assert box.call("relaunch", {"rationale": "x", "from_step": "load"})["status"] == "passed"


def test_relaunch_budget(clean_h5ad, tmp_path):
    box = _box("species_mislabel", clean_h5ad, tmp_path)
    box.max_relaunches = 1
    box.call("set_config", {"rationale": "x", "key": "input_scale", "value": "log1p"})
    box.call("relaunch", {"rationale": "x", "from_step": "load"})
    box.call("set_config", {"rationale": "x", "key": "input_scale", "value": "counts"})
    out = box.call("relaunch", {"rationale": "x", "from_step": "load"})
    assert "budget exhausted" in out["error"]
