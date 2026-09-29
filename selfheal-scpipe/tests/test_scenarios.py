import pytest

from healpipe.agent import RulePlanner, grade, run_scenario
from healpipe.faults import SCENARIOS


@pytest.mark.parametrize("name", list(SCENARIOS))
def test_rule_planner_handles_every_scenario(name, clean_h5ad, tmp_path):
    outcome, pipe, trace = run_scenario(name, clean_h5ad, tmp_path, RulePlanner(), trace_dir=tmp_path)
    assert outcome.outcome == SCENARIOS[name].expected
    assert grade(name, outcome, pipe.cfg)
    assert (tmp_path / f"{name}.rules.trace.md").exists()


def test_every_fault_actually_breaks_the_run(clean_h5ad, tmp_path):
    from healpipe.faults import materialize
    from healpipe.pipeline import Pipeline

    for s in SCENARIOS.values():
        result = Pipeline(materialize(s, clean_h5ad, tmp_path), s.config).run()
        assert (result.status == "passed") == (s.expected == "clean"), s.name


def test_fixed_runs_recover_clean_annotation(clean_h5ad, tmp_path):
    _, clean, _ = run_scenario("clean", clean_h5ad, tmp_path, RulePlanner())
    _, fixed, _ = run_scenario("transposed_matrix", clean_h5ad, tmp_path, RulePlanner())
    assert clean.reference_agreement() > 0.9
    assert fixed.reference_agreement() == pytest.approx(clean.reference_agreement())
