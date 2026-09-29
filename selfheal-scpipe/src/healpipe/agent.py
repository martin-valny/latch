"""The recovery agent: a planner that picks tools until the run passes or is escalated.

Two planners share one tool surface and one trace format:

* ClaudePlanner: an LLM plans. Claude reads the failure report, forms
  hypotheses, calls the inspection tools, and then fixes, relaunches, or
  escalates.
* RulePlanner: a hand-written decision table. It is a deterministic baseline,
  and it lets tests and the eval run offline with no API key.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

from .config import RunConfig
from .faults import SCENARIOS, materialize
from .pipeline import Pipeline
from .tools import TOOL_SPECS, Toolbox
from .trace import Trace

SYSTEM_PROMPT = """\
You are the on-call recovery agent for a single-cell RNA-seq pipeline \
(load -> qc -> normalize -> annotate). A monitor check has just failed.

Your job: work out the root cause from evidence, apply the smallest correct \
fix, and relaunch. If there is no safe fix, escalate to a human.

How to work:
- Start from the failed check. Gather evidence with the read-only inspect tools \
before changing anything. State your hypothesis in each call's `rationale`.
- You can only change sample-sheet metadata (species, input_scale, orientation). \
QC thresholds belong to the analysis owner. If the data only "passes" with looser \
thresholds, the data is the problem: escalate.
- Fix the cause, not the symptom. One config change per hypothesis, then relaunch \
from the earliest step the change affects.
- Corrupted values (negative counts, NaNs) and data that is too shallow or low-quality \
are not config problems: escalate with the numbers that show it.
- If a relaunch surfaces a new, different failure, diagnose that one the same way.
- When the pipeline passes, or after you escalate, reply with a two-sentence summary \
and stop calling tools.
"""


@dataclass
class AgentOutcome:
    outcome: str  # "clean" | "fixed" | "escalated" | "gave_up"
    changes: list[dict]
    relaunches: int
    tool_calls: int
    detail: str = ""


def _finish(box: Toolbox, trace: Trace, gave_up_reason: str = "") -> AgentOutcome:
    n_calls = sum(e["kind"] == "tool_call" for e in trace.events)
    if box.escalation:
        o = AgentOutcome("escalated", box.changes, box.relaunches, n_calls, box.escalation["reason"])
    elif box.pipeline.last_result.status == "passed":
        o = AgentOutcome("fixed", box.changes, box.relaunches, n_calls, f"config changes: {box.changes}")
    else:
        o = AgentOutcome("gave_up", box.changes, box.relaunches, n_calls, gave_up_reason)
    trace.add("outcome", outcome=o.outcome, detail=o.detail)
    return o


# ------------------------------------------------------------------ Claude


class ClaudePlanner:
    name = "claude"

    def __init__(self, model: str = "claude-opus-5-5", effort: str = "medium", max_turns: int = 12):
        import anthropic

        self.client = anthropic.Anthropic()
        self.model, self.effort, self.max_turns = model, effort, max_turns

    def recover(self, box: Toolbox, trace: Trace) -> AgentOutcome:
        failure = box.pipeline.last_result.to_dict()
        messages = [
            {
                "role": "user",
                "content": "Pipeline run failed. Monitor report:\n"
                + json.dumps(failure, indent=2)
                + "\nDiagnose and recover.",
            }
        ]
        tools = [{**t, "strict": True} for t in TOOL_SPECS]

        for _ in range(self.max_turns):
            response = self.client.beta.messages.create(
                model=self.model,
                max_tokens=16000,
                system=SYSTEM_PROMPT,
                tools=tools,
                messages=messages,
                thinking={"type": "adaptive"},
                output_config={"effort": self.effort},
                # If the primary model declines a request, the API re-runs it
                # server-side on a fallback model and routes by category.
                betas=["server-side-fallback-2026-07-01"],
                fallbacks="default",
            )
            if response.stop_reason in ("refusal", "max_tokens"):
                box.call("escalate", {"rationale": "planner could not continue", "reason": f"planner stopped: {response.stop_reason}", "evidence": ""})
                return _finish(box, trace)

            # Append content unchanged (thinking blocks included) to keep history append-only.
            messages.append({"role": "assistant", "content": response.content})
            tool_uses = [b for b in response.content if b.type == "tool_use"]
            if not tool_uses:
                text = " ".join(b.text for b in response.content if b.type == "text").strip()
                if text:
                    trace.add("note", text=f"agent summary: {text}")
                break

            results = []
            for tu in tool_uses:
                if box.done:
                    out = {"error": "run is already resolved; stop calling tools"}
                else:
                    out = box.call(tu.name, dict(tu.input))
                results.append(
                    {
                        "type": "tool_result",
                        "tool_use_id": tu.id,
                        "content": json.dumps(out, default=str),
                        "is_error": "error" in out,
                    }
                )
            messages.append({"role": "user", "content": results})
        return _finish(box, trace, gave_up_reason=f"no resolution within {self.max_turns} turns")


# ------------------------------------------------------------------ rules


class RulePlanner:
    """Hand-written triage. It encodes the same playbook the LLM is told to follow."""

    name = "rules"

    def recover(self, box: Toolbox, trace: Trace) -> AgentOutcome:
        for _ in range(box.max_relaunches + 1):
            if box.done:
                break
            status = box.call("get_run_status", {"rationale": "start from the failed check"})
            check = status["last_run"]["failed_checks"][0]
            key = f"{check['step']}.{check['name']}"
            fix = self._diagnose(box, key, check)
            if fix is None:
                break
            (k, v), from_step = fix
            box.call("set_config", {"rationale": f"root cause of {key}", "key": k, "value": v})
            box.call("relaunch", {"rationale": f"verify the {k} fix", "from_step": from_step})
        if not box.done:
            box.call("escalate", {"rationale": "out of playbook", "reason": "unresolved after playbook", "evidence": ""})
        return _finish(box, trace)

    def _diagnose(self, box: Toolbox, key: str, check: dict):
        """Return ((key, value), relaunch_step) for a fix, or None after escalating."""

        def esc(reason: str, evidence: dict) -> None:
            box.call("escalate", {"rationale": reason, "reason": reason, "evidence": json.dumps(evidence)})

        if key == "load.finite_nonnegative":
            m = box.call("inspect_matrix", {"rationale": "check how widespread the invalid values are"})
            return esc("matrix contains invalid values; upstream correction is broken", {k: m[k] for k in ("n_nan", "n_negative")})

        if key == "load.orientation":
            m = box.call("inspect_matrix", {"rationale": "observations are not barcodes; maybe the matrix is transposed"})
            if m["var_barcode_like_fraction"] > 0.5 > m["obs_barcode_like_fraction"]:
                return ("orientation", "genes_x_cells"), "load"
            return esc("cannot identify the cell axis", m)

        if key == "load.value_scale":
            m = box.call("inspect_matrix", {"rationale": "non-integer values: normalized input or corruption?"})
            if m["n_negative"] == 0 and (m["max"] or 0) < 20 and m.get("expm1_per_cell_total_cv", 1) < 0.05:
                return ("input_scale", "log1p"), "load"
            return esc("non-integer values that do not look like log-normalized data", m)

        if key in ("qc.mito_genes_detected", "annotate.marker_overlap"):
            g = box.call("inspect_gene_names", {"rationale": "maybe the gene nomenclature does not match the configured species"})
            inferred = "human" if g["fraction_all_uppercase"] > 0.8 and g["n_prefixed_MT-"] > 0 else (
                "mouse" if g["fraction_titlecase"] > 0.6 and g["n_prefixed_mt-"] > 0 else None
            )
            if inferred and inferred != box.pipeline.cfg.species:
                return ("species", inferred), "qc"
            return esc("gene naming does not explain the failure", g)

        if key == "qc.cells_retained":
            q = box.call("inspect_qc_distribution", {"rationale": "most cells fail QC: bad config or bad library?"})
            return esc("library too shallow: most cells are below the QC floor, and thresholds are not agent-editable", q)

        return esc(f"no playbook entry for {key}", check)


# ------------------------------------------------------------------ driver


def run_scenario(name: str, clean_path: Path, workdir: Path, planner, trace_dir: Path | None = None) -> tuple[AgentOutcome, Pipeline, Trace]:
    scenario = SCENARIOS[name]
    path = materialize(scenario, clean_path, workdir)
    pipe = Pipeline(path, scenario.config)
    trace = Trace(scenario=name, planner=planner.name)
    trace.add("note", text=scenario.description)

    first = pipe.run()
    trace.add("pipeline_run", from_step="load", **first.to_dict())
    if first.status == "passed":
        outcome = AgentOutcome("clean", [], 0, 0, "no failure; agent not invoked")
        trace.add("outcome", outcome="clean", detail=outcome.detail)
    else:
        box = Toolbox(pipe, trace)
        outcome = planner.recover(box, trace)

    agreement = pipe.reference_agreement()
    if agreement is not None:
        trace.add("note", text=f"agreement with published labels: {agreement:.1%}")
    if trace_dir:
        trace.write(trace_dir)
    return outcome, pipe, trace


def grade(name: str, outcome: AgentOutcome, cfg: RunConfig) -> bool:
    s = SCENARIOS[name]
    if s.expected != outcome.outcome:
        return False
    return all(getattr(cfg, k) == v for k, v in s.expected_fix.items())
