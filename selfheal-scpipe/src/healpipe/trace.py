"""Structured run log: symptom -> hypothesis -> action -> outcome."""

from __future__ import annotations

import json
import time
from dataclasses import dataclass, field
from pathlib import Path


@dataclass
class Trace:
    scenario: str
    planner: str
    events: list[dict] = field(default_factory=list)
    _t0: float = field(default_factory=time.monotonic)

    def add(self, kind: str, **data) -> None:
        self.events.append({"t": round(time.monotonic() - self._t0, 2), "kind": kind, **data})

    def write(self, out_dir: Path) -> tuple[Path, Path]:
        out_dir.mkdir(parents=True, exist_ok=True)
        stem = f"{self.scenario}.{self.planner}"
        j, m = out_dir / f"{stem}.trace.json", out_dir / f"{stem}.trace.md"
        j.write_text(json.dumps({"scenario": self.scenario, "planner": self.planner, "events": self.events}, indent=2, default=str))
        m.write_text(self.to_markdown())
        return j, m

    def to_markdown(self) -> str:
        lines = [f"# Run trace: `{self.scenario}` (planner: {self.planner})", ""]
        for e in self.events:
            k = e["kind"]
            if k == "pipeline_run":
                lines.append(f"**[{e['t']:>6}s] pipeline** from `{e['from_step']}` -> **{e['status'].upper()}**")
                for c in e.get("failed_checks", []):
                    lines.append(f"  - symptom: `{c['step']}.{c['name']}`: {c['message']}")
            elif k == "tool_call":
                args = {a: v for a, v in e["input"].items() if a != "rationale"}
                lines.append(f"**[{e['t']:>6}s] agent** -> `{e['tool']}({json.dumps(args)[1:-1]})`")
                if e["input"].get("rationale"):
                    lines.append(f"  - hypothesis: {e['input']['rationale']}")
            elif k == "tool_result":
                summary = e.get("summary") or json.dumps(e["result"], default=str)
                lines.append(f"  - result: {summary[:400]}")
            elif k == "outcome":
                lines += ["", f"**Outcome: {e['outcome'].upper()}**: {e.get('detail', '')}"]
            elif k == "note":
                lines.append(f"_{e['text']}_")
            lines.append("")
        return "\n".join(lines)
