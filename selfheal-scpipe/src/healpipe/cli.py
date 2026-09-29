"""healpipe CLI.

  healpipe scenarios                          list fault scenarios
  healpipe run --scenario prenormalized_input [--planner claude]
  healpipe eval [--planner claude]            run every scenario, print a scorecard
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from .agent import ClaudePlanner, RulePlanner, grade, run_scenario
from .faults import SCENARIOS


def _planner(args):
    if args.planner == "claude":
        try:
            return ClaudePlanner(model=args.model, effort=args.effort)
        except TypeError as e:  # raised by the SDK when no credentials resolve
            sys.exit(f"Claude planner needs Anthropic credentials (e.g. export ANTHROPIC_API_KEY=...): {e}")
    return RulePlanner()


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(prog="healpipe")
    sub = p.add_subparsers(dest="cmd", required=True)
    sub.add_parser("scenarios")
    for name in ("run", "eval"):
        sp = sub.add_parser(name)
        sp.add_argument("--data", default="data/pbmc3k_counts.h5ad")
        sp.add_argument("--planner", choices=["rules", "claude"], default="rules")
        sp.add_argument("--model", default="claude-opus-5-5")
        sp.add_argument("--effort", default="medium", choices=["low", "medium", "high", "xhigh", "max"])
        sp.add_argument("--out", default="runs")
        if name == "run":
            sp.add_argument("--scenario", required=True, choices=sorted(SCENARIOS))
    args = p.parse_args(argv)

    if args.cmd == "scenarios":
        for s in SCENARIOS.values():
            print(f"{s.name:<22} expect={s.expected:<9} {s.description}")
        return 0

    data = Path(args.data)
    if not data.exists():
        print(f"{data} not found. Run: python scripts/fetch_data.py", file=sys.stderr)
        return 2
    out = Path(args.out)
    planner = _planner(args)
    names = [args.scenario] if args.cmd == "run" else list(SCENARIOS)

    rows = []
    for name in names:
        outcome, pipe, trace = run_scenario(name, data, out / "inputs", planner, trace_dir=out / "traces")
        ok = grade(name, outcome, pipe.cfg)
        agree = pipe.reference_agreement()
        rows.append((name, SCENARIOS[name].expected, outcome, ok, agree))
        if args.cmd == "run":
            print(trace.to_markdown())

    lines = [
        f"| scenario | expected | outcome | correct | relaunches | tool calls | label agreement |",
        f"|---|---|---|---|---|---|---|",
    ]
    for name, exp, o, ok, agree in rows:
        a = f"{agree:.0%}" if agree is not None else "-"
        lines.append(f"| {name} | {exp} | {o.outcome} | {'yes' if ok else '**no**'} | {o.relaunches} | {o.tool_calls} | {a} |")
    n_ok = sum(r[3] for r in rows)
    lines.append(f"\n**{n_ok}/{len(rows)} correct** (planner: {planner.name})")
    table = "\n".join(lines)
    print(table)
    if args.cmd == "eval":
        (out / f"eval.{planner.name}.md").write_text(table + "\n")
    return 0 if n_ok == len(rows) else 1


if __name__ == "__main__":
    sys.exit(main())
