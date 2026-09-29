"""Exercise ClaudePlanner's tool loop against a scripted stand-in for the API client."""

from types import SimpleNamespace as NS

from healpipe.agent import ClaudePlanner, run_scenario


def _tool_use(i, name, **inp):
    return NS(type="tool_use", id=f"tu_{i}", name=name, input={"rationale": "test", **inp})


class ScriptedClient:
    def __init__(self, turns):
        self.turns, self.requests = list(turns), []
        self.beta = NS(messages=NS(create=self._create))

    def _create(self, **kw):
        self.requests.append({**kw, "messages": list(kw["messages"])})  # snapshot: the planner keeps appending
        content = self.turns.pop(0)
        stop = "tool_use" if any(b.type == "tool_use" for b in content) else "end_turn"
        return NS(content=content, stop_reason=stop)


def _planner(turns):
    p = ClaudePlanner.__new__(ClaudePlanner)
    p.client, p.model, p.effort, p.max_turns = ScriptedClient(turns), "claude-opus-5-5", "medium", 12
    return p


def test_claude_loop_fixes_species(clean_h5ad, tmp_path):
    planner = _planner(
        [
            [_tool_use(1, "inspect_gene_names")],
            [_tool_use(2, "set_config", key="species", value="human")],
            [_tool_use(3, "relaunch", from_step="qc")],
            [NS(type="text", text="Species was mislabeled; fixed and relaunched.")],
        ]
    )
    outcome, pipe, trace = run_scenario("species_mislabel", clean_h5ad, tmp_path, planner)
    assert outcome.outcome == "fixed" and pipe.cfg.species == "human"

    reqs = planner.client.requests
    assert all(t["strict"] for t in reqs[0]["tools"])
    # Each tool result answers the preceding tool_use id, in a user turn.
    last = reqs[-1]["messages"][-1]
    assert last["role"] == "user" and last["content"][0]["tool_use_id"] == "tu_3"
    assert any(e["kind"] == "note" and "agent summary" in e["text"] for e in trace.events)


def test_claude_loop_blocks_threshold_edit_then_escalates(clean_h5ad, tmp_path):
    planner = _planner(
        [
            [_tool_use(1, "set_config", key="min_genes", value="10")],
            [_tool_use(2, "escalate", reason="library too shallow", evidence="median 73 genes/cell")],
            [NS(type="text", text="Escalated.")],
        ]
    )
    outcome, _, trace = run_scenario("shallow_sequencing", clean_h5ad, tmp_path, planner)
    assert outcome.outcome == "escalated"
    first_result = planner.client.requests[1]["messages"][-1]["content"][0]
    assert first_result["is_error"] and "not agent-editable" in first_result["content"]
