| scenario | expected | outcome | correct | relaunches | tool calls | label agreement |
|---|---|---|---|---|---|---|
| clean | clean | clean | yes | 0 | 0 | 79% |
| prenormalized_input | fixed | fixed | yes | 1 | 4 | 79% |
| transposed_matrix | fixed | fixed | yes | 1 | 4 | 79% |
| species_mislabel | fixed | fixed | yes | 1 | 4 | 79% |
| shallow_sequencing | escalated | escalated | yes | 0 | 3 | - |
| negative_values | escalated | escalated | yes | 0 | 3 | - |

**6/6 correct** (planner: rules)
