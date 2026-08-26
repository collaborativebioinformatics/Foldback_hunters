# Table schemas

Frozen on day-1. Do not change columns without telling the group.

## Truth table

One row per simulated read. Written by the simulator (Task A) only. Read by everyone else.

| Column | Type | Meaning |
|---|---|---|
| `read_id` | string | The read's name in the FASTQ. Must match the `@` line exactly (no leading `@`, no trailing whitespace). |
| `is_foldback` | boolean | `True` if we spiked this read as a foldback, `False` for clean control reads. |
| `fold_position` | string | Where along the read the fold sits: `middle`, `off_center`, or `near_end`. Empty for clean reads (`is_foldback=False`). |
| `adapter_present` | boolean | `True` if we embedded an ONT adapter sequence at the fold junction, `False` otherwise. Empty for clean reads. |
| `source_locus` | string | Genomic region the read was simulated from, e.g. `chr6:12345678-12360000`. |

**File format:** TSV, one file per simulated condition, header row required. Filename convention: `truth_<condition>.tsv`.

## Per-read call table

One row per read per method. Each baseline and each detector writes its own file in this format.

| Column | Type | Meaning |
|---|---|---|
| `read_id` | string | Must match `read_id` in the truth table exactly. This is the join key. |
| `method` | string | Name of the tool or detector. Use one of: `breakinator`, `yacrd`, `duplex_tools`, `sniffles2_no_filter`, `sniffles2_with_filter`, `read_level_detector`. Add new methods to this list here before using them. |
| `flagged` | boolean or float | Did this method flag the read as a foldback? Either `True`/`False`, or a numeric score that can be thresholded downstream. Be consistent within a single file. |

**File format:** TSV, one file per method per condition, header row required. Filename convention: `calls_<method>_<condition>.tsv`.

## Notes

- All strings are case-sensitive.
- `read_id` values must be identical across the truth table and every call table; otherwise the join in the benchmarking step drops rows silently.
- If a method returns a score rather than a boolean, document the threshold used downstream in the benchmark script, not here.
- To add a column: edit this file first, commit, tell the group, then update code. Never the other way around.
