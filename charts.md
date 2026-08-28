# Methods Overview

```mermaid
flowchart TD
    B1[HG002] --> B2[simulated reads +<br/>ground truth labels]

    B2 --> Bhead["BASELINE"]
    B2 --> Chead["DETECTORS"]

    subgraph Bblock [" "]
        Bhead --> yacrd[yacrd]
    end

    subgraph Cblock [" "]
        Chead --> RB["reference-based:<br/>breakinator"]
        Chead --> RF["reference-free:<br/>foldback-hunter (ours)"]

        RB --> RBdesc["aligned to reference (minimap2);<br/>→ misses fold near read's end + reference bias"]
        RBdesc --> RBcalls[breakinator calls]

        RF --> RFdesc["2 modes: full / probe<br/><br/>full: SW local align</br> (read, rev_compl(read)) <br/><br/>probe: search read for its own tail, reverse-complemented probes;<br/><br/>score >= 99% → foldback call"]
    end

    yacrd --> perread["per-read calls<br/>(per method, per condition)"]
    RBcalls --> perread
    RFdesc --> perread

    perread --> D["BENCHMARKING<br/>compare against ground truth, matched read by read<br/>→ recall, precision, FPR per method<br/>facet: fold position × % of folds"]

    D --> results["results<br/>precision & recall by method, broken down by<br/>fold position and adapter presence<br/>recall per method, FPR on the clean (non-foldback)<br/>control, annotated<br/><br/>expected shape:<br/>yacrd ~0<br/>foldback-hunter (ours): high across all fold positions<br/>breakinator: high near the middle of the read → low near the end"]
```

Detail diagrams: [A]-detail, [C]-detail (below).

---

# [A]-detail — Simulation of foldback reads

```mermaid
flowchart TD
    G["HG002 diploid assembly, chr6 only"] --> P["PBSIM3 (Ono et al. 2022)<br/>ONT high-quality long-read error model<br/>substitution:insertion:deletion ratio = 39:24:36 (ONT-recommended)<br/>coverage: 50x per haplotype<br/>read length: mean 15kb, sd 8kb"]
    P --> CL["clean simulated long reads<br/>(both haplotypes)"]

    CL --> RS["reservoir sampling, fixed seed<br/>~230,000 reads (~20x coverage)"]
    CL --> FC["foldback conversion, fraction f of reads<br/>new read = forward segment joined with the<br/>reverse complement of the segment<br/>immediately preceding the fold point<br/>(input read length preserved)"]

    RS --> CB["clean-read baseline<br/>(identical across every condition in the sweep)<br/>used as the non-foldback control"]

    FC --> FPC{"fold-position category<br/>(which fold points are eligible)"}
    FPC --> MID["middle:<br/>fold at 45-55% of read length"]
    FPC --> OFF["off-center:<br/>fold at 70-80% of read length"]
    FPC --> NEAR["near-end:<br/>fold within last 50-500bp"]

    MID --> ADAPT
    OFF --> ADAPT
    NEAR --> ADAPT["embed ONT adapter at fold junction, 20% of foldback reads<br/>→ adapter-bridged and sequence-only foldback variants,<br/>in physiological proportions"]

    ADAPT --> SWEEP["sweep = 3 foldback fractions (1%, 5%, 10%)<br/>× 3 fold positions (middle, off-center, near-end)<br/>= 9 conditions, adapter presence stratified within each"]

    SWEEP --> OUT["simulated reads + ground truth (per condition)<br/>each simulated read is labeled with: whether it's a<br/>foldback, where the fold occurs, whether an adapter is<br/>present, and the genomic region it came from"]
```

---

# [C]-detail — foldback-hunter (ours): probe mode and full mode

Reference-free, read-level detector: looks for a read's own self-complementarity
instead of aligning to the genome or examining supplementary alignments.
Two modes considered: probe and full.
Probe mode is the primary mode — fully run across the dataset, detailed below.
Full mode was not run on the full dataset (time constraints during the
hackathon); treat its results as a partial spot-check, not a like-for-like
comparison against probe mode's coverage.

## Probe mode (primary)

`foldback_score.py: score_read_probe`

```mermaid
flowchart TD
    START(["read_id, seq, probe_lengths=(50,150,500,1500), min_len=1000"])
    START --> Q1{"len(seq) < min_len?"}
    Q1 -->|yes| R1(["return: status='too_short'"])
    Q1 -->|no| INIT["seq_rc = revcomp(seq)<br/>best_score, best_k, best_hit = None, None, None"]

    INIT --> LOOP{"for k in probe_lengths<br/>(50 / 150 / 500 / 1500bp)"}
    LOOP -->|next k| CALC["k = min(k, len(seq))<br/>probe = seq_rc[0:k]<br/>(edit_distance, locations) = edlib.align(probe, seq, mode=HW)"]
    CALC --> Q2{"edit_distance < 0?"}
    Q2 -->|yes| LOOP
    Q2 -->|no| SCORE["score = 1 - edit_distance / k"]
    SCORE --> Q3{"score > best_score OR<br/>(score == best_score AND k > best_k)?"}
    Q3 -->|yes| UPDATE["best_score, best_k, best_hit<br/>= score, k, locations[0]"]
    Q3 -->|"no (ties favor larger k)"| LOOP
    UPDATE --> LOOP

    LOOP -->|loop done| Q4{"best_score is None or<br/>best_hit is None?"}
    Q4 -->|yes| R2(["return: score=0.0, status='no_match'"])
    Q4 -->|no| POS["hit_start, hit_end = best_hit<br/>fold_position_bp = (hit_start + len(seq)) // 2"]
    POS --> R3(["return: best_score, fold_position_bp, status='ok'"])
```

> Note: score >= 99% similarity is applied downstream (benchmarking), not
> inside this function — `raw_score` is written out as-is (schema.md:
> "flagged is the raw float score, not a threshold").

## Full mode

```mermaid
flowchart TD
    A[read sequence] --> B["reverse-complement the entire read"]
    B --> C["align whole read against its own reverse complement"]
    C --> D["score = percent identity of aligned region<br/>fold position = midpoint of the aligned region"]
    D --> E["same rule as probe mode:<br/>score >= 99% similarity → foldback call"]
    E --> F["partial dataset only — not the full 9-condition sweep —<br/>due to time constraints; use as a directional check<br/>on full mode's behavior, not a final benchmark"]
```
