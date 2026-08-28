# Methods Overview

```mermaid
flowchart TD
    B1[HG002] --> B2[Simulated reads +<br/>ground truth]

    B2 --> BASE[yacrd<br/><i>baseline</i>]
    B2 --> RB["breakinator<br/><i>reference-based</i>"]
    B2 --> RF["foldback-hunter (ours)<br/><i>reference-free</i>"]

    BASE --> CALLS[Per-read calls]
    RB --> CALLS
    RF --> CALLS

    CALLS --> D["Benchmarking<br/>vs. ground truth →<br/>recall, precision, FPR"]
```

Detail diagrams: [A]-detail, [C]-detail (below).

---

# [A]-detail — Simulation of foldback reads

```mermaid
flowchart TD
    G[HG002 assembly] --> P[PBSIM3<br/>simulated ONT reads]
    P --> CL[Clean simulated reads]

    CL --> RS[Reservoir sampling] --> CB[Clean-read baseline]
    CL --> FC[Foldback conversion]

    FC --> FPC{Fold position}
    FPC --> MID[Middle]
    FPC --> OFF[Off-center]
    FPC --> NEAR[Near-end]

    MID --> ADAPT[± adapter at<br/>fold junction]
    OFF --> ADAPT
    NEAR --> ADAPT

    ADAPT --> SWEEP["Sweep: 3 fractions × 3 positions<br/>= 9 conditions"]
    SWEEP --> OUT[Simulated reads +<br/>ground truth labels]
```

---

# [C]-detail — foldback-hunter (ours)

Reference-free: detects a read's self-complementarity directly, no alignment to genome.
**Probe mode** is primary (ran on full dataset); **full mode** is a partial spot-check only.

## Probe mode

```mermaid
flowchart TD
    START([Read]) --> Q1{Too short?}
    Q1 -->|yes| R1([Skip])
    Q1 -->|no| LOOP["Try increasing probe lengths<br/>(short → long)"]
    LOOP --> CALC[Align rev-comp probe<br/>of that length against read]
    CALC --> KEEP[Keep best-scoring length so far]
    KEEP -->|next length| LOOP
    KEEP -->|done| SCORE[Best score across all lengths]
    SCORE --> Q3{Score high enough?}
    Q3 -->|yes| R3([Foldback call +<br/>fold position])
    Q3 -->|no| R2([No match])
```

## Full mode

```mermaid
flowchart TD
    A[Read] --> B[Reverse complement]
    B --> C[Align read vs.<br/>its reverse complement]
    C --> D{High similarity?}
    D -->|yes| E([Foldback call])
```
