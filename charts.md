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

    classDef detailA fill:#DCEEFB,stroke:#4A90D9,stroke-width:2px,color:#1a1a1a;
    classDef detailC fill:#FDECD8,stroke:#E8963C,stroke-width:2px,color:#1a1a1a;
    class B2 detailA;
    class RF detailC;
```

Detail diagrams: [A]-detail, [C]-detail (below).

---

# Simulation of foldback reads

```mermaid
flowchart TD
    G[HG002 assembly] --> P[PBSIM3<br/>simulated ONT reads]
    P --> CL["Clean simulated reads (50x)"]

    CL --> RS[Reservoir sampling] --> CB["Clean-read baseline (20x)"]
    CL --> FC["Foldback conversion<br/>S → S·revcomp(S)"]

    FC --> FRC{Foldback fraction}
    FRC --> F1["1%"]
    FRC --> F5["5%"]
    FRC --> F10["10%"]

    FC --> FPC{Fold position}
    FPC --> MID["Middle<br/>(45–55% of read)"]
    FPC --> OFF["Off-center<br/>(70–80%)"]
    FPC --> NEAR["Near-end<br/>(last 50–500 bp)"]

    F1 --> SWEEP["Sweep: 3 fractions × 3 positions<br/>= 9 conditions"]
    F5 --> SWEEP
    F10 --> SWEEP
    MID --> SWEEP
    OFF --> SWEEP
    NEAR --> SWEEP

    SWEEP --> ADAPT["± adapter at fold junction, p=0.2<br/>(stratified within each condition)"]
    ADAPT --> OUT[Simulated reads +<br/>ground truth labels]

    classDef detailA fill:#DCEEFB,stroke:#4A90D9,stroke-width:2px,color:#1a1a1a;
    class G,P,CL,RS,CB,FC,FRC,F1,F5,F10,FPC,MID,OFF,NEAR,SWEEP,ADAPT,OUT detailA;
```

---

# foldback-hunter (ours)

Reference-free: detects a read's self-complementarity directly, no alignment to genome.


**Probe mode** is primary (ran on full dataset); 

**full mode** is a partial spot-check only.

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

    classDef detailC fill:#FDECD8,stroke:#E8963C,stroke-width:2px,color:#1a1a1a;
    class START,Q1,R1,LOOP,CALC,KEEP,SCORE,Q3,R3,R2 detailC;
```

## Full mode

```mermaid
flowchart TD
    A[Read] --> B[Reverse complement]
    B --> C[Align read vs.<br/>its reverse complement]
    C --> D{High similarity?}
    D -->|yes| E([Foldback call])

    classDef detailC fill:#FDECD8,stroke:#E8963C,stroke-width:2px,color:#1a1a1a;
    class A,B,C,D,E detailC;
```
