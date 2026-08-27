#!/usr/bin/env python3
"""
compute_framing_a.py

Merges every per-read call table against the truth tables on read_id, then
computes recall, precision, and false positive rate per method per
condition, faceted by fold_position (primary) and adapter_present
(secondary), and plots the Framing A figure.

Directory layout expected:
    ~/foldback/sim_data/truth_<condition>.tsv        (9 spiked conditions)
    ~/foldback/sim_data/truth_clean_control.tsv       (optional -- FPR is
                                                        skipped with a
                                                        warning if missing)
    ~/foldback/results/**/calls_<method>_<condition>.tsv
    ~/foldback/results/**/calls_<method>_clean_control.tsv   (optional)

Call tables are found by exact filename match against the known METHODS x
CONDITIONS list (not by splitting on "_", since both method names and
condition names contain underscores).

Run:
    python compute_framing_a.py
"""

import glob
import os

import matplotlib.pyplot as plt
import pandas as pd

TRUTH_DIR = os.path.expanduser("~/foldback/sim_data")
RESULTS_DIR = os.path.expanduser("~/foldback/results")
OUT_DIR = os.path.expanduser("~/foldback/figures")

CONDITIONS = [
    "1pct_quarter", "1pct_middle", "1pct_near_end",
    "5pct_quarter", "5pct_middle", "5pct_near_end",
    "10pct_quarter", "10pct_middle", "10pct_near_end",
]
# Add methods here as they become available (e.g. duplex_tools,
# read_level_detector, alignment_level_detector).
METHODS = ["yacrd", "breakinator", "sniffles2_with_filter", "sniffles2_no_filter"]

SCORE_THRESHOLD = 0.5  # for methods whose `flagged` column is a float score

FOLD_ORDER = ["quarter", "middle", "near_end"]
FOLD_LABELS = {"quarter": "Quarter", "middle": "Middle", "near_end": "Near end"}


# ---------------------------------------------------------------------------
# File finding (exact-match against known method/condition names -- avoids
# ambiguous underscore-splitting since both contain underscores)
# ---------------------------------------------------------------------------

def find_call_file(method, condition):
    pattern = os.path.join(RESULTS_DIR, "**", f"calls_{method}_{condition}.tsv")
    matches = glob.glob(pattern, recursive=True)
    return matches[0] if matches else None


def find_clean_call_file(method):
    pattern = os.path.join(RESULTS_DIR, "**", f"calls_{method}_clean_control.tsv")
    matches = glob.glob(pattern, recursive=True)
    return matches[0] if matches else None


# ---------------------------------------------------------------------------
# Loading
# ---------------------------------------------------------------------------

def load_truth_tables():
    frames = []
    for condition in CONDITIONS:
        path = os.path.join(TRUTH_DIR, f"truth_{condition}.tsv")
        if not os.path.exists(path):
            print(f"WARNING: missing {path}, skipping this condition.")
            continue
        df = pd.read_csv(path, sep="\t")
        df["condition"] = condition
        df["is_foldback"] = df["is_foldback"].astype(str).str.strip() == "True"
        frames.append(df)
    return pd.concat(frames, ignore_index=True)


def _to_bool(value):
    if isinstance(value, bool):
        return value
    s = str(value).strip()
    if s in ("True", "False"):
        return s == "True"
    return float(s) >= SCORE_THRESHOLD


def load_call_table(path):
    df = pd.read_csv(path, sep="\t")
    df["flagged"] = df["flagged"].apply(_to_bool)
    return df[["read_id", "flagged"]]


# ---------------------------------------------------------------------------
# Metrics
# ---------------------------------------------------------------------------

def compute_recall_table(truth, methods):
    """Recall per (method, fold_position, adapter_present), pooled across
    the three foldback_fraction conditions -- recall is a per-read property
    and fraction isn't a facet in the spec."""
    foldback_truth = truth[truth["is_foldback"]].copy()

    rows = []
    for method in methods:
        calls = []
        for condition in CONDITIONS:
            path = find_call_file(method, condition)
            if path is None:
                continue
            c = load_call_table(path)
            c["condition"] = condition
            calls.append(c)
        if not calls:
            print(f"WARNING: no call tables found for method '{method}', skipping.")
            continue
        calls = pd.concat(calls, ignore_index=True)

        merged = foldback_truth.merge(calls, on=["read_id", "condition"], how="left")
        merged["flagged"] = merged["flagged"].fillna(False)

        for (fold_position, adapter_present), g in merged.groupby(
            ["fold_position", "adapter_present"]
        ):
            n = len(g)
            tp = g["flagged"].sum()
            rows.append({
                "method": method,
                "fold_position": fold_position,
                "adapter_present": adapter_present,
                "n_foldback_reads": n,
                "recall": tp / n if n else float("nan"),
            })
    return pd.DataFrame(rows)


def compute_precision_table(truth, methods):
    """Precision per method, pooled across all spiked conditions (uses each
    condition's own mix of foldback + clean reads)."""
    rows = []
    for method in methods:
        tp = fp = 0
        for condition in CONDITIONS:
            path = find_call_file(method, condition)
            if path is None:
                continue
            c = load_call_table(path)
            t = truth[truth["condition"] == condition]
            merged = t.merge(c, on="read_id", how="left")
            merged["flagged"] = merged["flagged"].fillna(False)
            flagged_rows = merged[merged["flagged"]]
            tp += flagged_rows["is_foldback"].sum()
            fp += (~flagged_rows["is_foldback"]).sum()
        precision = tp / (tp + fp) if (tp + fp) else float("nan")
        rows.append({"method": method, "tp": tp, "fp": fp, "precision": precision})
    return pd.DataFrame(rows)


def compute_fpr_table(methods):
    """FPR per method, measured on the dedicated clean control. Skips (with
    a warning) any method whose clean-control call table doesn't exist yet."""
    clean_truth_path = os.path.join(TRUTH_DIR, "truth_clean_control.tsv")
    if not os.path.exists(clean_truth_path):
        print("WARNING: truth_clean_control.tsv not found -- FPR will be NaN "
              "for all methods until the clean control run is available.")
        return pd.DataFrame([{"method": m, "fpr": float("nan")} for m in methods])

    clean_truth = pd.read_csv(clean_truth_path, sep="\t")
    n_clean = len(clean_truth)

    rows = []
    for method in methods:
        path = find_clean_call_file(method)
        if path is None:
            print(f"WARNING: no clean-control call table found for '{method}' yet -- FPR set to NaN.")
            rows.append({"method": method, "fpr": float("nan")})
            continue
        c = load_call_table(path)
        merged = clean_truth.merge(c, on="read_id", how="left")
        merged["flagged"] = merged["flagged"].fillna(False)
        fp = merged["flagged"].sum()
        rows.append({"method": method, "n_clean_reads": n_clean, "fpr": fp / n_clean})
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Plot: Framing A figure
# ---------------------------------------------------------------------------

def plot_framing_a(recall_table, fpr_table, methods, out_path):
    fpr_by_method = dict(zip(fpr_table["method"], fpr_table["fpr"]))

    fig, axes = plt.subplots(len(FOLD_ORDER), 1, figsize=(9, 4.5 * len(FOLD_ORDER)), sharey=True)
    bar_width = 0.35
    x = range(len(methods))

    for ax, fold_position in zip(axes, FOLD_ORDER):
        sub = recall_table[recall_table["fold_position"] == fold_position]

        present = [
            sub.query("method == @m and adapter_present == True")["recall"].mean() * 100
            if not sub.query("method == @m and adapter_present == True").empty else 0
            for m in methods
        ]
        absent = [
            sub.query("method == @m and adapter_present == False")["recall"].mean() * 100
            if not sub.query("method == @m and adapter_present == False").empty else 0
            for m in methods
        ]

        ax.bar([i - bar_width / 2 for i in x], present, bar_width, label="Adapter present", color="#2a78d6")
        ax.bar([i + bar_width / 2 for i in x], absent, bar_width, label="Adapter absent", color="#eb6834")

        for i, m in enumerate(methods):
            fpr = fpr_by_method.get(m, float("nan"))
            label = f"FPR {fpr * 100:.1f}%" if pd.notna(fpr) else "FPR n/a"
            bar_top = max(present[i], absent[i])
            ax.text(i, bar_top + 4, label, ha="center", fontsize=12, color="dimgray")

        ax.set_title(FOLD_LABELS[fold_position], loc="left", fontweight="bold", fontsize=16)
        ax.set_ylim(0, 112)
        ax.set_ylabel("Recall (%)", fontsize=14)
        ax.set_xticks(list(x))
        ax.set_xticklabels(methods, rotation=0, ha="center", fontsize=13)
        ax.tick_params(axis="y", labelsize=12)
        ax.spines[["top", "right"]].set_visible(False)

    axes[0].legend(loc="upper left", bbox_to_anchor=(0, 1.3), ncol=2, frameon=False, fontsize=13)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    print(f"Wrote {out_path}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    truth = load_truth_tables()

    recall_table = compute_recall_table(truth, METHODS)
    precision_table = compute_precision_table(truth, METHODS)
    fpr_table = compute_fpr_table(METHODS)

    recall_table.to_csv(os.path.join(OUT_DIR, "recall_table.tsv"), sep="\t", index=False)
    precision_table.to_csv(os.path.join(OUT_DIR, "precision_table.tsv"), sep="\t", index=False)
    fpr_table.to_csv(os.path.join(OUT_DIR, "fpr_table.tsv"), sep="\t", index=False)

    print("\n=== Recall (by method, fold_position, adapter_present) ===")
    print(recall_table.to_string(index=False))
    print("\n=== Precision (by method) ===")
    print(precision_table.to_string(index=False))
    print("\n=== FPR (by method) ===")
    print(fpr_table.to_string(index=False))

    plot_framing_a(recall_table, fpr_table, METHODS,
                    out_path=os.path.join(OUT_DIR, "framing_a.png"))


if __name__ == "__main__":
    main()
