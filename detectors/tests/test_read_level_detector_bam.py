"""
Unit tests for BAM input support in detectors/read_level_detector.py:
iter_bam's primary-only filtering, reverse-strand revcomp correctness,
FASTQ/BAM output equivalence, and an end-to-end smoke test on an
unaligned BAM.
"""

import random
import sys

import pysam
import pytest
import read_level_detector as rld
from foldback_score import revcomp
from test_read_level_detector import make_foldback_read


def make_aligned_bam(path):
    """One @SQ header, 3 records: forward primary, reverse primary,
    and a supplementary sharing the reverse record's QNAME."""
    header = pysam.AlignmentHeader.from_dict({
        "HD": {"VN": "1.6", "SO": "unsorted"},
        "SQ": [{"SN": "chr1", "LN": 10000}],
    })

    fwd_seq = "".join(random.Random(1).choice("ACGT") for _ in range(500))
    rev_seq = "".join(random.Random(2).choice("ACGT") for _ in range(500))

    with pysam.AlignmentFile(path, "wb", header=header) as af:
        fwd = pysam.AlignedSegment(header)
        fwd.query_name = "fwd_read"
        fwd.query_sequence = fwd_seq
        fwd.flag = 0
        fwd.reference_id = 0
        fwd.reference_start = 0
        fwd.mapping_quality = 60
        fwd.cigarstring = f"{len(fwd_seq)}M"
        af.write(fwd)

        rev = pysam.AlignedSegment(header)
        rev.query_name = "rev_read"
        rev.query_sequence = rev_seq
        rev.flag = 16  # reverse strand, primary
        rev.reference_id = 0
        rev.reference_start = 100
        rev.mapping_quality = 60
        rev.cigarstring = f"{len(rev_seq)}M"
        af.write(rev)

        supp = pysam.AlignedSegment(header)
        supp.query_name = "rev_read"
        supp.query_sequence = rev_seq[:200]
        supp.flag = 16 | 2048  # reverse strand, supplementary
        supp.reference_id = 0
        supp.reference_start = 5000
        supp.mapping_quality = 60
        supp.cigarstring = f"{len(rev_seq[:200])}M"
        af.write(supp)

    return rev_seq


def make_unaligned_bam(path):
    """A single uBAM-style unmapped record."""
    header = pysam.AlignmentHeader.from_dict({"HD": {"VN": "1.6"}})
    seq = "".join(random.Random(3).choice("ACGT") for _ in range(2000))

    with pysam.AlignmentFile(path, "wb", header=header) as af:
        rec = pysam.AlignedSegment(header)
        rec.query_name = "ubam_read"
        rec.query_sequence = seq
        rec.flag = 4  # unmapped
        af.write(rec)

    return seq


def test_iter_bam_skips_supplementary(tmp_path):
    bam_path = tmp_path / "aligned.bam"
    make_aligned_bam(bam_path)

    read_ids = [read_id for read_id, _seq in rld.iter_bam(str(bam_path))]

    assert read_ids == ["fwd_read", "rev_read"]


def test_iter_bam_revcomps_reverse_strand(tmp_path):
    bam_path = tmp_path / "aligned.bam"
    rev_seq_in_bam = make_aligned_bam(bam_path)

    reads = dict(rld.iter_bam(str(bam_path)))

    assert reads["rev_read"] == revcomp(rev_seq_in_bam).upper()


def test_fastq_and_bam_inputs_produce_identical_output(tmp_path, monkeypatch):
    """Same reads, same original orientation, fed as FASTQ vs. as an
    aligned BAM (one record stored reverse-complemented with flag=16,
    as a real aligner would) must score identically end-to-end."""
    reads = [
        ("r1", make_foldback_read(2000, 1000, seed=1)),
        ("r2", make_foldback_read(2000, 1500, seed=2)),
        ("r3", "".join(random.Random(3).choice("ACGT") for _ in range(2000))),
    ]

    fastq_path = tmp_path / "in.fastq"
    with open(fastq_path, "w") as fh:
        for read_id, seq in reads:
            fh.write(f"@{read_id}\n{seq}\n+\n{'I' * len(seq)}\n")

    header = pysam.AlignmentHeader.from_dict({
        "HD": {"VN": "1.6", "SO": "unsorted"},
        "SQ": [{"SN": "chr1", "LN": 10000}],
    })
    bam_path = tmp_path / "in.bam"
    with pysam.AlignmentFile(bam_path, "wb", header=header) as af:
        for i, (read_id, seq) in enumerate(reads):
            rec = pysam.AlignedSegment(header)
            rec.query_name = read_id
            if i == 1:
                rec.query_sequence = revcomp(seq)
                rec.flag = 16
            else:
                rec.query_sequence = seq
                rec.flag = 0
            rec.reference_id = 0
            rec.reference_start = i * 100
            rec.mapping_quality = 60
            rec.cigarstring = f"{len(seq)}M"
            af.write(rec)

    fastq_outdir = tmp_path / "out_fastq"
    monkeypatch.setattr(sys, "argv", [
        "read_level_detector.py", "--input", str(fastq_path),
        "--outdir", str(fastq_outdir), "--mode", "probe",
    ])
    rld.main()

    bam_outdir = tmp_path / "out_bam"
    monkeypatch.setattr(sys, "argv", [
        "read_level_detector.py", "--input", str(bam_path),
        "--outdir", str(bam_outdir), "--mode", "probe",
    ])
    rld.main()

    fastq_calls = (fastq_outdir / "calls_foldback_hunter_in.tsv").read_text()
    bam_calls = (bam_outdir / "calls_foldback_hunter_in.tsv").read_text()
    fastq_scores = (fastq_outdir / "scores_foldback_hunter_in.tsv").read_text()
    bam_scores = (bam_outdir / "scores_foldback_hunter_in.tsv").read_text()

    assert fastq_calls == bam_calls
    assert fastq_scores == bam_scores


def test_main_smoke_test_on_unaligned_bam(tmp_path, monkeypatch):
    bam_path = tmp_path / "in.bam"
    make_unaligned_bam(bam_path)
    outdir = tmp_path / "out"

    argv = [
        "read_level_detector.py",
        "--input", str(bam_path),
        "--outdir", str(outdir),
        "--mode", "probe",
    ]
    monkeypatch.setattr(sys, "argv", argv)
    rld.main()

    assert (outdir / "calls_foldback_hunter_in.tsv").exists()
    assert (outdir / "scores_foldback_hunter_in.tsv").exists()
    assert (outdir / "summary_foldback_hunter_in.json").exists()
