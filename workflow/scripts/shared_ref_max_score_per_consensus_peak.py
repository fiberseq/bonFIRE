#!/usr/bin/env python
"""
Shared-reference mode: picks one pileup row per consensus peak for one run.

Reads, on stdin, pileup rows intersected with the consensus peaks (bedtools intersect -wa -wb): chrom,
start, end, coverage, fire_coverage, score, coverage_H1, fire_coverage_H1, score_H1, coverage_H2,
fire_coverage_H2, score_H2 (NA where the pileup has no haplotype columns), then the consensus peak's
chrom, start, end and id.

For each consensus peak the row kept is the one with the highest score, then the highest
fire_coverage / coverage, then the highest fire_coverage, then the highest coverage; on a full tie the
first row (the leftmost) stays. This is the same choice as in graph mode
(260217_find_max_fire_score_per_consenses_peak.py). The haplotype values come from that same row, so
both haplotypes are measured at the same position.

Consensus peaks that no pileup row overlaps have no reads in this run: coverage and fire_coverage 0 and
score -1 (FIRE's score where there are no FIRE reads), and the same for each haplotype if the pileup has
haplotype columns (NA if it has none).

usage: ... | shared_ref_max_score_per_consensus_peak.py consensus_peaks.bed sample_id output.tsv
"""
import sys

VALUE_COLS = ["coverage", "fire_coverage", "score",
              "coverage_H1", "fire_coverage_H1", "score_H1",
              "coverage_H2", "fire_coverage_H2", "score_H2"]


def rank(cov, fire_cov, score):
    # a position with no reads has no FIRE fraction; count it as 0 so it never beats one with reads
    fire_frac = fire_cov / cov if cov > 0 else 0.0
    return (score, fire_frac, fire_cov, cov)


def main(consensus_file, sample_id, output_file):
    best = {}
    # whether the pileup has haplotype columns (they're NA when it doesn't)
    has_haplotypes = False
    for line in sys.stdin:
        f = line.rstrip("\n").split("\t")
        values = f[3:12]
        if values[3] != "NA":
            has_haplotypes = True
        peak_id = f[-1]
        r = rank(int(values[0]), int(values[1]), float(values[2]))
        # strictly better only, so on a full tie the first (leftmost) row stays
        if peak_id not in best or r > best[peak_id][0]:
            best[peak_id] = (r, values)

    with open(output_file, "w") as out:
        out.write("\t".join(["sample_id", "consensus_peak_id"] + VALUE_COLS) + "\n")
        with open(consensus_file) as cons:
            for line in cons:
                if line.startswith("#"):
                    continue
                peak_id = line.rstrip("\n").split("\t")[-1]
                if peak_id in best:
                    values = best[peak_id][1]
                else:
                    values = ["0", "0", "-1"] + (["0", "0", "-1"] * 2 if has_haplotypes else ["NA"] * 6)
                out.write("\t".join([sample_id, peak_id] + values) + "\n")


if __name__ == "__main__":
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    main(*sys.argv[1:])
