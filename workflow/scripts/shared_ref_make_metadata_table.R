## Shared-reference mode: builds the final tables from the per-run pileup values, the consensus peaks
## and the original peaks, all on the one reference.
##
## args: output folder (with trailing /), prefix, consensus peaks bed (chrom, start, end, id; header),
##       per-run pileup values (shared_ref_max_score_per_consensus_peak.py, all runs, with header),
##       original peaks intersected with the consensus peaks (bedtools -wa -wb: chrom, start, end, run,
##       peak_id, then the consensus peak's chrom, start, end, id), consensus peak sequences (bedtools
##       getfasta -bedOut -tab -name: chrom, start, end, id, seq), run names in input order (one per line),
##       threads
## out: [prefix]_final_tbl_w_metadata.tsv.gz (one row per run and consensus peak) and
##      [prefix]_consensus_peak_summary.tsv (one row per consensus peak)
suppressMessages(library(data.table))

args = commandArgs(trailingOnly = TRUE)
## the original-peak counting is shared with the graph-mode metadata script; it's in this script's folder
script_dir = dirname(sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)))
source(file.path(script_dir, "orig_peak_overlap_counts.R"))

output_dir = args[1]
prefix = args[2]
consensus_file = args[3]
values_file = args[4]
orig_intersect_file = args[5]
seq_file = args[6]
runs_file = args[7]
threads = if(length(args) >= 8) as.integer(args[8]) else 1L
setDTthreads(threads)

progress = function(msg){
  print(paste(format(Sys.time(), "%H:%M:%S"), msg))
}

progress("reading inputs")
cons = fread(consensus_file, select = 1:4, col.names = c("chr", "start", "end", "consensus_peak_id"))
runs = readLines(runs_file)
tbl = fread(values_file, colClasses = list(character = c("sample_id", "consensus_peak_id")))
ov = fread(orig_intersect_file, header = FALSE, select = 1:9,
           col.names = c("chr_orig", "start_orig", "end_orig", "run", "peak_id_orig", "chr_cons", "start_cons", "end_cons", "consensus_peak_id"))
seqs = fread(seq_file, header = FALSE, select = 4:5, col.names = c("consensus_peak_id", "seq"))

## every run has a row for every consensus peak
stopifnot(nrow(tbl) == length(runs) * nrow(cons), all(tbl$sample_id %in% runs))

progress("adding coordinates and accessibility")
m = match(tbl$consensus_peak_id, cons$consensus_peak_id)
tbl[, `:=`(cons_chr = cons$chr[m], cons_start = cons$start[m], cons_end = cons$end[m])]
tbl[, fire_cov_OV_cov := fire_coverage / coverage]
tbl[, fire_cov_OV_cov_H1 := fire_coverage_H1 / coverage_H1]
tbl[, fire_cov_OV_cov_H2 := fire_coverage_H2 / coverage_H2]

progress("adding original peak overlaps")
## each run's own original peaks: does one overlap the consensus peak, and how many separate ones
tbl[, key := (match(sample_id, runs) - 1) * nrow(cons) + m]
ov[, key := (match(run, runs) - 1) * nrow(cons) + match(consensus_peak_id, cons$consensus_peak_id)]
stopifnot(!anyNA(ov$key))
tbl[, overlaps_original_called_peak := key %in% ov$key]

ov[, group := key]
counts = count_orig_peaks_per_consensus_peak(as.data.frame(ov))
tbl[, num_overlapping_orig_peaks := 0L]
tbl[key %in% ov$key, num_overlapping_orig_peaks := 1L]
tbl[key %in% counts$pair_keys, num_overlapping_orig_peaks := 2L]
i = match(tbl$key, counts$gr2_keys)
tbl[!is.na(i), num_overlapping_orig_peaks := as.integer(counts$gr2_num_peaks[i[!is.na(i)]])]
tbl[, orig_peaks_overlapping_consensus := counts$summary$all_orig_overlapping_peaks[match(key, counts$summary$key)]]
rm(ov, counts, i)

progress("writing the table")
## runs in input order, consensus peaks along each chromosome
setorder(tbl, key)
tbl[, key := NULL]
setcolorder(tbl, c("sample_id", "consensus_peak_id", "coverage", "fire_coverage", "score", "cons_chr", "cons_start", "cons_end",
                   "fire_cov_OV_cov",
                   "coverage_H1", "fire_coverage_H1", "score_H1", "fire_cov_OV_cov_H1",
                   "coverage_H2", "fire_coverage_H2", "score_H2", "fire_cov_OV_cov_H2",
                   "overlaps_original_called_peak", "num_overlapping_orig_peaks", "orig_peaks_overlapping_consensus"))
fwrite(tbl, paste0(output_dir, prefix, "_final_tbl_w_metadata.tsv.gz"), sep = "\t", quote = FALSE, na = "NA")

progress("summarizing each consensus peak")
## chr_order numbers the consensus peaks 1..N along each chromosome
summary = copy(cons)
setorder(summary, chr, start, end)
summary[, chr_order := seq_len(.N), by = chr]
summary[, seq := seqs$seq[match(consensus_peak_id, seqs$consensus_peak_id)]]
setcolorder(summary, c("consensus_peak_id", "chr", "start", "end", "chr_order", "seq"))
fwrite(summary, paste0(output_dir, prefix, "_consensus_peak_summary.tsv"), sep = "\t", quote = FALSE, na = "NA")
progress("done")
