## Checks the placement of non-reference consensus peaks on CHM13 against every assembly the peak
## transferred to. The placement itself comes from the assembly the peak was called on
## (261030_determine_previous_peak_for_nonRef_chrs.R); this script only reports whether the other
## assemblies agree, and whether the peak looks like a copy of a CHM13 peak.
##
## In each assembly, the CHM13 peaks next to a non-reference peak on its contig say where it lies on
## CHM13, as an interval of CHM13 peak ranks (lo, hi): the peak is after CHM13 peak lo and before hi.
##   - CHM13 peaks on both sides, same chr: lo and hi are their ranks (either order, so the contig's
##     direction doesn't matter)
##   - contig end, or the two sides are on different chrs: only the side on the contig's most common
##     chr counts, read with the contig's direction, which leaves one end open (-Inf or Inf)
##   - no CHM13 peak on the contig: no information
## Ranks come from the full CHM13 consensus peak list.
##
## Another assembly agrees with the called one if their intervals overlap. Intersecting the intervals
## of every assembly gives the tightest placement they all allow (all_asms_* columns), when they agree.
##
## A peak is a copy of a CHM13 peak in an assembly where the two overlap there by at least a set
## fraction of each one's length; the copy columns are given for 10% and for 50%.
##
## args: combined transferred peaks (asm_chr, start, end, consensus_peak_id; all assemblies, with a
## header), graph prefix file, CHM13 chromosome list, CHM13 consensus peaks bed, output prefix, threads
## out: <prefix>_previous_peak_checks.tsv (one row per non-reference peak) and
##      <prefix>_previous_peak_intervals_per_asm.tsv.gz (one row per peak and assembly it transferred to)
library(data.table)

args = commandArgs(trailingOnly = TRUE)
combined_file = args[1]
sample_prefix_file_name = args[2]
ref_chr_file = args[3]
chm13_peaks_file = args[4]
out_prefix = args[5]
threads = if(length(args) >= 6) as.integer(args[6]) else 1L
setDTthreads(threads)

## a peak is a copy of a CHM13 peak where they overlap by at least this percent of both peaks' lengths
## (one set of copy columns per threshold)
copy_thresholds = c(10, 50)

progress = function(msg){
  print(paste(format(Sys.time(), "%H:%M:%S"), msg))
}

## the input-file asm a contig belongs to: its "#" prefix, or for assemblies without one
## (e.g. HG002's "chr10_MATERNAL") the subtag the name contains. NA if none matches
contig_asm = function(chr, asms){
  out = ifelse(grepl("#", chr), sub("^([^#]*)#([^#]*).*$", "\\1#\\2#", chr), NA_character_)
  for(a in asms[!grepl("#", asms)]){
    out[is.na(out) & grepl(a, chr, fixed = TRUE)] = a
  }
  out
}

ref_chrs = readLines(ref_chr_file)
sample_prefix_file = read.table(sample_prefix_file_name, sep = "\t", header = F, stringsAsFactors = F, comment.char = "")
colnames(sample_prefix_file) = c("asm", "graph_asm", "path_to_fasta")

## every CHM13 consensus peak, ranked along its chromosome
chm13 = fread(chm13_peaks_file, select = 1:3)
setnames(chm13, c("chr", "start", "end"))
setorder(chm13, chr, start, end)
chm13[, id := paste(chr, start, end, sep = "_")]
chm13[, rank := seq_len(.N), by = chr]
chm13[, key := paste(chr, rank)]
rank_to_id = function(chr, rank){
  chm13$id[match(paste(chr, rank), chm13$key)]
}
rank_to_start = function(chr, rank){
  chm13$start[match(paste(chr, rank), chm13$key)]
}

progress("reading transferred peaks")
x = fread(combined_file, header = TRUE, select = 1:4)
setnames(x, c("asm_chr", "asm_start", "asm_end", "consensus_peak_id"))
ids = unique(x$consensus_peak_id)
id_chr = sub("_[0-9]+_[0-9]+$", "", ids)
x[, cons_chr := id_chr[match(consensus_peak_id, ids)]]
x[, is_ref := cons_chr %in% ref_chrs]
## only contigs with a non-reference peak matter
x = x[asm_chr %in% unique(x[is_ref == FALSE]$asm_chr)]
x[, asm := sample_prefix_file$graph_asm[match(contig_asm(asm_chr, sample_prefix_file$asm), sample_prefix_file$asm)]]
stopifnot(!anyNA(x$asm))
x[is_ref == TRUE, rank := chm13$rank[match(consensus_peak_id, chm13$id)]]
if(anyNA(x[is_ref == TRUE]$rank)){
  print(paste("CHM13 peaks missing from the consensus peak list, ignored:", sum(is.na(x[is_ref == TRUE]$rank))))
  x = x[is_ref == FALSE | !is.na(rank)]
}
progress(paste(nrow(x), "rows on", uniqueN(x$asm_chr), "contigs"))

## neighbouring CHM13 peaks on the same contig
setorder(x, asm_chr, asm_start, asm_end)
x[, row := .I]
x[, ref_row := fifelse(is_ref, row, NA_integer_)]
x[, `:=`(prev_row = nafill(ref_row, "locf"), next_row = nafill(ref_row, "nocb")), by = asm_chr]

## each contig's most common chr and its direction along it (as in the called-assembly script)
refs = x[is_ref == TRUE]
contig_chr = refs[, .(n = uniqueN(consensus_peak_id)), by = .(asm_chr, cons_chr)][order(asm_chr, -n)][, .SD[1], by = asm_chr]
setnames(contig_chr, "cons_chr", "contig_chr")
contig_dir = merge(refs, contig_chr[, .(asm_chr, contig_chr)], by = "asm_chr")[cons_chr == contig_chr][order(asm_chr, asm_start)][
  , .(n_up = sum(diff(rank) > 0), n_down = sum(diff(rank) < 0)), by = asm_chr]
contig_dir[, contig_orientation := fifelse(n_down > n_up, "R", "F")]

progress("one interval per peak and assembly")
o = x[is_ref == FALSE, .(row, consensus_peak_id, cons_chr, asm, asm_chr, asm_start, asm_end, prev_row, next_row)]
o[, `:=`(p_chr = x$cons_chr[prev_row], p_rank = x$rank[prev_row], p_end = x$asm_end[prev_row],
         n_chr = x$cons_chr[next_row], n_rank = x$rank[next_row], n_start = x$asm_start[next_row])]
o = merge(o, contig_chr[, .(asm_chr, contig_chr)], by = "asm_chr", all.x = TRUE)
o = merge(o, contig_dir[, .(asm_chr, contig_orientation)], by = "asm_chr", all.x = TRUE)
o[, `:=`(lo = NA_real_, hi = NA_real_, kind = "no CHM13 peak on contig", local_orientation = NA_character_, chr = NA_character_)]
## both sides on the same chr
b = which(!is.na(o$p_chr) & !is.na(o$n_chr) & o$p_chr == o$n_chr)
o[b, `:=`(chr = p_chr, lo = pmin(p_rank, n_rank), hi = pmax(p_rank, n_rank), kind = "between",
          local_orientation = fifelse(p_rank < n_rank, "F", fifelse(p_rank > n_rank, "R", NA_character_)))]
o[b[which(o$p_rank[b] == o$n_rank[b])], `:=`(kind = "same CHM13 peak on both sides", lo = NA, hi = NA, chr = NA)]
## otherwise one side, on the contig's chr, read with the contig's direction
one = which(o$kind == "no CHM13 peak on contig" & (!is.na(o$p_chr) | !is.na(o$n_chr)))
o[one, kind := fifelse(!is.na(p_chr) & !is.na(n_chr), "sides on different chrs", "contig end")]
pf = one[which(o$p_chr[one] == o$contig_chr[one])]
nf = setdiff(one[which(o$n_chr[one] == o$contig_chr[one])], pf)
o[pf, `:=`(chr = p_chr, lo = fifelse(contig_orientation == "F", p_rank, -Inf), hi = fifelse(contig_orientation == "F", Inf, p_rank))]
o[nf, `:=`(chr = n_chr, lo = fifelse(contig_orientation == "F", -Inf, n_rank), hi = fifelse(contig_orientation == "F", n_rank, Inf))]
o[setdiff(one, c(pf, nf)), kind := "neighbours not on the contig's chr"]
## the peak and a neighbouring CHM13 peak overlap on the assembly, so their order there is a guess
o[, overlaps_neighbour := (!is.na(p_end) & p_end > asm_start) | (!is.na(n_start) & n_start < asm_end)]

## which assembly each peak was called on
called = unique(o[, .(consensus_peak_id, cons_chr)])
called[, called_asm := sample_prefix_file$graph_asm[match(contig_asm(cons_chr, sample_prefix_file$asm), sample_prefix_file$asm)]]
o = merge(o, called[, .(consensus_peak_id, called_asm)], by = "consensus_peak_id")
## a peak in the called assembly sits on its own contig
o[, is_called_asm := asm == called_asm & asm_chr == cons_chr]

progress("overlap with CHM13 peaks on each assembly")
## bed intervals are half-open; foverlaps wants closed ones
nr = o[, .(row, asm_chr, start = asm_start, end = asm_end - 1)]
rf = x[is_ref == TRUE, .(asm_chr, start = asm_start, end = asm_end - 1, ref_id = consensus_peak_id)]
setkey(rf, asm_chr, start, end)
ov = foverlaps(nr, rf, nomatch = NULL)
ov[, overlap_bp := pmin(end, i.end) - pmax(start, i.start) + 1]
ov[, reciprocal := pmin(overlap_bp / (i.end - i.start + 1), overlap_bp / (end - start + 1))]
ov = ov[order(row, -reciprocal)][, .SD[1], by = row]
o = merge(o, ov[, .(row, overlapping_CHM13_peak = ref_id, reciprocal_overlap = reciprocal)], by = "row", all.x = TRUE)
for(t in copy_thresholds){
  o[, paste0("is_copy_", t) := !is.na(reciprocal_overlap) & reciprocal_overlap >= t / 100]
}

progress("comparing assemblies")
## the chr: the called assembly's if it has one, else the most common
inf = o[!is.na(lo)]
chr_votes = inf[, .(n = uniqueN(asm), called = any(is_called_asm)), by = .(consensus_peak_id, chr)]
chr_pick = chr_votes[order(consensus_peak_id, -called, -n)][, .(chr = chr[1], n_asms_other_chr = sum(n[-1])), by = consensus_peak_id]
on = merge(inf, chr_pick[, .(consensus_peak_id, chosen = chr)], by = "consensus_peak_id")[chr == chosen]

own = on[is_called_asm == TRUE, .(olo = max(lo), ohi = min(hi), o_kind = kind[1], o_local = local_orientation[1],
                                  o_contig = contig_orientation[1], o_ov = any(overlaps_neighbour)), by = consensus_peak_id]
oth = merge(on[is_called_asm == FALSE], own, by = "consensus_peak_id")
oth[, agrees := lo < ohi & olo < hi]
## distance between a disagreeing interval and the called one (0: one says after CHM13 peak X, the
## other before X)
oth[, gap_peaks := fifelse(agrees, NA_real_, fifelse(hi <= olo, olo - hi, lo - ohi))]
oth[, gap_bp := fifelse(agrees, NA_real_, fifelse(hi <= olo, rank_to_start(chr, olo) - rank_to_start(chr, hi), rank_to_start(chr, lo) - rank_to_start(chr, ohi)))]
oth[agrees == FALSE, cause := fcase(
  o_ov | overlaps_neighbour, "overlaps_CHM13_peak",
  kind == "between" & o_kind == "between" & local_orientation != o_local, "inversion",
  !is.finite(lo) | !is.finite(hi) | !is.finite(olo) | !is.finite(ohi), "contig_end",
  gap_peaks == 0, "adjacent",
  gap_peaks <= 5, "nearby",
  default = "distant")]

checks = oth[, .(n_asms_checked = uniqueN(asm), n_asms_disagree = uniqueN(asm[agrees == FALSE]),
                 conflict_type = if(any(!agrees)) names(sort(table(cause[!agrees]), decreasing = TRUE))[1] else NA_character_,
                 conflict_gap_peaks = if(any(!agrees)) max(gap_peaks, na.rm = TRUE) else NA_real_,
                 conflict_gap_bp = if(any(!agrees)) max(abs(gap_bp), na.rm = TRUE) else NA_real_), by = consensus_peak_id]
checks[, called_is_minority := n_asms_disagree > n_asms_checked / 2]
## minor: the disagreement is only about which side of one CHM13 peak, or comes from overlapping peaks
checks[, conflict_level := fifelse(n_asms_disagree == 0, "none",
                           fifelse(conflict_type %in% c("overlaps_CHM13_peak", "adjacent") | conflict_gap_peaks <= 1, "minor", "major"))]

## the tightest placement allowed by every assembly, and the called assembly's own interval
all_asms = on[, .(L = max(lo), U = min(hi)), by = consensus_peak_id]
all_asms[, all_asms_agree := L < U]
called_iv = own[, .(consensus_peak_id, olo, ohi)]

## copies of CHM13 peaks: the CHM13 peak it copies on the called assembly, the number of assemblies
## where it is a copy, and a likely copy if it copies the same CHM13 peak in at least half of the
## assemblies it is in
copies_at = function(t){
  is_copy = o[[paste0("is_copy_", t)]]
  d = o[, .(n_asms_present = uniqueN(asm),
            copy_in_called_asm = if(any(is_called_asm & is_copy[.I])) overlapping_CHM13_peak[is_called_asm & is_copy[.I]][1] else NA_character_,
            n_asms_copy = uniqueN(asm[is_copy[.I]]),
            most_common_copy = if(any(is_copy[.I])) names(sort(table(overlapping_CHM13_peak[is_copy[.I]]), decreasing = TRUE))[1] else NA_character_),
        by = consensus_peak_id]
  n_same = o[is_copy, .(n_asms_most_common_copy = uniqueN(asm)), by = .(consensus_peak_id, most_common_copy = overlapping_CHM13_peak)]
  d = merge(d, n_same, by = c("consensus_peak_id", "most_common_copy"), all.x = TRUE)
  d[is.na(n_asms_most_common_copy), n_asms_most_common_copy := 0L]
  d[, likely_copy_of_CHM13_peak := fifelse(n_asms_most_common_copy > 0 & n_asms_most_common_copy * 2 >= n_asms_present, most_common_copy, NA_character_)]
  d = d[, .(consensus_peak_id, n_asms_present, copy_in_called_asm, n_asms_copy, likely_copy_of_CHM13_peak)]
  setnames(d, c("copy_in_called_asm", "n_asms_copy", "likely_copy_of_CHM13_peak"), paste0(c("copy_in_called_asm_", "n_asms_copy_", "likely_copy_of_CHM13_peak_"), t, "pct"))
  d
}
copies = Reduce(function(a, b) merge(a, b, by = c("consensus_peak_id", "n_asms_present")), lapply(copy_thresholds, copies_at))
copy_cols = setdiff(names(copies), c("consensus_peak_id", "n_asms_present"))

out = Reduce(function(a, b) merge(a, b, by = "consensus_peak_id", all.x = TRUE),
             list(copies, chr_pick, called_iv, all_asms, checks))
interval_ids = function(chr, lo, hi){
  list(fifelse(is.finite(lo), rank_to_id(chr, lo), NA_character_),
       fifelse(is.finite(hi), rank_to_id(chr, hi), NA_character_),
       fifelse(is.finite(lo) & is.finite(hi), hi - lo - 1, NA_real_))
}
out[, c("called_after_peak", "called_before_peak", "called_peaks_between") := interval_ids(chr, olo, ohi)]
out[all_asms_agree == TRUE, c("all_asms_after_peak", "all_asms_before_peak", "all_asms_peaks_between") := interval_ids(chr, L, U)]
out[is.na(chr), chr := "AMB"]
## peaks with no other assembly to compare with
out[is.na(n_asms_checked), `:=`(n_asms_checked = 0L, n_asms_disagree = 0L)]
out[is.na(conflict_level), conflict_level := fifelse(is.na(olo), "not_placed", "not_checked")]
out = out[, c("consensus_peak_id", "chr", "n_asms_present", "called_after_peak", "called_before_peak", "called_peaks_between",
              "n_asms_checked", "n_asms_disagree", "conflict_level", "conflict_type", "conflict_gap_peaks", "conflict_gap_bp",
              "called_is_minority", "n_asms_other_chr", "all_asms_agree", "all_asms_after_peak", "all_asms_before_peak",
              "all_asms_peaks_between", copy_cols), with = FALSE]
fwrite(out, paste0(out_prefix, "_previous_peak_checks.tsv"), sep = "\t", na = "NA", quote = FALSE)

oth_cols = oth[, .(row, agrees_with_called = agrees, gap_peaks, cause)]
o = merge(o, oth_cols, by = "row", all.x = TRUE)
o[, `:=`(after_peak = fifelse(is.finite(lo), rank_to_id(chr, lo), NA_character_),
         before_peak = fifelse(is.finite(hi), rank_to_id(chr, hi), NA_character_))]
per_asm_cols = c("consensus_peak_id", "called_asm", "asm", "is_called_asm", "asm_chr", "asm_start", "asm_end", "kind", "chr",
                 "after_peak", "before_peak", "local_orientation", "contig_chr", "contig_orientation", "overlaps_neighbour",
                 "overlapping_CHM13_peak", "reciprocal_overlap", paste0("is_copy_", copy_thresholds), "agrees_with_called", "gap_peaks", "cause")
fwrite(o[order(consensus_peak_id, -is_called_asm, asm), per_asm_cols, with = FALSE],
       paste0(out_prefix, "_previous_peak_intervals_per_asm.tsv.gz"), sep = "\t", na = "NA", quote = FALSE)
print(out[, .N, by = conflict_level])
for(t in copy_thresholds){
  print(paste0("likely copies of a CHM13 peak at ", t, "% overlap: ", sum(!is.na(out[[paste0("likely_copy_of_CHM13_peak_", t, "pct")]]))))
}
progress("done")
