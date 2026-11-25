## Counts the original peaks on each consensus peak, per sample (used by both metadata scripts).
## Two or more original peaks on one consensus peak count as separate peaks only where they're far
## enough apart to be distinct (e.g. copies from a duplication) rather than pieces of one peak:
##   - two peaks: separate if the gap between them is at least the consensus peak's length, or they
##     are on different contigs
##   - three or more: each peak whose gap to the previous one is longer than the consensus peak's
##     length starts a new peak
##
## ov: one row per original peak overlapping a consensus peak, with columns
##   group: one value per sample and consensus peak
##   key: that sample and consensus peak's row key in the final table
##   chr_orig, start_orig, end_orig, peak_id_orig: the original peak
##   start_cons, end_cons: the consensus peak
## Returns a list:
##   pair_keys: keys with two original peaks that count as 2
##   gr2_keys, gr2_num_peaks: keys with three or more original peaks, and how many they count as
##   summary: key and all_orig_overlapping_peaks (the original peak ids, comma-separated, in the
##            order of ov), one row per group
count_orig_peaks_per_consensus_peak = function(ov){
  temp = as.data.frame(table(ov$group))
  temp_filt = temp[which(temp$Freq > 1),]
  temp_filt_2 = temp[which(temp$Freq == 2),]
  temp_filt_gr2 = temp[which(temp$Freq > 2),]

  mult_to_one = ov[which(ov$group %in% temp_filt$Var1),]
  mult_to_one = mult_to_one[order(as.numeric(mult_to_one$start_orig)),]
  mult_to_one = mult_to_one[order(mult_to_one$group),]

  mult_to_one_2 = mult_to_one[which(mult_to_one$group %in% temp_filt_2$Var1),]
  mult_to_one_gr2 = mult_to_one[which(mult_to_one$group %in% temp_filt_gr2$Var1),]

  ## pairs of original peaks on one consensus peak. The table is sorted by group and then start, so
  ## each pair is two neighboring rows, the lower start first
  pair_keys = c()
  if(nrow(mult_to_one_2) > 0){
    first = mult_to_one_2[seq(1, nrow(mult_to_one_2), by = 2),]
    second = mult_to_one_2[seq(2, nrow(mult_to_one_2), by = 2),]
    stopifnot(first$group == second$group)
    space_between_peaks = as.numeric(second$start_orig) - as.numeric(first$end_orig)
    share_orig_chr = !is.na(first$chr_orig == second$chr_orig) & first$chr_orig == second$chr_orig
    cons_peak_length = as.numeric(first$end_cons) - as.numeric(first$start_cons)
    dist_less_than_cons_peak = !is.na(space_between_peaks < cons_peak_length) & space_between_peaks < cons_peak_length
    pair_keys = first$key[which(dist_less_than_cons_peak == FALSE | share_orig_chr == FALSE)]
  }

  ## three or more original peaks on one consensus peak: count them as separate peaks only where the
  ## gap to the previous one is longer than the consensus peak
  gr2_keys = c()
  gr2_num_peaks = c()
  if(nrow(mult_to_one_gr2) > 0){
    gr2 = mult_to_one_gr2
    cons_peak_length = gr2$end_cons - gr2$start_cons
    group_start = !duplicated(gr2$group)
    group = cumsum(group_start)
    start_orig = as.numeric(gr2$start_orig)
    prev_end_orig = c(NA, as.numeric(gr2$end_orig)[-nrow(gr2)])
    new_peak = !group_start & (start_orig - prev_end_orig > cons_peak_length[group_start][group])
    gr2_keys = gr2$key[group_start]
    gr2_num_peaks = 1 + as.vector(rowsum(as.integer(new_peak), group, reorder = FALSE))
  }

  summary = as.data.frame(as.data.table(ov)[, .(key = key[1], all_orig_overlapping_peaks = paste(peak_id_orig, collapse = ",")), by = group])
  list(pair_keys = pair_keys, gr2_keys = gr2_keys, gr2_num_peaks = gr2_num_peaks, summary = summary)
}
