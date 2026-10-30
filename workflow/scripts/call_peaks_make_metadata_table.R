library(ggplot2)
library(dplyr)
library(tidyr)
library(data.table)

## These parse ids like chrom_start_end. Chrom names can contain "_", so start and end
## are the last two fields. They work on whole vectors, parsing each distinct id once.
id_chr = function(a){
  u = unique(a)
  sub("_[^_]*_[^_]*$", "", u)[match(a, u)]
}

id_start = function(a){
  u = unique(a)
  sub("^.*_([^_]*)_[^_]*$", "\\1", u)[match(a, u)]
}

id_end = function(a){
  u = unique(a)
  sub("^.*_", "", u)[match(a, u)]
}

## "SMHT004#1#haplotype1-0000001" -> "SMHT004#1#"
id_chr_w_hash = function(a){
  sub("^([^#]*)#([^#]*).*$", "\\1#\\2#", a)
}

## one number per (a, b) pair, so large tables can be matched without pasting strings
## together. a and b are indices from match(); the pair is NA if either is NA
pair_key = function(a, b, n_b){
  a * (n_b + 1) + b
}

progress = function(msg){
  print(paste(format(Sys.time(), "%H:%M:%S"), msg))
}

read_tsv = function(file, header, unique_rows = FALSE, ...){
  x = fread(file, header = header, sep = "\t", quote = "", na.strings = "NA", ...)
  if(unique_rows){
    x = unique(x)
  }
  as.data.frame(x)
}

args <- commandArgs(trailingOnly = TRUE)

input_dir=args[1]
output_dir=input_dir
prefix=args[2]
sample_pefix_file_name=args[3]
threads=if(length(args) >= 4) as.integer(args[4]) else 1
setDTthreads(threads)

progress("reading inputs")
tbl = read_tsv(paste0(input_dir, prefix, "_max_score_per_peak_all_samples.bed"), header = T)
## a donor's assembly can be listed once per run; only the first copy of a peak is ever used
trans_peaks = read_tsv(paste0(input_dir, prefix, "_combined_tranferred_peaks_all_samples.bed"), header = T, unique_rows = T)
sample_prefix_file = read.table(sample_pefix_file_name, sep = "\t", header = F, stringsAsFactors = F, comment.char = "")
colnames(sample_prefix_file) = c("asm", "graph_asm", "path_to_fasta")

hg38_file=read_tsv(paste0(input_dir, prefix, "_All_consensus_peaks_moved_to_GRCh38.bed"), header=F)
colnames(hg38_file) = c("hg38_chr", "hg38_start", "hg38_end", "consensus_peak_id")
## only the run and the transferred peak are used from this file
orig_peak_trans_peak_intersect=read_tsv(paste0(input_dir, prefix, "_orig_called_peaks_intersect_transferred_consensus_peaks.bed"), header=F, unique_rows = T, select = c(4, 6, 7, 8))
colnames(orig_peak_trans_peak_intersect)=c("sample_id", "trans_chr", "trans_start", "trans_end")
orig_peak_trans_peak_intersect$trans_peak_id = paste0(orig_peak_trans_peak_intersect$trans_chr, "_", orig_peak_trans_peak_intersect$trans_start, "_",orig_peak_trans_peak_intersect$trans_end)

surj_peak_cons_peak_intersect=read_tsv(paste0(input_dir,prefix, "_all_orig_surjected_peaks_intersect_consensus_peaks.bed"), header=F)
colnames(surj_peak_cons_peak_intersect) = c("chr_surj", "start_surj", "end_surj", "sample_chr_id", "peak_id_orig", "chr_cons", "start_cons", "end_cons", "peak_id_cons" )

progress("adding assembly coordinates")
tbl$cons_chr = id_chr(tbl$consensus_peak_id)
tbl$cons_start = as.numeric(id_start(tbl$consensus_peak_id))
tbl$cons_end = as.numeric(id_end(tbl$consensus_peak_id))

tbl$fire_cov_OV_cov = tbl$fire_coverage/tbl$coverage

## sample_id is run_graphAsm, and graph_asm is its last two "_" fields
unique_sample_id = unique(tbl$sample_id)
sample_graph_asm = paste0(id_start(unique_sample_id), "_", id_end(unique_sample_id))
sample_asm = sample_prefix_file$asm[match(sample_graph_asm, sample_prefix_file$graph_asm)]
sample_run = substr(unique_sample_id, 1, nchar(unique_sample_id) - nchar(sample_graph_asm) - 1)
tbl_sample_idx = match(tbl$sample_id, unique_sample_id)
tbl$asm = sample_asm[tbl_sample_idx]

## rows are grouped by sample (in input order); within a sample, transferred peaks come
## first and peaks that could not be transferred (score -2) last
transferred = tbl$score != -2
row_order = order(tbl_sample_idx, !transferred, seq_len(nrow(tbl)))
row_order = row_order[!is.na(transferred[row_order])]
updated_tbl = tbl[row_order, ]
rownames(updated_tbl) = NULL
rm(tbl)
updated_sample_idx = tbl_sample_idx[row_order]
transferred = transferred[row_order]
rm(tbl_sample_idx, row_order)

## each transferred consensus peak gets the coordinates of its first copy on that assembly
## (assemblies are matched to contig names the same way as before: grepl of the asm name)
updated_tbl$asm_chr = NA_character_
updated_tbl$asm_start = NA_integer_
updated_tbl$asm_end = NA_integer_
trans_contigs = unique(trans_peaks$chr)
for(temp_asm in unique(sample_asm)){
  print(temp_asm)
  temp_trans_peaks = trans_peaks[which(trans_peaks$chr %in% trans_contigs[grepl(temp_asm, trans_contigs)]),]
  rows = which(transferred & updated_tbl$asm == temp_asm)
  m = match(updated_tbl$consensus_peak_id[rows], temp_trans_peaks$cons_peak_id)
  updated_tbl$asm_chr[rows] = temp_trans_peaks$chr[m]
  updated_tbl$asm_start[rows] = temp_trans_peaks$start[m]
  updated_tbl$asm_end[rows] = temp_trans_peaks$end[m]
}
rm(trans_peaks, transferred)

updated_tbl$asm_peak_id = paste0(updated_tbl$asm_chr, "_", updated_tbl$asm_start, "_", updated_tbl$asm_end)

### add hg38 info

progress("adding hg38 info")
hg38_file$cons_chr = id_chr(hg38_file$consensus_peak_id)
hg38_file$cons_start = as.numeric(id_start(hg38_file$consensus_peak_id))
hg38_file$cons_end = as.numeric(id_end(hg38_file$consensus_peak_id))

hg38_file$cons_peak_length = hg38_file$cons_end - hg38_file$cons_start
hg38_file$hg38_peak_length = hg38_file$hg38_end - hg38_file$hg38_start

hg38_file$length_diff = abs(hg38_file$hg38_peak_length-hg38_file$cons_peak_length)
hg38_file$percent_diff = hg38_file$length_diff/hg38_file$cons_peak_length*100

hg38_file_pass_filt = hg38_file[which(hg38_file$percent_diff < 20),]
hg38_file_pass_filt_for_merge = hg38_file_pass_filt[,c(1:4)]
if(anyDuplicated(hg38_file_pass_filt_for_merge$consensus_peak_id) == 0){
  ## one hg38 location per consensus peak, so a lookup gives the same rows as a join
  m = match(updated_tbl$consensus_peak_id, hg38_file_pass_filt_for_merge$consensus_peak_id)
  updated_tbl$hg38_chr = hg38_file_pass_filt_for_merge$hg38_chr[m]
  updated_tbl$hg38_start = hg38_file_pass_filt_for_merge$hg38_start[m]
  updated_tbl$hg38_end = hg38_file_pass_filt_for_merge$hg38_end[m]
} else {
  updated_tbl <- left_join(updated_tbl, hg38_file_pass_filt_for_merge, by = "consensus_peak_id")
  updated_sample_idx = match(updated_tbl$sample_id, unique_sample_id)
}
updated_tbl$exists_in_hg38 = FALSE
updated_tbl$exists_in_hg38[which(!is.na(updated_tbl$hg38_chr))] = TRUE

## add info about whether this run's own original peaks overlap the transferred peak
progress("adding original peak overlaps")
trans_peak_ids = unique(orig_peak_trans_peak_intersect$trans_peak_id)
orig_keys = pair_key(match(orig_peak_trans_peak_intersect$sample_id, sample_run),
                     match(orig_peak_trans_peak_intersect$trans_peak_id, trans_peak_ids), length(trans_peak_ids))
row_keys = pair_key(match(sample_run, sample_run)[updated_sample_idx],
                    match(updated_tbl$asm_peak_id, trans_peak_ids), length(trans_peak_ids))
updated_tbl$overlaps_original_called_peak = !is.na(row_keys) & row_keys %in% orig_keys
rm(orig_peak_trans_peak_intersect, orig_keys, row_keys)

## add multiple to one info

progress("counting overlapping original peaks")
## peak_id_orig is chrom_start_end-run, so drop "-run" before parsing the coordinates.
## The run is sample_chr_id, so this works even if run or chrom names contain "-"
run_suffix = paste0("-", surj_peak_cons_peak_intersect$sample_chr_id)
stopifnot(endsWith(surj_peak_cons_peak_intersect$peak_id_orig, run_suffix))
surj_peak_cons_peak_intersect$coords_orig = substr(surj_peak_cons_peak_intersect$peak_id_orig, 1, nchar(surj_peak_cons_peak_intersect$peak_id_orig) - nchar(run_suffix))
surj_peak_cons_peak_intersect$chr_orig = id_chr(surj_peak_cons_peak_intersect$coords_orig)
surj_peak_cons_peak_intersect$start_orig = id_start(surj_peak_cons_peak_intersect$coords_orig)
surj_peak_cons_peak_intersect$end_orig = id_end(surj_peak_cons_peak_intersect$coords_orig)

surj_peak_cons_peak_intersect$samp_chr_consPeak = paste0(surj_peak_cons_peak_intersect$sample_chr_id, "-", surj_peak_cons_peak_intersect$chr_orig, "-", surj_peak_cons_peak_intersect$peak_id_cons)

surj_peak_cons_peak_intersect_wHash = surj_peak_cons_peak_intersect[which(grepl("#",surj_peak_cons_peak_intersect$chr_orig)),]
surj_peak_cons_peak_intersect_NoHash = surj_peak_cons_peak_intersect[which(!(grepl("#",surj_peak_cons_peak_intersect$chr_orig))),]

if(nrow(surj_peak_cons_peak_intersect_wHash) > 0){
  surj_peak_cons_peak_intersect_wHash$asm = id_chr_w_hash(surj_peak_cons_peak_intersect_wHash$chr_orig)
}

if(nrow(surj_peak_cons_peak_intersect_NoHash) > 0){
  nonHash_asms = sample_prefix_file[which(!grepl("#",sample_prefix_file$asm)),]$asm

  surj_peak_cons_peak_intersect_NoHash$asm = "NONE"

  for(i in 1:length(nonHash_asms)){
    temp_asm = nonHash_asms[i]
    print(temp_asm)
    surj_peak_cons_peak_intersect_NoHash$asm[which(grepl(temp_asm,surj_peak_cons_peak_intersect_NoHash$chr_orig))] = temp_asm
  }
  if(nrow(surj_peak_cons_peak_intersect_wHash) > 0){
    surj_peak_cons_peak_intersect = as.data.frame(rbind(surj_peak_cons_peak_intersect_wHash, surj_peak_cons_peak_intersect_NoHash))
  } else {
    surj_peak_cons_peak_intersect = surj_peak_cons_peak_intersect_NoHash
  }
} else {
  surj_peak_cons_peak_intersect = surj_peak_cons_peak_intersect_wHash
}

sample_prefix_file_2col = sample_prefix_file[,c("asm", "graph_asm")]
surj_peak_cons_peak_intersect = left_join(surj_peak_cons_peak_intersect, sample_prefix_file_2col, by = "asm")

surj_peak_cons_peak_intersect$sample_graphHap_peak = paste0(surj_peak_cons_peak_intersect$sample_chr_id, "_", surj_peak_cons_peak_intersect$graph_asm, "-", surj_peak_cons_peak_intersect$peak_id_cons)

## the rows of updated_tbl are keyed by (sample_id, consensus_peak_id); sample_graphHap_peak
## is the same pair pasted together, so give both tables the same numeric key
unique_cons_peak_id = unique(updated_tbl$consensus_peak_id)
updated_tbl$sample_id_peak_id_TEMP_COL = pair_key(updated_sample_idx, match(updated_tbl$consensus_peak_id, unique_cons_peak_id), length(unique_cons_peak_id))
rm(updated_sample_idx)
surj_peak_cons_peak_intersect$key = pair_key(match(paste0(surj_peak_cons_peak_intersect$sample_chr_id, "_", surj_peak_cons_peak_intersect$graph_asm), unique_sample_id),
                                             match(surj_peak_cons_peak_intersect$peak_id_cons, unique_cons_peak_id), length(unique_cons_peak_id))

temp = as.data.frame(table(surj_peak_cons_peak_intersect$sample_graphHap_peak))
temp_filt = temp[which(temp$Freq > 1),]
temp_filt_2 = temp[which(temp$Freq == 2),]
temp_filt_gr2 = temp[which(temp$Freq >2),]

surj_peak_cons_peak_intersect_mult_to_one = surj_peak_cons_peak_intersect[which(surj_peak_cons_peak_intersect$sample_graphHap_peak %in% temp_filt$Var1),]
surj_peak_cons_peak_intersect_mult_to_one = surj_peak_cons_peak_intersect_mult_to_one[order(as.numeric(surj_peak_cons_peak_intersect_mult_to_one$start_orig)),]
surj_peak_cons_peak_intersect_mult_to_one = surj_peak_cons_peak_intersect_mult_to_one[order(surj_peak_cons_peak_intersect_mult_to_one$sample_graphHap_peak),]

surj_peak_cons_peak_intersect_mult_to_one_2 = surj_peak_cons_peak_intersect_mult_to_one[which(surj_peak_cons_peak_intersect_mult_to_one$sample_graphHap_peak %in% temp_filt_2$Var1),]
surj_peak_cons_peak_intersect_mult_to_one_gr2 = surj_peak_cons_peak_intersect_mult_to_one[which(surj_peak_cons_peak_intersect_mult_to_one$sample_graphHap_peak %in% temp_filt_gr2$Var1),]

## pairs of original peaks on one consensus peak. The table is sorted by sample_graphHap_peak
## and then start, so each pair is two neighboring rows, the lower start first
true_mult_to_one_keys = c()
if(nrow(surj_peak_cons_peak_intersect_mult_to_one_2) > 0){
  first = surj_peak_cons_peak_intersect_mult_to_one_2[seq(1, nrow(surj_peak_cons_peak_intersect_mult_to_one_2), by = 2),]
  second = surj_peak_cons_peak_intersect_mult_to_one_2[seq(2, nrow(surj_peak_cons_peak_intersect_mult_to_one_2), by = 2),]
  stopifnot(first$sample_graphHap_peak == second$sample_graphHap_peak)
  space_between_peaks = as.numeric(second$start_orig) - as.numeric(first$end_orig)
  share_orig_chr = !is.na(first$chr_orig == second$chr_orig) & first$chr_orig == second$chr_orig
  ## the surjected peak is the end of sample_graphHap_peak: ..._start_end
  surj_peak_length = as.numeric(id_end(first$sample_graphHap_peak)) - as.numeric(id_start(first$sample_graphHap_peak))
  dist_less_than_cons_peak = !is.na(space_between_peaks < surj_peak_length) & space_between_peaks < surj_peak_length
  true_mult_to_one_keys = first$key[which(dist_less_than_cons_peak == FALSE | share_orig_chr == FALSE)]
}

### greater than 2:

## three or more original peaks on one consensus peak: count them as separate peaks only
## where the gap to the previous one is longer than the consensus peak
gr2_keys = c()
gr2_num_peaks = c()
if(nrow(surj_peak_cons_peak_intersect_mult_to_one_gr2) > 0){
  gr2 = surj_peak_cons_peak_intersect_mult_to_one_gr2
  cons_peak_length = gr2$end_cons - gr2$start_cons
  group_start = !duplicated(gr2$sample_graphHap_peak)
  group = cumsum(group_start)
  start_orig = as.numeric(gr2$start_orig)
  prev_end_orig = c(NA, as.numeric(gr2$end_orig)[-nrow(gr2)])
  new_peak = !group_start & (start_orig - prev_end_orig > cons_peak_length[group_start][group])
  gr2_keys = gr2$key[group_start]
  gr2_num_peaks = 1 + as.vector(rowsum(as.integer(new_peak), group, reorder = FALSE))
}

all_peaks_summary = as.data.frame(as.data.table(surj_peak_cons_peak_intersect)[, .(key = key[1], all_orig_overlapping_peaks = paste(peak_id_orig, collapse = ",")), by = sample_graphHap_peak])

progress("putting the table together")
updated_tbl$num_overlapping_orig_peaks = 0
updated_tbl$num_overlapping_orig_peaks[which(updated_tbl$sample_id_peak_id_TEMP_COL %in% surj_peak_cons_peak_intersect$key)] = 1
updated_tbl$num_overlapping_orig_peaks[which(updated_tbl$sample_id_peak_id_TEMP_COL %in% true_mult_to_one_keys)] = 2

## rows with 3+ overlapping peaks get their count and move after the other rows
m = match(updated_tbl$sample_id_peak_id_TEMP_COL, gr2_keys)
in_gr2 = !is.na(m)
updated_tbl$num_overlapping_orig_peaks[in_gr2] = gr2_num_peaks[m[in_gr2]]

m = match(updated_tbl$sample_id_peak_id_TEMP_COL, all_peaks_summary$key)
non0 = updated_tbl$num_overlapping_orig_peaks > 0
updated_tbl$orig_peaks_overlapping_consensus = ifelse(non0, all_peaks_summary$all_orig_overlapping_peaks[m], NA)

## rows with overlapping peaks first, then the rest; each in the order above
final_order = c(which(!in_gr2), which(in_gr2))
final_order = c(final_order[non0[final_order]], final_order[!non0[final_order]])
rm(m, in_gr2, non0)
setDT(updated_tbl)
updated_tbl[, sample_id_peak_id_TEMP_COL := NULL]
updated_tbl = updated_tbl[final_order]

## sequences are in _transferred_cons_peaks_native_sequences_from_main_output.bed.gz, keyed by asm_peak_id
progress("writing the table")
## fwrite rounds the 15th digit differently from write.table, so write the ratio as R formats it
u = unique(updated_tbl$fire_cov_OV_cov)
updated_tbl[, fire_cov_OV_cov := as.character(u)[match(fire_cov_OV_cov, u)]]
## gzipped because of the .gz name; the table is large (~44 GB for 59 runs uncompressed)
fwrite(updated_tbl, paste0(output_dir, prefix, "_final_tbl_w_metadata.tsv.gz"), col.names = TRUE, row.names = FALSE, quote = F, sep = "\t", na = "NA")
progress("done")
