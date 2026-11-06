## this script takes in a combined peak file surjected into one assembly and attempts to place
# it in its correct spot if the peak was called in this asm
library(dplyr)
library(tidyr)
library(data.table)

## ids are chrom_start_end, and chrom names can contain "_", so start is the second-to-last field.
## x is NA when there is no ref peak before (or after) this one on the contig
return_start = function(x){
  splt = strsplit(x, "_")[[1]]
  if(length(splt) < 3){
    return(NA)
  }
  return(splt[length(splt)-1])
}

return_chr = function(x){
  splt = strsplit(x, "_")[[1]]
  if(length(splt) < 3){
    return(NA)
  }
  return(paste(splt[1:(length(splt)-2)], collapse = "_"))
}

args <- commandArgs(trailingOnly = TRUE)

input_bed=args[1]
sample_prefix_file_name=args[2]
asm=args[3]
ref_chr_file=args[4]
output_file=args[5]
## every CHM13 consensus peak (first three columns chrom, start, end; header allowed)
chm13_peaks_file=args[6]

#tbl = fread("/Users/minkinaa/Documents/StergachisLab/projects/HLA_peak_assignment/test_files/SMaHT_ElevenSampleConsPeakCalling_CLAUDETEST_ALL_global_peaks_surjected_onto_SMHT004_1_MANUAL_ADDITIONS.bed", header = F, sep = "\t")
#tbl = fread("/Users/minkinaa/Documents/StergachisLab/projects/HLA_peak_assignment/test_files/SMaHT_ElevenSampleConsPeakCalling_CLAUDETEST_ALL_global_peaks_surjected_onto_SMHT004_1.bed", header = F, sep = "\t")
tbl = as.data.frame(fread(input_bed, header = F, sep = "\t"))
colnames(tbl) = c("asm_chr", "asm_start", "asm_end", "consensus_peak_id")

## chromosome names in the reference (CHM13), one per line
ref_chrs = readLines(ref_chr_file)

## prefix of this assembly's contigs, e.g. SMHT004_1 -> "SMHT004#1#"
sample_prefix_file = read.table(sample_prefix_file_name, sep = "\t", header = F, stringsAsFactors = F, comment.char = "")
colnames(sample_prefix_file) = c("asm", "graph_asm", "path_to_fasta")
target_shared_prefix = sample_prefix_file$asm[which(sample_prefix_file$graph_asm == asm)]
stopifnot(length(target_shared_prefix) == 1)
print(paste("assembly", asm, "prefix", target_shared_prefix))

tbl$cons_chr = apply(tbl[,"consensus_peak_id",drop = F],1,return_chr)
tbl$cons_start = as.numeric(apply(tbl[,"consensus_peak_id",drop = F],1,return_start))
tbl$is_ref = tbl$cons_chr %in% ref_chrs

## table of just ref peaks in order by cons_chr and cons_start
tbl_ord_by_ref_coords <- tbl %>%
  arrange(cons_chr, cons_start)
tbl_ord_by_ref_coords = tbl_ord_by_ref_coords[which(tbl_ord_by_ref_coords$is_ref),]


### order by asm coordinates
tbl_ord_by_asm_coords <- tbl %>%
  arrange(asm_chr, asm_start)

### and we're only interested in arranging those consensus peaks that were called in this assembly.
## the subtag can be a prefix ("HG00438#1#CM089167.1") or appear later in the name ("chr10_MATERNAL"),
## so match it anywhere in the peak's chromosome, as the Snakefile's subtag filter does
tbl_ord_by_asm_coords <- tbl_ord_by_asm_coords[which(tbl_ord_by_asm_coords$is_ref | grepl(target_shared_prefix, tbl_ord_by_asm_coords$cons_chr, fixed = TRUE)),]

tbl_ord_by_asm_coords$order = seq_len(nrow(tbl_ord_by_asm_coords))

tbl_ord_by_asm_coords_NonRef = tbl_ord_by_asm_coords[which(!tbl_ord_by_asm_coords$is_ref),]
if(nrow(tbl_ord_by_asm_coords_NonRef) == 0){
  print("no peaks called on this assembly's own contigs, so nothing to place")
  write.table(data.frame(consensus_peak_id = character(), chr = character(), previous_peak = character(), orientation = character(),
                         previous_peak_full_list = character(), contig_position = character(), open_side = character(),
                         contig_anchor_peak = character()), output_file, col.names = TRUE, row.names = FALSE, quote = F, sep = "\t")
  quit(save = "no")
}
tbl_ord_by_asm_coords_Ref = tbl_ord_by_asm_coords[which(tbl_ord_by_asm_coords$is_ref),]
tbl_ord_by_asm_coords_Ref$prev_peak = tbl_ord_by_asm_coords_Ref$consensus_peak_id
tbl_ord_by_asm_coords_Ref$next_peak = tbl_ord_by_asm_coords_Ref$consensus_peak_id
tbl_ord_by_asm_coords_Ref$previous_index = tbl_ord_by_asm_coords_Ref$order
tbl_ord_by_asm_coords_Ref$next_index = tbl_ord_by_asm_coords_Ref$order


return_previous_ref_num = function(x,ref_num_list){
  return(max(ref_num_list[which(ref_num_list < x)]))
}

return_next_ref_num = function(x,ref_num_list){
  return(min(ref_num_list[which(ref_num_list > x)]))
}



## same as return_previous_ref_num / return_next_ref_num, but in one pass: the ref orders are
## sorted, so findInterval finds the last ref before each peak (NA if there is none)
ref_order = tbl_ord_by_asm_coords_Ref$order
n_refs_before = findInterval(tbl_ord_by_asm_coords_NonRef$order, ref_order)
tbl_ord_by_asm_coords_NonRef$previous_index = ifelse(n_refs_before == 0, NA, ref_order[pmax(n_refs_before, 1)])
tbl_ord_by_asm_coords_NonRef$next_index = ref_order[n_refs_before + 1]

tbl_ord_by_asm_coords_NonRef = left_join(tbl_ord_by_asm_coords_NonRef, tbl_ord_by_asm_coords_Ref[,c("prev_peak", "previous_index")], by = "previous_index")
tbl_ord_by_asm_coords_NonRef = left_join(tbl_ord_by_asm_coords_NonRef, tbl_ord_by_asm_coords_Ref[,c("next_peak", "next_index")], by = "next_index")

## only ref peaks on the peak's own contig count as its neighbours; past the end of the contig
## there is none, and the contig-end logic below places the peak instead
prev_contig = tbl_ord_by_asm_coords_Ref$asm_chr[match(tbl_ord_by_asm_coords_NonRef$previous_index, tbl_ord_by_asm_coords_Ref$order)]
next_contig = tbl_ord_by_asm_coords_Ref$asm_chr[match(tbl_ord_by_asm_coords_NonRef$next_index, tbl_ord_by_asm_coords_Ref$order)]
tbl_ord_by_asm_coords_NonRef$prev_peak[which(prev_contig != tbl_ord_by_asm_coords_NonRef$asm_chr)] = NA
tbl_ord_by_asm_coords_NonRef$next_peak[which(next_contig != tbl_ord_by_asm_coords_NonRef$asm_chr)] = NA

tbl_ord_by_asm_coords_NonRef$first_chr = apply(tbl_ord_by_asm_coords_NonRef[,"prev_peak",drop = F], 1, return_chr)
tbl_ord_by_asm_coords_NonRef$first_start = as.numeric(apply(tbl_ord_by_asm_coords_NonRef[,"prev_peak",drop = F], 1, return_start))

tbl_ord_by_asm_coords_NonRef$second_chr = apply(tbl_ord_by_asm_coords_NonRef[,"next_peak",drop = F], 1, return_chr)
tbl_ord_by_asm_coords_NonRef$second_start = as.numeric(apply(tbl_ord_by_asm_coords_NonRef[,"next_peak",drop = F], 1, return_start))

## (assigning into df[which(...),]$col, or one value into a column, fails when there are no rows,
## so these use df$col[which(...)] and rep(value, nrow(df)))
tbl_ord_by_asm_coords_NonRef$same_chr = rep(FALSE, nrow(tbl_ord_by_asm_coords_NonRef))
tbl_ord_by_asm_coords_NonRef$same_chr[which(tbl_ord_by_asm_coords_NonRef$first_chr == tbl_ord_by_asm_coords_NonRef$second_chr)] = TRUE

tbl_ord_by_asm_coords_NonRef_Noncomplex = tbl_ord_by_asm_coords_NonRef[which(tbl_ord_by_asm_coords_NonRef$same_chr == T),]
## AMB if the two ref peaks have the same start, so the direction can't be told
tbl_ord_by_asm_coords_NonRef_Noncomplex$True_Previous_Peak = rep("AMB", nrow(tbl_ord_by_asm_coords_NonRef_Noncomplex))
forward = which(tbl_ord_by_asm_coords_NonRef_Noncomplex$first_start < tbl_ord_by_asm_coords_NonRef_Noncomplex$second_start)
reverse = which(tbl_ord_by_asm_coords_NonRef_Noncomplex$first_start > tbl_ord_by_asm_coords_NonRef_Noncomplex$second_start)
tbl_ord_by_asm_coords_NonRef_Noncomplex$True_Previous_Peak[forward] = tbl_ord_by_asm_coords_NonRef_Noncomplex$prev_peak[forward]
tbl_ord_by_asm_coords_NonRef_Noncomplex$True_Previous_Peak[reverse] = tbl_ord_by_asm_coords_NonRef_Noncomplex$next_peak[reverse]

tbl_ord_by_asm_coords_NonRef_Noncomplex$TrueChr = tbl_ord_by_asm_coords_NonRef_Noncomplex$first_chr
## direction of the contig relative to the reference here, from the two ref peaks around the peak
tbl_ord_by_asm_coords_NonRef_Noncomplex$orientation = rep(NA_character_, nrow(tbl_ord_by_asm_coords_NonRef_Noncomplex))
tbl_ord_by_asm_coords_NonRef_Noncomplex$orientation[forward] = "F"
tbl_ord_by_asm_coords_NonRef_Noncomplex$orientation[reverse] = "R"

### ok, let's deal with contig ends

tbl_ord_by_asm_coords_NonRef_complex = tbl_ord_by_asm_coords_NonRef[which(tbl_ord_by_asm_coords_NonRef$same_chr == F),]

### let's assign contigs to chrs

contigs_to_ref_chrs = tbl_ord_by_asm_coords_Ref %>%
  group_by(asm_chr, cons_chr) %>%
  summarise(n = n_distinct(consensus_peak_id), .groups = "drop")

most_common_ref_chr <- contigs_to_ref_chrs %>%
  group_by(asm_chr) %>%
  slice_max(order_by = n, n = 1, with_ties = FALSE) %>%
  ungroup()

tbl_for_merging = most_common_ref_chr[,c("asm_chr", "cons_chr")]
colnames(tbl_for_merging) = c("asm_chr", "most_likely_chr")

## let's figure out orientation of each contig, from the ref peaks on its most likely chr:
## if cons_start goes down as asm_start goes up (more steps down than up), it is reverse
orientation_tbl_for_merging = tbl_ord_by_asm_coords_Ref %>%
  inner_join(tbl_for_merging, by = "asm_chr") %>%
  filter(cons_chr == most_likely_chr) %>%
  arrange(asm_chr, asm_start) %>%
  group_by(asm_chr) %>%
  summarise(n_up = sum(diff(cons_start) > 0), n_down = sum(diff(cons_start) < 0), .groups = "drop") %>%
  mutate(orientation = ifelse(n_down > n_up, "R", "F"))
print(as.data.frame(orientation_tbl_for_merging))
orientation_tbl_for_merging = orientation_tbl_for_merging[,c("asm_chr", "orientation")]

tbl_ord_by_asm_coords_NonRef_complex = left_join(tbl_ord_by_asm_coords_NonRef_complex, tbl_for_merging, by = "asm_chr")

tbl_ord_by_asm_coords_NonRef_complex_first_chr = tbl_ord_by_asm_coords_NonRef_complex[which(tbl_ord_by_asm_coords_NonRef_complex$most_likely_chr == tbl_ord_by_asm_coords_NonRef_complex$first_chr),]
tbl_ord_by_asm_coords_NonRef_complex_second_chr = tbl_ord_by_asm_coords_NonRef_complex[which(tbl_ord_by_asm_coords_NonRef_complex$most_likely_chr == tbl_ord_by_asm_coords_NonRef_complex$second_chr),]
tbl_ord_by_asm_coords_NonRef_complex_NonMatch_chr = tbl_ord_by_asm_coords_NonRef_complex[which(tbl_ord_by_asm_coords_NonRef_complex$most_likely_chr != tbl_ord_by_asm_coords_NonRef_complex$second_chr & tbl_ord_by_asm_coords_NonRef_complex$most_likely_chr != tbl_ord_by_asm_coords_NonRef_complex$first_chr),]


tbl_ord_by_asm_coords_NonRef_complex_first_chr = left_join(tbl_ord_by_asm_coords_NonRef_complex_first_chr, orientation_tbl_for_merging, by = "asm_chr")

tbl_ord_by_asm_coords_NonRef_complex_second_chr = left_join(tbl_ord_by_asm_coords_NonRef_complex_second_chr, orientation_tbl_for_merging, by = "asm_chr")


## a peak past the end of its contig's ref peaks sits next to the nearest one (the anchor): just
## after it, or just before it if the contig runs the other way. Just before the anchor means just
## after the CHM13 peak that comes before the anchor in the full CHM13 peak list (not just among
## the peaks in this assembly), or AT_START if the anchor is the chromosome's first peak
chm13_peaks = fread(chm13_peaks_file, select = 1:3)
setnames(chm13_peaks, c("chr", "start", "end"))
setorder(chm13_peaks, chr, start, end)
all_ref_ids = paste(chm13_peaks$chr, chm13_peaks$start, chm13_peaks$end, sep = "_")
all_ref_peak_before = c(NA, head(all_ref_ids, -1))
all_ref_peak_before[which(c(NA, head(chm13_peaks$chr, -1)) != chm13_peaks$chr)] = NA
peak_before = function(anchor){
  stopifnot(all(anchor %in% all_ref_ids))
  before = all_ref_peak_before[match(anchor, all_ref_ids)]
  ifelse(is.na(before), "AT_START", before)
}

## but if the anchor is this assembly's first CHM13 peak on its chromosome, the peak could be
## anywhere between the chromosome's start and the anchor, so it goes to AT_START; likewise past
## the assembly's last CHM13 peak on the chromosome it goes to AT_END
asm_ref_ids = all_ref_ids[all_ref_ids %in% tbl_ord_by_ref_coords$consensus_peak_id]
asm_ref_chr = sub("_[^_]*_[^_]*$", "", asm_ref_ids)
asm_first_ref = asm_ref_ids[!duplicated(asm_ref_chr)]
asm_last_ref = asm_ref_ids[!duplicated(asm_ref_chr, fromLast = TRUE)]

## the peak is after the anchor in contig order, so at the contig's end: after the anchor on CHM13
## if the contig runs forward, before it if reversed
place_at_contig_end = function(tbl, anchor, contig_position){
  n = nrow(tbl)
  tbl$True_Previous_Peak = rep("AMB", n)
  tbl$previous_peak_full_list = rep("AMB", n)
  tbl$contig_position = rep(contig_position, n)
  tbl$contig_anchor_peak = anchor
  tbl$open_side = rep(NA_character_, n)
  ## the anchor is before the peak on CHM13 when the contig end and its direction agree
  after = which((contig_position == "contig_end" & tbl$orientation == "F") | (contig_position == "contig_start" & tbl$orientation == "R"))
  before = which((contig_position == "contig_end" & tbl$orientation == "R") | (contig_position == "contig_start" & tbl$orientation == "F"))
  tbl$open_side[after] = "next_unknown"
  tbl$previous_peak_full_list[after] = anchor[after]
  tbl$True_Previous_Peak[after] = ifelse(anchor[after] %in% asm_last_ref, "AT_END", anchor[after])
  tbl$open_side[before] = "previous_unknown"
  tbl$previous_peak_full_list[before] = peak_before(anchor[before])
  tbl$True_Previous_Peak[before] = ifelse(anchor[before] %in% asm_first_ref, "AT_START", peak_before(anchor[before]))
  tbl
}

tbl_ord_by_asm_coords_NonRef_complex_first_chr = place_at_contig_end(tbl_ord_by_asm_coords_NonRef_complex_first_chr,
                                                                     tbl_ord_by_asm_coords_NonRef_complex_first_chr$prev_peak, "contig_end")


## the peak is before the anchor in contig order, so at the contig's start
tbl_ord_by_asm_coords_NonRef_complex_second_chr = place_at_contig_end(tbl_ord_by_asm_coords_NonRef_complex_second_chr,
                                                                      tbl_ord_by_asm_coords_NonRef_complex_second_chr$next_peak, "contig_start")

if(nrow(tbl_ord_by_asm_coords_NonRef_complex_NonMatch_chr) > 0){
  tbl_ord_by_asm_coords_NonRef_complex_NonMatch_chr$most_likely_chr = "AMB"
  tbl_ord_by_asm_coords_NonRef_complex_NonMatch_chr$True_Previous_Peak = "AMB"
  tbl_ord_by_asm_coords_NonRef_complex_NonMatch_chr = left_join(tbl_ord_by_asm_coords_NonRef_complex_NonMatch_chr, orientation_tbl_for_merging, by = "asm_chr")
  complex_amb = tbl_ord_by_asm_coords_NonRef_complex_NonMatch_chr[,c("consensus_peak_id", "most_likely_chr", "True_Previous_Peak", "orientation")]
  colnames(complex_amb) = c("consensus_peak_id", "chr", "previous_peak", "orientation")
  complex_amb$previous_peak_full_list = "AMB"
  complex_amb[, c("contig_position", "open_side", "contig_anchor_peak")] = NA_character_

}

### combine_all_outputs

## orientation: F if the contig runs the same direction as the reference here, R if reversed
## (for peaks at contig ends, the direction of the whole contig on its most likely chr)
## previous_peak_full_list: as previous_peak, but a peak at a contig end is placed beside its anchor
## even at the edges of the assembly's chromosome (where previous_peak is AT_START or AT_END).
## contig_position: contig_start or contig_end for peaks past the last CHM13 peak on their contig,
## with contig_anchor_peak the CHM13 peak they were placed by and open_side the side of the peak
## that isn't known on CHM13 (previous_unknown or next_unknown)
out_cols = c("consensus_peak_id", "chr", "previous_peak", "orientation", "previous_peak_full_list", "contig_position", "open_side", "contig_anchor_peak")
tbl_ord_by_asm_coords_NonRef_Noncomplex$previous_peak_full_list = tbl_ord_by_asm_coords_NonRef_Noncomplex$True_Previous_Peak
tbl_ord_by_asm_coords_NonRef_Noncomplex[, c("contig_position", "open_side", "contig_anchor_peak")] = NA_character_
non_complex = tbl_ord_by_asm_coords_NonRef_Noncomplex[,c("consensus_peak_id", "TrueChr","True_Previous_Peak", "orientation", out_cols[5:8])]
complex_first = tbl_ord_by_asm_coords_NonRef_complex_first_chr[,c("consensus_peak_id", "most_likely_chr", "True_Previous_Peak", "orientation", out_cols[5:8])]
complex_second = tbl_ord_by_asm_coords_NonRef_complex_second_chr[,c("consensus_peak_id", "most_likely_chr", "True_Previous_Peak", "orientation", out_cols[5:8])]

colnames(non_complex) = out_cols
colnames(complex_first) = out_cols
colnames(complex_second) = out_cols

combined_table = as.data.frame(rbind(non_complex, complex_first, complex_second))

if(nrow(tbl_ord_by_asm_coords_NonRef_complex_NonMatch_chr) > 0){
  combined_table = as.data.frame(rbind(combined_table, complex_amb))
}

## anything not placed above (e.g. peaks on contigs with no ref peaks) has no known chromosome
unplaced = setdiff(unique(tbl_ord_by_asm_coords_NonRef$consensus_peak_id), combined_table$consensus_peak_id)
if(length(unplaced) > 0){
  print(paste0("peaks with no chromosome (AMB): ", length(unplaced)))
  combined_table = as.data.frame(rbind(combined_table,
                                       data.frame(consensus_peak_id = unplaced, chr = "AMB", previous_peak = "AMB", orientation = NA_character_,
                                                  previous_peak_full_list = "AMB", contig_position = NA_character_, open_side = NA_character_,
                                                  contig_anchor_peak = NA_character_)))
}

write.table(combined_table, output_file, col.names = TRUE, row.names = FALSE, quote = F, sep = "\t")
