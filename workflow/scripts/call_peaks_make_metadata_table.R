library(ggplot2)
library(dplyr)
library(tidyr)

return_chr = function(a){
  num_elements = length(strsplit(a, "_")[[1]])
  return(paste(strsplit(a, "_")[[1]][1:(num_elements-2)], collapse = "_"))
}

return_start = function(a){
  num_elements = length(strsplit(a, "_")[[1]])
  return(strsplit(a, "_")[[1]][num_elements-1])
}

return_end = function(a){
  num_elements = length(strsplit(a, "_")[[1]])
  return(strsplit(a, "_")[[1]][num_elements])
}

return_surj_peak_length = function(a){
  splt=strsplit(a,"_")[[1]]
  return(as.numeric(splt[length(splt)])-as.numeric(splt[length(splt)-1]))
}

return_surj_chr = function(a){
  splt=strsplit(a,"-")[[1]][length(strsplit(a,"-")[[1]])]
  split2=paste(strsplit(splt, "_")[[1]][1:(length(strsplit(splt, "_")[[1]])-2)], collapse = "_")
  return(split2)
}

return_cons_peak = function(a){
  splt=strsplit(a,"-")[[1]][length(strsplit(a,"-")[[1]])]
  return(splt)
}

return_chr_w_hash = function(a){
  splt=strsplit(a, "#")[[1]]
  return(paste0(splt[1], "#", splt[2], "#"))
}

return_samp_hap_id = function(a, conversion_table){
  temp_sample = strsplit(a[1], "-")[[1]][1]
  for(i in 1:nrow(conversion_table)){
    temp_asm = conversion_table[i,"asm"]
    if(grepl(temp_asm, a[2])){
      temp_garph_asm=conversion_table[i,"graph_asm"]
      return(paste0(temp_sample, "_", temp_garph_asm))
    }
  }
}

return_samp_hap_id_2 = function(a, conversion_table){
  temp_sample = a[1]
  if(grepl("#", a[2])){
    temp_split=strsplit(a[2],"#")[[1]]
    temp_chr = paste0(temp_split[1], "#", temp_split[2], "#")
    temp_garph_asm=conversion_table[which(conversion_table$asm == temp_chr),"graph_asm"]
    return(paste0(temp_sample, "_", temp_garph_asm))
  }
  else {
    sub_conversion_table = conversion_table[which(!(grepl("#",conversion_table$asm))),]
    for(i in 1:nrow(sub_conversion_table)){
      temp_asm = sub_conversion_table[i,"asm"]
      if(grepl(temp_asm, a[2])){
        temp_garph_asm=sub_conversion_table[i,"graph_asm"]
        return(paste0(temp_sample, "_", temp_garph_asm))
      }
    }
  }
}

args <- commandArgs(trailingOnly = TRUE)

input_dir=args[1]
output_dir=input_dir
prefix=args[2]
sample_pefix_file_name=args[3]

tbl = read.table(paste0(input_dir, prefix, "_max_score_per_peak_all_samples.bed"), header = T, sep = "\t", stringsAsFactors = F, comment.char = "")
trans_peaks = read.table(paste0(input_dir, prefix, "_combined_tranferred_peaks_all_samples.bed"), header = T, sep = "\t", stringsAsFactors = F, comment.char = "")
sample_prefix_file = read.table(sample_pefix_file_name, sep = "\t", header = F, stringsAsFactors = F, comment.char = "")
colnames(sample_prefix_file) = c("asm", "graph_asm", "path_to_fasta")

hg38_file=read.table(paste0(input_dir, prefix, "_All_consensus_peaks_moved_to_GRCh38.bed"), sep="\t", header=F, stringsAsFactors = F, comment.char = "")
colnames(hg38_file) = c("hg38_chr", "hg38_start", "hg38_end", "consensus_peak_id")
orig_peak_trans_peak_intersect=read.table(paste0(input_dir, prefix, "_orig_called_peaks_intersect_transferred_consensus_peaks.bed"), sep="\t", header=F, stringsAsFactors = F, comment.char = "")
colnames(orig_peak_trans_peak_intersect)=c("orig_chr", "orig_start", "orig_end", "sample_id", "orig_peak_id", "trans_chr", "trans_start", "trans_end", "cons_peak_id")
orig_peak_trans_peak_intersect$trans_peak_id = paste0(orig_peak_trans_peak_intersect$trans_chr, "_", orig_peak_trans_peak_intersect$trans_start, "_",orig_peak_trans_peak_intersect$trans_end)

surj_peak_cons_peak_intersect=read.table(paste0(input_dir,prefix, "_all_orig_surjected_peaks_intersect_consensus_peaks.bed"), sep="\t", header=F, stringsAsFactors = F, comment.char = "")
colnames(surj_peak_cons_peak_intersect) = c("chr_surj", "start_surj", "end_surj", "sample_chr_id", "peak_id_orig", "chr_cons", "start_cons", "end_cons", "peak_id_cons" )

peak_sequence_file = read.table(paste0(input_dir,prefix,"_transferred_cons_peaks_native_sequences_from_main_output.bed"), sep="\t", header=F, stringsAsFactors = F, comment.char = "")
colnames(peak_sequence_file) = c("asm_chr", "asm_start", "asm_end", "cons_peak_id", "seq")
peak_sequence_file$asm_peak_id = paste0(peak_sequence_file$asm_chr, "_", peak_sequence_file$asm_start, "_", peak_sequence_file$asm_end)


if(1 == 0){
  input_dir="/mmfs1/gscratch/stergachislab/HPRC/fiber-seq-pilot/HLA_graphs/peak_calling_all_chrs/All_HPRC_peaks_NewGraph_w_HG002/"
  
  output_dir="/mmfs1/gscratch/stergachislab/HPRC/fiber-seq-pilot/HLA_graphs/peak_calling_all_chrs/All_HPRC_peaks_NewGraph_w_HG002/"
  prefix="All_HPRC_peaks_NewGraph_w_HG002"
  tbl = read.table(paste0(input_dir,"All_HPRC_peaks_NewGraph_w_HG002_max_score_per_peak_all_samples.bed"), header = T, sep = "\t", stringsAsFactors = F, comment.char = "")
  trans_peaks = read.table(paste0(input_dir,"All_HPRC_peaks_NewGraph_w_HG002_combined_tranferred_peaks_all_samples.bed"), header = T, sep = "\t", stringsAsFactors = F, comment.char = "")
  sample_prefix_file = read.table("/mmfs1/gscratch/stergachislab/HPRC/fiber-seq-pilot/HLA_graphs/alg_dev/sm_input_file.txt", sep = "\t", header = F, stringsAsFactors = F, comment.char = "")
  colnames(sample_prefix_file) = c("asm", "graph_asm", "path_to_fasta")
  
  hg38_file=read.table(paste0(input_dir,"All_HPRC_peaks_NewGraph_w_HG002_All_consensus_peaks_moved_to_GRCh38.bed"), sep="\t", header=F, stringsAsFactors = F, comment.char = "")
  colnames(hg38_file) = c("hg38_chr", "hg38_start", "hg38_end", "consensus_peak_id")
  orig_peak_trans_peak_intersect=read.table(paste0(input_dir, "All_HPRC_peaks_NewGraph_w_HG002_orig_called_peaks_intersect_transferred_consensus_peaks.bed"), sep="\t", header=F, stringsAsFactors = F, comment.char = "")
  colnames(orig_peak_trans_peak_intersect)=c("orig_chr", "orig_start", "orig_end", "sample_id", "orig_peak_id", "trans_chr", "trans_start", "trans_end", "cons_peak_id")
  orig_peak_trans_peak_intersect$trans_peak_id = paste0(orig_peak_trans_peak_intersect$trans_chr, "_", orig_peak_trans_peak_intersect$trans_start, "_",orig_peak_trans_peak_intersect$trans_end)
  
  surj_peak_cons_peak_intersect=read.table(paste0(input_dir,"All_HPRC_peaks_NewGraph_w_HG002_all_orig_surjected_peaks_intersect_consensus_peaks.bed"), sep="\t", header=F, stringsAsFactors = F, comment.char = "")
  colnames(surj_peak_cons_peak_intersect) = c("chr_surj", "start_surj", "end_surj", "sample_chr_id", "peak_id_orig", "chr_cons", "start_cons", "end_cons", "peak_id_cons" )
  
  peak_sequence_file = read.table(paste0(input_dir,"All_HPRC_peaks_NewGraph_w_HG002_transferred_cons_peaks_native_sequences_from_main_output.bed", sep="\t", header=F, stringsAsFactors = F, comment.char = ""))
  colnames(peak_sequence_file) = c("asm_chr", "asm_start", "asm_end", "cons_peak_id", "seq")
  peak_sequence_file$asm_peak_id = paste0(peak_sequence_file$asm_chr, "_", peak_sequence_file$asm_start, "_", peak_sequence_file$asm_end)
}

tbl$cons_chr = apply(tbl[,"consensus_peak_id",drop=F], 1, return_chr)
tbl$cons_start = as.numeric(apply(tbl[,"consensus_peak_id",drop=F], 1, return_start))
tbl$cons_end = as.numeric(apply(tbl[,"consensus_peak_id",drop=F], 1, return_end))

tbl$fire_cov_OV_cov = tbl$fire_coverage/tbl$coverage

return_asm=function(samp_id, samp_tbl){
  temp_split = strsplit(samp_id, "_")[[1]]
  temp_graph_asm = paste0(temp_split[length(temp_split)-1], "_", temp_split[length(temp_split)])
  return(samp_tbl[which(samp_tbl$graph_asm == temp_graph_asm),]$asm)
}

tbl$asm = apply(tbl[,"sample_id",drop=F],1,return_asm, samp_tbl = sample_prefix_file)

unique_sample_id = unique(tbl$sample_id)

updated_tbl = data.frame()

for(i in 1:length(unique_sample_id)){
  temp_sample_id = unique_sample_id[i]
  print(temp_sample_id)
  tbl_sub = tbl[which(tbl$sample_id == temp_sample_id),]
  temp_asm = unique(tbl_sub$asm)
  if(length(temp_asm) > 1){
    print("Move than one assembly!!")
  }
  temp_trans_peaks = trans_peaks[which(grepl(temp_asm, trans_peaks$chr)),]
  
  temp_tbl_sub_NA = tbl_sub[which(tbl_sub$score == -2),]
  temp_tbl_sub_NotNA = tbl_sub[which(tbl_sub$score != -2),]
  if(nrow(temp_trans_peaks) != nrow(temp_tbl_sub_NotNA)){
    print("Row numbers not matching")
  }
  temp_trans_peaks_ord = temp_trans_peaks[match(temp_tbl_sub_NotNA$consensus_peak_id,temp_trans_peaks$cons_peak_id),]
  temp_tbl_sub_NotNA$asm_chr = temp_trans_peaks_ord$chr
  temp_tbl_sub_NotNA$asm_start = temp_trans_peaks_ord$start
  temp_tbl_sub_NotNA$asm_end = temp_trans_peaks_ord$end
  temp_tbl_sub_NA$asm_chr = NA
  temp_tbl_sub_NA$asm_start = NA
  temp_tbl_sub_NA$asm_end = NA
  updated_tbl = as.data.frame(rbind(updated_tbl, temp_tbl_sub_NotNA, temp_tbl_sub_NA))
}

updated_tbl$asm_peak_id = paste0(updated_tbl$asm_chr, "_", updated_tbl$asm_start, "_", updated_tbl$asm_end)

### add hg38 info

### THIS IS A TEST FILE FROM A DIFFERENT DATASET
#hg38_file=read.table("/Users/minkinaa/Documents/StergachisLab/projects/HLA_peak_assignment/test_files/All_HPRC_peaks_NewGraph_w_HG002_All_consensus_peaks_moved_to_GRCh38.bed", sep="\t", header=F, stringsAsFactors = F, comment.char = "")
#colnames(hg38_file) = c("hg38_chr", "hg38_start", "hg38_end", "consensus_peak_id")

hg38_file$cons_chr = apply(hg38_file[,"consensus_peak_id",drop=F], 1, return_chr)
hg38_file$cons_start = as.numeric(apply(hg38_file[,"consensus_peak_id",drop=F], 1, return_start))
hg38_file$cons_end = as.numeric(apply(hg38_file[,"consensus_peak_id",drop=F], 1, return_end))

hg38_file$cons_peak_length = hg38_file$cons_end - hg38_file$cons_start
hg38_file$hg38_peak_length = hg38_file$hg38_end - hg38_file$hg38_start

hg38_file$length_diff = abs(hg38_file$hg38_peak_length-hg38_file$cons_peak_length)
hg38_file$percent_diff = hg38_file$length_diff/hg38_file$cons_peak_length*100

hg38_file_pass_filt = hg38_file[which(hg38_file$percent_diff < 20),]
hg38_file_pass_filt_for_merge = hg38_file_pass_filt[,c(1:4)]
updated_tbl <- left_join(updated_tbl, hg38_file_pass_filt_for_merge, by = "consensus_peak_id")
updated_tbl$exists_in_hg38 = FALSE
updated_tbl[which(!is.na(updated_tbl$hg38_chr)),]$exists_in_hg38 = TRUE

## add info about whether original peak existed
updated_tbl$overlaps_original_called_peak=FALSE
updated_tbl[which(updated_tbl$asm_peak_id %in% orig_peak_trans_peak_intersect$trans_peak_id),]$overlaps_original_called_peak=TRUE

## add multiple to one info

surj_peak_cons_peak_intersect$chr_orig = apply(surj_peak_cons_peak_intersect[,"peak_id_orig",drop = F], 1, return_chr)
surj_peak_cons_peak_intersect$start_orig = apply(surj_peak_cons_peak_intersect[,"peak_id_orig",drop = F], 1, return_start)
surj_peak_cons_peak_intersect$end_orig = apply(surj_peak_cons_peak_intersect[,"peak_id_orig",drop = F], 1, return_end)

surj_peak_cons_peak_intersect$samp_chr_consPeak = paste0(surj_peak_cons_peak_intersect$sample_chr_id, "-", surj_peak_cons_peak_intersect$chr_orig, "-", surj_peak_cons_peak_intersect$peak_id_cons)

surj_peak_cons_peak_intersect_wHash = surj_peak_cons_peak_intersect[which(grepl("#",surj_peak_cons_peak_intersect$chr_orig)),]
surj_peak_cons_peak_intersect_NoHash = surj_peak_cons_peak_intersect[which(!(grepl("#",surj_peak_cons_peak_intersect$chr_orig))),]

if(nrow(surj_peak_cons_peak_intersect_wHash) > 0){
  surj_peak_cons_peak_intersect_wHash$asm = apply(surj_peak_cons_peak_intersect_wHash[,"chr_orig", drop = FALSE], 1, return_chr_w_hash)
}

if(nrow(surj_peak_cons_peak_intersect_NoHash) > 0){
  nonHash_asms = sample_prefix_file[which(!grepl("#",sample_prefix_file$asm)),]$asm
  
  surj_peak_cons_peak_intersect_NoHash$asm = "NONE"
  
  for(i in 1:length(nonHash_asms)){
    temp_asm = nonHash_asms[i]
    print(temp_asm)
    surj_peak_cons_peak_intersect_NoHash[which(grepl(temp_asm,surj_peak_cons_peak_intersect_NoHash$chr_orig)),]$asm = temp_asm
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

temp = as.data.frame(table(surj_peak_cons_peak_intersect$sample_graphHap_peak))
temp_filt = temp[which(temp$Freq > 1),]
temp_filt_2 = temp[which(temp$Freq == 2),]
temp_filt_gr2 = temp[which(temp$Freq >2),]

surj_peak_cons_peak_intersect_mult_to_one = surj_peak_cons_peak_intersect[which(surj_peak_cons_peak_intersect$sample_graphHap_peak %in% temp_filt$Var1),]
surj_peak_cons_peak_intersect_mult_to_one = surj_peak_cons_peak_intersect_mult_to_one[order(as.numeric(surj_peak_cons_peak_intersect_mult_to_one$start_orig)),]
surj_peak_cons_peak_intersect_mult_to_one = surj_peak_cons_peak_intersect_mult_to_one[order(surj_peak_cons_peak_intersect_mult_to_one$sample_graphHap_peak),]

surj_peak_cons_peak_intersect_mult_to_one_2 = surj_peak_cons_peak_intersect_mult_to_one[which(surj_peak_cons_peak_intersect_mult_to_one$sample_graphHap_peak %in% temp_filt_2$Var1),]
surj_peak_cons_peak_intersect_mult_to_one_gr2 = surj_peak_cons_peak_intersect_mult_to_one[which(surj_peak_cons_peak_intersect_mult_to_one$sample_graphHap_peak %in% temp_filt_gr2$Var1),]

aggreg_by_cons_peak_2 = as.data.frame(surj_peak_cons_peak_intersect_mult_to_one_2 %>%
                                        group_by(sample_graphHap_peak) %>%
                                        summarise(vals = list(c(rbind(peak_id_orig, chr_orig, start_orig, end_orig))), .groups="drop") %>%
                                        unnest_wider(vals, names_sep = ""))

if(nrow(aggreg_by_cons_peak_2) > 0){
  colnames(aggreg_by_cons_peak_2) = c("sample_graphHap_peak", "peak_id_1", "orig_chr_1","orig_start_1", "orig_end_1", "peak_id_2", "orig_chr_2","orig_start_2", "orig_end_2")
  cols_to_convert <- c("orig_start_1", "orig_end_1", "orig_start_2", "orig_end_2")
  aggreg_by_cons_peak_2[cols_to_convert] <- lapply(aggreg_by_cons_peak_2[cols_to_convert], as.numeric)
  aggreg_by_cons_peak_2$start_diff =aggreg_by_cons_peak_2$orig_start_2-aggreg_by_cons_peak_2$orig_start_1
  aggreg_by_cons_peak_2$space_between_peaks = aggreg_by_cons_peak_2$orig_start_2-aggreg_by_cons_peak_2$orig_end_1
  aggreg_by_cons_peak_2$share_orig_chr = FALSE
  aggreg_by_cons_peak_2[which(aggreg_by_cons_peak_2$orig_chr_1 == aggreg_by_cons_peak_2$orig_chr_2),]$share_orig_chr = TRUE
  
  aggreg_by_cons_peak_2$surj_peak_length = apply(aggreg_by_cons_peak_2[,"sample_graphHap_peak",drop = F],1,return_surj_peak_length)
  
  aggreg_by_cons_peak_2$dist_less_than_cons_peak = FALSE
  aggreg_by_cons_peak_2[which(aggreg_by_cons_peak_2$space_between_peaks < aggreg_by_cons_peak_2$surj_peak_length),]$dist_less_than_cons_peak = TRUE
  
  aggreg_by_cons_peak_2_TRUE_MULT_TO_ONE = aggreg_by_cons_peak_2[which(aggreg_by_cons_peak_2$dist_less_than_cons_peak == FALSE | aggreg_by_cons_peak_2$share_orig_chr == FALSE),]
  aggreg_by_cons_peak_2_TRUE_MULT_TO_ONE$surj_chr = apply(aggreg_by_cons_peak_2_TRUE_MULT_TO_ONE[,"sample_graphHap_peak",drop = F], 1, return_surj_chr)
  aggreg_by_cons_peak_2_TRUE_MULT_TO_ONE$cons_peak = apply(aggreg_by_cons_peak_2_TRUE_MULT_TO_ONE[,"sample_graphHap_peak",drop = F], 1, return_cons_peak)
  
}

updated_tbl$sample_id_peak_id_TEMP_COL = paste0(updated_tbl$sample_id, "-", updated_tbl$consensus_peak_id)

### greater than 2:

if(nrow(aggreg_by_cons_peak_2) > 0){
  surj_peak_cons_peak_intersect_mult_to_one_gr2$cons_peak_length = surj_peak_cons_peak_intersect_mult_to_one_gr2$end_cons-surj_peak_cons_peak_intersect_mult_to_one_gr2$start_cons
  
  unique_sample_peaks = unique(surj_peak_cons_peak_intersect_mult_to_one_gr2$sample_graphHap_peak)
  
  num_peaks_col = c()
  peak_names = c()
  
  for(i in 1:length(unique_sample_peaks)){
    if(i%%100 == 0){
      print(i)
    }
    temp_samp_peak = unique_sample_peaks[i]
    temp_tbl = surj_peak_cons_peak_intersect_mult_to_one_gr2[which(surj_peak_cons_peak_intersect_mult_to_one_gr2$sample_graphHap_peak == temp_samp_peak),]
    temp_tbl_ord = temp_tbl[order(temp_tbl$start_orig),]
    cons_len = as.numeric(temp_tbl_ord[1,"cons_peak_length"])
    num_peaks = 1
    for(i in 2:nrow(temp_tbl)){
      if(as.numeric(temp_tbl[i,"start_orig"])-as.numeric(temp_tbl[i-1,"end_orig"]) > cons_len){
        num_peaks = num_peaks + 1
      }
    }
    all_peaks = paste(temp_tbl_ord$peak_id_orig, collapse = ",")
    peak_names = c(peak_names, all_peaks)
    num_peaks_col = c(num_peaks_col, num_peaks)
  }
  
  peak_id_and_num_peaks = as.data.frame(cbind(unique_sample_peaks, num_peaks_col, peak_names))
  colnames(peak_id_and_num_peaks) = c("sample_graphHap_peak", "num_ov_peaks", "all_overlapping_peaks")
} else {
  peak_id_and_num_peaks <- data.frame(
    sample_graphHap_peak = NULL,
    num_ov_peaks = NULL,
    all_overlapping_peaks = NULL
  )
}


all_peaks_summary = surj_peak_cons_peak_intersect %>%
  group_by(sample_graphHap_peak) %>%
  summarise(all_orig_overlapping_peaks = paste(peak_id_orig, collapse = ","))
all_peaks_summary = as.data.frame(all_peaks_summary)

updated_tbl$num_overlapping_orig_peaks = 0
updated_tbl[which(updated_tbl$sample_id_peak_id_TEMP_COL %in% surj_peak_cons_peak_intersect$sample_graphHap_peak),]$num_overlapping_orig_peaks = 1

if(nrow(aggreg_by_cons_peak_2) > 0){
  updated_tbl[which(updated_tbl$sample_id_peak_id_TEMP_COL %in% aggreg_by_cons_peak_2_TRUE_MULT_TO_ONE$sample_graphHap_peak),]$num_overlapping_orig_peaks = 2
}

updated_tbl_2_or_less = updated_tbl[which(!(updated_tbl$sample_id_peak_id_TEMP_COL %in% peak_id_and_num_peaks$sample_graphHap_peak)),]

if(nrow(aggreg_by_cons_peak_2) > 0){
  updated_tbl_gr_2 = updated_tbl[which(updated_tbl$sample_id_peak_id_TEMP_COL %in% peak_id_and_num_peaks$sample_graphHap_peak),]
  peak_id_and_num_peaks_ord = peak_id_and_num_peaks[match(updated_tbl_gr_2$sample_id_peak_id_TEMP_COL,peak_id_and_num_peaks$sample_graphHap_peak),]
  updated_tbl_gr_2$num_overlapping_orig_peaks = peak_id_and_num_peaks_ord$num_ov_peaks
  updated_tbl = as.data.frame(rbind(updated_tbl_2_or_less, updated_tbl_gr_2))
} else {
  updated_tbl = updated_tbl_2_or_less
}

#updated_tbl = as.data.frame(rbind(updated_tbl_2_or_less, updated_tbl_gr_2))

updated_tbl_non0 = updated_tbl[which(updated_tbl$num_overlapping_orig_peaks > 0),]
updated_tbl_0 = updated_tbl[which(updated_tbl$num_overlapping_orig_peaks == 0),]

nrow(updated_tbl_non0) == nrow(all_peaks_summary)

all_peaks_summary_ord = all_peaks_summary[match(updated_tbl_non0$sample_id_peak_id_TEMP_COL,all_peaks_summary$sample_graphHap_peak),]

updated_tbl_non0$orig_peaks_overlapping_consensus = all_peaks_summary_ord$all_orig_overlapping_peaks
updated_tbl_0$orig_peaks_overlapping_consensus = NA

updated_tbl = as.data.frame(rbind(updated_tbl_non0,updated_tbl_0))

#write.table(updated_tbl, paste0(output_dir, prefix, "_final_tbl_w_metadata.tsv"), col.names = TRUE, row.names = FALSE, quote = F, sep = "\t")
#updated_tbl = read.table("/mmfs1/gscratch/stergachislab/HPRC/fiber-seq-pilot/HLA_graphs/peak_calling_all_chrs/All_HPRC_peaks_NewGraph_w_HG002/All_HPRC_peaks_NewGraph_w_HG002_final_tbl_w_metadata.tsv", header = T, sep = "\t", stringsAsFactors = F, comment.char = "")
#sum(peak_sequence_file$asm_peak_id %in% updated_tbl$asm_peak_id)
#nrow(peak_sequence_file)

peak_sequence_file_2col = unique(peak_sequence_file[,c("asm_peak_id", "seq")])

updated_tbl = left_join(updated_tbl,peak_sequence_file_2col, by = "asm_peak_id")

updated_tbl <- updated_tbl %>%
  dplyr::select(-sample_id_peak_id_TEMP_COL)

write.table(updated_tbl, paste0(output_dir, prefix, "_final_tbl_w_metadata.tsv"), col.names = TRUE, row.names = FALSE, quote = F, sep = "\t")

