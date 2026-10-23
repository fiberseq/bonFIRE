#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

orig_bed_file=args[1]
surj_bed_file=args[2]
last_surj_asm=args[3]
output_all = args[4]
#output_surj = args[5]
#output_NOTsurj = args[6]
next_asm_file=args[5]

length_change_filter=as.numeric(args[6])
# graph names already surjected onto (one per line), and the asm -> graph name table
tried_asms_file=args[7]
asm_conversion_file=args[8]


#orig_bed_file = "/Users/minkinaa/Documents/StergachisLab/projects/HLA_peak_assignment/test_files/HLA-DQB1-10kb_wHapsAndPeakIDs_w_surject_info_NOT_SURJECTED_ONTO_CHM13.bed"
#surj_bed_file = "/Users/minkinaa/Documents/StergachisLab/projects/HLA_peak_assignment/test_files/HLA-DQB1-10kb_NonCHM13_peaks_fibertigs_surjected_on_newASM_processed.bed"
#last_surj_asm="HG02486#2#"

if(1 == 0){
  orig_bed_file = "/Users/minkinaa/Documents/StergachisLab/projects/HLA_peak_assignment/test_files/All_HPRC_peaks_NewGraph_wHapsAndPeakIDs_w_surject_info_NOT_SURJECTED_ONTO_CHM13.bed"
  surj_bed_file = "/Users/minkinaa/Documents/StergachisLab/projects/HLA_peak_assignment/test_files/All_HPRC_peaks_NewGraph_SurjLoop_TEMP_surjected_processed.bed"
  last_surj_asm= "HG02257#1#"
}

orig_bed = read.table(orig_bed_file, header = T, sep = "\t", stringsAsFactors = F, comment.char = "")
#colnames(orig_bed) = gsub("X.", "", colnames(orig_bed))
colnames(orig_bed)[1] <- "chrom"
# the surjected bed is empty when no peaks surjected onto this assembly
if(file.size(surj_bed_file) > 0){
  surj_bed = read.table(surj_bed_file, header = F, sep = "\t", stringsAsFactors = F, comment.char = "")
} else {
  surj_bed = orig_bed[0,]
}

colnames(surj_bed) = colnames(orig_bed)

orig_bed_previously_surjected = orig_bed[which(orig_bed$surjected_asm != "NONE"),]

orig_bed = orig_bed[which(orig_bed$surjected_asm == "NONE"),]
surj_bed = surj_bed[which(surj_bed$peak_id %in% orig_bed$peak_id),]

orig_bed_surj = orig_bed[which(orig_bed$peak_id %in% surj_bed$peak_id),]
orig_bed_NOT_surj = orig_bed[which(!(orig_bed$peak_id %in% surj_bed$peak_id)),]

surj_bed_ord = surj_bed[match(orig_bed_surj$peak_id,surj_bed$peak_id),]

orig_bed_surj$surjected_asm = rep(last_surj_asm, nrow(orig_bed_surj))
orig_bed_surj$surj_chr = surj_bed_ord$chrom
orig_bed_surj$surj_start = surj_bed_ord$start
orig_bed_surj$surj_end = surj_bed_ord$end

### new stuff
orig_bed_surj$surj_length = abs(orig_bed_surj$surj_end-orig_bed_surj$surj_start)
orig_bed_surj$orig_length = abs(orig_bed_surj$end-orig_bed_surj$start)
orig_bed_surj$abs_length_diff = abs(orig_bed_surj$surj_length-orig_bed_surj$orig_length)/orig_bed_surj$orig_length

if(nrow(orig_bed_surj[which(orig_bed_surj$abs_length_diff > length_change_filter),]) > 0){
  orig_bed_surj[which(orig_bed_surj$abs_length_diff > length_change_filter),]$surjected_asm = "NONE"
  orig_bed_surj[which(orig_bed_surj$abs_length_diff > length_change_filter),]$surj_chr = "NONE"
  orig_bed_surj[which(orig_bed_surj$abs_length_diff > length_change_filter),]$surj_start = "NONE"
  orig_bed_surj[which(orig_bed_surj$abs_length_diff > length_change_filter),]$surj_end = "NONE"
}

orig_bed_surj[,c("surj_length", "orig_length", "abs_length_diff")] <- list(NULL)
###

all = as.data.frame(rbind(orig_bed_previously_surjected, orig_bed_surj, orig_bed_NOT_surj))

write.table(all, output_all, sep = "\t", col.names = T, row.names = F, quote=F)

## this line is new
orig_bed_NOT_surj = all[which(all$surj_chr == "NONE"),]

## pick the assembly with the most unsurjected peaks, skipping ones already tried
tried_asms = readLines(tried_asms_file)
asm_input = read.table(asm_conversion_file, header = F, sep = "\t", stringsAsFactors = F, comment.char = "")
asm_to_graph = setNames(asm_input[,2], asm_input[,1])
candidate_asms = orig_bed_NOT_surj$asm[!(asm_to_graph[orig_bed_NOT_surj$asm] %in% c(tried_asms, NA))]

if(length(candidate_asms) > 0){
  next_asm_tbl = as.data.frame(table(candidate_asms))
  next_asm_tbl_ord = next_asm_tbl[order(-next_asm_tbl$Freq),]
  asm_to_return = as.character(next_asm_tbl_ord[1,1])
  writeLines(asm_to_return, next_asm_file)
} else {
  writeLines("NONE", next_asm_file)
}

