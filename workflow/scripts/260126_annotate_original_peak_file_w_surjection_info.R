#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

return_sample_name = function(a, asm_tbl){
  if(grepl("#", a)){
    return(paste0(strsplit(a, "#")[[1]][1], "#", strsplit(a, "#")[[1]][2], "#"))
  } else
    for(i in 1:nrow(asm_tbl)){
      temp_asm = asm_tbl[i,1]
      if(grepl(temp_asm, a)){
        return(temp_asm)
      }
    }
    #return(strsplit(a, "_")[[1]][2])
}

#orig_bed_file="/Users/minkinaa/Documents/StergachisLab/projects/HLA_peak_assignment/test_files/HLA-DQB1-10kb_wHapsAndPeakIDs.bed"
#surj_bed_file="/Users/minkinaa/Documents/StergachisLab/projects/HLA_peak_assignment/test_files/HLA-DQB1-10kb_all_original_peaks_fibertigs_surjected_on_CHM13_processed.bed"
#last_surj_asm="CHM13"
#output_dir="/Users/minkinaa/Documents/StergachisLab/projects/HLA_peak_assignment/test_files/"
#prefix="HLA-DQB1-10kb"
#output_all=paste0(output_dir, prefix, "_surjected_ALL_w_annotation.tsv")
#output_surj=paste0(output_dir, prefix, "_surjected_onto_CHM13_w_annotation.tsv")
#output_NOTsurj=paste0(output_dir, prefix, "_NOT_surjected_onto_CHM13_w_annotation.tsv")

orig_bed_file=args[1]
surj_bed_file=args[2]
last_surj_asm=args[3]
output_all = args[4]
output_surj = args[5]
output_NOTsurj = args[6]
next_asm_file=args[7]
asm_conversion_file=args[8]

#length_change_filter=.2
length_change_filter=args[9]

asm_input = read.table(asm_conversion_file, header = F, sep = "\t", stringsAsFactors = F, comment.char = "")
colnames(asm_input) = c("asm_orig_substring", "asm_graph", "fasta_input")
asm_input_sub = asm_input[!grepl("#", asm_input$asm_orig_substring),]

orig_bed = read.table(orig_bed_file, header = T, sep = "\t", stringsAsFactors = F, comment.char = "")
surj_bed = read.table(surj_bed_file, header = F, sep = "\t", stringsAsFactors = F, comment.char = "")
colnames(orig_bed) = gsub("X.", "", colnames(orig_bed))

colnames(surj_bed) = colnames(orig_bed)

orig_bed$asm=apply(orig_bed[,"chrom", drop = F], 1, return_sample_name, asm_tbl = asm_input_sub)

orig_bed$surjected_asm = "NONE"
orig_bed$surj_chr = "NONE"
orig_bed$surj_start = "NONE"
orig_bed$surj_end = "NONE"

orig_bed_surj = orig_bed[which(orig_bed$peak_id %in% surj_bed$peak_id),]
orig_bed_NOT_surj = orig_bed[which(!(orig_bed$peak_id %in% surj_bed$peak_id)),]

surj_bed_ord = surj_bed[match(orig_bed_surj$peak_id,surj_bed$peak_id),]

orig_bed_surj$surjected_asm= last_surj_asm
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

all = as.data.frame(rbind(orig_bed_surj, orig_bed_NOT_surj))

orig_bed_surj=all[which(all$surjected_asm != "NONE"),]
orig_bed_NOT_surj=all[which(all$surjected_asm == "NONE"),]

write.table(all, output_all, sep = "\t", col.names = T, row.names = F, quote=F)
write.table(orig_bed_surj, output_surj, sep = "\t", col.names = T, row.names = F, quote=F)
write.table(orig_bed_NOT_surj, output_NOTsurj, sep = "\t", col.names = T, row.names = F, quote=F)

print(head(orig_bed_NOT_surj))
next_asm_tbl = as.data.frame(table(orig_bed_NOT_surj$asm))
print(head(next_asm_tbl))
next_asm_tbl_ord = next_asm_tbl[order(-next_asm_tbl$Freq),]
asm_to_return = as.character(next_asm_tbl_ord[1,1])

writeLines(asm_to_return, next_asm_file)


