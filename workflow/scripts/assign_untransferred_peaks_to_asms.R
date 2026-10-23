#!/usr/bin/env Rscript
## Assign each peak that didn't surject onto CHM13 to an assembly.
##
## The peaks are surjected onto every assembly in parallel beforehand (one bed per
## assembly). This script then picks assemblies in the same greedy order the old
## surjection loop used: the assembly with the most unassigned peaks goes first, and
## its peaks that surjected (and passed the length filter) are assigned to it. Then
## the next assembly is picked from the peaks that are left, until none are left or
## every assembly has been tried. A peak surjects onto an assembly regardless of the
## other peaks, so this gives the same result as surjecting one assembly at a time.
##
## Ties (two assemblies with the same number of unassigned peaks) go to the one that
## sorts first alphabetically.
args = commandArgs(trailingOnly=TRUE)

orig_bed_file=args[1]
asm_conversion_file=args[2]
length_change_filter=as.numeric(args[3])
output_all=args[4]
output_tried=args[5]
# comma-separated graph assembly names, in the same order as the surjected beds that follow
surj_asms=strsplit(args[6], ",")[[1]]
surj_bed_files=args[7:length(args)]
stopifnot(length(surj_asms) == length(surj_bed_files))

orig_bed = read.table(orig_bed_file, header = T, sep = "\t", stringsAsFactors = F, comment.char = "")
colnames(orig_bed)[1] <- "chrom"

asm_input = read.table(asm_conversion_file, header = F, sep = "\t", stringsAsFactors = F, comment.char = "")
asm_to_graph = setNames(asm_input[,2], asm_input[,1])

# the surjected bed is empty when no peaks surjected onto this assembly
read_surj_bed = function(f){
  if(file.size(f) > 0){
    surj_bed = read.table(f, header = F, sep = "\t", stringsAsFactors = F, comment.char = "")
  } else {
    surj_bed = orig_bed[0,]
  }
  colnames(surj_bed) = colnames(orig_bed)
  return(surj_bed)
}
surj_beds = setNames(lapply(surj_bed_files, read_surj_bed), surj_asms)

all = orig_bed
tried_asms = c("CHM13")

repeat {
  orig_bed_NOT_surj = all[which(all$surj_chr == "NONE"),]
  candidate_asms = orig_bed_NOT_surj$asm[!(asm_to_graph[orig_bed_NOT_surj$asm] %in% c(tried_asms, NA))]
  if(length(candidate_asms) == 0){
    break
  }
  next_asm_tbl = as.data.frame(table(candidate_asms))
  next_asm_tbl_ord = next_asm_tbl[order(-next_asm_tbl$Freq),]
  last_surj_asm = as.character(next_asm_tbl_ord[1,1])
  asm_graph_name = asm_to_graph[[last_surj_asm]]
  if(!(asm_graph_name %in% surj_asms)){
    stop(paste0("No surjected peaks file for ", asm_graph_name, " (assembly ", last_surj_asm, ")"))
  }
  print(paste0("assigning peaks to ", last_surj_asm, " (", asm_graph_name, ")"))
  tried_asms = c(tried_asms, asm_graph_name)

  ## same steps as the old per-assembly annotation script
  surj_bed = surj_beds[[asm_graph_name]]

  orig_bed_previously_surjected = all[which(all$surjected_asm != "NONE"),]
  remaining = all[which(all$surjected_asm == "NONE"),]
  surj_bed = surj_bed[which(surj_bed$peak_id %in% remaining$peak_id),]

  orig_bed_surj = remaining[which(remaining$peak_id %in% surj_bed$peak_id),]
  orig_bed_NOT_surj = remaining[which(!(remaining$peak_id %in% surj_bed$peak_id)),]

  surj_bed_ord = surj_bed[match(orig_bed_surj$peak_id,surj_bed$peak_id),]

  orig_bed_surj$surjected_asm = rep(last_surj_asm, nrow(orig_bed_surj))
  orig_bed_surj$surj_chr = surj_bed_ord$chrom
  orig_bed_surj$surj_start = surj_bed_ord$start
  orig_bed_surj$surj_end = surj_bed_ord$end

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

  all = as.data.frame(rbind(orig_bed_previously_surjected, orig_bed_surj, orig_bed_NOT_surj))
}

write.table(all, output_all, sep = "\t", col.names = T, row.names = F, quote=F)
writeLines(tried_asms, output_tried)
