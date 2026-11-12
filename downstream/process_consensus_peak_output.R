## Downstream processing of a bonFIRE run: calls peaks in each sample haplotype, adds per-peak counts
## to the consensus peak summary, and tests each peak for haplotype selectivity.
##
## usage: Rscript process_consensus_peak_output.R input_dir output_dir prefix [exclude_runs]
##   input_dir: the pipeline's output folder, [output_dir]/[prefix]/
##   output_dir: where to write the outputs (created if needed)
##   prefix: the run prefix
##   exclude_runs: optional, comma-separated. Any sample_id containing one of these is left out, e.g.
##                 "HG002_PS01015_PacBio,HG01358_PS00999_ONT"
## out: [prefix]_consensus_peak_summary_processed.tsv and [prefix]_hap_selectivity_info.tsv
## The table is large for big runs (about 40 GB of memory for 78 HPRC assemblies), so run it as a job.

library(data.table)
library(Biostrings)

## for calling peaks
MIN_FIRE_COV = 4 #(minimum FIRE coverage to call a peak)
MIN_FRAC_ACC = 0.2 #(minimum fraction accessibility to call a peak)

## for calling haplotype-selective peaks
DIFF_THRESHOLD = 0.25 #(minimum H1-H2 accessibility difference for haplotype-selective)
MIN_HAP_COV = 10 #(minimum coverage on both haplotypes to test selectivity)
FISHER_P_CUTOFF = .05 #(maximum nominal Fisher's exact p value for haplotype-selective; no FDR correction)
## chromosomes left out of the haplotype-selectivity tests: one copy in males (X, Y), not diploid (M),
## or no chromosome (AMB)
NON_AUTOSOMAL = c("chrX", "chrY", "chrM", "AMB")

args = commandArgs(trailingOnly = TRUE)
if(length(args) < 3){
  stop("usage: Rscript process_consensus_peak_output.R input_dir output_dir prefix [exclude_runs]")
}
input_dir = args[1]
output_dir = args[2]
prefix = args[3]
exclude_runs = if(length(args) >= 4 && args[4] != "") strsplit(args[4], ",")[[1]] else character()

## make sure both directories end in exactly one "/", since file names are built by pasting onto them
input_dir = sub("/*$", "/", input_dir)
output_dir = sub("/*$", "/", output_dir)

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

print("reading in tables")
tbl = fread(paste0(input_dir, prefix, "_final_tbl_w_metadata.tsv.gz"), sep = "\t", header = T)
cons_peak_summary = fread(paste0(input_dir, prefix, "_consensus_peak_summary.tsv"), header = T, sep = "\t")

### leave out any runs you don't want in the analysis. Peaks only called in these runs are dropped below
excluded = unique(tbl$sample_id)[sapply(unique(tbl$sample_id), function(x) any(sapply(exclude_runs, grepl, x = x, fixed = TRUE)))]
print(paste0("leaving out ", length(excluded), " sample_ids: ", paste(excluded, collapse = ", ")))
tbl = tbl[!(sample_id %in% excluded)]

### call peaks (coverage columns can be NA where the peak didn't transfer, so NA counts as not a peak)
tbl[, is_peak := fcoalesce(fire_coverage >= MIN_FIRE_COV & fire_cov_OV_cov >= MIN_FRAC_ACC, FALSE)]

### remove any consensus_peak_ids not called a peak in any sample
cpid_peak_in_at_least_one_sample = unique(tbl[is_peak == TRUE, consensus_peak_id])
tbl = tbl[consensus_peak_id %in% cpid_peak_in_at_least_one_sample]
cons_peak_summary = cons_peak_summary[consensus_peak_id %in% cpid_peak_in_at_least_one_sample]

### for each consensus peak id, the number of sample haplotypes (sample_ids) and of assembly haplotypes (asm) in
### which it is a peak, and the number of assembly haplotypes it could be transferred to. The first two differ
### only when an assembly has more than one run
num_peaks_per_cons_peak_id_by_sample_hap = tbl[is_peak == TRUE, .(num_sample_haps_with_peak = uniqueN(sample_id)), by = consensus_peak_id]
num_peaks_per_cons_peak_id_by_asm = tbl[is_peak == TRUE, .(num_asms_with_peak = uniqueN(asm)), by = consensus_peak_id]
num_asm_per_cons_peak_id = tbl[!is.na(asm_chr), .(num_asms_with_genomic_location = uniqueN(asm)), by = consensus_peak_id]

cons_peak_summary = merge(cons_peak_summary, num_peaks_per_cons_peak_id_by_sample_hap, by = "consensus_peak_id", all.x = TRUE)
cons_peak_summary = merge(cons_peak_summary, num_peaks_per_cons_peak_id_by_asm, by = "consensus_peak_id", all.x = TRUE)
cons_peak_summary = merge(cons_peak_summary, num_asm_per_cons_peak_id, by = "consensus_peak_id", all.x = TRUE)

cons_peak_summary[, autosomal := !(chr %in% NON_AUTOSOMAL)]

fwrite(cons_peak_summary, paste0(output_dir, prefix, "_consensus_peak_summary_processed.tsv"), sep = "\t", col.names = T, row.names = F, quote = F, na = "NA")

### identify haplotype-selective peaks
### we will limit this analysis to autosomal peaks

autosomal_peaks = cons_peak_summary[autosomal == TRUE, consensus_peak_id]
tbl = tbl[consensus_peak_id %in% autosomal_peaks]

## sample_id ends in the haplotype, e.g. HG00438_PS00991_ONT_HG00438_1 -> sample HG00438_PS00991_ONT_HG00438, haplotype 1
tbl[, haplotype := sub("^.*_", "", sample_id)]
tbl[, sample := sub("_[^_]*$", "", sample_id)]

tbl_haplotype = dcast(tbl, consensus_peak_id + sample ~ haplotype,
                      value.var = c("coverage", "fire_coverage", "fire_cov_OV_cov"), sep = "_")

tbl_haplotype[, nonfire_1 := coverage_1 - fire_coverage_1]
tbl_haplotype[, nonfire_2 := coverage_2 - fire_coverage_2]
tbl_haplotype[, abs_hap_diff := abs(fire_cov_OV_cov_1 - fire_cov_OV_cov_2)]

## only rows with enough coverage on both haplotypes and a large enough difference can be
## haplotype-selective, so only those are tested (fewer tests, which also matters if an FDR correction is added)
tbl_haplotype[, tested := fcoalesce(coverage_1 >= MIN_HAP_COV & coverage_2 >= MIN_HAP_COV & abs_hap_diff >= DIFF_THRESHOLD, FALSE)]

print(paste0("calculating p values for ", sum(tbl_haplotype$tested), " of ", nrow(tbl_haplotype), " peak/sample pairs"))
tbl_haplotype[, fisher_p := NA_real_]
tbl_haplotype[tested == TRUE, fisher_p := mapply(function(a, b, c, d) {
    fisher.test(matrix(c(a, b, c, d), nrow = 2))$p.value
  }, fire_coverage_1, nonfire_1, fire_coverage_2, nonfire_2)]
print("done calculating p values")

## optional FDR correction over the tested pairs (not used for the calls below):
#tbl_haplotype[, fisher_fdr := NA_real_]
#tbl_haplotype[tested == TRUE, fisher_fdr := p.adjust(fisher_p, method = "BH")]

tbl_haplotype[, hap_selective := tested & fcoalesce(fisher_p <= FISHER_P_CUTOFF, FALSE)]

### determine whether sequences are shared

seqs_output = fread(paste0(input_dir, prefix, "_transferred_cons_peaks_native_sequences_from_main_output.bed.gz"), header = F, sep = "\t")
colnames(seqs_output) = c("asm_chr", "asm_start", "asm_end", "consensus_peak_id", "seq")

seq_output_w_sample_haplotype <- merge(
  tbl[, .(consensus_peak_id, sample, haplotype, asm_chr, asm_start, asm_end)],
  seqs_output,
  by = c("asm_chr", "asm_start", "asm_end", "consensus_peak_id")
)
rm(seqs_output)

seq_output_by_haplotype = dcast(seq_output_w_sample_haplotype, consensus_peak_id + sample ~ haplotype,
                                value.var = c("seq", "asm_chr", "asm_start", "asm_end"), sep = "_")
rm(seq_output_w_sample_haplotype)

## an assembly contig can run either way relative to the other haplotype's, so the sequences match if
## they are identical or one is the reverse complement of the other
seq_output_by_haplotype_has_both_seqs = seq_output_by_haplotype[!is.na(seq_1) & !is.na(seq_2)]
seq_output_by_haplotype_has_both_seqs[, seq_2_rev_comp := as.character(reverseComplement(DNAStringSet(seq_2)))]
seq_output_by_haplotype_has_both_seqs[, seqs_match := seq_1 == seq_2 | seq_1 == seq_2_rev_comp]

## all.x keeps peak/sample pairs without a sequence on both haplotypes (seqs_match is NA for them)
tbl_haplotype_w_seq_match = merge(tbl_haplotype, seq_output_by_haplotype_has_both_seqs[, .(consensus_peak_id, sample, seqs_match)],
                                  by = c("consensus_peak_id", "sample"), all.x = TRUE)
print(table(tbl_haplotype_w_seq_match[, .(hap_selective, seqs_match)], useNA = "ifany"))

fwrite(tbl_haplotype_w_seq_match, paste0(output_dir, prefix, "_hap_selectivity_info.tsv"), sep = "\t", col.names = T, row.names = F, quote = F, na = "NA")
