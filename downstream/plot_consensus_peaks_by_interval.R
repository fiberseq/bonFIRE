## Heatmaps of fire_cov_OV_cov (samples x consensus peaks) for intervals of interest, one plot per interval,
## plus a second plot with the samples clustered (_clustered.png) and a table of the plotted rows.
## Peaks are ordered along the chromosome by the chr_order column of [prefix]_consensus_peak_summary.tsv,
## so non-reference peaks appear between the CHM13 peaks they sit between. Both tables can be subset first
## (by sample, peak, chromosome...): a peak is plotted if it is in both.
##
## usage: Rscript plot_consensus_peaks_by_interval.R output_dir prefix final_tbl summary intervals.bed \
##          coverage_cutoff remove_single_sample_peaks [min_frac_actuation]
##   output_dir: where to write the plots (created if needed)
##   prefix: start of the output file names
##   final_tbl: [prefix]_final_tbl_w_metadata.tsv.gz, or a subset of it
##   summary: [prefix]_consensus_peak_summary.tsv (or the processed one, or a subset); needs
##            consensus_peak_id, chr, CHM13_start, CHM13_end, chr_slot, chr_order and at_end_of_chr
##   intervals.bed: tab-separated, no header: chr, start, end, name, TSS (CHM13 coordinates; TSS can be NA).
##                  The plot covers the CHM13 peaks overlapping the interval and every peak between them
##   coverage_cutoff: tiles with coverage below this are grey
##   remove_single_sample_peaks: TRUE to drop peaks transferred to only one assembly
##   min_frac_actuation: optional (default 0.10); peaks whose highest fire_cov_OV_cov is below this are dropped
## Tiles are black where the peak didn't transfer to that assembly, and * marks more than one original peak.
args <- commandArgs(trailingOnly = TRUE)
if(length(args) < 7){
  stop("usage: Rscript plot_consensus_peaks_by_interval.R output_dir prefix final_tbl summary intervals.bed coverage_cutoff remove_single_sample_peaks [min_frac_actuation]")
}

library(ggplot2)
library(dplyr)
library(tidyr)
library(data.table)

output_dir=args[1]
## make sure output_dir ends in exactly one "/", since file names are built by pasting onto it
output_dir = sub("/*$", "/", output_dir)
prefix=args[2]
peaks_input=args[3]
summary_input=args[4]
interval_input=args[5]
coverage_cutoff=as.numeric(args[6])
remove_single_sample_peaks = as.logical(args[7])
min_frac_actuation = if(length(args) >= 8) as.numeric(args[8]) else 0.10

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

print("reading_in_peaks_table")

## fread reads .gz directly (it needs the R.utils package for that)
tbl = as.data.frame(fread(peaks_input, header = TRUE, sep = "\t", na.strings = "NA"))
## tables from older pipeline runs still carry the peak sequences; they aren't used here
tbl$seq = NULL

colnames(tbl) <- gsub("^X\\.", "", colnames(tbl))

print("reading_in_peak_summary_table")

peak_summary = as.data.frame(fread(summary_input, header = TRUE, sep = "\t", na.strings = "NA"))
needed = c("consensus_peak_id", "chr", "CHM13_start", "CHM13_end", "chr_slot", "chr_order", "at_end_of_chr")
missing = setdiff(needed, colnames(peak_summary))
if(length(missing) > 0){
  stop("the peak summary is missing ", paste(missing, collapse = ", "), "; it needs a pipeline run that adds chr_order and at_end_of_chr")
}
## AMB peaks have no chromosome, so no order, and can't be plotted
print(paste0("peaks in the summary without a chromosome order (AMB): ", sum(is.na(peak_summary$chr_order))))
peak_summary = peak_summary[which(!is.na(peak_summary$chr_order)),]

if(remove_single_sample_peaks == TRUE){
  ## keep peaks transferred to more than one assembly
  tbl = as.data.frame(tbl %>%
                        group_by(consensus_peak_id) %>%
                        filter(n_distinct(asm[!is.na(asm_chr)]) > 1) %>%
                        ungroup())
}

print("reading_in_bed_intervals_table")

bed_intervals = read.table(interval_input, sep = "\t", header = F, stringsAsFactors = F)
temp_colnames = c("chr", "start", "end","gene", "tss")
colnames(bed_intervals) = temp_colnames[1:ncol(bed_intervals)]
if(!("tss" %in% colnames(bed_intervals))){
  bed_intervals$tss = NA
}


tbl$meets_cutoff = FALSE
tbl$meets_cutoff[which(tbl$coverage >= coverage_cutoff)] = TRUE

peaks_to_keep = unique(tbl[which(!is.na(tbl$fire_cov_OV_cov)),]$consensus_peak_id)

tbl = tbl[which(tbl$consensus_peak_id %in% peaks_to_keep),]

max_fire_score_per_peak = tbl %>%
  group_by(consensus_peak_id) %>%
  summarise(max_value = max(fire_cov_OV_cov, na.rm = TRUE))

peaks_to_remove = max_fire_score_per_peak[which(max_fire_score_per_peak$max_value < min_frac_actuation),]$consensus_peak_id

if(length(peaks_to_remove) > 0){
  tbl = tbl[which(!(tbl$consensus_peak_id %in% peaks_to_remove)),]
}

## only peaks in both tables can be plotted
peak_summary = peak_summary[which(peak_summary$consensus_peak_id %in% tbl$consensus_peak_id),]
tbl = tbl[which(tbl$consensus_peak_id %in% peak_summary$consensus_peak_id),]

### plotting

label_pt = 8

## heatmap of fire_cov_OV_cov, samples on y and peaks on x (in the order of the factor levels)
make_interval_plot = function(plot_tbl, title, tss_columns, legend_on_top){
  p = ggplot(plot_tbl, aes(x=consensus_peak_id, y=sample_id, fill = fire_cov_OV_cov)) + geom_tile() +
    theme(axis.text = element_text(size = label_pt),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    geom_tile(color = "black") +
    scale_fill_gradient(low = "white", high = "purple") +
    ggtitle(title) +
    geom_text(
      aes(label = ifelse(num_overlapping_orig_peaks > 1, "*", "")),
      color = "black",
      size = 6
    )

  if(nrow(subset(plot_tbl, meets_cutoff == FALSE)) > 0){
    p = p + geom_tile(data = subset(plot_tbl, meets_cutoff == FALSE),
                      fill = "grey50", color = "black")
  }

  p = p + geom_tile(data = subset(plot_tbl, score == -2),
                    fill = "black", color = "black")

  p = p + geom_tile(data = subset(plot_tbl, is.na(score)),
                    fill = "black", color = "black")

  if(legend_on_top){
    p = p + theme(legend.position = "top")
  }

  if(length(tss_columns) > 0){
    p = p + geom_tile(data = plot_tbl[which(plot_tbl$consensus_peak_id %in% tss_columns),],
                      fill = "NA", color = "black", linewidth = .8)
  }
  return(p)
}

## sample_ids ordered by hierarchical clustering on their fire_cov_OV_cov across the plotted peaks.
## Missing values (peak not transferred, or no score) are information too, so each peak gives two
## features: fire_cov_OV_cov (0 when missing) and whether it is missing (0/1). Two samples then
## differ by 1 for a peak that one has and the other doesn't, on top of any difference in value
cluster_sample_order = function(plot_tbl){
  wide = dcast(as.data.table(plot_tbl), sample_id ~ consensus_peak_id, value.var = "fire_cov_OV_cov",
               fun.aggregate = function(x) mean(x, na.rm = TRUE))
  mat = as.matrix(wide[, -1, with = FALSE])
  rownames(mat) = wide$sample_id
  if(nrow(mat) < 3){
    return(rownames(mat))
  }
  missing = is.na(mat)
  features = cbind(replace(mat, missing, 0), missing * 1)
  hc = hclust(dist(features), method = "average")
  return(rownames(mat)[hc$order])
}

### run through intervals

print("starting_to_run_through_intervals")

for(k in seq_len(nrow(bed_intervals))){
  temp_gene = bed_intervals[k,"gene"]
  interval_chr = bed_intervals[k,"chr"]
  interval_start = bed_intervals[k,"start"]
  interval_end = bed_intervals[k,"end"]
  tss = bed_intervals[k,"tss"]
  print(temp_gene)

  ## the interval is in CHM13 coordinates, so it is set by the CHM13 peaks in it: everything from the
  ## first one's chr_order to the last one's is plotted, non-reference peaks in between included.
  ## Non-reference peaks before the first or after the last CHM13 peak in the interval aren't, as they
  ## could be on either side of the interval's edge, except AT_START / AT_END peaks when the
  ## interval includes the chromosome's first / last CHM13 peak (chr_slot 2 / the one before AT_END)
  chr_summary = peak_summary[which(peak_summary$chr == interval_chr),]
  chm13_in_int = chr_summary[which(!is.na(chr_summary$CHM13_start) & chr_summary$CHM13_end > interval_start & chr_summary$CHM13_start < interval_end),]
  if(nrow(chm13_in_int) == 0){
    print("no CHM13 peaks in the interval, skipping")
    next
  }
  first_order = min(chm13_in_int$chr_order)
  last_order = max(chm13_in_int$chr_order)
  first_slot = min(chm13_in_int$chr_slot)
  last_slot = max(chm13_in_int$chr_slot)
  in_int = (chr_summary$chr_order >= first_order & chr_summary$chr_order <= last_order) |
    (first_slot == 2 & chr_summary$chr_slot == 1) |
    (chr_summary$at_end_of_chr %in% TRUE & chr_summary$chr_slot == last_slot + 1)
  int_summary = chr_summary[which(in_int),]
  int_summary = int_summary[order(int_summary$chr_order),]
  consensus_peak_vect_ord = int_summary$consensus_peak_id
  print(paste0("peaks: ", length(consensus_peak_vect_ord), " (non-reference: ", sum(is.na(int_summary$CHM13_start)), ")"))

  all_peaks_to_include_tbl = tbl[which(tbl$consensus_peak_id %in% consensus_peak_vect_ord),]
  print(paste0("rows_in_table: ", (nrow(all_peaks_to_include_tbl))))
  if(nrow(all_peaks_to_include_tbl) > 1){
    m = match(all_peaks_to_include_tbl$consensus_peak_id, int_summary$consensus_peak_id)
    all_peaks_to_include_tbl$chr_slot = int_summary$chr_slot[m]
    all_peaks_to_include_tbl$chr_order = int_summary$chr_order[m]
    all_peaks_to_include_tbl$consensus_peak_id = factor(all_peaks_to_include_tbl$consensus_peak_id, levels = consensus_peak_vect_ord)

    tss_columns = NULL
    if (!is.na(tss)){
      chm13_plotted = int_summary[which(!is.na(int_summary$CHM13_start)),]
      temp = chm13_plotted$consensus_peak_id[which(chm13_plotted$CHM13_start < tss & chm13_plotted$CHM13_end > tss)]
      if(length(temp) > 0){
        tss_columns = temp
      } else {
        ## the CHM13 peaks on either side of the TSS (only one side if the TSS is past the end)
        lower <- which(chm13_plotted$CHM13_end < tss)
        upper <- which(chm13_plotted$CHM13_start > tss)
        tss_rows = c(if(length(lower) > 0) max(lower), if(length(upper) > 0) min(upper))
        tss_columns = chm13_plotted$consensus_peak_id[tss_rows]
      }
    }

    print("plotting...")

    ## size the plot from what is in it, so tiles and labels look the same for 2 samples or 100+:
    ## a fixed size per tile, plus room for the longest label on each axis (x labels are
    ## rotated, so the longest consensus_peak_id sets how much height they take)
    inch_per_char = label_pt * 0.6 / 72   # about how wide one character is
    tile_height = 0.18                    # inches per sample row
    tile_width = 0.2                      # inches per peak column
    n_samples = length(unique(all_peaks_to_include_tbl$sample_id))
    n_peaks = length(consensus_peak_vect_ord)
    x_labels_height = max(nchar(consensus_peak_vect_ord)) * inch_per_char + 0.4
    y_labels_width = max(nchar(as.character(all_peaks_to_include_tbl$sample_id))) * inch_per_char + 0.4
    ## the legend needs about 1.5 inches of height beside the tiles; with fewer rows than that it
    ## goes above them instead, so the tiles keep their size
    legend_on_top = n_samples * tile_height < 1.5
    plot_height = n_samples * tile_height + x_labels_height + 0.8 + ifelse(legend_on_top, 0.8, 0)  # title and margins
    plot_width = n_peaks * tile_width + y_labels_width + ifelse(legend_on_top, 0.3, 1.5)            # legend and margins
    plot_width = max(plot_width, 6)

    plot_title = paste0(temp_gene, " ",  interval_chr, ":",interval_start,"_", interval_end)
    file_base = paste0(output_dir, prefix, "_", interval_chr, "_",interval_start,"-", interval_end, "_", temp_gene)

    full_int_plot = make_interval_plot(all_peaks_to_include_tbl, plot_title, tss_columns, legend_on_top)
    ggsave(
      filename = paste0(file_base, ".png"),
      plot = full_int_plot,
      width = plot_width,
      height = plot_height,
      dpi = 300,
      limitsize = FALSE
    )

    ## same plot with the samples in clustered order
    clustered_tbl = all_peaks_to_include_tbl
    clustered_tbl$sample_id = factor(clustered_tbl$sample_id, levels = cluster_sample_order(clustered_tbl))
    clustered_plot = make_interval_plot(clustered_tbl, paste0(plot_title, " (samples clustered)"), tss_columns, legend_on_top)
    ggsave(
      filename = paste0(file_base, "_clustered.png"),
      plot = clustered_plot,
      width = plot_width,
      height = plot_height,
      dpi = 300,
      limitsize = FALSE
    )

    write.table(all_peaks_to_include_tbl[order(all_peaks_to_include_tbl$consensus_peak_id),], paste0(output_dir, prefix, "_", interval_chr, "_",interval_start,"-", interval_end, "_", temp_gene, "_peaks_file.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)
  }
}
