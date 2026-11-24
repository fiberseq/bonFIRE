# Example data

A small subset of a bonFIRE run (multiple reference assembly mode, CHM13 as the main reference and
GRCh38 as the secondary reference) on 45 HPRC Fiber-seq runs: 100 consensus peaks around HLA-DRB1 and
HLA-DQB1 on chr6, for every sample haplotype (90). It's enough to try the downstream scripts and render
the vignette without running the pipeline.

| File | Contents |
|---|---|
| HPRC_HLA_example_final_tbl_w_metadata.tsv.gz | the main table for these consensus peaks (9,000 rows) |
| HPRC_HLA_example_consensus_peak_summary.tsv | the consensus peak summary for them |
| HPRC_HLA_example_transferred_cons_peaks_native_sequences_from_main_output.bed.gz | their sequences on each assembly |
| HPRC_HLA_example_intervals.bed | two intervals to plot (CHM13 coordinates), with each gene's TSS |

The prefix is `HPRC_HLA_example`. From the repo folder:

```bash
# calling peaks and haplotype-selective peaks
Rscript downstream/process_consensus_peak_output.R example_data my_output_folder HPRC_HLA_example

# plotting the two intervals
Rscript downstream/plot_consensus_peaks_by_interval.R my_plots HPRC_HLA_example \
  example_data/HPRC_HLA_example_final_tbl_w_metadata.tsv.gz \
  example_data/HPRC_HLA_example_consensus_peak_summary.tsv \
  example_data/HPRC_HLA_example_intervals.bed 6 FALSE
```

To render the vignette (in R; rendering needs pandoc, which RStudio includes). The vignette's code runs
from its own folder, so `input_dir` is relative to docs/vignettes/ (or give a full path):

```r
rmarkdown::render("docs/vignettes/vignette_process_main_output_and_consensus_peak_summary.Rmd",
                  output_dir = "my_output_folder",
                  params = list(input_dir = "../../example_data/", prefix = "HPRC_HLA_example"))
```
