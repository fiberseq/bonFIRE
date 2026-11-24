# bonFIRE

bonFIRE is a toolkit designed to simplify comparative analysis of Fiber-seq data across samples as well as enable **direct comparison of Fiber-seq data processed against multiple distinct genome assemblies.**

The input to bonFIRE is one or several [FIRE](https://fiberseq.github.io/fire/run.html) runs. See **Shared reference assembly (simple) mode: inputs and outputs** and **Multiple reference assembly mode: inputs and outputs** for more information on generating these.

### Some common use cases for bonFIRE are:

(1) Compare the chromatin landscape between samples processed against a **common haploid reference** (e.g. GRCh38).

(2) Compare the chromatin landscape between two haplotypes in a single sample, processed against its **own native diploid assembly**. We can perform this analysis for a single sample, or for a group of samples all processed against a common diploid assembly. These samples may be different tissues from one individual, different experimental conditions, etc. 

(3) Compare the chromatin landscape between samples which are processed in different assemblies. This enables us, for example, to evaluate the chromatin landscape of each sample in its own native donor assembly, and then combine them into a single dataset which **does not rely on a single reference**, and thus enabling comparison at loci that are only shared between a subset of assemblies. A combination of 2 & 3, where subset of samples share a reference is also allowed.

In order to perform inter-assembly comparisons, examples 2 and 3 make use of **pangenome graphs**, which inherently define homologies between assembly coordinates, and can be used to transfer chromatin accessibility peak calls between assemblies. This is described further below. 

See section titled **Consensus Peak Calling Approach** for an in depth overview of the method.

## Installation 

### 1. Install Pixi

Please start by installing [pixi](https://pixi.sh/latest/) which handles the environment of this Snakemake workflow.

### 2. Clone the repository

Clone the repository and move into the project directory:

```bash
git clone https://github.com/fiberseq/bonFIRE.git
cd bonFIRE
```

### 3. Install dependencies

You can then install the `pixi` environment by running the following in this repository:

```bash
pixi install
```

### 4. Build fibertools from GitHub


There are several bug fixes on the fibertools `main` branch ([3ba4545](https://github.com/fiberseq/fibertools-rs/commit/3ba4545)) that are not yet in the released version on bioconda. It is safest to build `ft` from GitHub and give its path as `ft_path` in the config. This needs Rust (install it with [rustup](https://rustup.rs) if `cargo` isn't available). (This build requires a recent version of Rust; if the build fails with 'upgrade rustc or select compatible dependency versions', run rustup update stable and try again.)

```bash
cargo install --git https://github.com/fiberseq/fibertools-rs --locked --root /path/to/fibertools_main
```

This puts the binary at `/path/to/fibertools_main/bin/ft`. Include this line in your config.yaml to use this version instead of the bioconda one.

```yaml
ft_path: /path/to/fibertools_main/bin/ft
```

The build still reports its version as 0.13.0, so the version number doesn't tell you whether you have the fix. If the build fails, see the fibertools [install instructions](https://fiberseq.github.io/fibertools/install.html). Once a fixed release is on bioconda, this step won't be needed.


## Usage

To run from inside the repo:

```bash
pixi run snakemake \
--configfile path/to/your/config.yaml \
-k \
...
```
To execute across a cluster (modify profiles/slurm-executor as needed):

```bash
pixi run snakemake \
--configfile path/to/your/config.yaml \
--profile path-to-repo/workflow/profiles/slurm-executor \
-k \
...
```

To run from another directory, specify repo via --manifest-path: 

```bash
pixi run \
--manifest-path path-to-repo/pixi.toml snakemake \
--configfile path/to/your/config.yaml \
-k \
...
```

And in place of `...` use all the normal Snakemake arguments for your workflow.

## Consensus Peak Calling Approach

bonFIRE is a computational tool designed to enable direct comparative analysis between Fiber-seq datasets (processed with FIRE). Its primary output is a matrix (and accompanying metadata), with sample-specific chromatin actuation values at a set of "consensus peaks" - genomic regions of interest where at least one input sample was accessible, and thus which may be of interest in a comparative analysis. Importantly, chromatin actuation need NOT be measured in a common genome reference space across input samples. Rather, a consensus peak can represent a set of homologous regions across multiple assemblies (or assembly haplotypes, for diploid inputs), identified via a pangenome graph-based approach described below. Notably, the set of consensus peaks defined here need not all exist in a single genome assembly, enabling comparative analysis at structurally complex regions which may be present in only a subset of input assemblies. 

<img width="582" height="414" alt="github_figs-05" src="https://github.com/user-attachments/assets/76b5283d-0b78-47e6-a9a5-742cbcb47b7c" />

See "**bonFIRE with shared reference assembly**" below for further details on running this tool with Fiber-seq datasets measured in a common (usually haploid) reference assembly, and "**bonFIRE with multiple reference assemblies (pangenome-based)**" for running bonFIRE on samples processed against different assemblies. 

Note, for a sample (or set of samples) processed against one diploid reference assembly, the pangenome-based mode allows us to compare not just between samples, but across homologous regions between the two haplotypes in a single sample. Thus, it is recommended that the multiple reference mode be used for this application.

### bonFIRE with shared reference assembly

The simplest case of bonFIRE takes as input Fiber-seq data from multiple samples which is mapped against a single common reference assembly (e.g. GRCh38). When each sample is processed via the FIRE pipeline, two relevant outputs are produced: a [sample].peaks.bed.gz file with accessibility peak coordinates indicated, and a [sample].pileup.bed.gz file with base level raw accessibility information. These two files (one per sample) are used to: <br>
(1) Identify a set of consensus peaks: genomic regions where a peak was called in at least one sample (using [sample].peaks.bed.gz) <br>
(2) Evaluate chromatin accessibility at these consensus peaks in each sample, whether or not a peak was originally called in that sample. (using [sample].pileup.bed.gz)

<img width="402" height="460" alt="github_figs-03" src="https://github.com/user-attachments/assets/4b05afe0-2880-483c-9277-66d9ebd2b4ef" />

The output contains raw counts from the pileup, and thus a chromatin actuation score can be recalculated with filters specific to the whole group. This approach is thus effectively normalizes for read depth differences between samples (where a peak may be called in one sample and not the other simply due to read depth differences), and is preferable to simply looking for overlaps of called peaks in the sample-specific peak outputs. 

### bonFIRE with multiple reference assemblies (pangenome-based)

Mapping samples to a common reference assembly results in loss and/or distortion of information at structurally complex regions, whereas analyzing samples in their own native genome context precludes inter-sample comparison, as each assembly has its own genomic coordinate system. bonFIRE was designed to address this challenge, enabling us to compare chromatin signal at homologous regions between assemblies. 

Note: This mode currently assumes that samples are processed against one or multiple diploid assemblies, but can be hacked to work with haploid assemblies as well. (More official haploid/haploid+diploid reference mode coming soon.)

The multiple reference assembly mode conceptually parallels the common reference assembly approach described above, where a set of consensus peaks are first defined, and then raw chromatin accessibility values are pulled from each sample's native FIRE output at those positions. However, a consensus peak no longer represents a static genomic interval in one assembly space; rather, it is a set of homologous loci across assemblies, associated with a single consensus peak label. To facilitate the identification of homologous regions, we first construct a pangenome graph of all input assemblies as well as one or two 'reference' assemblies for easy downstream interpretation of data against annotations in those reference spaces. We can use these graphs to move features (e.g. FIRE peak coordinates from a single sample) between any pair of assemblies in the graph. Defining consensus peaks involves progressively moving the full set to natively-called FIRE peaks across all samples to shared assembly spaces, and defining consensus peaks within that shared coordinate system. Each consensus peak is ultimately transferred back through the graph to each native assembly, and natively measured chromatin accessibility data is pulled from that region. This process retains genomic loci that only exist in one or a subset of assemblies, as well as those ubiquitously present across the population. 

<img width="1121" height="1633" alt="github_figs-06" src="https://github.com/user-attachments/assets/de343e71-d678-45d9-b5e0-96ca4bab883b" />



## Shared reference assembly (simple) mode: inputs and outputs

The main inputs to both bonFIRE modes are the single sample outputs of [FIRE](https://fiberseq.github.io/fire/run.html), generated using 'pixi run fire' (see **3.1 Running and installing** under the 3. FIRE on the left-hand side menu). Specifically, it will look for the *.peaks.bed.gz file and *.pileup.bed.gz file in the FIRE output folder indicated in the input file (described below). **Be sure that only one file in the folder provided has the suffixes above.**

Both modes above take in an input.txt file with sample information and a config.yaml. See sections below for mode-specific structure.

### Shared reference assembly inputs

#### config.yaml
An example config can be found under **config/sample_configs/config_shared_ref.yaml**.
The **required** fields are listed below. Other fields (see Multiple reference assemblies) should in theory be ignored (though this is not fully tested).

**output_dir**: directory to write to, no trailing / <br>
**prefix**: run prefix. A folder within output_dir will be generated with this name and all outputs will be written there. <br>
**input_file**: Full path to inputs.txt (different name ok), described below <br>
**ft_path:** path to the `ft` (fibertools) binary to use. If removed, the `ft` in the pipeline's environment is used. **Set this to a fibertools build from GitHub until the next release (beyond v0.13.0; see Installation, step 4)** <br>
**shared_ref**: True <br>
**main_reference_fasta**: reference against which samples are processed (gzipped or not). If it has no .fai index (or is gzipped), a copy or link is indexed in the output folder, so nothing is written next to it.

#### inputs.txt

A sample input file for the single reference mode can be found here: <br>
config/sample_configs/input_shared_ref.txt

It should include one line per FIRE run, with a unique 'sample_name' (first column) for each run. These will appear in the output file, so information that would be useful for downstream analysis can be helpful to include (e.g. if multiple FIRE runs on different tissues in the same individual are included, the sample name can be \[donor\]_\[tissue\].)

Input file for this mode has just two columns (tab delimited): <br>
**sample_name**: must be unique for each line. Can't contain '|' or ';' <br>
**path_to_fire_output_folder**: full path to FIRE output folder for this sample, without a trailing /. It must contain exactly one *peaks.bed.gz and one *pileup.bed.gz file directly inside it (the pipeline stops and lists any folder that doesn't).

### Shared reference assembly outputs (simple)
(in [output_dir]/[prefix]/)

**[prefix]_final_tbl_w_metadata.tsv.gz**: one row per run and consensus peak.
- **sample_id**: the run's sample_name
- **consensus_peak_id**, **cons_chr**, **cons_start**, **cons_end**: the consensus peak and its coordinates
- **coverage**: number of reads overlapping this consensus peak region at the point of maximum accessibility 
- **fire_coverage**: number of reads with FIRE elements overlapping this consensus peak region at the point of maximum accessibility.
- **score**:	an accessibility score at this point calculated during original FIRE peak calling and defined [here](https://fiberseq.github.io/fire/methods/aggregation.html). I rarely use this column in downstream analyses, instead focusing on the raw fire_coverage/coverage calculation (in the fire_cov_OV_cov).
- **fire_cov_OV_cov**: raw accessibility value calculated as fire_coverage/coverage. Where the run's pileup has no data at the consensus peak, coverage and fire_coverage are 0, score is -1 and fire_cov_OV_cov is NA.
- **coverage_H1**, **fire_coverage_H1**, **score_H1**, **fire_cov_OV_cov_H1**, and the same for **H2**: the values for each haplotype at the same pileup position. NA for runs whose pileup has no haplotype columns (e.g. unphased runs).
- **overlaps_original_called_peak**: TRUE if there was at least one peak overlapping this consensus peak location in the original FIRE dataset. This information should be used with caution, as differences in coverage between samples can impact whether a peak rose above the FDR threshold in the original FIRE calling scheme.
- **num_overlapping_orig_peaks**: (Probably not useful) The number of originally called FIRE peaks that overlap this consensus peak region. A relic from multi-reference mode and probably not useful here.
- **orig_peaks_overlapping_consensus**: (Probably not useful) Positions of originally called peaks which overlap this consensus peak region, each as chrom_start_end-[sample_name] and separated by commas. *This set is NOT filtered for very nearby peaks as is done for calculating num_overlapping_orig_peaks.
  
**[prefix]_consensus_peak_summary.tsv**: one row per consensus peak.
- **consensus_peak_id**, **chr**, **start**, **end**
- **chr_order**: the peak's position along its chromosome, 1 to the number of peaks on it
- **seq**: the peak's sequence in the reference

Also kept: [prefix]_Ref_shared_peaks.bed and _6col.bed (the consensus peaks), [prefix]_wPeakIDs.bed.gz (the runs' peaks), [prefix]_max_score_per_peak_all_samples.tsv and [prefix]_orig_called_peaks_intersect_consensus_peaks.bed (inputs to the final tables), [prefix]_Ref_chromosomes.txt and [prefix]_runs.txt, and logs/.

The plotting script (see **Downstream analysis**) also reads these outputs; remove_single_sample_peaks has no effect for them, since every peak is in every run.


## Multiple reference assembly mode: inputs and outputs

### Multiple reference assemblies inputs (pangenome-based)

Graph construction is the most time intensive part of this pipeline and **does not** need to be performed each time if a graph has already been constructed. One can use this pipeline to: <br>
- Construct graph only (make_graph_only=True in the config file) <br>
- Construct graph and perform consensus peak calling <br>
- Perform consensus peak calling with provided graph (This is safest to do with a graph generated via this pipeline!) <br>
  
As with the common reference mode above, as input, option 2 and 3 require the outputs of [FIRE](https://fiberseq.github.io/fire/run.html), generated using 'pixi run fire' (see **3.1 Running and installing** under the 3. FIRE on the left-hand side menu). Specifically, it will look for the *.peaks.bed.gz file and *.pileup.bed.gz file in the FIRE output folder indicated in the input file (described below). **Be sure that only one file in the folder provided has the suffixes above.**

bonFIRE requires an **input file** and a **config file**. 

Sample **input files** can be found here. They can be identical, except for an NA in the second column instead of a FIRE output. For graph_only mode, you can also include just one line per assembly, but this isn't strictly necessary. Lines with different sample_name but identical values in all other columns will only produce two graph paths (one per haplotype). sample_name is not used in graph construction -- haplotypes will be denoted by values in hap1_in_graph/hap2_in_graph: <br>
config/sample_configs/input_multiref.txt 
config/sample_configs/input_multiref_graph_only.txt

For this mode, bonFIRE expects that mapping and FIRE has been run against a **diploid** assembly. Downstream, each haplotype will be treated as a separate entity, but the pipeline will take care of splitting the peaks associated with each haplotype for you, based on haplotype-specific subtags you provide below.

The input file for this mode include one line per FIRE run, with a unique 'sample_name' (first column) for each run. These will appear in the output file, so information that would be useful for downstream analysis can be helpful to include (e.g. if multiple FIRE runs on different tissues in the same individual are included, the sample name can be \[donor\]_\[tissue\].)

Input file columns (tab delimited): <br>
**sample_name**: must be unique for each line. Can't contain '|' or ';'. It may be useful for downstream analysis to include any associated info (e.g. tissue, condition, etc.) in this name for parsing in downstream steps. This name will not be used in the graph construction step. <br>
**path_to_fire_output_folder**: full path to FIRE output folder for this sample, without a trailing /. It must contain exactly one *peaks.bed.gz and one *pileup.bed.gz file directly inside it (the pipeline stops and lists any folder that doesn't). If generating graph only, put NA here. <br>
**hap1_orig_subtag** & **hap2_orig_subtag**: This is any substring that can be used to differentiate hap1 chromosome names from hap2 chromosome names. These can be found anywhere in the chromosome name. (e.g. h1/h2, #1#/#2#, MATERNAL/PATERNAL, etc.) <br>
**hap1_in_graph** & **hap2_in_graph**: These are the prefixes used for each haplotype in the graph. I commonly use \[donor-id\]_1 and \[donor-id\]_2. ** If multiple lines of this input file are associated with the same assembly (e.g. multiple tissue samples per donor), these values should be identical in those lines. This way, each unique assembly will only appear in the graph once. <br>
**hap1_fasta_path** & **hap2_fasta_path**: path to hap1 and hap2 assembly fasta files. Make sure unassigned contigs are not included in these files, and all contigs that are included can be distinguished by **hap1_orig_subtag** & **hap2_orig_subtag**. Lines that share an assembly must have the exact same path (not just an identical copy of the file); otherwise that assembly is added to the graph more than once.

Sample **config.yaml** files can be found here: <br>
config/sample_configs/config_multiref_graph_only.yaml <br>
config/sample_configs/config_multiref_graph_provided.yaml <br>
config/sample_configs/config_multiref_build_graph_and_run_consensus_peak_calling.yaml


Notes:
- if generating graph only, set **make_graph_only: True** and remove the **graph:** line
- If using existing graph generated via a previous run, add full path to **graph:**, using the version ending in *.full_all_chrs.gbz
- **make_graph_only: False** plus setting other graph parameters will run both graph building and consensus peak calling
- By default, the graph includes one or two haploid reference genomes (main and secondary below). The main reference is used as illustrated in the figure above (as CHM13 is there). If an optional secondary reference is added, the coordinates for all consensus peaks in that reference space will be returned as additional columns in the consensus peak summary output. 

Fields in config.yaml: <br>
**output_dir**: directory to write to, no trailing / <br>
**prefix**: run prefix. A folder within output_dir will be generated with this name and all outputs will be written there. <br>
**input_file**: Full path to inputs.txt (different name ok), described below <br>
**ft_path:** path to the `ft` (fibertools) binary to use. If removed, the `ft` in the pipeline's environment is used. **Set this to a fibertools build from GitHub until the next release (beyond v0.13.0; see Installation, step 4)** <br>                                       
**graph**: path to graph if one has already been generated. If not, remove this line, and a graph will be generated for you. A provided graph must contain every contig in the input fastas, including those which do not align during graph building. This pipeline adds them back into the graph file to generate [output_dir]/[prefix]/graph/[prefix].full_all_chrs.gbz; use this file rather than [prefix].full.gbz. Else, downstream steps may fail (e.g.`vg inject` fails with "alignment references path not present in graph".) <br>
**main_reference_fasta**: Path to the main reference fasta (haploid; gzipped or not), e.g. CHM13. This reference will serve as the graph scaffold. <br>
**main_reference_name**: e.g. CHM13. This name will be used as a path name in the pangenome graph, and cannot be a substring of hap1_in_graph, hap2_in_graph, hap1_orig_subtag, hap2_orig_subtag from the input file above. Letters, numbers and _ only. <br>
**secondary_reference_fasta**: Optional. Path to a second haploid reference fasta, e.g. GRCh38. It isn't used in consensus peak calling; the pipeline outputs each consensus peak's coordinates on it, which can be useful for annotations defined in that reference. Remove this and secondary_reference_name to skip it. <br>
**secondary_reference_name**: e.g. GRCh38. Same rules as main_reference_name. <br>
(Older configs with CHM13_fasta_path / hg38_fasta_path still work, as references named CHM13 and GRCh38.) <br>
**ft_param**: Default 1000. Peaks are 'strung together' into longer sequences to aid surjection onto the correct path. Keep as is for most cases. <br>
**Peak_length_change_param**: 1 (Keep as is). If a peak length changes by more than 100% after surjection, it will be considered not surjected to that assembly. <br>
**make_graph_only**: Default False. Make True if you'd like to make the graph only, without consensus peak calling. <br>
**split_graph**: If graph is being generated, should it be split by chromosome? Default is False, but may need to split for larger graphs. <br>
**graph_mem**: 300000 (Can increase for large graphs) <br>
**graph_runtime**: 1440 (In minutes. Can increase for large graphs.) <br>
**vg_mem_mb**: Default 80000. Memory for the vg inject/surject steps; depends on the graph (about 11 GB for a 4-haplotype graph, 60 GB for an 11-sample graph, 180-250 GB for an 82-sample graph). Jobs that fail are retried twice with double the memory each time. <br>
**vg_threads**: Default 8. CPUs for the vg inject/surject steps.

If you would like the pipeline to generate a graph, remove the "graph:" line from this config file and set the other graph parameters. make_graph_only can be True or False depending on whether you want your pipeline to continue with downstream consensus peak calling. If a graph is provided, the other graph parameters won't be used.

## Outputs

All outputs are written to **[output_dir]/[prefix]/**, and file names start with **[prefix]_**.

### Main outputs

| File | Contents |
|---|---|
| **[prefix]_final_tbl_w_metadata.tsv.gz** | The main output table: one row per consensus peak per sample haplotype (Output further described below). Chromatin accessibility data is found here. |
| **[prefix]_consensus_peak_summary.tsv** | One row per consensus peak region with associated metadata: reference chromosome for both reference and non-reference peaks, main and secondary reference coordinates where available, relative order of reference and non-reference peaks on a chromosome. (Output further described below). |
| **[prefix]_transferred_cons_peaks_native_sequences_from_main_output.bed.gz** | This file returns the native sequence on each assembly associated with each consensus peak region which could be transferred to that assembly. Useful to determine whether haplotype-selective chromatin might be explained by underlying sequence variation, for example. See downstream vignette for halotype-selective analysis. (Output further described below). |
| **logs/** | One log per rule (and per sample or assembly where a rule runs for each). If a job fails, its log is printed at the end of the run. |

The `.gz` files are gzipped: read them with `zless`, or directly with `data.table::fread()` in R.

### Intermediate files that are kept

These are kept so the final tables can be regenerated, or to help with troubleshooting. Most other intermediate files are deleted when the run finishes. The last step of the pipeline can be rerun with --forcerun make_metadata_table in the snakemake call; all inputs for this step are retained here.

- **Inputs to the final metadata step**:
  - [prefix]_max_score_per_peak_all_samples.bed
  - [prefix]_combined_tranferred_peaks_all_samples.bed
  - [prefix]_All_consensus_peaks_moved_to_SecondaryRef.bed (only with a secondary reference)
  - [prefix]_orig_called_peaks_intersect_transferred_consensus_peaks.bed
  - [prefix]_all_orig_surjected_peaks_intersect_consensus_peaks.bed
  - [prefix]_graph_prefixed_input.tsv
  - [prefix]_sample_ids.tsv (each sample_id with its sample_name and graph name)
  - [prefix]_Ref_chromosomes.txt
  - temp_previous_peak_files/[prefix]_previous_peaks_[assembly].tsv, one per assembly
- **Consensus peaks**: [prefix]_Ref_shared_peaks.bed and _6col.bed, [prefix]_NonRef_shared_peaks.bed and _6col.bed, and [prefix]_ALL_consensus_peaks_combined.bed
- **Original peaks and where they were transferred**: [prefix]_wHapsAndPeakIDs.bed.gz, the [prefix]_wHapsAndPeakIDs_w_surject_info*.bed files, and [prefix]_list_of_surj_assemblies_in_order.txt (the order assemblies were tried in)
- **Graph** (only if the pipeline built it): graph/[prefix].full.gbz, graph/[prefix].full_all_chrs.gbz (the graph used, with the contigs Minigraph-Cactus left out added back). Various other Minigraph-Cactus outputs also remain in graph folder.

### [prefix]_final_tbl_w_metadata.tsv.gz

All columns are described below, with the most commonly used columns starred (*).

\***sample_id**: **sample_name** from the input file + "_" + **hap1/2_in_graph** from input file. This sample_id describes a sample haplotype pair. e.g. SMHT004_liver_SMHT004_1 and SMHT004_liver_SMHT004_2. Each line in the input file produces two sample_id's. 

\***consensus_peak_id**: consensus peak region identifier. This is built from the consensus peak coordinates, and thus indicates the assembly on which that consensus peak region was called. Additional metadata for each consensus_peak_id, such as the assembly it was called on, reference chromosome it is associated with, etc., is recorded in **[prefix]_consensus_peak_summary.tsv**. **All possible consensus_peak_id - sample_id combinations are represented in this output file, whether or not accessibility data is available. See downstream analysis vignettes for filtering suggestions.**

\***coverage**: number of reads overlapping this consensus peak region at the point of maximum accessibility **in native assembly context**. NA indicates that this consensus_peak_id has no homologous region in this assembly.  

\***fire_coverage**: number of reads with FIRE elements overlapping this consensus peak region at the point of maximum accessibility. NA indicates that this consensus_peak_id has no homologous region in this assembly.

**score**:	an accessibility score at this point calculated during original FIRE peak calling and defined [here](https://fiberseq.github.io/fire/methods/aggregation.html). (Where a consensus peak was not able to be transferred to the assembly because there is no homologous region, the score column will have a value of -2.) I rarely use this column in downstream analyses, instead focusing on the raw fire_coverage/coverage calculation (in the fire_cov_OV_cov), though the -2 designation is used in downstream plotting.

**cons_chr**, **cons_start**, **cons_end**: consensus_peak_id, split into coordinates.

\***fire_cov_OV_cov**: raw accessibility value calculated as fire_coverage/coverage

**asm**: the assembly haplotype of this sample_id, as its hap subtag from the input file (e.g. SMHT004#1#)

**asm_chr**, **asm_start**, **asm_end**: coordinates of the consensus peak transferred to its native assembly space. This is the genomic region from which we pull the 'coverage' and 'fire_coverage' values. NA if the consensus peak could not be transferred to this assembly.

**asm_peak_id**: the three previous columns, combined

**secondary_ref_chr_[name]**, **secondary_ref_start_[name]**, **secondary_ref_end_[name]**: consensus coordinates transferred to the secondary reference, e.g. secondary_ref_chr_GRCh38. Useful for overlapping with annotations defined in that reference. NA means this consensus region has no homology in the secondary reference, or its length changed by more than Peak_length_change_param (1 = 100%) when transferred to this reference. These columns are left out if no secondary reference is given.

**exists_in_secondary_ref_[name]**: TRUE/FALSE. TRUE if values in the previous 3 columns are not NA

**overlaps_original_called_peak**: TRUE if there was at least one peak overlapping this consensus peak location in the original natively-called FIRE dataset. This information should be used with caution, as differences in coverage between samples can impact whether a peak rose above the FDR threshold in the original FIRE calling scheme.

**num_overlapping_orig_peaks**: The number of originally called FIRE peaks that overlap this consensus peak region. This can be informative for identifying structural differences that resulted in multiple distinct peaks being transferred to a single location on the assembly where the consensus peak was defined. This can happen, for example, when a duplication event results in multiple regions in one assembly being associated with a single region in another. Because the main purpose of this column is to flag these events, this number is corrected to only count as distinct peaks that are far enough apart from one another to plausibly arise from this type of event. I would use information in this column for flagging potential regions that require further analysis, but not as a clear indicator of the presence/absence of duplicated event.

**orig_peaks_overlapping_consensus**: Positions of originally called peaks which overlap this consensus peak region, each as chrom_start_end-[sample_name] and separated by commas. This can be useful for troubleshooting/identifying regions where structural differences may complicate direct comparison. *This set is NOT filtered for very nearby peaks as is done for calculating '**num_overlapping_orig_peaks**'; thus there are occasionally more peaks listed here than the number found in that column.


### [prefix]_consensus_peak_summary.tsv

One row per consensus peak.

**consensus_peak_id**: as above

**consensus_peak_asm**: the assembly the consensus peak was called on: the main reference (main_reference_name), or the assembly's name in the graph (hap1_in_graph or hap2_in_graph from the input file, e.g. SMHT004_1)

**chr**: the main reference chromosome the peak is on. Consensus peaks which were not called in main reference assembly space can in most cases still be associated with a reference chromosome based on the chromosome designation of nearby reference peaks. "AMB" indicates a reference chromosome could not be determined, and occurs when a contig does not align to the main reference at all during graph building.

**exists_in_main_ref_[name]**: TRUE if the consensus peak was called in main reference space. e.g. exists_in_main_ref_CHM13

**main_ref_chr_[name]**, **main_ref_start_[name]**, **main_ref_end_[name]**: its main reference coordinates (NA if not called on the main reference)

**exists_in_secondary_ref_[name]**: TRUE if the consensus peak could be transferred to the secondary reference (only with a secondary reference)

**secondary_ref_chr_[name]**, **secondary_ref_start_[name]**, **secondary_ref_end_[name]**: its secondary reference coordinates (NA if it couldn't be transferred; only with a secondary reference)

\***chr_order**: the peak's position along its chromosome, numbered 1 to the number of peaks on that chromosome, with no two peaks sharing a number. Main reference peaks are in reference order, and each non-reference peak is placed between its closest neighboring reference peaks. Non-reference peaks in the same gap are ordered by assembly, then contig, then position along the contig (reversed where the contig is oriented opposite to reference in the native assembly). NA for AMB peaks.

**chr_slot**: the gap chr_order is built from. The k-th reference peak on the chromosome has chr_slot 2k; non-reference peaks after it have 2k+1, those before the first reference peak have 1. Non-reference peaks in the same gap share a chr_slot, and odd numbers with no peaks are skipped by chr_order, so chr_slot is usually larger than chr_order. NA for AMB peaks.

**at_start_of_chr**, **at_end_of_chr**: TRUE for non-reference peaks placed before the chromosome's first reference peak, or after its last. NA for AMB peaks.

### [prefix]_transferred_cons_peaks_native_sequences_from_main_output.bed.gz

The underlying sequence of each consensus peak region on each assembly. bgzipped, with no header; columns: asm_chr, asm_start, asm_end, consensus_peak_id, seq. Join it to the final table on asm_chr, asm_start, asm_end (and consensus_peak_id).

This can be useful in order to determine whether differences in accessibility are mediated by sequence changes. Keep in mind that the way assemblies are constructed may result in unsynchronized sequence complementarity. Thus, one must consider the reverse complement of sequences in order to compare between them. See downstream analysis vignette for example code.

## Downstream analysis

(example_data/ has a small subset of a multiple reference assembly run (CHM13 main reference, GRCh38 secondary reference) on 45 HPRC Fiber-seq runs: 100 consensus peaks around HLA-DRB1 and HLA-DQB1, for all 90 sample haplotypes, with the prefix HPRC_HLA_example. Use it to try the scripts and vignette below without running the pipeline; example_data/README.md has the commands.)

There are several common next steps we may want to perform with these outputs:
- Call chromatin accessibility peaks
- Identify haplotype-selective peaks
- Generate a matrix plot of chromatin accessibility for all peaks within a genomic region for all samples, similar to the one above.

See sections below for tools to do these. The scripts in **downstream/** start from the pipeline outputs. They need R with data.table, Biostrings and ggplot2 (plus dplyr and tidyr for plotting); the pipeline's R environment (workflow/envs/R_env.yml) has them.

### Calling peaks and haplotype-selective peaks

**This script is written for processing multi-reference runs and will not work for the shared reference mode as is. However, peak calling and hap selectivity filters from this script (and described below) can be applied to shared reference runs.**

[docs/vignettes/vignette_process_main_output_and_consensus_peak_summary.Rmd](docs/vignettes/vignette_process_main_output_and_consensus_peak_summary.Rmd) walks through the steps of process_consensus_peak_output.R with explanations, with each cutoff as a parameter. Render it with your own outputs:                                          

```r
rmarkdown::render("docs/vignettes/vignette_process_main_output_and_consensus_peak_summary.Rmd",
      params = list(input_dir = "/path/to/output_dir/prefix/", prefix = "prefix"))
```

To try it on the example data use `params = list(input_dir = "../../example_data/", prefix = "HPRC_HLA_example")`. The vignette's code runs from docs/vignettes/, so a relative input_dir is relative to that folder (or give a full path). Rendering needs pandoc, which RStudio includes.

This pipeline can be run from the command line as shown below, **but I do not recommend this for new datasets** as parameter choices are dataset specific. Instead, use the vignette above to run script step by step, and determine appropriate parameters based on the explanations there.

```bash
Rscript downstream/process_consensus_peak_output.R input_dir output_dir prefix [exclude_runs]
```
**input_dir** is the pipeline's [output_dir]/[prefix]/. **exclude_runs** is optional: a comma-separated list, and any sample_id containing one of them is left out. The cutoffs are set at the top of the script:

- By default, a consensus peak is a **peak** (accessible chromatin) in a sample haplotype if fire_coverage ≥ 4 and fire_cov_OV_cov ≥ 0.2. Consensus peaks that are a peak in no sample haplotype are dropped.
- A peak is **tested** for haplotype selectivity in a run if both haplotypes have coverage ≥ 10 and their fire_cov_OV_cov differ by ≥ 0.25. Tested peaks are compared with a Fisher's exact test (FIRE vs. non-FIRE reads on each haplotype) and are called **haplotype-selective** via this nominally significant value. Only autosomal peaks are tested (not chrX, chrY, chrM or AMB).
- **seqs_match** is TRUE if the peak's sequence is the same on both haplotypes (or the reverse complement), so a haplotype-selective peak may be driven by the sequence change if FALSE.

### Plotting chromatin actuation matrix

This script takes in a genomic interval or set of genomic intervals and generates plots with all intervening peaks as columns, all samples haplotypes as rows, with box color denoting relative fire_cov_OV_cov. Because we have ordered non-reference peaks relative to reference peaks in the consensus peak summary file, non-reference peaks will also appear on this plot if they exist in the interval provided. It outputs two plots per interval, one with samples in name order and one with samples clustered by their accessibility and missing peaks (_clustered.png)

#### Inputs:

- **output_dir**: whatever you want, with trailing slash
- **prefix**: also whatever you want. Only used for output file name
- **final_tbl_w_metadata**: main output from consensus peak calling pipeline. Can used processed version from the downstream analysis script above.
- **consensus_peak_summary**: other main output from consensus peak calling pipeline. Can used processed version from the downstream analysis script above.
- **intervals.bed**: denotes the intervals to plot **in main reference coordinates**. See sample input file example_data/HPRC_HLA_example_intervals.bed. This is a tab-delimited file with columns:
chr, start, end, name, and TSS (currently no header, but should change this)
name: a name for your interval (e.g. gene name); these are useful when multiple plots are being generated of different intervals.
TSS: Optional (NA if not providing). This is a single genomic coordinate within the region provided (e.g. TSS position) or comma-separate list of coordinates (no spaces). If this coordinate falls within in consensus peak, that peak will be marked with a heavy black outline.
- **coverage_cutoff**: I default to 6 here. Tiles with a coverage below this are gray, to indicate that the actuation information at that position is likely noisy.
- **remove_single_sample_peaks**: (TRUE/FALSE) If a **genomic region** where a consensus peak is called is only present in one sample, all other boxes in that column will be black. You can choose to remove these columns. (Note, this is NOT the same as the case where a peak is called in just one sample, but homologous genomic loci without peaks exists in other sample; these will not be removed).
- **min_frac_actuation**: Optional. I default to 0.1. This removes any column where the highest chromatin actuation is below this value.

```bash
Rscript downstream/plot_consensus_peaks_by_interval.R output_dir prefix final_tbl_w_metadata consensus_peak_summary intervals.bed coverage_cutoff remove_single_sample_peaks [min_frac_actuation]
```

#### Outputs:
Two plots per input interval (one clustered, one not) are generated via this script, as well as the data subset from which the plot was generated. 

Interpreting the plot:
- **Purple**: chromatin accessibility. Darker = more accessible.
- **Gray**: Not enough data to evaluate (Total coverage at this site is below the coverage_cutoff parameter)
- **Black** consensus peak cannot be transferred to this assembly. This means a homologous genomic region doesn't exist in this assembly.
- **Black outline around column**: User provided coordinate falls within this peak. See TSS column explanation above
- **Star in box**: This indicates that more than one distinct original peak in this sample overlaps this consensus peak when transferred to the consensus assembly (there is a number >1 in the **num_overlapping_orig_peaks** column in the main output -- see how these are filtered above). This **may** indicate copy number differences. These are returned purely to alert us to locations where we need to take a closer look, as chromatin actuation values may not actually represent consistent copies between samples.



