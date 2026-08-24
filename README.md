# bonFIRE

bonFIRE is a toolkit designed to simplify comparative analysis of Fiber-seq data across samples as well as enable **direct comparison of Fiber-seq data processed against multiple distinct genome assemblies.**

The input to bonFIRE is one or several [FIRE](https://fiberseq.github.io/fire/run.html) runs. See **Input files** for more information on generating these.

### Some common use cases for bonFIRE are:

(1) Compare the chromatin landscape between samples processed against a **common haploid reference** (e.g. GRCh38).

(2) Compare the chromatin landscape between two haplotypes in a single sample, processed against its **own native diploid assembly**. We can perform this analysis for a single sample, or for a group of samples all processed against a common diploid assembly. These samples may be different tissues from one individual, different experimental conditions, etc. 

(3) Compare the chromatin landscape between samples which are processed in different assemblies. This enables us, for example, to evaluate the chromatin landscape of each sample in its own native donor assembly, and then combine them into a single dataset which **does not rely on a single reference**, and thus enabling comparison at loci that are only shared between a subset of assemblies. 

In order to perform inter-assembly comparisons, examples 2 and 3 make use of **pangenome graphs**, which inherently define homologies between assembly coordinates, and can be used to transfer accessibility peak calls between assemblies. This is described further below. 

See section titled **Consensus Peak Calling Approach** for an in depth overview of the method.

## Installation 

### 1. Clone the repository

Clone the repository and move into the project directory:

```bash
git clone https://github.com/minkinaa/bonFIRE.git
cd bonFIRE
```

### 2. Install Pixi

Please start by installing [pixi](https://pixi.sh/latest/) which handles the environment of this Snakemake workflow.

### 3. Install dependencies

You can then install the `pixi` environment by running the following in this repository:

```bash
pixi install
```

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

## Input files

A sample config can be found in 

## Consensus Peak Calling Approach

We begin by describing the simplest case, where all data is processed against a common haploid reference (see **Some common use cases for bonFIRE** above). The multi-assembly, graph-dependent approach builds upon the concepts illustrated below, while moving through different assembly spaces. 

### Consensus peak calling in a common haploid reference space

**This simplest case is a variation on the bonFIRE approach but exists as a separate repo**: [fire_consensus_pipeline](https://github.com/StergachisLab/fire_consensus_pipeline)  

The consensus FIRE peak calling approach enables comparative analysis of Fiber-seq chromatin accessibility between Fiber-seq datasets, and is the recommended approach for comparing individual Fiber-seq FIRE datasets.

The consensus peak calling approach has several steps. First, a set of 'consensus peaks' are generated using the individual Fiber-seq dataset FIRE peak calls. Second, the chromatin accessibility value at these 'consensus peaks' is quantified in each of the Fiber-seq datasets, enabling the direct comparison of chromatin accessibility at each position of interest. 

Individual Fiber-seq dataset FIRE peak calls serve as inputs for consensus peak calling, where a 'consensus peak region' is called at any position where a FIRE peak was called in at least one sample. Where multiple samples have overlapping peaks, a single consensus peak region is defined, using the median start and end positions of the underlying peaks. 

<img width="1021" height="377" alt="consensus_peak_calling_image1" src="https://github.com/user-attachments/assets/d9ccd9fc-473a-46aa-93bc-fd0d47886db2" />

Once this set of consensus peak regions is defined, we pull raw coverage (number of reads) and FIRE coverage (number of reads with FIRE elements) at each consensus peak region for each sample, reporting the position with the maximum FIRE coverage if it differs across the length of the region.  

The output is a table containing values shown below (see **Output** below for all column descriptions), containing these values for each sample at each consensus peak region. We can thus calculate a chromatin accessibility score, and define custom filters for samples and/or consensus peak regions based on minimum coverage, for example. Thus, peaks not called in some samples due to relatively low coverage, can be rescued with this approach, and chromatin accessibility can be compared quantitatively even at regions where a peak was not originally called in all samples.

<img width="774" height="284" alt="consensus_peak_calling_image2" src="https://github.com/user-attachments/assets/59320439-e7a8-47f2-bef6-4daf17227cc3" />

### Consensus peak calling across multiple genome assemblies

As with the above approach, the input for this pipeline are the peak outputs from the FIRE pipeline run separately on each sample. Even a single sample processed in its native diploid assembly space contains chromatin information for two haplotypes, and thus this approach can be used to compare the chromatin landscape between them. 

Conceptually, this approach mirrors the single assembly approach, but contains additional steps to first move peaks called across different assemblies into a common assembly space via a pangenome graph, followed by calling 'consensus peaks' in that reference space. While the majority of peaks called across assembly spaces are able to be transferred to a common reference (in this case, CHM13), the subset will not. These are peaks in genomic regions which have no homology to CHM13 and thus cannot be represented there. However, these regions may be homologous across multiple other assemblies, and thus we don't want them dropping out of our comparative analysis. Thus, the remaining peaks are transferred to a new assembly, and 'consensus peaks' are defined as above from all peaks that were able to be transferred to this assembly. Because some peaks will still remain as untransferrable to either CHM13 or the new assembly, we iteratively move remaining peaks to other assemblies present in the graph, and call consensus peaks, until every starting peak has contributed to a 'consensus peak'. This process culminates with a set of consensus peaks regions -- a set of genomic chr:start-end intervals defined across multiple assembly spaces.

As in the single reference approach above, we ultimately want to report raw coverage and FIRE coverage values, enabling us to calculate chromatin accessibility at each of these consensus peak regions for each sample. To accomplish this, we again use the graph to move the entire set of consensus peak regions back to each assembly. Using these assembly-specific coordinates, we can pull raw coverage and fire coverage data from the pileup.bed file generated during the original, sample-specific FIRE peak calling step. We merge these into a table containing chromatin accessibility value and metadata for each [sample haplotype]-[consensus peak region] pair, enabling direct chromatin comparison at any region where a peak was originally called in any sample.

<img width="1057" height="495" alt="consensus_peak_calling_image3" src="https://github.com/user-attachments/assets/24d56462-88b8-48ca-ac8d-40cf38f131a1" />

Many metadata columns are appended to this output to facilitate downstream analysis (see **Outputs**).

## Outputs

The main output from this pipeline is called **[sample]_final_tbl_w_metadata.tsv** . All columns are described below, with the most commonly used columns starred (*).

\***sample_id**: A sample id that distinguished between two haplotypes (e.g. TestSample_1, TestSample_2)

\***consensus_peak_id**: An id defining a consensus peak region. These are genomic coordinates of the consensus peak region, from which one can infer the assembly on which they were initially defined. 	

\***coverage**: number of reads overlapping this consensus peak region at the point of maximum accessibility

\***fire_coverage**: number of reads with FIRE elements overlapping this consensus peak region at the point of maximum accessibility

**score**:	an accessibility score at this point calculated during original FIRE peak calling and defined [here](https://fiberseq.github.io/fire/methods/aggregation.html). (Where a consensus peak was not able to be transferred to the assembly because there is no homologous region, the score column will have a value of -2.) I rarely use this column in downstream analyses, instead focusing on the raw fire_coverage/coverage calculation (in the fire_cov_OV_cov).

**cons_chr**, **cons_start**, **cons_end**: consensus_peak_id, split into coordinates.

\***fire_cov_OV_cov**: raw accessibility value calculated as fire_coverage/coverage

**asm**: assembly associated with this sample 	

**asm_chr**, asm_start, asm_end: coordinates of the consensus peak in the native assembly space. This is the genomic region from which we pull the 'coverage' and 'fire_coverage' values.

**asm_peak_id**: the three previous columns, combined 

**hg38_chr**, **hg38_start**, **hg38_end**: consensus coordinates transferred to GRCh38. This information is useful for overlapping with annotations defined in GRCh38 coordinates. NAs in these columns mean that this consensus region has no homology in GRCh38.	

**exists_in_hg38**: TRUE/FALSE. TRUE is values in the previous 3 columns are not NA

**overlaps_original_called_peak**: TRUE if there was at least one peak overlapping this consensus peak location in the original FIRE dataset for this sample. This information should be used with caution, as differences in coverage between samples can impact whether a peak rose above the FDR threshold in the original FIRE calling scheme. 

**num_overlapping_orig_peaks**: The number of originally called FIRE peaks that	overlap this consensus peak region. This can be informative for identifying structural differences that resulted in multiple distinct peaks being transferred to a single location on the assembly where the consensus peak was defined. This can happen, for example, when a duplication event results in multiple regions in one assembly being associated with a single region in another. Because the main purpose of this column is to flag these events, this number is corrected to only count as distinct peaks that are far enough apart from one another to plausibly arise from this type of event. 

**orig_peaks_overlapping_consensus**: Positions of originally called peaks which overlap this consensus peak region. This can be useful for troubleshooting/identifying regions where structural differences may complicate direct comparison. *This set is NOT filtered for very nearby peaks as is done for calculating '**num_overlapping_orig_peaks**'; thus there are occasionally more peaks listed here than the number found in that column. 

**seq**: The underlying sequence of this consensus peak region. This can be useful in order to determine whether differences in accessibility are mediated by sequence changes. Keep in mind that the way assemblies are constructed may result in unsynchronized sequence complementarity. Thus, one must consider the reverse complement of sequences in order to compare between them. For an example of how to do this see this vignette (**TO DO**).

