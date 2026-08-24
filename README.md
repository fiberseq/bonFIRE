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






## 
