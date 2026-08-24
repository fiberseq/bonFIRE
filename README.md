# bonFIRE

bonFIRE is a toolkit designed to simplify comparative analysis of Fiber-seq data across samples as well as enable **direct comparison of Fiber-seq data processed against multiple distinct genome assemblies.**

The input to bonFIRE is one or several [FIRE](https://fiberseq.github.io/fire/run.html) runs. See **Input files** for more information on generating these.

Some common use cases for bonFIRE are:

(1) Compare the chromatin landscape between samples processed against a common haploid reference (e.g. GRCh38).

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


## 
