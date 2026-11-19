## Shared-reference mode: every run was processed against the same haploid reference (main_reference_fasta),
## so nothing is moved between assemblies. Consensus peaks are called on that reference from all runs'
## peaks, and each run's coverage and FIRE coverage (overall and per haplotype) are read from its pileup at
## each consensus peak. Included by the Snakefile when shared_ref is set; uses the common setup there.

graph_only_keys = ["graph", "ft_param", "Peak_length_change_param", "main_reference_name", "secondary_reference_fasta",
                   "secondary_reference_name", "hg38_fasta_path", "split_graph", "graph_mem", "graph_runtime",
                   "vg_mem_mb", "vg_threads", "threads"]
ignored_keys = [k for k in graph_only_keys if k in config]
if ignored_keys:
    print(f"shared_ref: ignoring these config keys, which only apply to runs on their own assemblies: {ignored_keys}", file=sys.stderr)

runs = input_tbl["sample_name"].tolist()
duplicated_runs = sorted({r for r in runs if runs.count(r) > 1})
if duplicated_runs:
    sys.exit(f"sample_name must be unique; these appear more than once: {duplicated_runs}")

wildcard_constraints:
    run="|".join(re.escape(r) for r in runs)

# run names in input order, for the metadata step
runs_file = f"{global_output_dir}/{sample}/{sample}_runs.txt"
write_tsv_if_changed(pd.DataFrame({"run": runs}), runs_file)

All_Targets = [f"{global_output_dir}/{sample}/{sample}_final_tbl_w_metadata.tsv.gz",
               f"{global_output_dir}/{sample}/{sample}_consensus_peak_summary.tsv"]

rule list_Ref_chromosomes:
    input:
        ref_fasta=MAIN_REF_FASTA
    output:
        chroms=f"{global_output_dir}/{sample}/{sample}_Ref_chromosomes.txt"
    log:
        f"{log_dir}/list_Ref_chromosomes.log"
    shell:
        """
        exec > {log} 2>&1
        # the name is the first word of each header line; zcat -f reads gzipped and plain fastas
        zcat -f {input.ref_fasta} | grep '^>' | sed -e 's/^>//' -e 's/[[:space:]].*//' > {output.chroms}
        """

# the original peaks, with peak_id chrom_start_end-run as in graph mode. Peaks on contigs that aren't in
# the reference are dropped (and counted in the log); a run with no peaks on the reference stops the run,
# since it was most likely processed against a different reference
rule make_bed_with_peakIDs:
    input:
        input_bed=rules.combine_peaks.output.combined_peaks_bed,
        ref_chroms=rules.list_Ref_chromosomes.output.chroms
    output:
        output_bed=f"{global_output_dir}/{sample}/{sample}_wPeakIDs.bed.gz",
        output_bed_unzipped=temp(f"{global_output_dir}/{sample}/{sample}_wPeakIDs.bed")
    conda:
        "../envs/env.yml"
    log:
        f"{log_dir}/make_bed_with_peakIDs.log"
    shell:
        """
        exec > {log} 2>&1
        gzip -cd {input.input_bed} | awk -F'\\t' -v OFS='\\t' '
            NR == FNR {{ ref[$1] = 1; next }}
            FNR == 1 {{ print "#chrom", "start", "end", "sample", "peak_id" }}
            {{
                runs[$4] = 1
                if (!($1 in ref)) {{ off_ref[$4 " on " $1]++; next }}
                on_ref[$4]++
                peak_id = $1"_"$2"_"$3"-"$4
                if (seen[peak_id]++) {{ n_dup++; next }}
                print $1, $2, $3, $4, peak_id
            }}
            END {{
                if (n_dup) print "Removed " n_dup " duplicate peaks (same run and coordinates)" > "/dev/stderr"
                for (k in off_ref) print "Removed " off_ref[k] " peaks of " k ", which is not in the reference" > "/dev/stderr"
                bad = 0
                for (r in runs) if (!(r in on_ref)) {{ print "Run " r " has no peaks on the reference contigs; was it processed against this reference?" > "/dev/stderr"; bad = 1 }}
                if (bad) exit 1
            }}' {input.ref_chroms} - > {output.output_bed_unzipped}
        bgzip -c {output.output_bed_unzipped} > {output.output_bed}
        """

# consensus peaks from all runs' peaks, the same way as in graph mode: each original peak becomes a mock
# FIRE read and ft call-peaks finds where they pile up
rule call_consensus_peaks:
    input:
        bed=rules.make_bed_with_peakIDs.output.output_bed
    output:
        bed_4col=temp(f"{global_output_dir}/{sample}/{sample}_wPeakIDs_4col.bed"),
        mock_bam=temp(f"{global_output_dir}/{sample}/{sample}_wPeakIDs_4col_mockbam.bam"),
        ref_peaks=f"{global_output_dir}/{sample}/{sample}_Ref_shared_peaks.bed",
        ref_peaks_5col=temp(f"{global_output_dir}/{sample}/{sample}_Ref_shared_peaks_5col.bed"),
        ref_peaks_5col_w_peak_id=f"{global_output_dir}/{sample}/{sample}_Ref_shared_peaks_6col.bed",
        ref_peaks_5col_w_peak_id_zipped=temp(f"{global_output_dir}/{sample}/{sample}_Ref_shared_peaks_6col.bed.gz")
    conda:
        "../envs/env.yml"
    threads: 4
    resources:
        mem_mb=mem_for(32000),  # 13 GB for 10 runs on GRCh38
        runtime=60*12
    log:
        f"{log_dir}/call_consensus_peaks.log"
    shell:
        """
        exec > {log} 2>&1
        rm -f {global_output_dir}/{sample}/{sample}_wPeakIDs_4col_mockbam.bam.tmp*
        gzip -cd {input.bed} | awk 'NR==1 {{print "#chrom\\tstart\\tend\\tsample"; next}} {{print $1, $2, $3, $5}}' OFS="\\t" > {output.bed_4col}
        {PATH_TO_FT} mock-fire -t {threads} {output.bed_4col} | samtools sort -@ {threads} -m 2G -o {output.mock_bam} --write-index
        {PATH_TO_FT} call-peaks -t {threads} {output.mock_bam} --min-fire-coverage 1 --min-fire-frac 0.01 > {output.ref_peaks}
        cat {output.ref_peaks} | cut -f1,2,3 > {output.ref_peaks_5col}

        bash {WORKFLOW_DIR}/scripts/add_peak_id_col_to_bed_file.sh \
            {output.ref_peaks_5col} \
            {output.ref_peaks_5col_w_peak_id} \
            {output.ref_peaks_5col_w_peak_id_zipped}
        """

# each run's values at each consensus peak: the pileup columns are found by name (the haplotype ones are
# NA if the pileup has none), the rows overlapping a consensus peak kept, and one row per peak chosen
rule pileup_values_at_consensus_peaks:
    input:
        pileup=get_pileup,
        consensus_peaks=rules.call_consensus_peaks.output.ref_peaks_5col_w_peak_id
    output:
        values=temp(f"{global_output_dir}/{sample}/run_pileup_counts/{{run}}/{{run}}_max_score_per_peak.tsv")
    conda:
        "../envs/env.yml"
    resources:
        mem_mb=mem_for(4000),  # under 1 GB; 35-60 min for a 6 GB pileup
        runtime=runtime_for(180)
    log:
        f"{log_dir}/pileup_values_at_consensus_peaks_{{run}}.log"
    shell:
        """
        exec > {log} 2>&1
        set -o pipefail
        zcat {input.pileup} | awk -F'\\t' -v OFS='\\t' '
            NR == 1 {{
                if ($1 !~ /^#/) {{ print "the pileup has no header line" > "/dev/stderr"; exit 1 }}
                for (i = 1; i <= NF; i++) col[$i] = i
                n = split("coverage fire_coverage score coverage_H1 fire_coverage_H1 score_H1 coverage_H2 fire_coverage_H2 score_H2", want, " ")
                for (k = 1; k <= 3; k++) if (!(want[k] in col)) {{ print "the pileup has no " want[k] " column" > "/dev/stderr"; exit 1 }}
                next
            }}
            {{
                out = $1 OFS $2 OFS $3
                for (k = 1; k <= n; k++) out = out OFS ((want[k] in col) ? $(col[want[k]]) : "NA")
                print out
            }}' \
        | bedtools intersect -a stdin -b {input.consensus_peaks} -wa -wb \
        | python3 {WORKFLOW_DIR}/scripts/shared_ref_max_score_per_consensus_peak.py {input.consensus_peaks} {wildcards.run} {output.values}
        """

rule combine_max_score_outputs:
    input:
        expand(f"{global_output_dir}/{sample}/run_pileup_counts/{{run}}/{{run}}_max_score_per_peak.tsv", run=runs)
    output:
        combined_out=f"{global_output_dir}/{sample}/{sample}_max_score_per_peak_all_samples.tsv"
    log:
        f"{log_dir}/combine_max_score_outputs.log"
    shell:
        """
        exec > {log} 2>&1
        head -n 1 {input[0]} > {output.combined_out}
        for f in {input}; do tail -n +2 $f; done >> {output.combined_out}
        """

# each run's original peaks on each consensus peak
rule intersect_original_and_consensus_peaks:
    input:
        orig_bed=rules.make_bed_with_peakIDs.output.output_bed,
        consensus_peaks=rules.call_consensus_peaks.output.ref_peaks_5col_w_peak_id
    output:
        intersect_out=f"{global_output_dir}/{sample}/{sample}_orig_called_peaks_intersect_consensus_peaks.bed"
    conda:
        "../envs/env.yml"
    resources:
        mem_mb=mem_for(16000),
        runtime=60*2
    log:
        f"{log_dir}/intersect_original_and_consensus_peaks.log"
    shell:
        """
        exec > {log} 2>&1
        bedtools intersect -a <(zcat {input.orig_bed} | tail -n +2) -b <(tail -n +2 {input.consensus_peaks}) -wa -wb > {output.intersect_out}
        """

rule get_consensus_peak_sequences:
    input:
        consensus_peaks=rules.call_consensus_peaks.output.ref_peaks_5col_w_peak_id,
        fasta=MAIN_REF_FASTA_PLAIN,
        # a gzipped reference's decompressed copy is indexed when it's made; a plain one is indexed by bedtools if needed
        fai=[f"{MAIN_REF_FASTA_PLAIN}.fai"] if MAIN_REF_FASTA.endswith(".gz") else []
    output:
        seqs=temp(f"{global_output_dir}/{sample}/{sample}_consensus_peak_sequences.tsv")
    conda:
        "../envs/env.yml"
    resources:
        mem_mb=mem_for(16000),
        runtime=60
    log:
        f"{log_dir}/get_consensus_peak_sequences.log"
    shell:
        """
        exec > {log} 2>&1
        bedtools getfasta -bedOut -tab -name -fi {input.fasta} -bed <(tail -n +2 {input.consensus_peaks}) > {output.seqs}
        """

rule make_metadata_table:
    input:
        consensus_peaks=rules.call_consensus_peaks.output.ref_peaks_5col_w_peak_id,
        values=rules.combine_max_score_outputs.output.combined_out,
        orig_intersect=rules.intersect_original_and_consensus_peaks.output.intersect_out,
        seqs=rules.get_consensus_peak_sequences.output.seqs,
        runs=runs_file
    output:
        meta_tbl=f"{global_output_dir}/{sample}/{sample}_final_tbl_w_metadata.tsv.gz",
        peak_summary=f"{global_output_dir}/{sample}/{sample}_consensus_peak_summary.tsv"
    conda:
        "../envs/R_env.yml"
    threads: 4
    resources:
        mem_mb=mem_for(16000),  # 3 GB for 10 runs and 431,000 consensus peaks
        runtime=runtime_for(60*2)
    log:
        f"{log_dir}/make_metadata_table.log"
    shell:
        """
        exec > {log} 2>&1
        Rscript {WORKFLOW_DIR}/scripts/shared_ref_make_metadata_table.R \
            {global_output_dir}/{sample}/ \
            {sample} \
            {input.consensus_peaks} \
            {input.values} \
            {input.orig_intersect} \
            {input.seqs} \
            {input.runs} \
            {threads}
        """
