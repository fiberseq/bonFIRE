## Graph mode: each run is on its own assembly; peaks move between assemblies through a pangenome graph.
## Included by the Snakefile unless shared_ref is set; uses the common setup there.

graphs_dir=f"{global_output_dir}/{sample}/graph"
graph=config.get("graph", "NA")
FT_PARAM = int(config.get("ft_param", 1000))
# The main reference: peaks are moved onto it first, and consensus peaks not on it are placed relative
# to it. The optional secondary reference only gets each consensus peak's coordinates. Both are haploid,
# and their names are their path names in the graph. The older CHM13_fasta_path / hg38_fasta_path keys
# still work, as references named CHM13 and GRCh38.
MAIN_REF_NAME = config.get("main_reference_name", None if "main_reference_fasta" in config else "CHM13")
SECONDARY_REF_FASTA = config.get("secondary_reference_fasta", config.get("hg38_fasta_path"))
if SECONDARY_REF_FASTA in ("", "NA"):
    SECONDARY_REF_FASTA = None
SECONDARY_REF_NAME = config.get("secondary_reference_name", None if "secondary_reference_fasta" in config else "GRCh38")
HAS_SECONDARY_REF = SECONDARY_REF_FASTA is not None
if not HAS_SECONDARY_REF:
    SECONDARY_REF_NAME = None
if not MAIN_REF_NAME:
    sys.exit("Set main_reference_name in the config")
if HAS_SECONDARY_REF and not SECONDARY_REF_NAME:
    sys.exit("Set secondary_reference_name in the config (or remove secondary_reference_fasta)")
REF_NAMES = [MAIN_REF_NAME] + ([SECONDARY_REF_NAME] if HAS_SECONDARY_REF else [])
# for shell commands
REF_NAMES_STR = " ".join(REF_NAMES)
length_change_param=float(config.get("Peak_length_change_param", 1))
TEMP_THREADS=int(config.get("threads", 32))
SPLIT_GRAPH=config.get("split_graph", False)
GRAPH_MEM=config.get("graph_mem", 300000)
GRAPH_RUNTIME=config.get("graph_runtime", 1440)
# memory for the vg steps; depends on the graph (~11 GB for SMHT043, ~180 GB for an 82-sample HPRC graph)
VG_MEM_MB=int(config.get("vg_mem_mb", 80000))
VG_THREADS=int(config.get("vg_threads", 8))

graph_cols = ["hap1_orig_subtag", "hap2_orig_subtag", "hap1_in_graph", "hap2_in_graph", "hap1_fasta_path", "hap2_fasta_path"]
missing_cols = [c for c in graph_cols if c not in input_tbl.columns]
if missing_cols:
    sys.exit(f"The input file is missing these columns, which runs on their own assemblies need: {missing_cols} (set shared_ref: True if every run is on the same reference)")

## make combined fasta and put and graph_prefixes_tbl
hap1_fasta_tbl = input_tbl[
    ["hap1_orig_subtag","hap1_in_graph", "hap1_fasta_path"]
].drop_duplicates()

hap2_fasta_tbl = input_tbl[
    ["hap2_orig_subtag","hap2_in_graph", "hap2_fasta_path"]
].drop_duplicates()

hap1_fasta_tbl.columns = ["hap_subtag","graph_prefix", "fasta_path"]
hap2_fasta_tbl.columns = ["hap_subtag","graph_prefix", "fasta_path"]

combined_fasta_tbl = pd.concat(
    [hap1_fasta_tbl, hap2_fasta_tbl],
    ignore_index=True
).drop_duplicates()

ALL_FASTAS=combined_fasta_tbl["fasta_path"].tolist()

## check the reference names: they become graph path names, column names and file names, and contigs and
## peaks are told apart from the assemblies' by name (often by substring), so they must not overlap
asm_names = set(combined_fasta_tbl["graph_prefix"]) | set(combined_fasta_tbl["hap_subtag"])
for ref_name in REF_NAMES:
    if not re.fullmatch(r"[A-Za-z0-9_]+", ref_name):
        sys.exit(f"Reference name '{ref_name}' can only contain letters, numbers and '_'")
    if ref_name in ("NONE", "AMB", "AT_START", "AT_END", "NA", "Ref", "NonRef", "SecondaryRef"):
        sys.exit(f"Reference name '{ref_name}' is a word the pipeline uses for something else; choose another")
    clashes = [a for a in asm_names if ref_name in a or a in ref_name]
    if clashes:
        sys.exit(f"Reference name '{ref_name}' overlaps these assembly names or hap subtags: {clashes}")
if HAS_SECONDARY_REF and SECONDARY_REF_NAME == MAIN_REF_NAME:
    sys.exit("main_reference_name and secondary_reference_name must differ")

## graph input: the main reference first, then the assemblies, then the secondary reference if there is one
combined_fasta_tbl_w_refs = pd.concat(
    [pd.DataFrame({'graph_prefix': [MAIN_REF_NAME], 'fasta_path': [MAIN_REF_FASTA]}),
     combined_fasta_tbl[["graph_prefix", "fasta_path"]]] +
    ([pd.DataFrame({'graph_prefix': [SECONDARY_REF_NAME], 'fasta_path': [SECONDARY_REF_FASTA]})] if HAS_SECONDARY_REF else []),
    ignore_index=True
)

graph_input_file=f"{global_output_dir}/{sample}/{sample}_graph_input.tsv"
write_tsv_if_changed(combined_fasta_tbl_w_refs, graph_input_file)

graph_prefixes_input = f"{global_output_dir}/{sample}/{sample}_graph_prefixed_input.tsv"
write_tsv_if_changed(combined_fasta_tbl, graph_prefixes_input)

## wildcard definitions
run_wildcards=input_tbl["sample_name"].values.tolist()
asm_wildcards = input_tbl[["hap1_in_graph", "hap2_in_graph"]].values.flatten().tolist()

wildcard_constraints:
    asm="|".join(asm_wildcards),
    run="|".join(run_wildcards)

# each assembly once, in input order (a donor with several runs is listed more than once)
unique_asms = list(dict.fromkeys(asm_wildcards))

assemblies_surjected_onto_file=f"{global_output_dir}/{sample}/{sample}_list_of_surj_assemblies_in_order.txt"

pairs_1_3 = list(zip(input_tbl["sample_name"], input_tbl["hap1_in_graph"]))
pairs_1_4 = list(zip(input_tbl["sample_name"], input_tbl["hap2_in_graph"]))

all_sample_pileup_pairs = pairs_1_3 + pairs_1_4

# each run's sample_id for each of its assembly haplotypes is run_graphAsm (see process_intersection_file);
# the metadata step looks the run, assembly and haplotype (1 for hap1_in_graph, 2 for hap2_in_graph) up here
# instead of splitting the name, so run and graph names can be anything
sample_ids_tbl = pd.DataFrame(
    [(f"{run}_{asm}", run, asm, hap) for pairs, hap in ((pairs_1_3, 1), (pairs_1_4, 2)) for run, asm in pairs],
    columns=["sample_id", "run", "graph_asm", "haplotype"]
).drop_duplicates()
ambiguous = sample_ids_tbl[sample_ids_tbl["sample_id"].duplicated(keep=False)]
if len(ambiguous) > 0:
    sys.exit("These run and graph name pairs give the same sample_id (run_graphName), so they can't be told apart:\n" + ambiguous.to_string(index=False))
sample_ids_file = f"{global_output_dir}/{sample}/{sample}_sample_ids.tsv"
write_tsv_if_changed(sample_ids_tbl, sample_ids_file)

## Old stuff

Ref_global_peaks_surjected_onto_all_asm = [f"{global_output_dir}/{sample}/Ref_global_peaks_surjected_onto_asms/{sample}_Ref_global_peaks_surjected_onto_{asm}.bed" for asm in unique_asms]
NonRef_global_peaks_surjected_onto_all_asm = [f"{global_output_dir}/{sample}/NonRef_global_peaks_surjected_onto_asms/{sample}_NonRef_global_peaks_surjected_onto_{asm}.bed" for asm in unique_asms]
combined_global_peaks_surjected_onto_all_asm = [f"{global_output_dir}/{sample}/Combined_global_peaks_surjected_onto_asms/{sample}_ALL_global_peaks_surjected_onto_{asm}.bed" for asm in unique_asms]

run_modified_pileup=[f"{global_output_dir}/{sample}/run_pileup_counts/{run}_{asm}/{run}_{asm}_modified_pileup.bed" for run, asm in all_sample_pileup_pairs]

overlap_outs=[f"{global_output_dir}/{sample}/run_pileup_counts/{run}_{asm}/{run}_{asm}_pileup_intervals_overlapping_peaks.bed" for run, asm in all_sample_pileup_pairs]

max_score_per_consensus_peak=[f"{global_output_dir}/{sample}/run_pileup_counts/{run}_{asm}/{run}_{asm}_max_score_per_peak.bed" for run, asm in all_sample_pileup_pairs]
#print(max_score_per_consensus_peak)

Target=[Ref_global_peaks_surjected_onto_all_asm, NonRef_global_peaks_surjected_onto_all_asm]

if graph=="NA":
    graph = f"{graphs_dir}/{sample}.full_all_chrs.gbz"

if GRAPH_ONLY == True:
    All_Targets = f"{graphs_dir}/{sample}.full_all_chrs.gbz"
else:
    All_Targets = [f"{global_output_dir}/{sample}/{sample}_final_tbl_w_metadata.tsv.gz",
                   f"{global_output_dir}/{sample}/additional_outputs/{sample}_previous_peaks_nonRef_chrs.tsv",
                   f"{global_output_dir}/{sample}/additional_outputs/{sample}_previous_peak_checks.tsv",
                   f"{global_output_dir}/{sample}/{sample}_consensus_peak_summary.tsv",
                   f"{global_output_dir}/{sample}/{sample}_transferred_cons_peaks_native_sequences_from_main_output.bed.gz"]

rule make_bed_with_haps_and_peakIDs:
    input:
        input_bed=rules.combine_peaks.output.combined_peaks_bed,
        graph_prefixes_input=graph_prefixes_input
    output:
        output_bed=f"{global_output_dir}/{sample}/{sample}_wHapsAndPeakIDs.bed.gz",
        output_bed_unzipped=temp(f"{global_output_dir}/{sample}/{sample}_wHapsAndPeakIDs.bed")
    conda:
        "../envs/env.yml"
    log:
        f"{log_dir}/make_bed_with_haps_and_peakIDs.log"
    shell:
        """
        exec > {log} 2>&1
        (
            gzip -cd "{input.input_bed}" 2>/dev/null || cat "{input.input_bed}"
        ) | awk -F'\t' -v prefixes="{input.graph_prefixes_input}" '
            BEGIN {{
                OFS="\t"
                # hap subtags from the input file (first column); a contig is kept only if its
                # name contains one, since only those contigs get a graph path name
                while ((getline line < prefixes) > 0) {{
                    split(line, f, "\t")
                    subtags[f[1]] = 1
                }}
            }}
            NR == 1 && $2 !~ /^[0-9]+$/ {{
                print $0, "peak_id"
                next
            }}
            NR == 1 {{
                header = "#chrom\tstart\tend"
                for (i = 4; i <= NF; i++) {{
                    header = header "\tcol" i
                }}
                print header, "peak_id"
            }}
            /unassigned/ {{ next }}
            {{
                if (!($1 in has_subtag)) {{
                    has_subtag[$1] = 0
                    for (t in subtags) if (index($1, t)) {{ has_subtag[$1] = 1; break }}
                }}
                if (!has_subtag[$1]) {{
                    n_no_subtag[$1]++
                    next
                }}
                # peak_id is chrom_start_end-run, so runs sharing an assembly keep their own peaks
                peak_id = $1"_"$2"_"$3"-"$4
                if (seen[peak_id]++) {{
                    n_dup++
                    next
                }}
                print $0, peak_id
            }}
            END {{
                if (n_dup) print "Removed " n_dup " duplicate peaks (same run and coordinates)" > "/dev/stderr"
                for (c in n_no_subtag) print "Removed " n_no_subtag[c] " peaks on " c ", which matches no hap subtag in the input file" > "/dev/stderr"
            }}' > {output.output_bed_unzipped}

        bgzip -c {output.output_bed_unzipped} > {output.output_bed}
        """

rule combine_fastas:
    input:
        ALL_FASTAS
    output:
        combined_fastas = temp(f"{global_output_dir}/{sample}/{sample}_combined_fastas.fa")
    log:
        f"{log_dir}/combine_fastas.log"
    shell:
        r"""
        exec > {log} 2>&1
        for f in {input}; do
            if [[ "$f" == *.gz ]]; then
                zcat "$f"
            else
                cat "$f"
            fi
        done > {output.combined_fastas}
        """

rule make_graph:
    input:
        graph_input = graph_input_file
        #hap1_fasta_path = rules.make_fasta_files.output.hap1_fasta,
        #hap2_fasta_path = rules.make_fasta_files.output.hap2_fasta,
    output:
        #cp_graph_input=f"{graphs_dir}/{{sample}}/{{sample}}_graph_input.txt",
        cp_graph=f"{graphs_dir}/{sample}.full.gbz"
    conda:
        "../envs/env.yml"
    #params:
    #    hap1_new_tag=get_hap1_new_tag,
    #    hap2_new_tag=get_hap2_new_tag
    threads: TEMP_THREADS
    resources:
        mem_mb=GRAPH_MEM,
        runtime=GRAPH_RUNTIME
    log:
        f"{log_dir}/make_graph.log"
    shell:
        """
        exec > {log} 2>&1
        if [ -d "{graphs_dir}/js" ]; then
            rm -r {graphs_dir}/js
        fi
        #rm -r {graphs_dir}/js

        cat {input.graph_input} > {graphs_dir}/inputs.txt

        sample_list=$(cat {input.graph_input} | cut -f1) 

        if [[ "{SPLIT_GRAPH}" == "False" ]]; then
            cactus-pangenome \
            {graphs_dir}/js \
            {graphs_dir}/inputs.txt \
            --outDir {graphs_dir}/ \
            --outName {sample} \
            --reference $sample_list \
            --vcfReference {REF_NAMES_STR} \
            --defaultCores 32 \
            --gfa full \
            --vcf full \
            --gbz full \
            --filter 0 \
            --noSplit \
            --binariesMode singularity
        else
            cactus-pangenome \
            {graphs_dir}/js \
            {graphs_dir}/inputs.txt \
            --outDir {graphs_dir}/ \
            --outName {sample} \
            --reference $sample_list \
            --vcfReference {REF_NAMES_STR} \
            --defaultCores 32 \
            --gfa full \
            --vcf full \
            --gbz full \
            --filter 0 \
            --binariesMode singularity
        fi

        """

rule make_list_of_chrs_in_graph_for_modifying_graph:
    input: 
        #graph=get_graph
        graph=rules.make_graph.output.cp_graph
    output:
        chrs_to_keep_list=temp(f"{graphs_dir}/{sample}-chrs_list_for_graph_mod.txt")
    #params:
    #    common_chr_tag = get_common_chr_tag,
        #vg=f"{vg}"
    conda:
        "../envs/env.yml"
    threads: 1
    resources:
        mem_mb=10240,
        runtime=30
    log:
        f"{log_dir}/make_list_of_chrs_in_graph_for_modifying_graph.log"
    shell:
        """
        exec > {log} 2>&1
        #vg='singularity exec /mmfs1/gscratch/stergachislab/minkinaa/vg_v1.65.0.sif vg'
        # every path except the references' (SAMPLE, column 3, is the reference name)
        vg paths -M -x {input.graph} | awk -F'\t' -v refs="{REF_NAMES_STR}" '
            BEGIN {{ n = split(refs, r, " "); for (i = 1; i <= n; i++) is_ref[r[i]] = 1 }}
            NR > 1 && !($3 in is_ref) {{ print $1 }}' | uniq > {output.chrs_to_keep_list}
        """

rule make_fasta_of_missing_chrs:
    input:
        fasta=rules.combine_fastas.output.combined_fastas,
        chrs_in_graph=rules.make_list_of_chrs_in_graph_for_modifying_graph.output.chrs_to_keep_list,
        graph_prefixes_input=graph_prefixes_input
    output:
        fasta_of_missing_chrs=temp(f"{graphs_dir}/{sample}-missing_chrs.fa")
    #params:
    #    hap1_orig_subtag=get_hap1_subtag,
    #    hap2_orig_subtag=get_hap2_subtag,
    #    hap1_new_tag=get_hap1_new_tag,
    #    hap2_new_tag=get_hap2_new_tag
    conda:
        "../envs/env.yml"
    log:
        f"{log_dir}/make_fasta_of_missing_chrs.log"
    shell:
        """
        exec > {log} 2>&1
        python {WORKFLOW_DIR}/scripts/make_fasta_of_missing_contigs_multiple_asms.py \
            {input.chrs_in_graph} \
            {input.graph_prefixes_input} \
            {input.fasta} \
            {output.fasta_of_missing_chrs} \
        """

rule make_graph_with_missing_chrs:
    input:
        missing_chr_fasta=rules.make_fasta_of_missing_chrs.output.fasta_of_missing_chrs,
        orig_graph=rules.make_graph.output.cp_graph
    output:
        graph_of_missing_contigs=temp(f"{graphs_dir}/{sample}.missing_contigs_graph.vg"),
        orig_graph_vg=temp(f"{graphs_dir}/{sample}.full.vg"),
        full_graph_vg=temp(f"{graphs_dir}/{sample}.full_graph_w_all_contigs.vg"),
        full_graph_gfa=temp(f"{graphs_dir}/{sample}.full_graph_w_all_contigs.gfa"),
        full_graph=f"{graphs_dir}/{sample}.full_all_chrs.gbz"
    conda:
        "../envs/env.yml"
    resources:
        mem_mb=50000,
        runtime=60*2
    log:
        f"{log_dir}/make_graph_with_missing_chrs.log"
    shell:
        """
        exec > {log} 2>&1
        LINE_COUNT=$(wc -l < {input.missing_chr_fasta})
        echo $LINE_COUNT

        if [ $LINE_COUNT -gt 1 ]; then
            vg construct -r {input.missing_chr_fasta} > {output.graph_of_missing_contigs}
            vg convert --packed-out --no-translation {input.orig_graph} > {output.orig_graph_vg}
            vg combine {output.orig_graph_vg} {output.graph_of_missing_contigs} > {output.full_graph_vg}
            vg convert -f {output.full_graph_vg} > {output.full_graph_gfa}
            vg gbwt -G {output.full_graph_gfa} --gbz-format -g {output.full_graph}
        else
            #vg convert --packed-out --no-translation {input.orig_graph} > {output.orig_graph_vg}
            # make files that will be deleted so rule doesn't crash
            touch {output.graph_of_missing_contigs}
            touch {output.orig_graph_vg}
            touch {output.full_graph_vg}
            touch {output.full_graph_gfa}
            cp {input.orig_graph} {output.full_graph}
        fi

        """

rule move_peaks_to_Ref:
    input:
        graph=graph,
        fasta=rules.combine_fastas.output.combined_fastas,
        bed=rules.make_bed_with_haps_and_peakIDs.output.output_bed,
        graph_prefixes_input=graph_prefixes_input
    output:
        fibertigs_out=temp(f"{global_output_dir}/{sample}/{sample}_all_original_peaks_fibertigs.bam"),
        fibertigs_for_inject_out=temp(f"{global_output_dir}/{sample}/{sample}_all_original_peaks_fibertigs_ToInject.bam"),
        gam=temp(f"{global_output_dir}/{sample}/{sample}_all_original_peaks_fibertigs_ToInject.gam"),
        surjected=temp(f"{global_output_dir}/{sample}/{sample}_all_original_peaks_fibertigs_surjected_on_Ref.bam"),
        surj_processed=temp(f"{global_output_dir}/{sample}/{sample}_all_original_peaks_fibertigs_surjected_on_Ref_processed.bam"),
        bed_final=temp(f"{global_output_dir}/{sample}/{sample}_all_original_peaks_fibertigs_surjected_on_Ref_processed.bed")
    params:
        asm_to_surject_onto=MAIN_REF_NAME
    conda:
        "../envs/env.yml"
    threads: VG_THREADS
    resources:
        mem_mb=mem_for(VG_MEM_MB),
        runtime=60*12
    log:
        f"{log_dir}/move_peaks_to_Ref.log"
    shell:
        """
        exec > {log} 2>&1
        sh {WORKFLOW_DIR}/scripts/inject_surject_bed_w_fibertigs.txt \
            {input.bed} \
            {input.graph} \
            {input.fasta} \
            {PATH_TO_FT} \
            {FT_PARAM} \
            {input.graph_prefixes_input} \
            {output.fibertigs_out} \
            {output.fibertigs_for_inject_out} \
            {output.gam} \
            {output.surjected} \
            {output.surj_processed} \
            {output.bed_final} \
            {WORKFLOW_DIR}/scripts \
            {params.asm_to_surject_onto} \
            {threads}
        
        """

rule make_untransferred_peaks_list:
    input:
        orig_bed = rules.make_bed_with_haps_and_peakIDs.output.output_bed,
        surj_bed = rules.move_peaks_to_Ref.output.bed_final,
        graph_prefixes_input=graph_prefixes_input
    output:
        orig_bed_annotated = f"{global_output_dir}/{sample}/{sample}_wHapsAndPeakIDs_w_surject_info.bed",
        orig_bed_annotated_surjected = f"{global_output_dir}/{sample}/{sample}_wHapsAndPeakIDs_w_surject_info_SURJECTED_ONTO_Ref.bed",
        orig_bed_annotated_not_surjected = f"{global_output_dir}/{sample}/{sample}_wHapsAndPeakIDs_w_surject_info_NOT_SURJECTED_ONTO_Ref.bed",
        # assign_untransferred_peaks_to_asms picks its own assemblies, so this is only kept until the job ends
        next_asm_after_ref = temp(f"{global_output_dir}/{sample}/{sample}_next_asm_after_Ref.txt"),
    params:
        last_surjected_asm=MAIN_REF_NAME
    conda:
        "../envs/R_env.yml"
    threads: 4
    resources:
        mem_mb=mem_for(16000),
        runtime=60*5
    log:
        f"{log_dir}/make_untransferred_peaks_list.log"
    shell:
        """
        exec > {log} 2>&1
        Rscript {WORKFLOW_DIR}/scripts/260126_annotate_original_peak_file_w_surjection_info.R \
            {input.orig_bed} \
            {input.surj_bed} \
            {params.last_surjected_asm} \
            {output.orig_bed_annotated} \
            {output.orig_bed_annotated_surjected} \
            {output.orig_bed_annotated_not_surjected} \
            {output.next_asm_after_ref} \
            {input.graph_prefixes_input} \
            {length_change_param}

        """

rule make_gam_of_untransferred_peaks:
    input:
        graph=graph,
        fasta=rules.combine_fastas.output.combined_fastas,
        bed=rules.make_untransferred_peaks_list.output.orig_bed_annotated_not_surjected,
        graph_prefixes_input=graph_prefixes_input
    output:
        bed_for_inject=temp(f"{global_output_dir}/{sample}/{sample}_untransferred_TEMP_peaks_for_inject.bed"),
        fibertigs_out=temp(f"{global_output_dir}/{sample}/{sample}_untransferred_TEMP_peaks_fibertigs.bam"),
        fibertigs_for_inject_out=temp(f"{global_output_dir}/{sample}/{sample}_untransferred_TEMP_peaks_fibertigs_ToInject.bam"),
        gam=temp(f"{global_output_dir}/{sample}/{sample}_NonRef_peaks.gam")
    conda:
        "../envs/env.yml"
    threads: VG_THREADS
    resources:
        mem_mb=mem_for(VG_MEM_MB),
        runtime=60*12
    log:
        f"{log_dir}/make_gam_of_untransferred_peaks.log"
    shell:
        """
        exec > {log} 2>&1
        # pg-inject needs the header commented out; edit a copy so the input (and its mtime) is untouched
        sed '1s/^#*/#/' {input.bed} > {output.bed_for_inject}

        sh {WORKFLOW_DIR}/scripts/260127_inject_w_fibertigs.txt \
            {output.bed_for_inject} \
            {input.graph} \
            {input.fasta} \
            {PATH_TO_FT} \
            {FT_PARAM} \
            {input.graph_prefixes_input} \
            {output.fibertigs_out} \
            {output.fibertigs_for_inject_out} \
            {output.gam} \
            {WORKFLOW_DIR}/scripts \
            {threads}

        """

# Surject the peaks that didn't go to the main reference onto every input assembly, one job per assembly.
# assign_untransferred_peaks_to_asms then picks assemblies in greedy order from these results.
rule surject_untransferred_peaks_onto_asm:
    input:
        graph=graph,
        gam=rules.make_gam_of_untransferred_peaks.output.gam,
        graph_prefixes_input=graph_prefixes_input
    output:
        surjected=temp(f"{global_output_dir}/{sample}/untransferred_peaks_surjected_onto_asms/{sample}_TEMP_surjected_onto_{{asm}}.bam"),
        surj_processed=temp(f"{global_output_dir}/{sample}/untransferred_peaks_surjected_onto_asms/{sample}_TEMP_surjected_onto_{{asm}}_processed.bam"),
        bed_final=temp(f"{global_output_dir}/{sample}/untransferred_peaks_surjected_onto_asms/{sample}_untransferred_peaks_surjected_onto_{{asm}}.bed")
    conda:
        "../envs/env.yml"
    threads: VG_THREADS
    resources:
        mem_mb=mem_for(VG_MEM_MB),
        runtime=60*12
    log:
        f"{log_dir}/surject_untransferred_peaks_onto_asm_{{asm}}.log"
    shell:
        """
        exec > {log} 2>&1
        sh {WORKFLOW_DIR}/scripts/260127_surject_w_fibertigs.txt \
            {input.graph} \
            {PATH_TO_FT} \
            {input.graph_prefixes_input} \
            {input.gam} \
            {output.surjected} \
            {output.surj_processed} \
            {output.bed_final} \
            {WORKFLOW_DIR}/scripts \
            {wildcards.asm} \
            {threads}
        """

rule assign_untransferred_peaks_to_asms:
    input:
        annotated_bed_noRef=rules.make_untransferred_peaks_list.output.orig_bed_annotated_not_surjected,
        graph_prefixes_input=graph_prefixes_input,
        surjected_beds=expand(rules.surject_untransferred_peaks_onto_asm.output.bed_final, asm=unique_asms)
    output:
        fully_annotated_bed=f"{global_output_dir}/{sample}/{sample}_wHapsAndPeakIDs_w_surject_info_NOT_SURJECTED_ONTO_Ref_w_FullyAnnotated.bed",
        # assemblies tried, in order
        assemblies_tried=assemblies_surjected_onto_file
    params:
        asm_names=",".join(unique_asms)
    conda:
        "../envs/R_env.yml"
    threads: 1
    resources:
        mem_mb=mem_for(16000),
        runtime=60*2
    log:
        f"{log_dir}/assign_untransferred_peaks_to_asms.log"
    shell:
        """
        exec > {log} 2>&1
        Rscript {WORKFLOW_DIR}/scripts/assign_untransferred_peaks_to_asms.R \
            {input.annotated_bed_noRef} \
            {input.graph_prefixes_input} \
            {length_change_param} \
            {output.fully_annotated_bed} \
            {output.assemblies_tried} \
            {MAIN_REF_NAME} \
            {params.asm_names} \
            {input.surjected_beds}
        """

rule call_global_peaks_on_Ref:
    input:
        bed=rules.make_untransferred_peaks_list.output.orig_bed_annotated_surjected
    output:
        bed_4col=temp(f"{global_output_dir}/{sample}/{sample}_wHapsAndPeakIDs_w_surject_info_SURJECTED_ONTO_Ref_4col.bed"),
        mock_bam=temp(f"{global_output_dir}/{sample}/{sample}_wHapsAndPeakIDs_w_surject_info_SURJECTED_ONTO_Ref_4col_mockbam.bam"),
        ref_peaks=f"{global_output_dir}/{sample}/{sample}_Ref_shared_peaks.bed",
        ref_peaks_5col=temp(f"{global_output_dir}/{sample}/{sample}_Ref_shared_peaks_5col.bed"),
        ref_peaks_5col_w_peak_id=f"{global_output_dir}/{sample}/{sample}_Ref_shared_peaks_6col.bed",
        ref_peaks_5col_w_peak_id_zipped=temp(f"{global_output_dir}/{sample}/{sample}_Ref_shared_peaks_6col.bed.gz")
    conda:
        "../envs/env.yml"
    threads: 4
    resources:
        mem_mb=mem_for(64000),
        runtime=60*12
    log:
        f"{log_dir}/call_global_peaks_on_Ref.log"
    shell:
        """
        exec > {log} 2>&1
        rm -f {global_output_dir}/{sample}/{sample}*SURJECTED_ONTO_Ref_4col_mockbam.bam.tmp*
        awk 'NR==1 {{print "#chrom\tstart\tend\tsample"; next}} {{print $(NF-2), $(NF-1), $NF, $4"_"$(NF-2)}}' OFS="\t" {input.bed} > {output.bed_4col}
        {PATH_TO_FT} mock-fire -t {threads} {output.bed_4col} | samtools sort -@ {threads} -m 2G -o {output.mock_bam} --write-index
        {PATH_TO_FT} call-peaks -t {threads} {output.mock_bam} --min-fire-coverage 1 --min-fire-frac 0.01 > {output.ref_peaks}
        cat {output.ref_peaks} | cut -f1,2,3 > {output.ref_peaks_5col}

        bash {WORKFLOW_DIR}/scripts/add_peak_id_col_to_bed_file.sh \
            {output.ref_peaks_5col} \
            {output.ref_peaks_5col_w_peak_id} \
            {output.ref_peaks_5col_w_peak_id_zipped}
        """

rule call_global_peaks_on_NonRef:
    input:
        bed=rules.assign_untransferred_peaks_to_asms.output.fully_annotated_bed
    output:
        bed_4col=temp(f"{global_output_dir}/{sample}/{sample}_wHapsAndPeakIDs_w_surject_info_NOT_SURJECTED_ONTO_Ref_w_FullyAnnotated_4col.bed"),
        mock_bam=temp(f"{global_output_dir}/{sample}/{sample}_wHapsAndPeakIDs_w_surject_info_NOT_SURJECTED_ONTO_Ref_w_FullyAnnotated_4col_mockbam.bam"),
        nonref_peaks=f"{global_output_dir}/{sample}/{sample}_NonRef_shared_peaks.bed",
        nonref_peaks_5col=temp(f"{global_output_dir}/{sample}/{sample}_NonRef_shared_peaks_5col.bed"),
        nonref_peaks_5col_w_peak_id=f"{global_output_dir}/{sample}/{sample}_NonRef_shared_peaks_6col.bed",
        nonref_peaks_5col_w_peak_id_zipped=temp(f"{global_output_dir}/{sample}/{sample}_NonRef_shared_peaks_6col.bed.gz")
    conda:
        "../envs/env.yml"
    threads: 4
    resources:
        mem_mb=mem_for(64000),
        runtime=60*12
    log:
        f"{log_dir}/call_global_peaks_on_NonRef.log"
    shell:
        """
        exec > {log} 2>&1
        rm -f {global_output_dir}/{sample}/{sample}*NOT_SURJECTED_ONTO_Ref_w_FullyAnnotated_4col_mockbam.bam.tmp*

        cat {input.bed} | awk '$7!="NONE"' | awk 'NR==1 {{print "#chrom\tstart\tend\tsample"; next}} {{print $(NF-2), $(NF-1), $NF, $4"_"$(NF-2)}}' OFS="\t" > {output.bed_4col}
        {PATH_TO_FT} mock-fire -t {threads} {output.bed_4col} | samtools sort -@ {threads} -m 2G -o {output.mock_bam} --write-index
        {PATH_TO_FT} call-peaks -t {threads} {output.mock_bam} --min-fire-coverage 1 --min-fire-frac 0.01 > {output.nonref_peaks}
        cat {output.nonref_peaks} | cut -f1,2,3 > {output.nonref_peaks_5col}

        bash {WORKFLOW_DIR}/scripts/add_peak_id_col_to_bed_file.sh \
            {output.nonref_peaks_5col} \
            {output.nonref_peaks_5col_w_peak_id} \
            {output.nonref_peaks_5col_w_peak_id_zipped}
        """

rule make_Ref_global_peaks_gam:
    input:
        bed=rules.call_global_peaks_on_Ref.output.ref_peaks_5col_w_peak_id,
        graph=graph,
        fasta=MAIN_REF_FASTA_PLAIN
        #graph_prefixes_input=graph_prefixes_input
    output:
        fibertigs_out=temp(f"{global_output_dir}/{sample}/{sample}_Ref_TEMP_peaks_fibertigs.bam"),
        fibertigs_for_inject_out=temp(f"{global_output_dir}/{sample}/{sample}_Ref_TEMP_peaks_fibertigs_ToInject.bam"),
        gam=temp(f"{global_output_dir}/{sample}/{sample}_Ref_global_peaks.gam")
    conda:
        "../envs/env.yml"
    threads: VG_THREADS
    resources:
        mem_mb=mem_for(VG_MEM_MB),
        runtime=60*12
    log:
        f"{log_dir}/make_Ref_global_peaks_gam.log"
    shell:
        """
        exec > {log} 2>&1
        sh {WORKFLOW_DIR}/scripts/260127_inject_w_fibertigs_single_asm.txt \
            {input.bed} \
            {input.graph} \
            {input.fasta} \
            {PATH_TO_FT} \
            {FT_PARAM} \
            {MAIN_REF_NAME} \
            {output.fibertigs_out} \
            {output.fibertigs_for_inject_out} \
            {output.gam} \
            {WORKFLOW_DIR}/scripts \
            {threads}
        """

rule make_NonRef_global_peaks_gam:
    input:
        bed=rules.call_global_peaks_on_NonRef.output.nonref_peaks_5col_w_peak_id,
        graph=graph,
        fasta=rules.combine_fastas.output.combined_fastas,
        graph_prefixes_input=graph_prefixes_input
    output:
        fibertigs_out=temp(f"{global_output_dir}/{sample}/{sample}_NonRef_TEMP_peaks_fibertigs.bam"),
        fibertigs_for_inject_out=temp(f"{global_output_dir}/{sample}/{sample}_NonRef_TEMP_peaks_fibertigs_ToInject.bam"),
        gam=temp(f"{global_output_dir}/{sample}/{sample}_NonRef_global_peaks.gam")
        #peaks_on_asm=f"{global_output_dir}/{sample}/{peaks_on_all_asm}/{{asm}}/{sample}_nonref_peaks_on_{{asm}}.bed"
    conda:
        "../envs/env.yml"
    threads: VG_THREADS
    resources:
        mem_mb=mem_for(VG_MEM_MB),
        runtime=60*12
    log:
        f"{log_dir}/make_NonRef_global_peaks_gam.log"
    shell:
        """
        exec > {log} 2>&1
        sh {WORKFLOW_DIR}/scripts/260127_inject_w_fibertigs.txt \
            {input.bed} \
            {input.graph} \
            {input.fasta} \
            {PATH_TO_FT} \
            {FT_PARAM} \
            {input.graph_prefixes_input} \
            {output.fibertigs_out} \
            {output.fibertigs_for_inject_out} \
            {output.gam} \
            {WORKFLOW_DIR}/scripts \
            {threads}
        """

rule surject_global_Ref_peaks_onto_all_assemblies:
    input:
        gam=rules.make_Ref_global_peaks_gam.output.gam,
        graph=graph,
        graph_prefixes_input=graph_prefixes_input
    output:
        surjected=temp(f"{global_output_dir}/{sample}/Ref_global_peaks_surjected_onto_asms/{sample}_TEMP_surjected_onto_{{asm}}.bam"),
        surj_processed=temp(f"{global_output_dir}/{sample}/Ref_global_peaks_surjected_onto_asms/{sample}_TEMP_surjected_onto_{{asm}}_processed.bam"),
        bed_final=temp(f"{global_output_dir}/{sample}/Ref_global_peaks_surjected_onto_asms/{sample}_Ref_global_peaks_surjected_onto_{{asm}}.bed")
    conda:
        "../envs/env.yml"
    threads: VG_THREADS
    resources:
        mem_mb=mem_for(VG_MEM_MB),
        runtime=60*12
    log:
        f"{log_dir}/surject_global_Ref_peaks_onto_all_assemblies_{{asm}}.log"
    shell:
        """
        exec > {log} 2>&1
            sh {WORKFLOW_DIR}/scripts/260127_surject_w_fibertigs.txt \
                {input.graph} \
                {PATH_TO_FT} \
                {input.graph_prefixes_input} \
                {input.gam} \
                {output.surjected} \
                {output.surj_processed} \
                {output.bed_final} \
                {WORKFLOW_DIR}/scripts \
                {wildcards.asm} \
                {threads}
        """

rule surject_global_NonRef_peaks_onto_all_assemblies:
    input:
        gam=rules.make_NonRef_global_peaks_gam.output.gam,
        graph=graph,
        graph_prefixes_input=graph_prefixes_input
    output:
        surjected=temp(f"{global_output_dir}/{sample}/NonRef_global_peaks_surjected_onto_asms/{sample}_TEMP_surjected_onto_{{asm}}.bam"),
        surj_processed=temp(f"{global_output_dir}/{sample}/NonRef_global_peaks_surjected_onto_asms/{sample}_TEMP_surjected_onto_{{asm}}_processed.bam"),
        bed_final=temp(f"{global_output_dir}/{sample}/NonRef_global_peaks_surjected_onto_asms/{sample}_NonRef_global_peaks_surjected_onto_{{asm}}.bed"),
    conda:
        "../envs/env.yml"
    threads: VG_THREADS
    resources:
        mem_mb=mem_for(VG_MEM_MB),
        runtime=60*12
    log:
        f"{log_dir}/surject_global_NonRef_peaks_onto_all_assemblies_{{asm}}.log"
    shell:
        """
        exec > {log} 2>&1
            sh {WORKFLOW_DIR}/scripts/260127_surject_w_fibertigs.txt \
                {input.graph} \
                {PATH_TO_FT} \
                {input.graph_prefixes_input} \
                {input.gam} \
                {output.surjected} \
                {output.surj_processed} \
                {output.bed_final} \
                {WORKFLOW_DIR}/scripts \
                {wildcards.asm} \
                {threads}
        """

rule combine_all_consensus_peaks:
    input:
        ref_consensus_peaks=rules.call_global_peaks_on_Ref.output.ref_peaks_5col_w_peak_id,
        nonref_consensus_peaks=rules.call_global_peaks_on_NonRef.output.nonref_peaks_5col_w_peak_id
    output:
        consensus_bed_combined=f"{global_output_dir}/{sample}/{sample}_ALL_consensus_peaks_combined.bed"
    conda:
        "../envs/env.yml"
    log:
        f"{log_dir}/combine_all_consensus_peaks.log"
    shell:
        """
        exec > {log} 2>&1
        #cat {input.ref_consensus_peaks} {input.nonref_consensus_peaks} > {output.consensus_bed_combined}
        cat {input.ref_consensus_peaks} <( grep -v "#chrom" {input.nonref_consensus_peaks}) > {output.consensus_bed_combined}
        """

rule combined_Ref_and_NonRef_transferred_peak_files:
    input:
        ref_peaks_bed=rules.surject_global_Ref_peaks_onto_all_assemblies.output.bed_final,
        nonref_peaks_bed=rules.surject_global_NonRef_peaks_onto_all_assemblies.output.bed_final
    output:
        bed_combined=temp(f"{global_output_dir}/{sample}/Combined_global_peaks_surjected_onto_asms/{sample}_ALL_global_peaks_surjected_onto_{{asm}}.bed")
    conda:
        "../envs/env.yml"
    log:
        f"{log_dir}/combined_Ref_and_NonRef_transferred_peak_files_{{asm}}.log"
    shell:
        """
        exec > {log} 2>&1
        cat {input.ref_peaks_bed} {input.nonref_peaks_bed} > {output.bed_combined}
        """

rule list_Ref_chromosomes:
    input:
        ref_fasta=MAIN_REF_FASTA,
        graph=graph,
        graph_prefixes_input=graph_prefixes_input
    output:
        chroms=f"{global_output_dir}/{sample}/{sample}_Ref_chromosomes.txt"
    conda:
        "../envs/env.yml"
    resources:
        mem_mb=mem_for(16000),
        runtime=60
    log:
        f"{log_dir}/list_Ref_chromosomes.log"
    shell:
        """
        exec > {log} 2>&1
        # the name is the first word of each header line; zcat -f reads gzipped and plain fastas
        zcat -f {input.ref_fasta} | grep '^>' | sed -e 's/^>//' -e 's/[[:space:]].*//' > {output.chroms}

        # each reference must be a path in the graph, under its name
        graph_samples=$(vg paths -M -x {input.graph} | awk -F'\t' 'NR > 1 {{ print $3 }}' | sort -u)
        for ref in {REF_NAMES_STR}; do
            if ! echo "$graph_samples" | grep -qxF "$ref"; then
                echo "Reference $ref is not a path in {input.graph}; its name must match the graph's"
                exit 1
            fi
        done

        # reference contigs are told apart from assembly contigs by the assemblies' hap subtags
        bad=$(cut -f1 {input.graph_prefixes_input} | grep -F -f - {output.chroms} || true)
        if [ -n "$bad" ]; then
            echo "These main reference contig names contain an assembly's hap subtag:"
            echo "$bad" | head -20
            exit 1
        fi
        """

rule determine_previous_peak_for_nonRef_chrs:
    input:
        bed=rules.combined_Ref_and_NonRef_transferred_peak_files.output.bed_combined,
        graph_prefixes_input=graph_prefixes_input,
        ref_chroms=rules.list_Ref_chromosomes.output.chroms,
        ref_peaks=rules.call_global_peaks_on_Ref.output.ref_peaks_5col_w_peak_id
    output:
        previous_peaks=f"{global_output_dir}/{sample}/temp_previous_peak_files/{sample}_previous_peaks_{{asm}}.tsv"
    conda:
        "../envs/R_env.yml"
    resources:
        mem_mb=mem_for(16000),
        runtime=60*2
    log:
        f"{log_dir}/determine_previous_peak_for_nonRef_chrs_{{asm}}.log"
    shell:
        """
        exec > {log} 2>&1
        Rscript {WORKFLOW_DIR}/scripts/261030_determine_previous_peak_for_nonRef_chrs.R \
            {input.bed} \
            {input.graph_prefixes_input} \
            {wildcards.asm} \
            {input.ref_chroms} \
            {output.previous_peaks} \
            {input.ref_peaks}
        """

rule subset_pileup_file:
    input:
        pileup=get_pileup,
        graph_prefixes_input=graph_prefixes_input
    output:
        mod_pileup=temp(f"{global_output_dir}/{sample}/run_pileup_counts/{{run}}_{{asm}}/{{run}}_{{asm}}_modified_pileup.bed")
    conda:
        "../envs/env.yml"
    resources:
        mem_mb=mem_for(16000),
        runtime=60
    log:
        f"{log_dir}/subset_pileup_file_{{run}}_{{asm}}.log"
    shell:
        """
        exec > {log} 2>&1
        sample_hap_name=$(awk -F'\t' -v a="{wildcards.asm}" '$2 == a {{ print $1 }}' {input.graph_prefixes_input})
        if [[ $(echo "$sample_hap_name" | wc -l) -ne 1 || -z "$sample_hap_name" ]]; then
            echo "Expected one entry for {wildcards.asm} in {input.graph_prefixes_input}, found: $sample_hap_name" >&2
            exit 1
        fi
        zcat {input.pileup} | awk -F'\t' -v p="$sample_hap_name" 'index($1, p)' | cut -f1-6 > {output.mod_pileup}
        """

rule intersect_peaks_and_pileup:
    input:
        mod_pileup=rules.subset_pileup_file.output.mod_pileup,
        peaks=rules.combined_Ref_and_NonRef_transferred_peak_files.output.bed_combined
    output:
        overlaps=temp(f"{global_output_dir}/{sample}/run_pileup_counts/{{run}}_{{asm}}/{{run}}_{{asm}}_pileup_intervals_overlapping_peaks.bed")
    conda:
        "../envs/env.yml"
    resources:
        mem_mb=mem_for(16000),
        runtime=60
    log:
        f"{log_dir}/intersect_peaks_and_pileup_{{run}}_{{asm}}.log"
    shell:
        """
        exec > {log} 2>&1
        bedtools intersect -a {input.mod_pileup} -b {input.peaks} -wa -wb > {output.overlaps}
        """    

rule process_intersection_file:
    input:
        overlaps=rules.intersect_peaks_and_pileup.output.overlaps,
        consensus_peaks=rules.combine_all_consensus_peaks.output.consensus_bed_combined,
        consensus_peaks_transferred=rules.combined_Ref_and_NonRef_transferred_peak_files.output.bed_combined
    output:
        max_score_per_consensus_peak=temp(f"{global_output_dir}/{sample}/run_pileup_counts/{{run}}_{{asm}}/{{run}}_{{asm}}_max_score_per_peak.bed")
    conda:
        "../envs/env.yml"
    resources:
        mem_mb=mem_for(16000),
        runtime=60
    log:
        f"{log_dir}/process_intersection_file_{{run}}_{{asm}}.log"
    shell:
        """
        exec > {log} 2>&1
        python3 {WORKFLOW_DIR}/scripts/260217_find_max_fire_score_per_consenses_peak.py \
            {input.overlaps} \
            {input.consensus_peaks} \
            {input.consensus_peaks_transferred} \
            {wildcards.run}_{wildcards.asm} \
            {output.max_score_per_consensus_peak}
        """

rule combine_max_score_outputs:
    input:
        ind_input=max_score_per_consensus_peak
    output:
        combined_out=f"{global_output_dir}/{sample}/{sample}_max_score_per_peak_all_samples.bed"
    conda:
        "../envs/env.yml"
    resources:
        mem_mb=mem_for(16000),
        runtime=60
    log:
        f"{log_dir}/combine_max_score_outputs.log"
    shell:
        """
        exec > {log} 2>&1
        echo "sample_id\tconsensus_peak_id\tcoverage\tfire_coverage\tscore" > {output.combined_out}
        cat <(grep -v -h "sample_id" {input.ind_input}) >> {output.combined_out}
        """

rule combine_all_transferred_peak_files:
    input:
        combined_global_peaks_surjected_onto_all_asm
    output:
        combined_out=f"{global_output_dir}/{sample}/{sample}_combined_tranferred_peaks_all_samples.bed"
    conda:
        "../envs/env.yml"
    resources:
        mem_mb=mem_for(16000),
        runtime=60
    log:
        f"{log_dir}/combine_all_transferred_peak_files.log"
    shell:
        """
        exec > {log} 2>&1
        echo "chr\tstart\tend\tcons_peak_id" > {output.combined_out}
        cat {input} >> {output.combined_out}
        """

# checks each non-reference peak's placement against every assembly it transferred to, and flags
# peaks that look like copies of a main reference peak (the placement itself still comes from the assembly
# the peak was called on)
rule check_previous_peaks_across_assemblies:
    input:
        combined_transferred_peaks=rules.combine_all_transferred_peak_files.output.combined_out,
        graph_prefixes_input=graph_prefixes_input,
        ref_chroms=rules.list_Ref_chromosomes.output.chroms,
        ref_peaks=rules.call_global_peaks_on_Ref.output.ref_peaks_5col_w_peak_id
    output:
        checks=f"{global_output_dir}/{sample}/additional_outputs/{sample}_previous_peak_checks.tsv",
        per_asm=f"{global_output_dir}/{sample}/additional_outputs/{sample}_previous_peak_intervals_per_asm.tsv.gz"
    conda:
        "../envs/R_env.yml"
    threads: 4
    resources:
        mem_mb=mem_for(40000),  # 26 GB for HPRC (78 assemblies)
        runtime=runtime_for(60)
    log:
        f"{log_dir}/check_previous_peaks_across_assemblies.log"
    shell:
        """
        exec > {log} 2>&1
        Rscript {WORKFLOW_DIR}/scripts/261105_check_previous_peaks_across_assemblies.R \
            {input.combined_transferred_peaks} \
            {input.graph_prefixes_input} \
            {input.ref_chroms} \
            {input.ref_peaks} \
            {global_output_dir}/{sample}/additional_outputs/{sample} \
            {threads}
        """

rule surject_Ref_consensus_peaks_onto_SecondaryRef:
    input:
        gam=rules.make_Ref_global_peaks_gam.output.gam,
        graph=graph,
        graph_prefixes_input=graph_prefixes_input
    output:
        surjected=temp(f"{global_output_dir}/{sample}/{sample}_Ref_consensus_peaks_moved_to_SecondaryRef_raw_surjected.bam"),
        surj_processed=temp(f"{global_output_dir}/{sample}/{sample}_Ref_consensus_peaks_moved_to_SecondaryRef_processed.bam"),
        bed_final=temp(f"{global_output_dir}/{sample}/{sample}_Ref_consensus_peaks_moved_to_SecondaryRef.bed")
    threads: VG_THREADS
    resources:
        mem_mb=mem_for(VG_MEM_MB),
        runtime=60*12
    params:
        asm_to_surject_onto=SECONDARY_REF_NAME
    conda:
        "../envs/env.yml"
    log:
        f"{log_dir}/surject_Ref_consensus_peaks_onto_SecondaryRef.log"
    shell:
        """
        exec > {log} 2>&1
        sh {WORKFLOW_DIR}/scripts/surject_gam_w_fibertigs.txt \
            {input.graph} \
            {PATH_TO_FT} \
            {input.graph_prefixes_input} \
            {input.gam} \
            {output.surjected} \
            {output.surj_processed} \
            {output.bed_final} \
            {WORKFLOW_DIR}/scripts \
            {params.asm_to_surject_onto} \
            {threads}
        """

rule surject_NonRef_consensus_peaks_onto_SecondaryRef:
    input:
        gam=rules.make_NonRef_global_peaks_gam.output.gam,
        graph=graph,
        graph_prefixes_input=graph_prefixes_input
    output:
        surjected=temp(f"{global_output_dir}/{sample}/{sample}_NonRef_consensus_peaks_moved_to_SecondaryRef_raw_surjected.bam"),
        surj_processed=temp(f"{global_output_dir}/{sample}/{sample}_NonRef_consensus_peaks_moved_to_SecondaryRef_processed.bam"),
        bed_final=temp(f"{global_output_dir}/{sample}/{sample}_NonRef_consensus_peaks_moved_to_SecondaryRef.bed")
    threads: VG_THREADS
    resources:
        mem_mb=mem_for(VG_MEM_MB),
        runtime=60*12
    params:
        asm_to_surject_onto=SECONDARY_REF_NAME
    conda:
        "../envs/env.yml"
    log:
        f"{log_dir}/surject_NonRef_consensus_peaks_onto_SecondaryRef.log"
    shell:
        """
        exec > {log} 2>&1
        sh {WORKFLOW_DIR}/scripts/surject_gam_w_fibertigs.txt \
            {input.graph} \
            {PATH_TO_FT} \
            {input.graph_prefixes_input} \
            {input.gam} \
            {output.surjected} \
            {output.surj_processed} \
            {output.bed_final} \
            {WORKFLOW_DIR}/scripts \
            {params.asm_to_surject_onto} \
            {threads}
        """

rule combine_SecondaryRef_transferred_peaks:
    input:
        ref_peaks=rules.surject_Ref_consensus_peaks_onto_SecondaryRef.output.bed_final,
        nonref_peaks=rules.surject_NonRef_consensus_peaks_onto_SecondaryRef.output.bed_final 
    output:
        combined_out=f"{global_output_dir}/{sample}/{sample}_All_consensus_peaks_moved_to_SecondaryRef.bed"
    resources:
        mem_mb=mem_for(16000),
        runtime=60*5
    conda:
        "../envs/env.yml"
    log:
        f"{log_dir}/combine_SecondaryRef_transferred_peaks.log"
    shell:
        """
        exec > {log} 2>&1
        cat {input.ref_peaks} {input.nonref_peaks} > {output.combined_out}
        """

rule intersect_original_and_transferred_consensus_peaks:
    input: 
        orig_bed=rules.make_bed_with_haps_and_peakIDs.output.output_bed,
        transferred_bed=rules.combine_all_transferred_peak_files.output.combined_out
    output:
        intersect_out=f"{global_output_dir}/{sample}/{sample}_orig_called_peaks_intersect_transferred_consensus_peaks.bed"
    threads: 4
    resources:
        mem_mb=mem_for(16000),
        runtime=60*5
    conda:
        "../envs/env.yml"
    log:
        f"{log_dir}/intersect_original_and_transferred_consensus_peaks.log"
    shell:
        """
        exec > {log} 2>&1
        # with -sorted, bedtools streams both inputs instead of loading -b into memory
        export LC_ALL=C
        sort_bed="sort -k1,1 -k2,2n -S 4G --parallel={threads} -T ${{TMPDIR:-/tmp}}"
        bedtools intersect -sorted -a <(zcat {input.orig_bed} | tail -n +2 | $sort_bed) -b <(tail -n +2 {input.transferred_bed} | $sort_bed) -wa -wb > {output.intersect_out}
        """
    
rule all_peaks_surjected_onto_consensus_assemblies:
    input:
        ref=rules.make_untransferred_peaks_list.output.orig_bed_annotated_surjected,
        nonref=rules.assign_untransferred_peaks_to_asms.output.fully_annotated_bed
    output: 
        out=temp(f"{global_output_dir}/{sample}/{sample}_all_peaks_surjected_onto_consensus_assemblies_5col_for_intersect.bed")
    resources:
        mem_mb=mem_for(16000),
        runtime=60*2
    conda:
        "../envs/env.yml"
    log:
        f"{log_dir}/all_peaks_surjected_onto_consensus_assemblies.log"
    shell:
        """
        exec > {log} 2>&1
        cat {input.ref} | tail -n +2 | awk  'BEGIN {{OFS="\t"}} {{print($8,$9,$10,$4,$5)}}' | grep -v NONE > {output.out}
        cat {input.nonref} | tail -n +2 | awk  'BEGIN {{OFS="\t"}} {{print($8,$9,$10,$4,$5)}}' | grep -v NONE >> {output.out}
        """

rule intersect_cons_surjected_peaks_w_consensus_peaks:
    input:
        cons_surj_peaks=rules.all_peaks_surjected_onto_consensus_assemblies.output.out,
        cons_peaks=rules.combine_all_consensus_peaks.output.consensus_bed_combined
    output:
        out=f"{global_output_dir}/{sample}/{sample}_all_orig_surjected_peaks_intersect_consensus_peaks.bed"
    resources:
        mem_mb=mem_for(16000),
        runtime=60*5
    conda:
        "../envs/env.yml"
    log:
        f"{log_dir}/intersect_cons_surjected_peaks_w_consensus_peaks.log"
    shell:
        """
        exec > {log} 2>&1
        bedtools intersect -a {input.cons_surj_peaks} -b {input.cons_peaks} -wa -wb > {output.out}
        """

rule get_native_peak_sequences:
    input:
        combined_transferred_peaks_bed=rules.combine_all_transferred_peak_files.output.combined_out,
        fasta=rules.combine_fastas.output.combined_fastas
    output:
        seq_out=f"{global_output_dir}/{sample}/{sample}_transferred_cons_peaks_native_sequences_from_main_output.bed.gz"
    threads: 4
    resources:
        mem_mb=mem_for(16000),
        runtime=60*5
    conda:
        "../envs/env.yml"
    log:
        f"{log_dir}/get_native_peak_sequences.log"
    shell:
        """
        exec > {log} 2>&1
        # GNU sort spills to disk; bedtools sort loads the whole file
        export LC_ALL=C
        bedtools getfasta -bedOut -tab -name -fi {input.fasta} -bed <(tail -n +2 {input.combined_transferred_peaks_bed} | sort -k1,1 -k2,2n -S 4G --parallel={threads} -T ${{TMPDIR:-/tmp}}) | bgzip -@ {threads} > {output.seq_out}
        """
    
rule make_metadata_table:
    input:
        max_score_tbl=rules.combine_max_score_outputs.output.combined_out,
        combined_transferred_peaks=rules.combine_all_transferred_peak_files.output.combined_out,
        # only with a secondary reference
        moved_to_secondary_ref=rules.combine_SecondaryRef_transferred_peaks.output.combined_out if HAS_SECONDARY_REF else [],
        orig_peaks_intersect_cons=rules.intersect_original_and_transferred_consensus_peaks.output.intersect_out,
        surj_peaks_intersect_cons=rules.intersect_cons_surjected_peaks_w_consensus_peaks.output.out,
        graph_prefixes_input=graph_prefixes_input,
        sample_ids=sample_ids_file,
        ref_chroms=rules.list_Ref_chromosomes.output.chroms,
        previous_peaks=expand(rules.determine_previous_peak_for_nonRef_chrs.output.previous_peaks, asm=unique_asms)
    output:
        meta_tbl=f"{global_output_dir}/{sample}/{sample}_final_tbl_w_metadata.tsv.gz",
        previous_peaks=f"{global_output_dir}/{sample}/additional_outputs/{sample}_previous_peaks_nonRef_chrs.tsv",
        peak_summary=f"{global_output_dir}/{sample}/{sample}_consensus_peak_summary.tsv"
    params:
        secondary_ref_name=SECONDARY_REF_NAME if HAS_SECONDARY_REF else "NA"
    conda:
        "../envs/R_env.yml"
    threads: 4
    resources:
        mem_mb=mem_for(128000),  # peaked at 85 GB for 59 runs
        runtime=runtime_for(60*5)
    log:
        f"{log_dir}/make_metadata_table.log"
    shell:
        """
        exec > {log} 2>&1
        Rscript {WORKFLOW_DIR}/scripts/call_peaks_make_metadata_table.R \
            {global_output_dir}/{sample}/ \
            {sample} \
            {input.graph_prefixes_input} \
            {threads} \
            {input.ref_chroms} \
            {MAIN_REF_NAME} \
            {params.secondary_ref_name} \
            {length_change_param} \
            {input.sample_ids} \
            {input.previous_peaks}
        """

