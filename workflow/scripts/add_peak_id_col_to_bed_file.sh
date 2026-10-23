#!/bin/bash

input_bed=$1
unzipped_output=$2
zipped_output=$3

(
    gzip -cd "${input_bed}" 2>/dev/null || cat "${input_bed}"
) | awk -F'\t' '
    BEGIN { OFS="\t" }
    NR == 1 && $2 !~ /^[0-9]+$/ {
        print $0, "peak_id"
        next
    }
    NR == 1 {
        header = "#chrom\tstart\tend"
        for (i = 4; i <= NF; i++) {
            header = header "\tcol" i
        }
        print header, "peak_id"
    }
    {
        peak_id = $1"_"$2"_"$3
        # keep the first copy of each peak
        if (seen[peak_id]++) {
            n_dup++
            next
        }
        print $0, peak_id
    }
    END {
        if (n_dup) print "Removed " n_dup " duplicate peaks" > "/dev/stderr"
    }' > "${unzipped_output}"

bgzip -c "${unzipped_output}" > "${zipped_output}"
