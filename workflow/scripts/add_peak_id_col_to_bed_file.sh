#!/bin/bash   

input_bed=$1
unzipped_output=$2
zipped_output=$3

(
    gzip -cd "${input_bed}" 2>/dev/null || cat "${input_bed}"
) | awk -F'\t' '
    BEGIN { OFS="\t" }
    NR == 1 {
        if ($2 ~ /^[0-9]+$/) {
            header = "#chrom\tstart\tend"
            for (i = 4; i <= NF; i++) {
                header = header "\tcol" i
            }
            print header, "peak_id"
        } else {
            print $0, "peak_id"
            next
        }
    }
    {
        peak_id = $1"_"$2"_"$3
        print $0, peak_id
    }' > ${unzipped_output}

    colname="peak_id"
    col_num=$(awk -F '\t' -v colname="$colname" ' 
    NR==1 {
        for (i=1; i <= NF; i++)
            if ($i == colname) 
                print i
    }' "${unzipped_output}")

    duplicates=$(tail -n +2 "${unzipped_output}" \
    | cut -f"$col_num" \
    | sort \
    | uniq -c \
    | awk '$1 > 1 {print $2}')

    dup_array=($duplicates)

    if [[ -n "$duplicates" ]]; then
        echo "Duplicates found:"
        grep -v -F -f <(printf "%s\n" "${dup_array[@]}") ${unzipped_output} |  bgzip > ${zipped_output}
    else
        bgzip -c ${unzipped_output} > ${zipped_output}
    fi