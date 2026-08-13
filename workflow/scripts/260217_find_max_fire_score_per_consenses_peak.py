#!/usr/bin/env python

import sys
from pathlib import Path
from typing import Optional

import defopt

def main(
    overlaps_file: Path,
    consenseus_peak_file: Path,
    cons_peaks_transferred: Path,
    sample_id: str,
    output_file: Path,
    #infile: Optional[Path] = None,
    #outfile: Optional[Path] = None,
    *,
    threads: int = 32,
):
    """
    Author Anna Minkina

    :param infile: Input file, stdin by default
    :param outfile: Output file, stdout by default
    """

    out_dict = {}

    with open(overlaps_file, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue
            else:
                vect = line.rstrip("\n").split("\t")
                cons_peak_id = vect[-1]
                score = float(vect[5])
                cov = int(vect[3])
                fire_cov = int(vect[4])
                if cons_peak_id not in out_dict:
                    out_dict[cons_peak_id] = [cov, fire_cov, score]
                else:
                    curr_score = float(out_dict[cons_peak_id][2])
                    if score > curr_score:
                        out_dict[cons_peak_id] = [cov, fire_cov, score]
                    elif score == curr_score:
                        curr_fire_frac=float(out_dict[cons_peak_id][1])/float(out_dict[cons_peak_id][0])
                        new_fire_frac=float(fire_cov)/float(cov)
                        if new_fire_frac > curr_fire_frac:
                            out_dict[cons_peak_id] = [cov, fire_cov, score]
                        elif new_fire_frac == curr_fire_frac:
                            curr_fire_cov = float(out_dict[cons_peak_id][1])
                            new_fire_cov = float(fire_cov)
                            if new_fire_cov > curr_fire_cov:
                                out_dict[cons_peak_id] = [cov, fire_cov, score]
                            elif new_fire_cov == curr_fire_cov:
                                curr_cov = float(out_dict[cons_peak_id][0])
                                new_cov = float(cov)
                                if new_cov > curr_cov:
                                    out_dict[cons_peak_id] = [cov, fire_cov, score]

    with open(cons_peaks_transferred, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue
            else:
                vect = line.rstrip("\n").split("\t")
                cons_peak_id = vect[-1]
                if cons_peak_id not in out_dict:
                    out_dict[cons_peak_id] = ["NA", "NA", -1]
    
    with open(consenseus_peak_file, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue
            else:
                vect = line.rstrip("\n").split("\t")
                cons_peak_id = vect[-1]
                if cons_peak_id not in out_dict:
                    out_dict[cons_peak_id] = ["NA", "NA", -2]

    with open(output_file, 'w') as file:
        file.write(f"sample_id\tcoverage\tfire_coverage\tscore\n")
        for key in out_dict:
            to_write = f"{sample_id}\t{key}\t{out_dict[key][0]}\t{out_dict[key][1]}\t{out_dict[key][2]}\n"
            file.write(to_write)

    return 0

if __name__ == "__main__":
    defopt.run(main, show_types=True, version="0.0.1")