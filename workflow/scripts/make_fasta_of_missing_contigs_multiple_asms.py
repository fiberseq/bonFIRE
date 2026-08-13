
import sys
from pathlib import Path
from typing import Optional

import defopt
from Bio import SeqIO
import gzip

def main(
    chrs_in_graph: Path,
    prefix_file: Path,
    input_fasta: Path,
    output_fasta: Path
    #*,
    #hap1_tag: str = "haplotype1",
    #hap2_tag: str = "haplotype2",
    #hap1_new_prefix: str = None,
    #hap2_new_prefix: str = None,
):
    """
    Author:Anna Minkina

    :param chrs_in_graph: List of chromosomes in graph
    :param fasta: full diploid fasta
    """

    with open(chrs_in_graph, "r") as chrs_in_graph_file:
        chrs_in_graph_list = [line.strip() for line in chrs_in_graph_file]

    prefix_dictionary = {}

    print("making dictionary")

    with open(prefix_file) as f:
        for line in f:
            #print(line)
            key, value, _ = line.rstrip("\n").split("\t")
            prefix_dictionary[key] = value

    print("done making dict")

    # Open the output file for writing
    with open(output_fasta, "w") as out_fasta:
        # Iterate through each record in the FASTA
        with open(input_fasta, "rt") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                new_record_name="NA"
                for key in prefix_dictionary:
                    if key in record.name:
                        if "#1#" in record.name:
                            new_record_name = prefix_dictionary[key] + "#0#" + record.name.split("#1#")[1]
                        elif "#2#" in record.name:
                            new_record_name = prefix_dictionary[key] + "#0#" + record.name.split("#2#")[1]
                        else: 
                            new_record_name = prefix_dictionary[key] + "#0#" + record.name
                        break
                if new_record_name not in chrs_in_graph_list and new_record_name != "NA":
                    record.id = new_record_name
                    record.description = ""
                    SeqIO.write(record, out_fasta, "fasta")

if __name__ == "__main__":
    defopt.run(main, show_types=True, version="0.0.1")

