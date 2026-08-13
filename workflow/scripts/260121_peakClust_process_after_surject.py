#!/usr/bin/env python

import sys
from pathlib import Path
from typing import Optional

import defopt
import pysam

## in this script we'll
### modify header
## add original tags
## move flag tag
## add secondary alignment tag (if applicable?)

def create_modified_header_from_surjected_bam(old_header, prefix_dict):
    for line in old_header['SQ']:
        old_chr_name=line['SN']
        new_chr_name="NA"

        for key in prefix_dict:
            if key in old_chr_name:
                if "#1#" in prefix_dict[key] or "#2#" in prefix_dict[key]:
                    new_chr_name = prefix_dict[key] + old_chr_name.split("#0#")[1]
                else:
                    new_chr_name = old_chr_name.split("#0#")[1]
                break
            else:
                new_chr_name = old_chr_name.split("#0#")[1]

        if new_chr_name != "NA":
            line['SN']=new_chr_name
        # remove M5
        line.pop('M5', None)
    return old_header

def write_bam(bam,obam):
    for rec in bam.fetch(until_eof=True):
        obam.write(rec)


def main(
    prefix_file: Path,
    infile: Optional[Path] = None,
    outfile: Optional[Path] = None,
    *,
    threads: int = 32,
):
    """
    Author Anna Minkina

    :param infile: Input file, stdin by default
    :param outfile: Output file, stdout by default
    """

    if infile is None:
        infile = sys.stdin
    if outfile is None:
        outfile = sys.stdout.buffer

    if str(infile)[-4:] == "cram":
        bam = pysam.AlignmentFile(infile, "rc", threads=threads)
    elif str(infile)[-3:] == "bam":
        bam = pysam.AlignmentFile(infile, "rb", threads=threads)

    prefix_dictionary = {}

    with open(prefix_file) as f:
        for line in f:
            value, key, _ = line.rstrip("\n").split("\t")
            prefix_dictionary[key] = value


    bam_header_dict = bam.header.to_dict()
    new_header_dict = create_modified_header_from_surjected_bam(bam_header_dict, prefix_dictionary) 

    obam = pysam.AlignmentFile(outfile, "wb", header=new_header_dict, threads=threads)

    #obam_with_secondary_flag = pysam.AlignmentFile(outfile_w_secondary_tags, "wb", header=new_header_dict, threads=threads)

    #reassign_supp_flag_and_add_secondary_flag_and_add_original_tags(bam,obam,obam_with_secondary_flag,hap1_new_prefix, hap2_new_prefix, bam_with_original_tags)
    write_bam(bam,obam)

    obam.close()
    bam.close()

    return 0

if __name__ == "__main__":
    defopt.run(main, show_types=True, version="0.0.1")

    