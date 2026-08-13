import sys
from pathlib import Path
from typing import Optional

import defopt
import pysam

def create_modified_header_from_dict(old_header, prefix_dict):
    for line in old_header['SQ']:
        old_chr_name=line['SN']
        new_chr_name="NA"
        for key in prefix_dict:
            if key in old_chr_name:
                if "#1#" in old_chr_name:
                    new_chr_name = prefix_dict[key] + "#0#" + old_chr_name.split("#1#")[1]
                elif "#2#" in old_chr_name:
                    new_chr_name = prefix_dict[key] + "#0#" + old_chr_name.split("#2#")[1]
                else:
                    new_chr_name=prefix_dict[key] + "#0#" + old_chr_name
                break
        if new_chr_name != "NA":
            line['SN']=new_chr_name
    return old_header

def write_bam(bam,obam):
    for rec in bam.fetch(until_eof=True):
        obam.write(rec)

def main(
    #tags_output_file: Path,
    #tags_out_bam: Path,
    prefix_file: Path,
    infile: Optional[Path] = None,
    outfile: Optional[Path] = None,
    *,
    threads: int = 32,
    min_mapq: int = 0
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
            key, value, _ = line.rstrip("\n").split("\t")
            prefix_dictionary[key] = value

    bam_header_dict = bam.header.to_dict()
    new_header_dict = create_modified_header_from_dict(bam_header_dict, prefix_dictionary)

    obam = pysam.AlignmentFile(outfile, "wb", header=new_header_dict, threads=threads)

    write_bam(bam,obam)

    obam.close()
    bam.close()
    #tags_obam.close()

    return 0

if __name__ == "__main__":
    defopt.run(main, show_types=True, version="0.0.1")