import sys
from pathlib import Path
from typing import Optional

import defopt
import pysam

def create_modified_header_from_dict_single_asm(old_header, asm):
    for line in old_header['SQ']:
        old_chr_name=line['SN']
        new_chr_name=asm + "#0#" + old_chr_name
        if new_chr_name != "NA":
            line['SN']=new_chr_name
    return old_header

def write_bam(bam,obam):
    for rec in bam.fetch(until_eof=True):
        obam.write(rec)

def main(
    #tags_output_file: Path,
    #tags_out_bam: Path,
    asm_to_surject_onto: str,
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

    bam_header_dict = bam.header.to_dict()
    new_header_dict = create_modified_header_from_dict_single_asm(bam_header_dict, asm_to_surject_onto)

    obam = pysam.AlignmentFile(outfile, "wb", header=new_header_dict, threads=threads)

    write_bam(bam,obam)

    obam.close()
    bam.close()
    #tags_obam.close()

    return 0

if __name__ == "__main__":
    defopt.run(main, show_types=True, version="0.0.1")