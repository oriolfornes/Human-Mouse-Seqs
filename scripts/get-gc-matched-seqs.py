#!/usr/bin/env python

from Bio import SeqIO
from Bio.SeqUtils import GC
import click
import pandas as pd
import random

CONTEXT_SETTINGS = {
    "help_option_names": ["-h", "--help"],
}

@click.command(no_args_is_help=True, context_settings=CONTEXT_SETTINGS)
@click.argument(
    "human_fasta",
    type=click.Path(exists=True, resolve_path=True, allow_dash=True),
)
@click.argument(
    "mouse_fasta",
    type=click.Path(exists=True, resolve_path=True, allow_dash=True),
)
@click.option(
    "-f", "--filter-masked",
    help="Filter masked DNA sequences.",
    is_flag=True
)
@click.option(
    "-o", "--output-file",
    help="Output file.  [default: STDOUT]",
    type=click.Path(writable=True, readable=False, resolve_path=True,
        allow_dash=True),
    default="-",
)

def cli(**params):

    # Group sequences by %GC content
    fasta_files = ["human_fasta", "mouse_fasta"]
    gc_groups = {}
    for i in range(len(fasta_files)):
        for record in SeqIO.parse(params[fasta_files[i]], "fasta"):
            if params["filter_masked"]:
                if "N" in record.seq:
                    continue
                gc = round(GC(record.seq))
                gc_groups.setdefault(gc, [[], []])
                gc_groups[gc][i].append(record)

    # Downsampling
    sampled = []
    random_seed = 123
    for i in sorted(gc_groups):
        random.Random(random_seed).shuffle(gc_groups[i][0])
        random.Random(random_seed).shuffle(gc_groups[i][1])
        for j in range(min([len(gc_groups[i][0]), len(gc_groups[i][1])])):
            h = gc_groups[i][0][j]
            m = gc_groups[i][1][j]
            sampled.append([i, h.id, str(h.seq), m.id, str(m.seq)])

    # Write
    column_names = ["%GC", "HumanSeqId", "HumanSeq", "MouseSeqId", "MouseSeq"]
    df = pd.DataFrame(sampled, columns=column_names)
    with click.open_file(params["output_file"], "wt") as f:
        df.to_csv(f, sep="\t")

if __name__ == "__main__":
    cli()