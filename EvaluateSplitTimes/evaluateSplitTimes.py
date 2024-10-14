#!/usr/bin/env python3
import shutil
import subprocess
from multiprocessing import Pool
from pathlib import Path

import pandas as pd
from Bio import SeqIO

import paths

cutoffs = range(100000, 1600001, 100000)[::2]

aa_counts = pd.read_csv(paths.gene_family_counts, index_col=0)
uniref_path = Path("uniref90_sequences.fasta")
template_path = Path("template")
humann_input_fastq = Path("input.fastq")


def delete_template_files(split_path):
    # for file_to_delete in template_path.rglob("*"):
    #    file_to_delete = split_path / file_to_delete.relative_to(template_path)
    #    if file_to_delete.is_file():
    #        print(f"deleting {file_to_delete}")
    #        file_to_delete.unlink()
    # shutil.rmtree(split_path / "db")

    # delete everything except for the log containing the timestamps
    for file in split_path.rglob("*"):
        if not file.name.endswith(".log"):
            file.unlink()


def evaluate_cutoff(cutoff):
    small_db_whitelist = set(aa_counts.index[:cutoff])

    uniref_input_db = uniref_path
    current_folder = Path(__file__).parent
    split_path = current_folder / "split_results" / str(cutoff)
    if split_path.exists():
        delete_template_files(split_path)
        return

    shutil.copytree(template_path, split_path)

    db_path = split_path / "db"
    uniref_small_output_db = db_path / "1_small_201901b.fasta"
    uniref_large_output_db = db_path / "2_large_201901b.fasta"

    db_path.mkdir(exist_ok=True, parents=True)

    print("Writing split database")
    small_db_records = (r for r in SeqIO.parse(uniref_input_db, "fasta") if
                        r.id.split("|")[0] in small_db_whitelist)
    count = SeqIO.write(small_db_records, uniref_small_output_db, "fasta-2line")
    print(f"Wrote {count} sequences to small db")

    large_db_records = (r for r in SeqIO.parse(uniref_input_db, "fasta") if
                        r.id.split("|")[0] not in small_db_whitelist)
    count = SeqIO.write(large_db_records, uniref_large_output_db, "fasta-2line")
    print(f"Wrote {count} sequences to large db")

    cmd = ["diamond",
           "makedb",
           "--in", str(uniref_small_output_db),
           "-d", str(uniref_small_output_db.with_suffix(".dmnd")),
           ]
    subprocess.check_call(cmd)
    uniref_small_output_db.unlink()

    cmd = ["diamond",
           "makedb",
           "--in", str(uniref_large_output_db),
           "-d", str(uniref_large_output_db.with_suffix(".dmnd")),
           ]
    subprocess.check_call(cmd)
    uniref_large_output_db.unlink()

    print("Running humann")
    cmd = ["humann",
           "-i", str(humann_input_fastq.resolve()),
           "-o", str(split_path.resolve()),
           "--threads", "20",
           "--translated-alignment", "diamond_multistage",
           "--protein-database", str(db_path.resolve()),
           "--resume"
           ]
    subprocess.check_call(cmd)
    delete_template_files(split_path)


if __name__ == '__main__':
    with Pool(processes=5) as pool:
        pool.map(evaluate_cutoff, cutoffs)
