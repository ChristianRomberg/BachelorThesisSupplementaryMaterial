import argparse
import multiprocessing
import re
import shutil
import tarfile
import traceback
from collections import Counter, defaultdict
from io import TextIOWrapper
from pathlib import Path
import json
import lzma

import paths

input_folder = Path("humann_results")

output_folder = paths.cohort1_per_sample_results


def process_log(log_stream):
    log = log_stream.read()

    # b'200000 reads; of these:
    total_match = re.search(r"(\d+) reads; of these:", log)
    total_reads = int(total_match.group(1))

    # Unaligned reads after nucleotide alignment: 74.7759751156 %
    nucl_unaligned_percentage_match = re.search(r"Unaligned reads after nucleotide alignment: (\d+\.\d+) %", log)
    nucleotide_unaligned_reads = int(total_reads * float(nucl_unaligned_percentage_match.group(1)) / 100)

    # Unaligned reads after translated alignment: 98.9000000000 %
    transl_unaligned_percentage_match = re.search(r"Unaligned reads after translated alignment: (\d+\.\d+) %", log)
    translated_unaligned_reads = int(total_reads * float(transl_unaligned_percentage_match.group(1)) / 100)

    return total_reads, nucleotide_unaligned_reads, translated_unaligned_reads


def process_diamond_file(diamond_stream):
    gene_families = Counter(line.split("\t")[1].split("|")[0] for line in diamond_stream)
    return gene_families


def extract_counts(results_folder: Path):
    print(f"Parsing data for accession {results_folder.name}")
    item_output_folder = output_folder / results_folder.name
    try:
        log_path = results_folder / "input.log.xz"
        diamond_aligned_path = results_folder / "diamond_aligned.tsv.xz"
        with lzma.open(log_path, "rt") as log_reader:
            total_reads, nucleotide_unaligned_reads, translated_unaligned_reads = process_log(log_reader)
        with lzma.open(diamond_aligned_path, "rt") as diamond_reader:
            gene_family_counter = process_diamond_file(diamond_reader)

        item_output_folder.mkdir(exist_ok=True, parents=True)

        diamond_condensed_path = item_output_folder / "diamond_condensed.tsv"
        with diamond_condensed_path.open("w") as diamond_output:
            lines_iterator = (f"{gene_family}\t{gf_count}\n" for gene_family, gf_count in gene_family_counter.items())
            diamond_output.writelines(lines_iterator)

        counts_path = item_output_folder / "counts.tsv"
        with counts_path.open("w") as counts_output:
            lines = ["type\tcount\n",
                     f"total\t{total_reads}\n",
                     f"nucleotide_unaligned\t{nucleotide_unaligned_reads}\n",
                     f"translated_unaligned\t{translated_unaligned_reads}\n", ]
            counts_output.writelines(lines)

        return gene_family_counter, (total_reads, nucleotide_unaligned_reads, translated_unaligned_reads)

    except Exception as e:
        print(f"Failed for folder {results_folder}")
        traceback.print_exc()
        shutil.rmtree(item_output_folder, ignore_errors=True)
        return dict(), (0,0,0)


total_gene_families_count = defaultdict(int)
total_reads_count = 0
total_nucleotide_unaligned_count = 0
total_translated_unaligned_count = 0

pool = multiprocessing.Pool(10)

for gene_family_counter, (total_reads, nucleotide_unaligned_reads, translated_unaligned_reads) \
    in pool.imap_unordered(extract_counts, input_folder.iterdir()):

    for gene_family, count in gene_family_counter.items():
        total_gene_families_count[gene_family] += count

    total_reads_count += total_reads
    total_nucleotide_unaligned_count += nucleotide_unaligned_reads
    total_translated_unaligned_count += translated_unaligned_reads

total_gene_families_count = list(total_gene_families_count.items())
total_gene_families_count.sort(reverse=True, key=lambda x: x[1])

with paths.gene_family_counts.open("w") as f:
    f.write("gene_family\tcount\n")
    for amino_acid, count in total_gene_families_count:
        f.write(f"{amino_acid}\t{count}\n")

with paths.counts_summary.open("w") as f:
    lines = ["type\tcount\n",
             f"total\t{total_reads_count}\n",
             f"nucleotide_unaligned\t{total_nucleotide_unaligned_count}\n",
             f"translated_unaligned\t{total_translated_unaligned_count}\n", ]
    f.writelines(lines)

