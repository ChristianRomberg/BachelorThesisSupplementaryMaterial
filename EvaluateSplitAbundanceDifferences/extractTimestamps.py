#!/usr/bin/env python3
import lzma
import re
from datetime import datetime
from pathlib import Path

import pandas as pd

import paths

slow_path = Path(__file__).parent / "humann_results_slow"
fast_path = Path(__file__).parent / "humann_results_fast"

input_paths = {
    "slow": slow_path,
    "fast": fast_path,
}

def extract_timestamps(accession_directory):
    print("Reading logs from", accession_directory.name)
    log_file = next(accession_directory.rglob("*.log.xz"))

    with lzma.open(log_file, "rt") as f:
        log = f.read()
        timestamps = re.finditer("INFO: TIMESTAMP: Completed\s+([A-Za-z -]+)\s+:\s+(\d+)\s+seconds", log)
        data = {name.strip(): int(seconds) for name, seconds in (ts.groups() for ts in timestamps)}
        data["accession"] = accession_directory.name
        try:
            small_starttime, small_endtime, large_endtime = [
                datetime.strptime(s, "%m/%d/%Y %I:%M:%S %p")
                for s in re.findall(
                    r"([0-9]+/[0-9]+/[0-9]+ [0-9]+:[0-9]+:[0-9]+ [AP]M) - (?:humann\.search\.translated - INFO: Aligning to reference database|humann\.utilities - DEBUG: Using software: /usr/bin/cat)",
                    log
                )
            ]
            data["translated1"] = (small_endtime - small_starttime).total_seconds()
            data["translated2"] = (large_endtime - small_endtime).total_seconds()
        except ValueError:
            pass
        return data

def extract_all_timestamps():
    all_data = []
    for mode, path in input_paths.items():
        for accession_path in path.iterdir():
            data = extract_timestamps(accession_path)
            data["mode"] = mode
            all_data.append(data)

    df = pd.DataFrame(all_data)
    df.set_index(["accession", "mode"], inplace=True, drop=True)
    return df

extract_all_timestamps().to_csv(paths.evaluation_times_by_sample)
