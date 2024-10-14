from pathlib import Path

import paths

split_results_folder = Path(__file__).parent / "split_results"

#!/usr/bin/env python3
import re
import traceback
from pathlib import Path
from datetime import datetime

import pandas as pd


def extract_timestamps(benchmark_directory):
    print("Reading logs from", benchmark_directory.name)
    log_file = next(benchmark_directory.rglob("*.log"))

    with log_file.open() as f:
        log = f.read()
        timestamps = re.finditer("INFO: TIMESTAMP: Completed\s+([A-Za-z ]+)\s+:\s+(\d+)\s+seconds", log)
        data = {name.strip(): int(seconds) for name, seconds in (ts.groups() for ts in timestamps)}
        data["split"] = int(benchmark_directory.name)

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

if __name__ == '__main__':
    current_folder = Path(__file__).parent.resolve()
    results_dicts = []
    for benchmark_directory in split_results_folder.iterdir():
        try:
            results_dicts.append(extract_timestamps(benchmark_directory))
        except Exception:
            traceback.print_exc()

    print(results_dicts)
    timestamps = pd.DataFrame.from_records(results_dicts)
    print(timestamps.head())
    timestamps.set_index("split", inplace=True, drop=True)
    timestamps.sort_index(inplace=True)
    timestamps.to_csv(paths.evaluation_times_by_split)