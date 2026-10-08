#!/usr/bin/env python3

import argparse
import json
import re
import statistics
import sys
import os
from datetime import datetime
from pathlib import Path

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np


DURATION_PATTERN = re.compile(
    r"Done script, total duration\s+([\d.]+)\s+seconds"
)


def collect_durations(directory: Path) -> list[dict]:
    """Recursively extract durations from all .out files."""
    records = []

    for filename in sorted(directory.rglob("*.out")):
        try:
            with filename.open(encoding="utf-8", errors="ignore") as f:
                for line_number, line in enumerate(f, start=1):
                    match = DURATION_PATTERN.search(line)

                    if match:
                        records.append(
                            {
                                "duration_seconds": float(match.group(1)),
                                "file": str(filename),
                                "line": line_number,
                            }
                        )

        except OSError as error:
            print(f"Warning: could not read {filename}: {error}")

    return records


def calculate_statistics(durations: list[float]) -> dict:
    """Calculate summary statistics."""
    result = {
        "number_of_jobs": len(durations),
        "mean_seconds": statistics.mean(durations),
        "median_seconds": statistics.median(durations),
        "minimum_seconds": min(durations),
        "maximum_seconds": max(durations),
        "standard_deviation_seconds": None,
    }

    if len(durations) > 1:
        result["standard_deviation_seconds"] = statistics.stdev(durations)

    return result


def print_statistics(stats: dict) -> None:
    """Print summary statistics."""
    print(f"Number of jobs : {stats['number_of_jobs']}")
    print(f"Average        : {stats['mean_seconds']:.2f} seconds")
    print(f"Median         : {stats['median_seconds']:.2f} seconds")
    print(f"Minimum        : {stats['minimum_seconds']:.2f} seconds")
    print(f"Maximum        : {stats['maximum_seconds']:.2f} seconds")

    std = stats["standard_deviation_seconds"]

    if std is not None:
        print(f"Std. deviation : {std:.2f} seconds")
    else:
        print("Std. deviation : n/a")


def write_json(
    records: list[dict],
    output_dir: str,
) -> None:
    durations = [record["duration_seconds"] for record in records]
    stats = calculate_statistics(durations)

    data = {
        "statistics": stats,
        "records": records,
    }

    with open(f"{output_dir}/time_distribution.json", "w") as f:
        json.dump(data, f)


def read_json(json_file: Path) -> tuple[list[dict], dict]:
    """Read durations and metadata from an existing JSON file."""
    try:
        with json_file.open(encoding="utf-8") as f:
            data = json.load(f)
    except OSError as error:
        print(f"Error reading '{json_file}': {error}")
        sys.exit(1)
    except json.JSONDecodeError as error:
        print(f"Invalid JSON in '{json_file}': {error}")
        sys.exit(1)

    records = data.get("records")

    if not isinstance(records, list):
        print(f"Error: '{json_file}' does not contain a valid records list.")
        sys.exit(1)

    valid_records = []

    for record in records:
        try:
            duration = float(record["duration_seconds"])
        except (KeyError, TypeError, ValueError):
            print(f"Warning: skipping invalid JSON record: {record}")
            continue

        valid_record = dict(record)
        valid_record["duration_seconds"] = duration
        valid_records.append(valid_record)

    if not valid_records:
        print(f"No valid durations found in '{json_file}'.")
        sys.exit(0)

    return valid_records, data


def plot_durations(
    durations: list[float],
    output_dir: str | None = None,
) -> None:
    """Plot the runtime distribution."""
    values = np.asarray(durations, dtype=float)

    mean = np.mean(values)
    median = np.median(values)
    std = np.std(values, ddof=1) if len(values) > 1 else None

    plt.figure(figsize=(8, 6))

    # Freedman–Diaconis automatic binning.
    # Fall back to a single bin if all values are identical.
    bins = "fd" if np.ptp(values) > 0 else 1

    plt.hist(
        values,
        bins=bins,
        edgecolor="black",
        alpha=0.75,
    )

    plt.axvline(
        mean,
        linestyle="--",
        linewidth=2,
        label=f"Mean = {mean:.1f} s",
    )

    plt.axvline(
        median,
        linestyle=":",
        linewidth=2,
        label=f"Median = {median:.1f} s",
    )

    title = f"Job runtime distribution\nN = {len(values)}"

    if std is not None:
        title += f", standard deviation = {std:.1f} s"

    plt.xlabel("Total duration [seconds]")
    plt.ylabel("Number of jobs")
    plt.title(title)
    plt.legend()
    plt.grid(axis="y", alpha=0.3)
    plt.tight_layout()

    plt.savefig(f"{output_dir}/time_distribution.png", dpi=150)
    plt.close()


def main() -> None:

    input_dir = "logdir/ipc_primary_grid_studies/FCCee_Z_GHC_V25p1/CFG_GRIDD1_256_256_256"

    input_dir = input_dir.rstrip('\\').rstrip('/').replace("//", "/")
    parameter_set = input_dir.split("/")[-1]
    accelerator = input_dir.split("/")[-2]
    campaign =  input_dir.split("/")[-3]


    output_dir = f"/home/submit/jaeyserm/public_html/fccee/guineapig/validation/{campaign}/{accelerator}/{parameter_set}/"
    os.system(f"cp /home/submit/jaeyserm/public_html/fccee/guineapig/validation/index.php {output_dir}")
    os.system(f"mkdir -p {output_dir}")

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--fromJson",
        action="store_true",
        help="Take values from JSON file",
    )
    args = parser.parse_args()

    if args.fromJson:
        records, json_data = read_json(args.fromJson)

        source = json_data.get("source_directory", "unknown")
        print(f"Loaded JSON     : {args.from_json}")
        print(f"Original source : {source}")

    else:
        input_dir_path = Path(input_dir)
        records = collect_durations(input_dir_path)

        if not records:
            print(
                f"No matching durations found under "
                f"'{input_dir}'."
            )
            sys.exit(0)

        print(f"Directory       : {input_dir}")

        write_json(
            records=records,
            output_dir=output_dir,
        )

    durations = [
        record["duration_seconds"]
        for record in records
    ]

    stats = calculate_statistics(durations)
    print_statistics(stats)
    plot_durations(durations, output_dir)


if __name__ == "__main__":
    main()