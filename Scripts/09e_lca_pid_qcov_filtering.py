#!/usr/bin/env python3

import sys
from collections import defaultdict

# Rank order from finest to coarsest, with depth = number of lineage
# components to keep and default PID / qcov thresholds.
RANKS = [
    {"name": "species",      "depth": 8, "pidThreshold": 95, "qcovThreshold": 90},
    {"name": "genus",        "depth": 7, "pidThreshold": 92, "qcovThreshold": 0},
    {"name": "family",       "depth": 6, "pidThreshold": 80, "qcovThreshold": 0},
    {"name": "order",        "depth": 5, "pidThreshold": 70, "qcovThreshold": 0},
    {"name": "class",        "depth": 4, "pidThreshold": 70, "qcovThreshold": 0},
    {"name": "phylum",       "depth": 3, "pidThreshold": 70, "qcovThreshold": 0},
    {"name": "kingdom",      "depth": 2, "pidThreshold": 0,  "qcovThreshold": 0},
    {"name": "superkingdom", "depth": 1, "pidThreshold": 0,  "qcovThreshold": 0},
]

# Rank prefixes from depth 1 (superkingdom) to depth 8 (species).
PREFIXES = ["sk__", "k__", "p__", "c__", "o__", "f__", "g__", "s__"]


def parse_threshold_arg(arg):
    """
    Parse:
        'species=95,genus=92,family=85'
    into:
        {'species': 95, 'genus': 92, 'family': 85}
    """
    overrides = {}

    if not arg:
        return overrides

    for pair in arg.split(","):
        parts = pair.split("=")
        if len(parts) != 2:
            continue

        name = parts[0].strip().lower()
        value = parts[1].strip()

        try:
            overrides[name] = float(value)
        except ValueError:
            pass

    return overrides


def build_ranks(pid_override_str=None, qcov_override_str=None):
    pid_overrides = parse_threshold_arg(pid_override_str)
    qcov_overrides = parse_threshold_arg(qcov_override_str)

    ranks = []

    for rank in RANKS:
        new_rank = rank.copy()

        if rank["name"] in pid_overrides:
            new_rank["pidThreshold"] = pid_overrides[rank["name"]]

        if rank["name"] in qcov_overrides:
            new_rank["qcovThreshold"] = qcov_overrides[rank["name"]]

        ranks.append(new_rank)

    return ranks


def parse_lineage(lineage_str):
    return [
        x.strip()
        for x in lineage_str.split(";")
        if x.strip()
    ]


def lca(paths):
    """
    Lowest Common Ancestor across lineage paths.
    """
    if not paths:
        return []

    base = paths[0]
    common_len = len(base)

    for current in paths[1:]:
        j = 0

        while (
            j < common_len
            and j < len(current)
            and base[j] == current[j]
        ):
            j += 1

        common_len = j

        if common_len == 0:
            break

    return base[:common_len]


def pad_lineage(lineage):
    """
    Pad lineage to 8 ranks using rank prefixes.
    """
    padded = list(lineage)

    for i in range(len(padded), len(PREFIXES)):
        padded.append(PREFIXES[i])

    return padded


def ranked_lca(hits, ranks):
    """
    Walk ranks finest -> coarsest.

    At the first rank where some hits clear BOTH the PID and qcov
    thresholds, compute the LCA of those hits and truncate to that
    rank depth.
    """
    for rank in ranks:
        passing = [
            hit for hit in hits
            if hit["pid"] >= rank["pidThreshold"]
            and hit["qcov"] >= rank["qcovThreshold"]
        ]

        if not passing:
            continue

        common = lca([hit["lineage"] for hit in passing])
        return common[:rank["depth"]]

    return []


def main(input_file, pid_threshold_arg=None, qcov_threshold_arg=None):
    ranks = build_ranks(pid_threshold_arg, qcov_threshold_arg)

    groups = defaultdict(list)

    with open(input_file, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()

            if not line:
                continue

            fields = line.split("\t")

            if len(fields) < 2:
                continue

            identifier = fields[0]
            lineage = fields[1]

            pid = float(fields[2]) if len(fields) > 2 and fields[2] else 0.0
            qcov = float(fields[3]) if len(fields) > 3 and fields[3] else 0.0

            groups[identifier].append({
                "lineage": parse_lineage(lineage),
                "pid": pid,
                "qcov": qcov
            })

    for identifier, hits in groups.items():
        common = ranked_lca(hits, ranks)
        padded = pad_lineage(common)

        print(f"{identifier}\t{';'.join(padded)};")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit(
            "Usage: python lca.py input.tsv "
            "'species=95,genus=92,family=85' "
            "'species=90,genus=85,family=80'"
        )

    input_file = sys.argv[1]
    pid_arg = sys.argv[2] if len(sys.argv) > 2 else None
    qcov_arg = sys.argv[3] if len(sys.argv) > 3 else None

    main(input_file, pid_arg, qcov_arg)
