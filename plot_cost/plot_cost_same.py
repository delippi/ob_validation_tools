#!/usr/bin/env python3

import argparse
import glob
import os
import re
import sys

import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt

series1_color = "tab:red"
series2_color = "tab:blue"
diff_color = "black"


def outer_loop_lengths(records):
    lengths = {}
    for r in records:
        o = r["outer"]
        lengths[o] = lengths.get(o, 0) + 1
    return lengths


def resolve_one_path(pattern):
    matches = sorted(glob.glob(pattern))
    if len(matches) == 0:
        raise RuntimeError(f"No files match: {pattern}")
    if len(matches) > 1:
        msg = [f"Multiple files match: {pattern}"]
        msg.extend(matches)
        raise RuntimeError("\n".join(msg))
    return matches[0]


def parse_jedi_log(path):
    """
    Parse JEDI DRPCG cost-function output.

    Expected lines like:
      Quadratic cost function: J   ( 0) = ...
      Quadratic cost function: Jb  ( 0) = ...
      Quadratic cost function: JoJc( 0) = ...

    Returns records:
      {
        "outer": int,
        "inner": int,
        "cumiter": int,
        "J": float,
        "Jb": float or None,
        "JoJc": float or None,
      }

    Outer loop is inferred by detecting resets/decreases in inner iteration.
    """
    pat_j = re.compile(
        r"Quadratic cost function:\s*J\s*\(\s*(\d+)\s*\)\s*=\s*"
        r"([+-]?\d+(?:\.\d*)?(?:[Ee][+-]?\d+)?)"
    )
    pat_jb = re.compile(
        r"Quadratic cost function:\s*Jb\s*\(\s*(\d+)\s*\)\s*=\s*"
        r"([+-]?\d+(?:\.\d*)?(?:[Ee][+-]?\d+)?)"
    )
    pat_jojc = re.compile(
        r"Quadratic cost function:\s*JoJc\s*\(\s*(\d+)\s*\)\s*=\s*"
        r"([+-]?\d+(?:\.\d*)?(?:[Ee][+-]?\d+)?)"
    )

    with open(path, "r", errors="ignore") as f:
        lines = f.readlines()

    seq = []
    jb_pending = {}
    jojc_pending = {}

    for lineno, line in enumerate(lines):
        m = pat_j.search(line)
        if m:
            seq.append(
                {
                    "lineno": lineno,
                    "inner": int(m.group(1)),
                    "J": float(m.group(2)),
                    "Jb": None,
                    "JoJc": None,
                }
            )
            continue

        m = pat_jb.search(line)
        if m:
            inner = int(m.group(1))
            val = float(m.group(2))
            jb_pending.setdefault(inner, []).append((lineno, val))
            continue

        m = pat_jojc.search(line)
        if m:
            inner = int(m.group(1))
            val = float(m.group(2))
            jojc_pending.setdefault(inner, []).append((lineno, val))
            continue

    for rec in seq:
        inner = rec["inner"]

        if inner in jb_pending:
            for i, (lineno, val) in enumerate(jb_pending[inner]):
                if lineno >= rec["lineno"]:
                    rec["Jb"] = val
                    del jb_pending[inner][i]
                    break

        if inner in jojc_pending:
            for i, (lineno, val) in enumerate(jojc_pending[inner]):
                if lineno >= rec["lineno"]:
                    rec["JoJc"] = val
                    del jojc_pending[inner][i]
                    break

    out = []
    outer = 1
    last_inner = None

    for rec in seq:
        inner = rec["inner"]

        if last_inner is not None and inner <= last_inner:
            outer += 1

        out.append(
            {
                "outer": outer,
                "inner": inner,
                "cumiter": len(out),
                "J": rec["J"],
                "Jb": rec["Jb"],
                "JoJc": rec["JoJc"],
            }
        )
        last_inner = inner

    return out


def get_outer_loops(records):
    return sorted(set(r["outer"] for r in records))


def require_outer_loop(records, label, required_outer=2):
    outers = get_outer_loops(records)
    if required_outer not in outers:
        raise RuntimeError(
            f"{label} log does not contain outer loop {required_outer}. "
            f"Found outer loops: {outers}"
        )


def build_outer_boundaries(records):
    if not records:
        return []

    bounds = []
    prev_outer = records[0]["outer"]
    for rec in records[1:]:
        if rec["outer"] != prev_outer:
            bounds.append(rec["cumiter"] - 0.5)
            prev_outer = rec["outer"]
    return bounds


def default_label(path):
    return os.path.basename(path)


def plot_one_series(ax, records, color, label):
    xvals = [r["cumiter"] for r in records]
    yvals = [r["J"] for r in records]
    line, = ax.plot(
        xvals, yvals, marker="o", markersize=3, linewidth=1.5, color=color, label=label
    )
    for xb in build_outer_boundaries(records):
        ax.axvline(x=xb, linestyle="--", linewidth=0.8, alpha=0.5, color=color)
    return line


def plot_difference(ax2, rec1, rec2, label1, label2):
    n = min(len(rec1), len(rec2))
    if n == 0:
        return None

    xvals = [rec1[i]["cumiter"] for i in range(n)]
    ydiff = [rec1[i]["J"] - rec2[i]["J"] for i in range(n)]

    line, = ax2.plot(
        xvals,
        ydiff,
        marker="s",
        markersize=2.5,
        linewidth=1.2,
        color=diff_color,
        label=f"{label1} - {label2}",
    )
    max_abs = max(abs(v) for v in ydiff) if ydiff else 1.0
    max_abs = max_abs*1.2
    if max_abs == 0:
        max_abs = 1.0

    ax2.set_ylim(-max_abs, max_abs)
    ax2.axhline(0.0, color="black", linewidth=1.0, linestyle="--", alpha=0.6)
    ax2.set_ylabel("Difference in J")
    return line


def print_summary(label, records):
    yvals = [r["J"] for r in records]
    print(f"{label} first/last J: {yvals[0]}  {yvals[-1]}")
    lengths = outer_loop_lengths(records)
    for o in sorted(lengths):
        print(f"{label} outer loop {o} length: {lengths[o]}")


def main():
    parser = argparse.ArgumentParser(
        description="Parse one or two JEDI cost-function logs and plot raw J"
    )
    parser.add_argument(
        "--log1",
        required=True,
        help="Path or glob to first JEDI log file",
    )
    parser.add_argument(
        "--log2",
        help="Path or glob to second JEDI log file",
    )
    parser.add_argument(
        "--label1",
        help="Legend label for first log (default: basename of file)",
    )
    parser.add_argument(
        "--label2",
        help="Legend label for second log (default: basename of file)",
    )
    parser.add_argument(
        "--cdate",
        required=True,
        help="Cycle date (e.g. 2024050606)",
    )
    parser.add_argument(
        "-o",
        "--output",
        required=True,
        help="Output PNG filename",
    )
    parser.add_argument(
        "--title",
        default="JEDI raw cost function",
        help="Plot title",
    )
    parser.add_argument(
        "--require_outer2",
        action="store_true",
        help="Fail unless all provided logs contain outer loop 2",
    )
    args = parser.parse_args()

    try:
        log1_path = resolve_one_path(args.log1)
        log2_path = resolve_one_path(args.log2) if args.log2 else None

        rec1 = parse_jedi_log(log1_path)
        if not rec1:
            raise RuntimeError(f"No JEDI cost-function records found in {log1_path}")

        rec2 = None
        if log2_path is not None:
            rec2 = parse_jedi_log(log2_path)
            if not rec2:
                raise RuntimeError(f"No JEDI cost-function records found in {log2_path}")

        label1 = args.label1 if args.label1 else default_label(log1_path)
        label2 = None
        if log2_path is not None:
            label2 = args.label2 if args.label2 else default_label(log2_path)

        print(f"JEDI file 1: \n{log1_path}")
        print(f"JEDI 1 outer loops found: {get_outer_loops(rec1)}")
        print("")

        if log2_path is not None:
            print(f"JEDI file 2: \n{log2_path}")
            print(f"JEDI 2 outer loops found: {get_outer_loops(rec2)}")
            print("")

        if args.require_outer2:
            require_outer_loop(rec1, label1, required_outer=2)
            if rec2 is not None:
                require_outer_loop(rec2, label2, required_outer=2)

        fig = plt.figure(figsize=(10, 6))
        ax = fig.add_subplot(1, 1, 1)

        line1 = plot_one_series(ax, rec1, series1_color, label1)
        lines = [line1]

        if rec2 is not None:
            line2 = plot_one_series(ax, rec2, series2_color, label2)
            lines.append(line2)

            ax2 = ax.twinx()
            diff_line = plot_difference(ax2, rec1, rec2, label1, label2)
            if diff_line is not None:
                lines.append(diff_line)

        ax.set_title(f"{args.title} ({args.cdate})")
        print(f"{args.title} ({args.cdate})")
        print("")
        ax.set_xlabel("Cumulative iteration")
        ax.set_ylabel("Raw cost J")
        ax.grid(True, alpha=0.3)

        labels = [line.get_label() for line in lines]
        ax.legend(lines, labels, loc="best")

        fig.tight_layout()
        fig.savefig(args.output, dpi=200, bbox_inches="tight")

        print(f"Wrote: {args.output}")
        print("")

        print_summary(label1, rec1)
        if rec2 is not None:
            print("")
            print_summary(label2, rec2)

            n = min(len(rec1), len(rec2))
            if n > 0:
                first_diff = rec1[0]["J"] - rec2[0]["J"]
                last_diff = rec1[n - 1]["J"] - rec2[n - 1]["J"]
                print("")
                print(f"Difference first/last J ({label1} - {label2}): {first_diff}  {last_diff}")

    except Exception as e:
        print(f"ERROR: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
