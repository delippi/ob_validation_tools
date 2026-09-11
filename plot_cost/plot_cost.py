#!/usr/bin/env python3

import argparse
import glob
import re
import sys

import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt

gsi_color = "tab:blue"
jedi_color = "tab:red"
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


def parse_gsi_log(path, divide_by_2=False):
    """
    Parse GSI minimization output.

    Expected lines like:
      cost,grad,step,b,step? = outer inner cost grad step bstep

    Returns records:
      {
        "outer": int,
        "inner": int,
        "cumiter": int,
        "J": float,
        "grad": float,
        "step": float,
        "bstep": float,
      }
    """
    pat = re.compile(
        r"cost,grad,step,b,step\?\s*=\s*"
        r"(\d+)\s+(\d+)\s+"
        r"([+-]?\d+(?:\.\d*)?(?:[Ee][+-]?\d+)?)\s+"
        r"([+-]?\d+(?:\.\d*)?(?:[Ee][+-]?\d+)?)\s+"
        r"([+-]?\d+(?:\.\d*)?(?:[Ee][+-]?\d+)?)\s+"
        r"([+-]?\d+(?:\.\d*)?(?:[Ee][+-]?\d+)?)"
    )

    out = []
    with open(path, "r", errors="ignore") as f:
        for line in f:
            m = pat.search(line)
            if not m:
                continue

            outer = int(m.group(1))
            inner = int(m.group(2))
            cost = float(m.group(3))
            grad = float(m.group(4))
            step = float(m.group(5))
            bstep = float(m.group(6))

            if divide_by_2:
                cost *= 0.5

            out.append(
                {
                    "outer": outer,
                    "inner": inner,
                    "cumiter": len(out),
                    "J": cost,
                    "grad": grad,
                    "step": step,
                    "bstep": bstep,
                }
            )

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


def plot_difference(ax2, jedi, gsi):
    n = min(len(jedi), len(gsi))
    if n == 0:
        return None

    xvals = [jedi[i]["cumiter"] for i in range(n)]
    ydiff = [jedi[i]["J"] - gsi[i]["J"] for i in range(n)]

    line, = ax2.plot(
        xvals,
        ydiff,
        marker="^",
        markersize=2.5,
        linewidth=1.2,
        color=diff_color,
        label="JEDI - GSI",
    )
    max_abs = max(abs(v) for v in ydiff) if ydiff else 1.0
    max_abs = max_abs*1.2
    if max_abs == 0:
        max_abs = 1.0

    ax2.set_ylim(-max_abs, max_abs)
    ax2.axhline(0.0, color="black", linewidth=1.0, linestyle="--", alpha=0.6)
    ax2.set_ylabel("Difference in J")
    return line


def main():
    parser = argparse.ArgumentParser(
        description="Parse JEDI and GSI cost-function logs and plot raw J"
    )
    parser.add_argument(
        "--jedi",
        required=True,
        help="Path or glob to JEDI log file",
    )
    parser.add_argument(
        "--gsi",
        required=True,
        help="Path or glob to GSI log file",
    )
    parser.add_argument(
        "--cdate",
        required=True,
        help="Cycle date (e.g., 2024050606)",
    )
    parser.add_argument(
        "-o",
        "--output",
        required=True,
        help="Output PNG filename",
    )
    parser.add_argument(
        "--title",
        default="GSI vs JEDI raw cost function",
        help="Plot title",
    )
    parser.add_argument(
        "--require_outer2",
        action="store_true",
        help="Fail unless both logs contain outer loop 2",
    )
    parser.add_argument(
        "--gsi_half",
        action="store_true",
        help="Divide GSI printed cost by 2 before plotting",
    )
    args = parser.parse_args()

    try:
        jedi_path = resolve_one_path(args.jedi)
        gsi_path = resolve_one_path(args.gsi)

        jedi = parse_jedi_log(jedi_path)
        gsi = parse_gsi_log(gsi_path, divide_by_2=args.gsi_half)

        if not jedi:
            raise RuntimeError(f"No JEDI cost-function records found in {jedi_path}")
        if not gsi:
            raise RuntimeError(f"No GSI cost-function records found in {gsi_path}")

        print(f"JEDI file: \n{jedi_path}")
        print(f"GSI  file: \n{gsi_path}")
        print("")
        print(f"JEDI outer loops found: {get_outer_loops(jedi)}")
        print(f"GSI  outer loops found: {get_outer_loops(gsi)}")
        print("")

        if args.gsi_half:
            print("Applied GSI cost scaling: plotted GSI cost as 0.5 * printed cost")
            print("")

        if args.require_outer2:
            require_outer_loop(jedi, "JEDI", required_outer=2)
            require_outer_loop(gsi, "GSI", required_outer=2)

        xj = [r["cumiter"] for r in jedi]
        yj = [r["J"] for r in jedi]

        xg = [r["cumiter"] for r in gsi]
        yg = [r["J"] for r in gsi]

        fig = plt.figure(figsize=(10, 6))
        ax = fig.add_subplot(1, 1, 1)

        jedi_line, = ax.plot(
            xj,
            yj,
            marker="o",
            markersize=3,
            linewidth=1.5,
            color=jedi_color,
            label="JEDI J",
        )
        gsi_line, = ax.plot(
            xg,
            yg,
            marker="s",
            markersize=3,
            linewidth=1.5,
            color=gsi_color,
            label="GSI J",
        )

        for xb in build_outer_boundaries(jedi):
            ax.axvline(x=xb, linestyle="--", linewidth=0.8, alpha=0.5, color=jedi_color)

        for xb in build_outer_boundaries(gsi):
            ax.axvline(x=xb, linestyle=":", linewidth=0.8, alpha=0.5, color=gsi_color)

        ax2 = ax.twinx()
        diff_line = plot_difference(ax2, jedi, gsi)

        ax.set_title(f"{args.title} ({args.cdate})")
        print(f"{args.title} ({args.cdate})")
        print("")
        ax.set_xlabel("Cumulative iteration")
        ax.set_ylabel("Raw cost J")
        ax.grid(True, alpha=0.3)

#        ax2 = ax.twinx()
#        diff_line = plot_difference(ax2, jedi, gsi)
#
#        # Show the actual blue-line values inside the saved image.
#        blue_y = gsi_line.get_ydata()
#        blue_summary = (f"GSI plotted J: {blue_y[0]:,.2f} to {blue_y[-1]:,.2f}")
#        ax.set_title(f"{args.title} ({args.cdate})\n{blue_summary}")
#        print(f"{args.title} ({args.cdate})")
#        print(blue_summary)
#        print("")
        ax.set_xlabel("Cumulative iteration")
        ax.set_ylabel("Raw cost J")
        ax.grid(True, alpha=0.3)

        lines = [jedi_line, gsi_line]
        if diff_line is not None:
            lines.append(diff_line)
        labels = [line.get_label() for line in lines]
        ax.legend(lines, labels, loc="best")

        fig.tight_layout()
        fig.savefig(args.output, dpi=200, bbox_inches="tight")

        print(f"Wrote: {args.output}")
        print("")

        print(f"JEDI first/last J: {yj[0]}  {yj[-1]}")
        jedi_lengths = outer_loop_lengths(jedi)
        for o in sorted(jedi_lengths):
            print(f"JEDI outer loop {o} length: {jedi_lengths[o]}")

        print("")

        print(f"GSI  first/last J: {yg[0]}  {yg[-1]}")
        gsi_lengths = outer_loop_lengths(gsi)
        for o in sorted(gsi_lengths):
            print(f"GSI  outer loop {o} length: {gsi_lengths[o]}")

        n = min(len(jedi), len(gsi))
        if n > 0:
            first_diff = jedi[0]["J"] - gsi[0]["J"]
            last_diff = jedi[n - 1]["J"] - gsi[n - 1]["J"]
            print("")
            print(f"Difference first/last J (JEDI - GSI): {first_diff}  {last_diff}")

    except Exception as e:
        print(f"ERROR: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
