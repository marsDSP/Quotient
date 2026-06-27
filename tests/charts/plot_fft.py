#!/usr/bin/env python3
"""Render FFT verification charts from the harness-exported CSV data.

The CSVs are produced by the FFT test harness:

    quotient_fft_tests --export-charts --chart-dir tests/charts/data

and this script turns them into PNG figures next to (by default) this file.

    python3 tests/charts/plot_fft.py --data tests/charts/data --out tests/charts
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import sys

import matplotlib

matplotlib.use("Agg")  # headless / CI friendly
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def read_csv(path: str) -> tuple[list[str], dict[str, np.ndarray]]:
    """Return (header, columns-dict) from a CSV file."""
    with open(path, newline="") as fh:
        reader = csv.reader(fh)
        header = next(reader)
        rows = [[float(x) for x in row] for row in reader if row]
    data = np.array(rows, dtype=float) if rows else np.empty((0, len(header)))
    columns = {name: data[:, i] for i, name in enumerate(header)}
    return header, columns


# ---------------------------------------------------------------------------
# Plot: magnitude spectrum of a pure tone
# ---------------------------------------------------------------------------

def plot_spectrum(path: str, out: str) -> str:
    header, cols = read_csv(path)
    bins = cols["bin"].astype(int)
    mag  = cols["magnitude_db"]
    N    = len(bins)

    # Only show the positive-frequency half (bins 0 … N/2)
    half   = N // 2 + 1
    b_half = bins[:half]
    m_half = mag[:half]

    fig, ax = plt.subplots(figsize=(10, 5))

    # Noise floor as a filled region
    noise_floor = np.percentile(m_half, 10)
    ax.fill_between(b_half, noise_floor - 5, m_half,
                    where=(m_half < noise_floor + 3),
                    color="tab:blue", alpha=0.15, label="noise floor")

    # Main spectrum line
    ax.plot(b_half, m_half, color="tab:blue", lw=0.8, alpha=0.7)

    # Highlight the peaks above noise floor + 20 dB
    threshold = noise_floor + 20.0
    peaks = b_half[m_half > threshold]
    ax.vlines(peaks, noise_floor - 5, m_half[m_half > threshold],
              color="tab:orange", lw=2.0, label=f"tone bin(s): {peaks.tolist()}")

    ax.set_xlabel("Bin index")
    ax.set_ylabel("Magnitude (dB)")
    ax.set_title(f"FFT spectrum — cosine tone, N={N}")
    ax.set_xlim(0, half - 1)
    ax.set_ylim(noise_floor - 10, 5)
    ax.grid(True, ls=":", alpha=0.5)
    ax.legend(fontsize=9)
    fig.tight_layout()
    fig.savefig(out, dpi=130)
    plt.close(fig)
    return out


# ---------------------------------------------------------------------------
# Plot: forward and round-trip error vs N  (log-log)
# ---------------------------------------------------------------------------

def plot_accuracy(path: str, out: str) -> str:
    header, cols = read_csv(path)
    N       = cols["N"]
    fwd_err = cols["fwd_max_error"]
    rt_err  = cols["roundtrip_max_error"]

    fig, ax = plt.subplots(figsize=(9, 5))

    ax.loglog(N, fwd_err,  "o-", color="tab:blue",   lw=1.6, label="forward max error")
    ax.loglog(N, rt_err,   "s-", color="tab:orange",  lw=1.6, label="round-trip max error")

    # Reference: k · N · log2(N) · ε  (machine epsilon for double ≈ 2.2e-16)
    eps = np.finfo(float).eps
    ref = N * np.log2(N) * eps
    ax.loglog(N, ref, "--", color="gray", lw=1.2, label=r"$N \cdot \log_2 N \cdot \varepsilon$  (ref)")

    ax.set_xlabel("Transform size N")
    ax.set_ylabel("Max absolute error")
    ax.set_title("FFT numerical accuracy vs N")
    ax.legend(fontsize=9)
    ax.grid(True, which="both", ls=":", alpha=0.5)
    # x-axis tick labels as powers of 2
    ax.set_xticks(N)
    ax.set_xticklabels([str(int(n)) for n in N], rotation=45, fontsize=8)
    fig.tight_layout()
    fig.savefig(out, dpi=130)
    plt.close(fig)
    return out


# ---------------------------------------------------------------------------
# Plot: Parseval relative error vs N  (log-log)
# ---------------------------------------------------------------------------

def plot_parseval(path: str, out: str) -> str:
    header, cols = read_csv(path)
    N        = cols["N"]
    rel_err  = cols["parseval_relative_error"]

    fig, ax = plt.subplots(figsize=(9, 5))
    ax.loglog(N, rel_err, "D-", color="tab:green", lw=1.6, label="Parseval relative error")

    eps = np.finfo(float).eps
    ref = np.log2(N) * eps
    ax.loglog(N, ref, "--", color="gray", lw=1.2, label=r"$\log_2 N \cdot \varepsilon$  (ref)")

    ax.set_xlabel("Transform size N")
    ax.set_ylabel("Relative energy error  |E_t − E_f| / E_t")
    ax.set_title("FFT Parseval energy conservation error vs N")
    ax.legend(fontsize=9)
    ax.grid(True, which="both", ls=":", alpha=0.5)
    ax.set_xticks(N)
    ax.set_xticklabels([str(int(n)) for n in N], rotation=45, fontsize=8)
    fig.tight_layout()
    fig.savefig(out, dpi=130)
    plt.close(fig)
    return out


# ---------------------------------------------------------------------------
# Figure registry
# ---------------------------------------------------------------------------

FIGURES = [
    ("fft_spectrum.csv",  plot_spectrum,  "fft_spectrum.png"),
    ("fft_accuracy.csv",  plot_accuracy,  "fft_accuracy.png"),
    ("fft_parseval.csv",  plot_parseval,  "fft_parseval.png"),
]


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main(argv=None):
    parser = argparse.ArgumentParser(description="Render FFT charts from CSV data.")
    here = os.path.dirname(os.path.abspath(__file__))
    parser.add_argument("--data", default=os.path.join(here, "data"),
                        help="directory containing the exported CSVs")
    parser.add_argument("--out", default=here,
                        help="directory to write PNG charts into")
    args = parser.parse_args(argv)

    if not os.path.isdir(args.data):
        print(f"error: data directory not found: {args.data}", file=sys.stderr)
        print("Run the harness first: quotient_fft_tests --export-charts", file=sys.stderr)
        return 1

    os.makedirs(args.out, exist_ok=True)

    written = []
    missing = []
    for csv_name, fn, png_name in FIGURES:
        csv_path = os.path.join(args.data, csv_name)
        if not os.path.isfile(csv_path):
            missing.append(csv_name)
            continue
        out_path = os.path.join(args.out, png_name)
        written.append(fn(csv_path, out_path))

    for path in written:
        print(f"wrote {path}")
    if missing:
        print(f"warning: missing CSVs (skipped): {', '.join(missing)}", file=sys.stderr)

    if not written:
        print("error: no charts produced", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
