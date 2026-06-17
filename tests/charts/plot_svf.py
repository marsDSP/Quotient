#!/usr/bin/env python3
"""Render SVF verification charts from the harness-exported CSV data.

The CSVs are produced by the filter test harness:

    quotient_filter_tests --export-charts --chart-dir tests/charts/data

and this script turns them into PNG figures next to (by default) this file.

    python3 tests/charts/plot_svf.py --data tests/charts/data --out tests/charts

It reproduces, from the *production* MarsDSP::Filters::TwoPoleSVF, the kinds of
figures shown in docs/papers/SvfLinearTrapOptimised2.pdf (magnitude sweeps,
all-pass phase, time-domain saw responses) plus explicit code-vs-paper overlays
that demonstrate the verification.
"""

from __future__ import annotations

import argparse
import csv
import os
import sys

import matplotlib

matplotlib.use("Agg")  # headless / CI friendly
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402


def read_csv(path):
    """Return (header: list[str], columns: dict[str, np.ndarray])."""
    with open(path, newline="") as handle:
        reader = csv.reader(handle)
        header = next(reader)
        rows = [[float(x) for x in row] for row in reader if row]
    data = np.array(rows, dtype=float) if rows else np.empty((0, len(header)))
    columns = {name: data[:, i] for i, name in enumerate(header)}
    return header, columns


def _style_freq_axis(ax, freqs):
    ax.set_xscale("log")
    ax.set_xlim(max(freqs.min(), 10.0), freqs.max())
    ax.set_xlabel("Frequency (Hz)")
    ax.grid(True, which="both", ls=":", alpha=0.5)


def plot_magnitude(path, out, title, ylabel="Magnitude (dB)"):
    header, cols = read_csv(path)
    freqs = cols[header[0]]
    fig, ax = plt.subplots(figsize=(9, 5))
    for name in header[1:]:
        ax.plot(freqs, cols[name], label=name, lw=1.6)
    _style_freq_axis(ax, freqs)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.legend(loc="best", fontsize=9)
    fig.tight_layout()
    fig.savefig(out, dpi=130)
    plt.close(fig)
    return out


def plot_verification(path, out, title):
    """Overlay code vs paper; show the difference on a twin axis."""
    header, cols = read_csv(path)
    freqs = cols[header[0]]
    code = cols[header[1]]
    paper = cols[header[2]]

    fig, (ax, axd) = plt.subplots(
        2, 1, figsize=(9, 6), sharex=True, gridspec_kw={"height_ratios": [3, 1]}
    )
    ax.plot(freqs, paper, label=header[2], lw=3.0, alpha=0.4, color="tab:orange")
    ax.plot(freqs, code, label=header[1], lw=1.4, color="tab:blue")
    ax.set_ylabel("Magnitude (dB)")
    ax.set_title(title)
    ax.legend(loc="best", fontsize=9)
    ax.grid(True, which="both", ls=":", alpha=0.5)

    axd.plot(freqs, code - paper, lw=1.2, color="tab:red")
    axd.axhline(0.0, color="k", lw=0.6)
    axd.set_ylabel("code - paper\n(dB)")
    _style_freq_axis(axd, freqs)
    fig.tight_layout()
    fig.savefig(out, dpi=130)
    plt.close(fig)
    return out


def plot_phase(path, out, title):
    header, cols = read_csv(path)
    freqs = cols[header[0]]
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.plot(freqs, cols[header[2]], label=header[2], lw=3.0, alpha=0.4, color="tab:orange")
    ax.plot(freqs, cols[header[1]], label=header[1], lw=1.4, color="tab:blue")
    _style_freq_axis(ax, freqs)
    ax.set_ylabel("Phase (degrees)")
    ax.set_title(title)
    ax.legend(loc="best", fontsize=9)
    fig.tight_layout()
    fig.savefig(out, dpi=130)
    plt.close(fig)
    return out


def plot_time(path, out):
    header, cols = read_csv(path)
    t_ms = cols["t"] * 1000.0
    series = [name for name in header if name not in ("t", "input")]

    n = len(series)
    fig, axes = plt.subplots(n + 1, 1, figsize=(9, 1.6 * (n + 1)), sharex=True)
    axes[0].plot(t_ms, cols["input"], color="k", lw=1.2)
    axes[0].set_ylabel("input", rotation=0, ha="right", va="center")
    axes[0].grid(True, ls=":", alpha=0.5)

    for ax, name in zip(axes[1:], series):
        ax.plot(t_ms, cols[name], lw=1.2)
        ax.set_ylabel(name, rotation=0, ha="right", va="center")
        ax.grid(True, ls=":", alpha=0.5)

    axes[-1].set_xlabel("Time (ms)")
    fig.suptitle("Time-domain response to a 500 Hz saw (fc = 1 kHz)")
    fig.tight_layout(rect=(0, 0, 1, 0.98))
    fig.savefig(out, dpi=130)
    plt.close(fig)
    return out


# (csv filename, plot function, output png, *extra args)
FIGURES = [
    ("mag_core.csv", plot_magnitude, "svf_magnitude_core.png",
     "SVF core responses (fc = 1 kHz, Q = 0.707)"),
    ("mag_lowpass_resonance.csv", plot_magnitude, "svf_lowpass_resonance.png",
     "Low-pass resonance vs Q (fc = 1 kHz)"),
    ("mag_bell.csv", plot_magnitude, "svf_bell.png",
     "Bell / peaking EQ (fc = 1 kHz, Q = 2)"),
    ("mag_shelf.csv", plot_magnitude, "svf_shelves.png",
     "Low/high shelving filters (fc = 1 kHz)"),
    ("verify_lowpass.csv", plot_verification, "svf_verify_lowpass.png",
     "Verification: low-pass magnitude() vs paper transfer function (Q = 2)"),
    ("verify_bandpass.csv", plot_verification, "svf_verify_bandpass.png",
     "Band-pass: 0 dB-normalised code vs paper peak-gain-Q band (Q = 4)"),
    ("phase_allpass.csv", plot_phase, "svf_allpass_phase.png",
     "All-pass phase: measured impulse response vs paper (fc = 1 kHz)"),
    ("time_saw.csv", plot_time, "svf_time_domain.png"),
]


def main(argv=None):
    parser = argparse.ArgumentParser(description="Render SVF charts from CSV data.")
    here = os.path.dirname(os.path.abspath(__file__))
    parser.add_argument("--data", default=os.path.join(here, "data"),
                        help="directory containing the exported CSVs")
    parser.add_argument("--out", default=here,
                        help="directory to write PNG charts into")
    args = parser.parse_args(argv)

    if not os.path.isdir(args.data):
        print(f"error: data directory not found: {args.data}", file=sys.stderr)
        print("Run the harness first: quotient_filter_tests --export-charts", file=sys.stderr)
        return 1

    os.makedirs(args.out, exist_ok=True)

    written = []
    missing = []
    for csv_name, fn, png_name, *extra in FIGURES:
        csv_path = os.path.join(args.data, csv_name)
        if not os.path.isfile(csv_path):
            missing.append(csv_name)
            continue
        out_path = os.path.join(args.out, png_name)
        written.append(fn(csv_path, out_path, *extra))

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
