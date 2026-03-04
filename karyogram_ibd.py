"""
karyogram_ibd.py
----------------
Reads a TRUFFLE .segments file and produces a karyogram-style plot
showing IBD1 and IBD2 segments across chromosome(s).

Usage:
    python karyogram_ibd.py truffle.segments
    python karyogram_ibd.py truffle.segments --out karyogram.png
"""

import sys
import argparse
import re
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch
from matplotlib.collections import PatchCollection
import matplotlib.patheffects as pe

# ── Chromosome lengths (GRCh38, Mbp) ──────────────────────────────────────────
CHR_LENGTHS_MBP = {
    1: 248.956, 2: 242.193, 3: 198.295, 4: 190.214, 5: 181.538,
    6: 170.806, 7: 159.345, 8: 145.138, 9: 138.394, 10: 133.797,
    11: 135.086, 12: 133.275, 13: 114.364, 14: 107.043, 15: 101.991,
    16:  90.338, 17:  83.257, 18:  80.373, 19:  58.617, 20:  64.444,
    21:  46.709, 22:  50.818, 23: 156.040, 24:  57.227,
}

# Approximate centromere positions (Mbp) — GRCh38
CENTROMERES_MBP = {
    1: 125.0, 2: 93.3,  3: 91.2,  4: 50.4,  5: 48.4,
    6: 61.0,  7: 59.9,  8: 45.6,  9: 49.0,  10: 40.2,
    11: 53.7, 12: 35.8, 13: 17.9, 14: 17.6, 15: 19.0,
    16: 36.8, 17: 25.1, 18: 18.5, 19: 26.2, 20: 28.1,
    21: 12.0, 22: 14.7, 23: 61.0, 24: 12.5,
}

# ── Colour palette ─────────────────────────────────────────────────────────────
BG_COLOR      = "#0f1117"
CHR_BASE      = "#2a2d3e"
CHR_EDGE      = "#4a4e69"
IBD1_COLOR    = "#4cc9f0"   # cyan-blue
IBD2_COLOR    = "#f72585"   # hot pink
CENTRO_COLOR  = "#ffd60a"   # amber
TEXT_COLOR    = "#e2e8f0"
GRID_COLOR    = "#1e2130"

LABEL_NAMES = {23: "X", 24: "Y"}


# ── Parsing ────────────────────────────────────────────────────────────────────

def parse_segments(filepath):
    """Return list of dicts with keys: type, id1, id2, chrom, pos_mbp, length_mbp."""
    segments = []
    with open(filepath) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("TYPE"):
                continue
            parts = line.split()
            # cols: TYPE ID1 ID2 CHROM VARSTART VAREND POS(Mbp) Mbp LENGTH(Mbp) Mbp NMARKERS
            ibd_type = parts[0]          # IBD1 or IBD2
            id1      = parts[1]
            id2      = parts[2]
            chrom    = int(parts[3])
            pos_mbp  = float(parts[6])   # start position in Mbp
            length   = float(parts[8])   # segment length in Mbp
            segments.append(dict(
                type=ibd_type, id1=id1, id2=id2,
                chrom=chrom, pos=pos_mbp, length=length
            ))
    return segments


# ── Drawing helpers ────────────────────────────────────────────────────────────

def rounded_bar(ax, x, y, width, height, radius=None, **kwargs):
    """Draw a rounded rectangle centred at (x, y) with given width/height."""
    if radius is None:
        radius = min(width, height) * 0.35
    rect = FancyBboxPatch(
        (x - width / 2, y - height / 2),
        width, height,
        boxstyle=f"round,pad=0,rounding_size={radius}",
        **kwargs
    )
    ax.add_patch(rect)
    return rect


def draw_chromosome(ax, chrom_num, x_center, chr_length, bar_width=0.55, bar_height_scale=1.0):
    """Draw the grey chromosome bar and centromere."""
    h = chr_length * bar_height_scale
    rounded_bar(
        ax, x_center, chr_length / 2,
        width=bar_width, height=h,
        facecolor=CHR_BASE, edgecolor=CHR_EDGE, linewidth=0.8, zorder=2
    )

    # centromere pinch
    centro = CENTROMERES_MBP.get(chrom_num, chr_length / 2)
    pinch_h = max(chr_length * 0.025, 1.5)
    rounded_bar(
        ax, x_center, centro,
        width=bar_width * 1.25, height=pinch_h,
        facecolor=CENTRO_COLOR, edgecolor=CENTRO_COLOR,
        alpha=0.85, zorder=4
    )


def draw_ibd_segment(ax, x_center, pos, length, ibd_type, bar_width=0.55):
    """Overlay a coloured IBD segment on the chromosome bar."""
    color = IBD2_COLOR if ibd_type == "IBD2" else IBD1_COLOR
    alpha = 0.88
    seg_w = bar_width * 0.95

    rounded_bar(
        ax, x_center, pos + length / 2,
        width=seg_w, height=length,
        facecolor=color, edgecolor="none",
        alpha=alpha, zorder=5
    )


# ── Main plot ──────────────────────────────────────────────────────────────────

def build_karyogram(segments, out_path="karyogram_ibd.png"):

    # Determine which chromosomes are present; fall back to all autosomes if needed
    chroms_in_data = sorted(set(s["chrom"] for s in segments))
    # For display, show all 22 autosomes + X, highlighting those with data
    all_chroms = list(range(1, 23)) #+ [23]
    display_chroms = all_chroms  # always show full karyotype layout

    n = len(display_chroms)
    x_positions = {c: i + 1 for i, c in enumerate(display_chroms)}

    max_len = max(CHR_LENGTHS_MBP[c] for c in display_chroms)

    # Figure dimensions
    fig_w = max(16, n * 0.75)
    fig_h = 10
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    fig.patch.set_facecolor(BG_COLOR)
    ax.set_facecolor(BG_COLOR)

    # ── Background grid lines ──
    for y in np.arange(0, max_len + 50, 50):
        ax.axhline(y, color=GRID_COLOR, linewidth=0.4, zorder=1)

    # ── Draw all chromosomes ──
    bar_width = 0.55
    for chrom in display_chroms:
        x = x_positions[chrom]
        length = CHR_LENGTHS_MBP[chrom]
        draw_chromosome(ax, chrom, x, length, bar_width=bar_width)

    # ── Overlay IBD segments — IBD1 first, IBD2 on top so overlaps are visible ──
    for ibd_type in ("IBD1", "IBD2"):
        for seg in segments:
            if seg["type"] != ibd_type:
                continue
            c = seg["chrom"]
            if c not in x_positions:
                continue
            x = x_positions[c]
            draw_ibd_segment(ax, x, seg["pos"], seg["length"], seg["type"], bar_width=bar_width)

    # ── Axes cosmetics ──
    ax.set_xlim(0.2, n + 0.8)
    ax.set_ylim(-8, max_len + 15)
    ax.invert_yaxis()

    ax.set_xticks([x_positions[c] for c in display_chroms])
    ax.set_xticklabels(
        [LABEL_NAMES.get(c, str(c)) for c in display_chroms],
        fontsize=7.5, color=TEXT_COLOR, fontfamily="monospace"
    )
    ax.tick_params(axis="x", length=0, pad=4)

    ax.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:.0f}"))
    ax.set_ylabel("Position (Mbp)", color=TEXT_COLOR, fontsize=9, labelpad=8)
    ax.tick_params(axis="y", colors=TEXT_COLOR, labelsize=8)

    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.tick_params(axis="x", bottom=False)

    # ── Title ──
    pair_label = ""
    if segments:
        pair_label = f"  —  {segments[0]['id1']} × {segments[0]['id2']}"
    ax.set_title(
        f"IBD Karyogram{pair_label}",
        color=TEXT_COLOR, fontsize=13, fontweight="bold",
        fontfamily="monospace", pad=14
    )

    # ── Legend ──
    legend_patches = [
        mpatches.Patch(facecolor=IBD1_COLOR, label="IBD1", alpha=0.9),
        mpatches.Patch(facecolor=IBD2_COLOR, label="IBD2", alpha=0.9),
        mpatches.Patch(facecolor=CENTRO_COLOR, label="Centromere", alpha=0.85),
    ]
    leg = ax.legend(
        handles=legend_patches,
        loc="lower right", frameon=True,
        framealpha=0.25, facecolor="#1a1d2e",
        edgecolor=CHR_EDGE, labelcolor=TEXT_COLOR,
        fontsize=9, handlelength=1.6, handleheight=1.1
    )

    # ── Segment count annotation ──
    ibd1_n = sum(1 for s in segments if s["type"] == "IBD1")
    ibd2_n = sum(1 for s in segments if s["type"] == "IBD2")
    ibd1_tot = sum(s["length"] for s in segments if s["type"] == "IBD1")
    ibd2_tot = sum(s["length"] for s in segments if s["type"] == "IBD2")
    info = (
        f"IBD1: {ibd1_n} segments  ({ibd1_tot:.1f} Mbp total)\n"
        f"IBD2: {ibd2_n} segments  ({ibd2_tot:.1f} Mbp total)"
    )
    ax.text(
        0.01, 0.01, info,
        transform=ax.transAxes,
        fontsize=8, color=TEXT_COLOR, fontfamily="monospace",
        va="bottom", ha="left",
        bbox=dict(facecolor="#1a1d2e", alpha=0.5, edgecolor=CHR_EDGE, boxstyle="round,pad=0.4")
    )

    plt.tight_layout(pad=1.5)
    fig.savefig(out_path, dpi=180, bbox_inches="tight", facecolor=BG_COLOR)
    print(f"Saved → {out_path}")
    plt.close(fig)


# ── Entry point ───────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description="Karyogram-style IBD plot from TRUFFLE .segments file")
    parser.add_argument("segments_file", help="Path to the .segments file")
    parser.add_argument("--out", default="karyogram_ibd.png", help="Output image path (default: karyogram_ibd.png)")
    args = parser.parse_args()

    segments = parse_segments(args.segments_file)
    if not segments:
        print("No segments found — check the input file format.", file=sys.stderr)
        sys.exit(1)

    print(f"Loaded {len(segments)} segments across chromosomes: "
          f"{sorted(set(s['chrom'] for s in segments))}")

    build_karyogram(segments, out_path=args.out)


if __name__ == "__main__":
    main()