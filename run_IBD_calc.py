import argparse
import os
import numpy as np
import pysam
import pandas as pd
import json
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


def is_ibs1(gt1, gt2):
  gt1_set = set(gt1)
  gt2_set = set(gt2)
  jacard = len(gt1_set.intersection(gt2_set)) / len(gt1_set.union(gt2_set))
  if jacard >= 0.5:
    return True
  return False

def is_ibs2(gt1, gt2):
  gt1_set = set(gt1)
  gt2_set = set(gt2)
  jaccard = len(gt1_set.intersection(gt2_set)) / len(gt1_set.union(gt2_set))
  if jaccard == 1:
    return True
  return False




def plot_ibd_chromosome(ibd1_segments, ibd2_segments, chrom_length=None, chrom_name="chr", ax=None):
    """
    Plot IBD0, IBD1, and IBD2 regions along a chromosome.

    Parameters:
        ibd1_segments: list of (start, end) tuples
        ibd2_segments: list of (start, end) tuples
        chrom_length:  total chromosome length in bp. If None, inferred from segments.
        chrom_name:    label for the y-axis
        ax:            optional matplotlib Axes to draw on
    """
    if chrom_length is None:
        all_ends = [e for _, e in ibd1_segments] + [e for _, e in ibd2_segments]
        chrom_length = max(all_ends) if all_ends else 1

    if ax is None:
        fig, ax = plt.subplots(figsize=(14, 2))

    bar_height = 0.6
    y = 0

    colors = {
        "IBD0": "#D3D3D3",
        "IBD1": "#4A90D9",
        "IBD2": "#D94A4A",
    }

    # draw full chromosome as IBD0 background
    ax.barh(y, chrom_length, height=bar_height, color=colors["IBD0"], edgecolor="black", linewidth=0.5)

    # overlay IBD1
    for start, end in ibd1_segments:
        ax.barh(y, end - start, left=start, height=bar_height, color=colors["IBD1"])

    # overlay IBD2 on top (takes priority)
    for start, end in ibd2_segments:
        ax.barh(y, end - start, left=start, height=bar_height, color=colors["IBD2"])

    # formatting
    ax.set_xlim(0, chrom_length)
    ax.set_ylim(-0.5, 0.5)
    ax.set_yticks([0])
    ax.set_yticklabels([chrom_name])
    ax.set_xlabel("Position (bp)")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)

    # scale x-axis to Mb if chromosome is large
    if chrom_length > 1_000_000:
        ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f"{x/1e6:.0f}"))
        ax.set_xlabel("Position (Mb)")

    legend_patches = [mpatches.Patch(color=c, label=l) for l, c in colors.items()]
    ax.legend(handles=legend_patches, loc="upper right", frameon=False, ncol=3)

    plt.tight_layout()
    return ax



def compute_ibd_prop_and_segments(
    vcf, samp1, samp2,
    min_af=0.1, A=1, B=10, min_markers=2000,
    r0=np.inf, r1=np.inf, r2=np.inf,
    ibd_function=is_ibs1
):
    L_shared = 0
    s1, s2, s3, e = 0, 0, 0, 0
    n_valid_markers = 0
    segments = []
    pos_start = None
    prev_pos = None
    prev_chrom = None

    def do_break():
        nonlocal s1, s2, s3, e, L_shared, pos_start, prev_pos
        length = s3
        if (s1 > r0) or (s2 > r0):
            length = s3 + s2 + s1
        if (s1 > r1) or (s2 > r2):
            s1 = s1 + s2
        else:
            s1 = s2
        if length > min_markers:
            L_shared += length
            segments.append((pos_start, prev_pos))
            s1, s2, s3 = 0, 0, 0
        s2 = s3
        e = 0
        s3 = 0

    for record in vcf.fetch():
        chrom = record.chrom
        pos = record.pos
        af = record.info.get("AF", [0])[0]

        if (af < min_af) or (af > (1 - min_af)):
            continue

        # chromosome change — flush and reset
        if chrom != prev_chrom:
            if prev_chrom is not None:
                do_break()
                s1, s2, s3, e = 0, 0, 0, 0
            prev_chrom = chrom
            pos_start = pos

        n_valid_markers += 1
        s3 += 1
        gt_samp1 = record.samples[samp1]["GT"]
        gt_samp2 = record.samples[samp2]["GT"]
        if ibd_function(gt_samp1, gt_samp2):
            e = max(e - 1, 0)
        else:
            e += A

        if (e + A) > B:
            do_break()
            pos_start = pos
        prev_pos = pos
    # final segment check
    if s3 > 0:
        do_break()

    return segments, L_shared, n_valid_markers


def reconcile_ibd_segments(ibd1_segments, ibd2_segments):
    ibd2 = sorted(ibd2_segments, key=lambda x: x[0])
    ibd1_out = []
    for start, end in sorted(ibd1_segments, key=lambda x: x[0]):
        # Subtract all IBD2 regions from this IBD1 segment
        remaining = [(start, end)]
        for s2_start, s2_end in ibd2:
            if s2_start > end:
                break
            new_remaining = []
            for r_start, r_end in remaining:
                if s2_end < r_start or s2_start > r_end:
                    # no overlap
                    new_remaining.append((r_start, r_end))
                else:
                    # left fragment
                    if r_start < s2_start:
                        new_remaining.append((r_start, s2_start - 1))
                    # right fragment
                    if r_end > s2_end:
                        new_remaining.append((s2_end + 1, r_end))
            remaining = new_remaining
        ibd1_out.extend(remaining)
    return ibd1_out, list(ibd2_segments)

def parse_args():
    parser = argparse.ArgumentParser(description="Compute IBD segments and proportions between two samples in a VCF.")
    parser.add_argument("--config", required=True, help="Config file with sample pairs and parameters.")
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    
    config_path = args.config
    config = json.load(open(config_path))
    
    vcf_path = config["vcf"]
    A = config.get("A", 1)
    B = config.get("B", 10)
    min_af = config.get("min_af", 0.1)
    min_markers = config.get("min_markers", 2000)
    r0 = config.get("r0", np.inf)
    r1 = config.get("r1", np.inf)
    r2 = config.get("r2", np.inf)
    samp1 = config["sample1"]
    samp2 = config["sample2"]
    
    
    vcf = pysam.VariantFile(vcf_path)
    
    ibd1_segments, L_shared, n_valid_markers = compute_ibd_prop_and_segments(
        vcf, samp1, samp2,
        min_af=min_af, A=A, B=B, min_markers=min_markers,
        r0=r0, r1=r1, r2=r2,
        ibd_function=is_ibs1
    )
    
    ibd2_segments, _, _ = compute_ibd_prop_and_segments(
        vcf, samp1, samp2,
        min_af=min_af, A=A, B=B, min_markers=min_markers,
        r0=r0, r1=r1, r2=r2,
        ibd_function=is_ibs2
    )
    
    ibd1_final, ibd2_final = reconcile_ibd_segments(ibd1_segments, ibd2_segments)
    
    out_p = config.get("output_prefix", "ibd_output")
    fig, ax = plt.subplots(figsize=(14, 2))
    plot_ibd_chromosome(ibd1_final, ibd2_final, chrom_name=f"{samp1} vs {samp2}",ax=ax)
    plt.savefig(f"{out_p}_ibd_plot.png", dpi=300)
    plt.close()
    
    
    





