import argparse
import os
import numpy as np
import pysam
import pandas as pd
import json
import time

# def is_ibs1(gt1, gt2):
#   gt1_set = set(gt1)
#   gt2_set = set(gt2)
#   jacard = len(gt1_set.intersection(gt2_set)) / len(gt1_set.union(gt2_set))
#   if jacard >= 0.5:
#     return True
#   return False

# def is_ibs2(gt1, gt2):
#   gt1_set = set(gt1)
#   gt2_set = set(gt2)
#   jaccard = len(gt1_set.intersection(gt2_set)) / len(gt1_set.union(gt2_set))
#   if jaccard == 1:
#     return True
#   return False


def is_ibs1(gt1, gt2):
    return abs(sum(gt1) - sum(gt2)) < 2

def is_ibs2(gt1, gt2):
    return sum(gt1) == sum(gt2)


def compute_ibd_prop_and_segments(
    vcf, samp1, samp2,
    min_af=0.1, A=1, B=3, min_sig_length=2000000, min_markers=1000,
    r0=np.inf,
    ibd_function=is_ibs1
):
    e = 0
    s1_pos = 0
    s2_pos = 0
    s3_pos = 0
    s0_pos = 0
    

    segments = []
    prev_pos = None
    prev_chrom = None
    
    def do_break():
        nonlocal s0_pos, s1_pos, s2_pos, s3_pos, e, r0, min_sig_length, segments, chrom
        length = s3_pos - s2_pos
        s2_length = s2_pos - s1_pos
        s1_length = s1_pos - s0_pos
    
        if s1_length > r0:
            length += (s1_length + s2_length)
        elif s2_length > r0:
            length += s2_length
            
        if length > min_sig_length:
            pos_start = s3_pos - length
            segments.append((chrom, pos_start, s3_pos))
            # segment added, reset segs
            s0_pos = s3_pos
            s1_pos = s3_pos
            s2_pos = s3_pos
            # s3_pos = 0
        else:
            s0_pos = s1_pos
            s1_pos = s2_pos
            s2_pos = s3_pos

        e = 0

    for record in vcf.fetch():

        af = record.info.get("AF", [0])[0]
        if (af < min_af) or (af > (1 - min_af)):
            continue

        new_chrom = record.chrom
        # chromosome change — flush and reset
        if new_chrom != prev_chrom:
            if prev_chrom is not None:
                do_break()
            prev_chrom = new_chrom
            chrom = new_chrom
            print(f"Processing chromosome {chrom}")
            s0_pos = 0
            s1_pos = 0
            s2_pos = 0
            s3_pos = 0

        s3_pos = record.pos

        gt_samp1 = record.samples[samp1]["GT"]
        gt_samp2 = record.samples[samp2]["GT"]
        if ibd_function(gt_samp1, gt_samp2):
            e = max(e - 1, 0)
        else:
            e += A
        if (e + A) > B:
            do_break()

    do_break()
    return segments

def reconcile_ibd_segments(ibd1_segments, ibd2_segments, chr_sizes_d):
    ibd2 = sorted(ibd2_segments, key=lambda x: (x[0], x[1]))
    ibd1_out = []
    ibd1_out_d = {}
    for chrom, start, end in sorted(ibd1_segments, key=lambda x: (x[0], x[1])):
        remaining = [(start, end)]
        for s2_chrom, s2_start, s2_end in ibd2:
            if s2_chrom != chrom:
                continue
            if s2_start > end:
                break
            new_remaining = []
            for r_start, r_end in remaining:
                if s2_end < r_start or s2_start > r_end:
                    new_remaining.append((r_start, r_end))
                else:
                    if r_start < s2_start:
                        new_remaining.append((r_start, s2_start - 1))
                    if r_end > s2_end:
                        new_remaining.append((s2_end + 1, r_end))
            remaining = new_remaining
        for r_start, r_end in remaining:
            ibd1_out.append((chrom, r_start, r_end))
            current_eles = ibd1_out_d.get(chrom, [])
            current_eles.append((r_start, r_end, "ibd1"))
            ibd1_out_d[chrom] = current_eles
    ibd2_out_d = {}
    for chrom, start, end in ibd2_segments:
        lis = ibd2_out_d.get(chrom, [])
        lis.append((start, end, "ibd2"))
        ibd2_out_d[chrom] = lis

    chroms = sorted(np.unique([ele for ele in list(ibd1_out_d.keys()) + list(ibd2_out_d.keys())]))
    data_d = {"chrom": [], "start": [], "end": [], "ibd_status": []}
    for chrom in chroms:
        ibd1_eles = ibd1_out_d.get(chrom, [])
        ibd2_eles = ibd2_out_d.get(chrom, [])
        comb_lis = sorted(ibd1_eles + ibd2_eles)

        chrom_len = chr_sizes_d.get(chrom, 0)
        cursor = 0

        for start, end, ibd_state in comb_lis:
            if start > cursor:
                data_d["chrom"].append(chrom)
                data_d["start"].append(cursor)
                data_d["end"].append(start)
                data_d["ibd_status"].append("ibd0")
            data_d["chrom"].append(chrom)
            data_d["start"].append(start)
            data_d["end"].append(end)
            data_d["ibd_status"].append(ibd_state)
            cursor = end

        if cursor < chrom_len:
            data_d["chrom"].append(chrom)
            data_d["start"].append(cursor)
            data_d["end"].append(chrom_len)
            data_d["ibd_status"].append("ibd0")

    df = pd.DataFrame(data_d)
    return df

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
    B = config.get("B", 3)
    min_af = config.get("min_af", 0.1)
    samp1 = config["sample1"]
    samp2 = config["sample2"]
    out_p = config.get("out_dir", f"./{samp1}_vs_{samp2}")
    
    if not os.path.exists(out_p):
        os.makedirs(out_p)
    
    min_sig_length_ibd1 = config.get("min_sig_length_ibd1", 5_000_000)
    r0_ibd1 = min_sig_length_ibd1 // 2
    
    min_sig_length_ibd2 = config.get("min_sig_length_ibd2", 2_000_000)
    r0_ibd2 = min_sig_length_ibd2 // 2
    
    start_time = time.time()
    vcf = pysam.VariantFile(vcf_path)
    
    ibd1_segments = compute_ibd_prop_and_segments(
        vcf, samp1, samp2,
        min_af=min_af, A=A, B=B, min_sig_length=min_sig_length_ibd1,
        r0=r0_ibd1,
        ibd_function=is_ibs1
    )
    
    # reopen vcf to reset iterator
    vcf = pysam.VariantFile(vcf_path) 
    ibd2_segments = compute_ibd_prop_and_segments(
        vcf, samp1, samp2,
        min_af=min_af, A=A, B=B, min_sig_length=min_sig_length_ibd2,
        r0=r0_ibd2,
        ibd_function=is_ibs2
    )
    
    # hard coded from the vcf being used
    chrom_size_d = {
        "1": 249250621,
        "2": 243199373,
        "3": 198022430,
        "4": 191154276,
        "5": 180915260,
        "6": 171115067,
        "7": 159138663,
        "8": 146364022,
        "9": 141213431,
        "10": 135534747,
        "11": 135006516,
        "12": 133851895,
        "13": 115169878,
        "14": 107349540,
        "15": 102531392,
        "16": 90354753,
        "17": 81195210,
        "18": 78077248,
        "19": 59128983,
        "20": 63025520,
        "21": 48129895,
        "22": 51304566,
    }
    
    out_d = {"ibd1_segments": ibd1_segments, "ibd2_segments": ibd2_segments}
    out_json_p = os.path.join(out_p, f"ibd_segments.json")
    with open(out_json_p, "w") as f:
        json.dump(out_d, f, indent=4)
    
    segment_df = reconcile_ibd_segments(ibd1_segments, ibd2_segments, chrom_size_d)
    out_p = os.path.join(out_p, f"ibd_segments.tsv")
    segment_df.to_csv(out_p, sep="\t", index=False)
    
    end_time = time.time()
    print(f"IBD calculation completed in {end_time - start_time:.2f} seconds")
    
    





