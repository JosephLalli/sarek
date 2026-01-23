#!/usr/bin/env python3

import argparse
import gzip
import os
import sys

def open_maybe_gz(path):
    if not path: return None
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")

def parse_info(info_str):
    info = {}
    for part in info_str.split(";"):
        if "=" in part:
            k, v = part.split("=", 1)
            info[k] = v
        else:
            info[part] = True
    return info

def safe_int(value):
    try:
        return int(value)
    except (TypeError, ValueError):
        return None

def parse_repcn(value):
    if not value or value in (".", ".,."):
        return None
    parts = value.split(",")
    repcn = []
    for part in parts:
        try:
            repcn.append(float(part))
        except ValueError:
            return None
    return repcn if repcn else None

def parse_vcf(path, tool):
    records = []
    header_lines = []
    sample_name = "SAMPLE"
    if not path:
        return records, header_lines, sample_name

    with open_maybe_gz(path) as f:
        for line in f:
            if line.startswith("##"):
                header_lines.append(line.rstrip("\n"))
                continue
            if line.startswith("#CHROM"):
                parts = line.strip().split("\t")
                if len(parts) > 9:
                    sample_name = parts[9]
                continue

            parts = line.strip().split("\t")
            if len(parts) < 10:
                continue
            chrom = parts[0]
            pos = safe_int(parts[1])
            locus_id = parts[2]
            ref = parts[3]
            alt = parts[4]
            qual = parts[5]
            filt = parts[6]
            info_str = parts[7]
            fmt_str = parts[8]
            sample_str = parts[9]

            if pos is None:
                continue

            info = parse_info(info_str)
            info_id = info.get("ID", locus_id)
            motif = info.get("RU", "N")
            period = safe_int(info.get("PERIOD"))
            if period is None and motif not in ("", "N"):
                period = len(motif)
            end = safe_int(info.get("END"))
            if end is None:
                end = pos + len(ref) - 1

            fmt_keys = fmt_str.split(":")
            fmt_vals = sample_str.split(":")
            fmt_map = dict(zip(fmt_keys, fmt_vals))
            repcn = parse_repcn(fmt_map.get("REPCN"))

            records.append({
                "tool": tool,
                "chrom": chrom,
                "pos": pos,
                "end": end,
                "id": locus_id,
                "info_id": info_id,
                "ref": ref,
                "alt": alt,
                "qual": qual,
                "filter": filt,
                "motif": motif,
                "period": period,
                "gt": fmt_map.get("GT", "./."),
                "repcn": repcn,
                "info": info
            })
    return records, header_lines, sample_name

def is_pass_filter(value):
    return value in ("PASS", ".", "", None)

def repcn_within_delta(repcn_a, repcn_b, max_delta):
    if repcn_a is None or repcn_b is None:
        return True
    if len(repcn_a) != len(repcn_b):
        return False
    a_sorted = sorted(repcn_a)
    b_sorted = sorted(repcn_b)
    diffs = [abs(a - b) for a, b in zip(a_sorted, b_sorted)]
    return max(diffs) <= max_delta

def overlaps_with_slop(start_a, end_a, start_b, end_b, max_bp_shift):
    return (start_a <= (end_b + max_bp_shift)) and (end_a >= (start_b - max_bp_shift))

def main():
    parser = argparse.ArgumentParser(description="Consensus merging for STR VCFs (EH, GS, SL).")
    parser.add_argument("--eh", help="ExpansionHunter VCF")
    parser.add_argument("--gs", help="GangSTR VCF")
    parser.add_argument("--sl", help="STRling standardized VCF")
    parser.add_argument("--output", required=True, help="Output consensus VCF")
    parser.add_argument("--motif-threshold", type=int, default=30, help="Motif length threshold for EH preference. Default: 30")
    parser.add_argument("--max-bp-shift", type=int, default=1, help="Max bp shift to allow when matching loci. Default: 1")
    parser.add_argument("--max-copy-delta", type=float, default=4.0, help="Max REPCN delta to consider calls equivalent. Default: 4")
    args = parser.parse_args()

    eh_recs, eh_header, sample = parse_vcf(args.eh, "EH")
    gs_recs, gs_header, _ = parse_vcf(args.gs, "GS")
    sl_recs, sl_header, _ = parse_vcf(args.sl, "SL")

    all_records = eh_recs + gs_recs + sl_recs
    clusters = []

    for rec in all_records:
        assigned = False
        for cluster in clusters:
            if rec["chrom"] != cluster["chrom"]:
                continue
            if not overlaps_with_slop(rec["pos"], rec["end"], cluster["start"], cluster["end"], args.max_bp_shift):
                continue
            if cluster["motif"] not in (None, "N") and rec["motif"] not in ("N", None) and rec["motif"] != cluster["motif"]:
                continue
            if cluster["period"] is not None and rec["period"] is not None and rec["period"] != cluster["period"]:
                continue
            if not repcn_within_delta(cluster["repcn"], rec["repcn"], args.max_copy_delta):
                continue

            cluster["records"][rec["tool"]] = rec
            cluster["start"] = min(cluster["start"], rec["pos"])
            cluster["end"] = max(cluster["end"], rec["end"])
            if cluster["motif"] in (None, "N") and rec["motif"] not in ("N", None):
                cluster["motif"] = rec["motif"]
            if cluster["period"] is None and rec["period"] is not None:
                cluster["period"] = rec["period"]
            if cluster["repcn"] is None and rec["repcn"] is not None:
                cluster["repcn"] = rec["repcn"]
            assigned = True
            break

        if not assigned:
            clusters.append({
                "chrom": rec["chrom"],
                "start": rec["pos"],
                "end": rec["end"],
                "motif": rec["motif"] if rec["motif"] not in ("N", None) else None,
                "period": rec["period"],
                "repcn": rec["repcn"],
                "records": {rec["tool"]: rec}
            })

    clusters.sort(key=lambda item: (item["chrom"], item["start"]))

    with open(args.output, "w") as out:
        out.write("##fileformat=VCFv4.2\n")
        out.write('##INFO=<ID=MERGE_SRC,Number=1,Type=String,Description="Tool that provided the consensus GT">\n')
        out.write('##INFO=<ID=MERGE_RULE,Number=1,Type=String,Description="Rule used for selection">\n')
        out.write('##INFO=<ID=RU,Number=1,Type=String,Description="Repeat Unit">\n')
        out.write('##INFO=<ID=PERIOD,Number=1,Type=Integer,Description="Repeat unit length">\n')
        out.write('##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the locus">\n')
        out.write('##INFO=<ID=ID,Number=1,Type=String,Description="Locus ID from catalog">\n')
        out.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Consensus Genotype">\n')
        out.write('##FORMAT=<ID=REPCN,Number=2,Type=Float,Description="Consensus Repeat Copy Number">\n')
        out.write('##FORMAT=<ID=EH_GT,Number=1,Type=String,Description="ExpansionHunter Genotype">\n')
        out.write('##FORMAT=<ID=GS_GT,Number=1,Type=String,Description="GangSTR Genotype">\n')
        out.write('##FORMAT=<ID=SL_GT,Number=1,Type=String,Description="STRling Genotype">\n')
        out.write('##FORMAT=<ID=EH_REPCN,Number=2,Type=Float,Description="ExpansionHunter REPCN">\n')
        out.write('##FORMAT=<ID=GS_REPCN,Number=2,Type=Float,Description="GangSTR REPCN">\n')
        out.write('##FORMAT=<ID=SL_REPCN,Number=2,Type=Float,Description="STRling REPCN">\n')
        out.write(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}\n")

        for cluster in clusters:
            eh = cluster["records"].get("EH")
            gs = cluster["records"].get("GS")
            sl = cluster["records"].get("SL")

            winner_src = "None"
            rule = "Default"
            winner = None

            motif = "N"
            if eh and eh["motif"] not in ("", "N"):
                motif = eh["motif"]
            elif gs and gs["motif"] not in ("", "N"):
                motif = gs["motif"]
            elif sl and sl["motif"] not in ("", "N"):
                motif = sl["motif"]
            elif cluster["motif"]:
                motif = cluster["motif"]
            motif_len = len(motif) if motif else 0

            if motif_len >= args.motif_threshold:
                if eh and is_pass_filter(eh["filter"]):
                    winner, winner_src, rule = eh, "EH", f"MotifLen>={args.motif_threshold}"
                elif gs and is_pass_filter(gs["filter"]):
                    winner, winner_src, rule = gs, "GS", "Fallback_GS"
                elif sl and is_pass_filter(sl["filter"]):
                    winner, winner_src, rule = sl, "SL", "Fallback_SL"
            else:
                if gs and is_pass_filter(gs["filter"]):
                    winner, winner_src, rule = gs, "GS", f"MotifLen<{args.motif_threshold}"
                elif eh and is_pass_filter(eh["filter"]):
                    winner, winner_src, rule = eh, "EH", "Fallback_EH"
                elif sl and is_pass_filter(sl["filter"]):
                    winner, winner_src, rule = sl, "SL", "Fallback_SL"

            if not winner:
                if eh:
                    winner, winner_src, rule = eh, "EH", "NoPass_Priority"
                elif gs:
                    winner, winner_src, rule = gs, "GS", "NoPass_Priority"
                elif sl:
                    winner, winner_src, rule = sl, "SL", "NoPass_Priority"

            if not winner:
                continue

            chrom = winner["chrom"]
            pos = winner["pos"]
            ref = winner["ref"]
            alt = winner["alt"]
            qual = "."
            filt = winner["filter"] if winner["filter"] else "PASS"
            rec_id = winner["info_id"]
            period = winner["period"] if winner["period"] is not None else (len(motif) if motif else 0)
            end = winner["end"]

            info = f"ID={rec_id};RU={motif};PERIOD={period};END={end};MERGE_SRC={winner_src};MERGE_RULE={rule}"

            fmt_keys = ["GT", "REPCN", "EH_GT", "GS_GT", "SL_GT", "EH_REPCN", "GS_REPCN", "SL_REPCN"]
            fmt_vals = [
                winner["gt"],
                ",".join([f"{val:.2f}" for val in winner["repcn"]]) if winner["repcn"] else ".,.",
                eh["gt"] if eh else "./.",
                gs["gt"] if gs else "./.",
                sl["gt"] if sl else "./.",
                ",".join([f"{val:.2f}" for val in eh["repcn"]]) if eh and eh["repcn"] else ".,.",
                ",".join([f"{val:.2f}" for val in gs["repcn"]]) if gs and gs["repcn"] else ".,.",
                ",".join([f"{val:.2f}" for val in sl["repcn"]]) if sl and sl["repcn"] else ".,."
            ]

            out.write(f"{chrom}\t{pos}\t{rec_id}\t{ref}\t{alt}\t{qual}\t{filt}\t{info}\t{':'.join(fmt_keys)}\t{':'.join(fmt_vals)}\n")

if __name__ == "__main__":
    main()
