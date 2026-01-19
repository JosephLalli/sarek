#!/usr/bin/env python3

import argparse
import gzip
import os
import sys

def open_maybe_gz(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")

def parse_tags(field4):
    tags = {}
    for part in field4.strip().split(";"):
        if "=" in part:
            k, v = part.split("=", 1)
            tags[k] = v
        else:
            tags[part] = True
    return tags

def main():
    parser = argparse.ArgumentParser(description="Convert STRling output to VCF using a shared catalog.")
    parser.add_argument("--strling", required=True, help="STRling genotype.txt file")
    parser.add_argument("--catalog", required=True, help="Catalog BED file (TRGT/ExpansionHunter style)")
    parser.add_argument("--output", required=True, help="Output VCF file")
    parser.add_argument("--sample", required=True, help="Sample name")
    args = parser.parse_args()

    # 1. Parse Catalog
    catalog = {}
    with open_maybe_gz(args.catalog) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split("\t")
            chrom, start, end, field4 = parts[0], int(parts[1]), int(parts[2]), parts[3]
            tags = parse_tags(field4)
            locus_id = tags.get("ID", f"{chrom}_{start}_{end}")
            motifs = tags.get("MOTIFS", "").split(",")
            motif = motifs[0] if motifs else "N"
            ref_len = end - start
            motif_len = len(motif)
            ref_copies = ref_len / motif_len if motif_len > 0 else 0
            
            # Key by (chrom, start, end)
            catalog[(chrom, start, end)] = {
                "id": locus_id,
                "motif": motif,
                "ref_copies": ref_copies
            }

    # 2. Write VCF Header
    with open(args.output, "w") as out:
        out.write("##fileformat=VCFv4.2\n")
        out.write('##FILTER=<ID=PASS,Description="All filters passed">\n')
        out.write('##INFO=<ID=ID,Number=1,Type=String,Description="Locus ID from catalog">\n')
        out.write('##INFO=<ID=RU,Number=1,Type=String,Description="Repeat Unit">\n')
        out.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        out.write('##FORMAT=<ID=REPCN,Number=2,Type=Float,Description="Repeat Copy Number">\n')
        out.write('##FORMAT=<ID=AD,Number=2,Type=Integer,Description="Anchored and Spanning read counts">\n')
        out.write('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Total Depth">\n')
        out.write(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{args.sample}\n")

        # 3. Parse STRling and Write Records
        with open_maybe_gz(args.strling) as f:
            header = f.readline()
            if not header:
                return
            col_map = {name: i for i, name in enumerate(header.strip().split("\t"))}
            
            for line in f:
                parts = line.strip().split("\t")
                chrom = parts[col_map["#chrom"]]
                left = int(parts[col_map["left"]])
                right = int(parts[col_map["right"]])
                
                # Match against catalog
                locus = catalog.get((chrom, left, right))
                if not locus:
                    continue
                
                al1_est = float(parts[col_map["allele1_est"]])
                al2_est = float(parts[col_map["allele2_est"]])
                anchored = int(parts[col_map["anchored_reads"]])
                spanning = int(parts[col_map["spanning_reads"]])
                depth = int(parts[col_map["depth"]])
                
                abs_al1 = locus["ref_copies"] + al1_est
                abs_al2 = locus["ref_copies"] + al2_est
                
                # GT logic: if estimates are 0, it's the reference allele.
                # In TR-VCF, we often just use 0/1 if they are different from ref, 
                # but with REPCN we can be more explicit.
                # For now, let's use 0/1 style or ./. if missing.
                gt = "0/1"
                if al1_est == 0 and al2_est == 0:
                    gt = "0/0"
                elif al1_est != 0 and al2_est != 0:
                    gt = "1/2" # Simplified
                
                info = f"ID={locus['id']};RU={locus['motif']}"
                fmt = "GT:REPCN:AD:DP"
                sample_data = f"{gt}:{abs_al1:.2f},{abs_al2:.2f}:{anchored},{spanning}:{depth}"
                
                # POS in VCF is 1-based. Catalog/STRling 'left' is 0-based start.
                pos = left + 1
                out.write(f"{chrom}\t{pos}\t{locus['id']}\tN\t<STR>\t.\tPASS\t{info}\t{fmt}\t{sample_data}\n")

if __name__ == "__main__":
    main()
