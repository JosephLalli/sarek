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

def safe_int(value, default=0):
    try:
        return int(value)
    except (TypeError, ValueError):
        return default

def safe_float(value, default=0.0):
    try:
        return float(value)
    except (TypeError, ValueError):
        return default

def parse_fai(path):
    lengths = {}
    order = []
    with open(path, "r") as f:
        for line in f:
            if not line.strip():
                continue
            parts = line.strip().split("\t")
            if len(parts) < 2:
                continue
            name = parts[0]
            try:
                length = int(parts[1])
            except ValueError:
                continue
            lengths[name] = length
            order.append(name)
    return lengths, order

def parse_fasta_lengths(path):
    lengths = {}
    order = []
    current = None
    with open_maybe_gz(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                current = line[1:].split()[0]
                if current not in lengths:
                    lengths[current] = 0
                    order.append(current)
            elif current:
                lengths[current] += len(line)
    return lengths, order

def collect_strling_contigs(path):
    contigs = set()
    with open_maybe_gz(path) as f:
        header = f.readline()
        if not header:
            return contigs
        col_map = {name: i for i, name in enumerate(header.strip().split("\t"))}
        chrom_idx = col_map.get("#chrom")
        if chrom_idx is None:
            return contigs
        for line in f:
            if not line.strip():
                continue
            parts = line.strip().split("\t")
            if len(parts) <= chrom_idx:
                continue
            contigs.add(parts[chrom_idx])
    return contigs

def build_repeat_sequence(motif, length):
    if length <= 0:
        return "N"
    if not motif or motif == "N":
        return "N" * length
    repeats = (motif * ((length // len(motif)) + 1))[:length]
    return repeats

def main():
    parser = argparse.ArgumentParser(description="Convert STRling output to VCF using a shared catalog.")
    parser.add_argument("--strling", required=True, help="STRling genotype.txt file")
    parser.add_argument("--catalog", required=True, help="Catalog BED file (TRGT/ExpansionHunter style)")
    parser.add_argument("--ref-fai", help="Reference FASTA index (.fai) for contig lengths")
    parser.add_argument("--ref", help="Reference FASTA (used to infer contig lengths if .fai not provided)")
    parser.add_argument("--output", required=True, help="Output VCF file")
    parser.add_argument("--sample", required=True, help="Sample name")
    args = parser.parse_args()

    # 1. Parse Catalog
    catalog = {}
    catalog_contigs = set()
    with open_maybe_gz(args.catalog) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split("\t")
            chrom, start, end, field4 = parts[0], int(parts[1]), int(parts[2]), parts[3]
            catalog_contigs.add(chrom)
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
                "ref_copies": ref_copies,
                "start": start,
                "end": end
            }

    # 2. Collect contig lengths
    contig_lengths = {}
    contig_order = []
    if args.ref_fai:
        contig_lengths, contig_order = parse_fai(args.ref_fai)
    elif args.ref:
        fai_path = args.ref + ".fai"
        if os.path.exists(fai_path):
            contig_lengths, contig_order = parse_fai(fai_path)
        else:
            contig_lengths, contig_order = parse_fasta_lengths(args.ref)

    strling_contigs = collect_strling_contigs(args.strling)
    contigs = catalog_contigs | strling_contigs

    if not contig_lengths:
        sys.exit("Missing contig lengths; provide --ref-fai or --ref to populate VCF contig headers.")

    missing = sorted([c for c in contigs if c not in contig_lengths])
    if missing:
        sys.exit(f"Missing contig lengths for: {', '.join(missing)}")

    ordered_contigs = [c for c in contig_order if c in contigs]
    for c in sorted(contigs):
        if c not in ordered_contigs:
            ordered_contigs.append(c)

    # 3. Write VCF Header
    with open(args.output, "w") as out:
        out.write("##fileformat=VCFv4.2\n")
        out.write("##command=hipstr\n")
        for contig in ordered_contigs:
            out.write(f"##contig=<ID={contig},length={contig_lengths[contig]}>\n")
        out.write('##FILTER=<ID=PASS,Description="All filters passed">\n')
        out.write('##INFO=<ID=ID,Number=1,Type=String,Description="Locus ID from catalog">\n')
        out.write('##INFO=<ID=RU,Number=1,Type=String,Description="Repeat Unit">\n')
        out.write('##INFO=<ID=START,Number=1,Type=Integer,Description="Start position of the locus">\n')
        out.write('##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the locus">\n')
        out.write('##INFO=<ID=PERIOD,Number=1,Type=Integer,Description="Repeat unit length">\n')
        out.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        out.write('##FORMAT=<ID=REPCN,Number=2,Type=Float,Description="Repeat Copy Number">\n')
        out.write('##FORMAT=<ID=Q,Number=1,Type=Integer,Description="Quality score derived from STRling evidence">\n')
        out.write('##FORMAT=<ID=AD,Number=2,Type=Integer,Description="Anchored and Spanning read counts">\n')
        out.write('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Total Depth">\n')
        out.write(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{args.sample}\n")

        # 4. Parse STRling and Write Records
        with open_maybe_gz(args.strling) as f:
            header = f.readline()
            if not header:
                return
            col_map = {name: i for i, name in enumerate(header.strip().split("\t"))}
            
            for line in f:
                parts = line.strip().split("\t")
                chrom = parts[col_map["#chrom"]]
                left = safe_int(parts[col_map["left"]])
                right = safe_int(parts[col_map["right"]])
                repeatunit = parts[col_map["repeatunit"]] if "repeatunit" in col_map else "N"
                repeatunit = repeatunit if repeatunit else "N"

                # Match against catalog by exact coordinates; otherwise treat as novel
                locus = catalog.get((chrom, left, right))
                if locus and repeatunit == "N":
                    repeatunit = locus["motif"] if locus["motif"] else "N"
                period = len(repeatunit) if repeatunit != "N" else 1
                if locus:
                    locus_id = locus["id"]
                    ref_len = locus["end"] - locus["start"]
                    pos = locus["start"] + 1
                    end_1based = locus["end"]
                else:
                    locus_id = f"STRling_Novel_{chrom}_{left}"
                    ref_len = right - left
                    pos = left + 1
                    end_1based = right

                ref_copies = ref_len / period if period > 0 else 0.0
                ref_copies_rounded = round(ref_copies, 2)

                al1_est = safe_float(parts[col_map["allele1_est"]])
                al2_est = safe_float(parts[col_map["allele2_est"]])
                anchored = safe_int(parts[col_map["anchored_reads"]])
                spanning = safe_int(parts[col_map["spanning_reads"]])
                depth = safe_int(parts[col_map["depth"]])

                abs_al1 = ref_copies + al1_est
                abs_al2 = ref_copies + al2_est
                abs_al1_rounded = round(abs_al1, 2)
                abs_al2_rounded = round(abs_al2, 2)

                abs_al1_len = int(round(abs_al1 * period))
                abs_al2_len = int(round(abs_al2 * period))

                alt_lengths = []
                for allele_len in (abs_al1_len, abs_al2_len):
                    if allele_len != ref_len and allele_len not in alt_lengths:
                        alt_lengths.append(allele_len)
                alt_lengths.sort()

                alt_alleles = ",".join(
                    [build_repeat_sequence(repeatunit, allele_len) for allele_len in alt_lengths]
                ) if alt_lengths else "."

                allele_to_index = {val: idx + 1 for idx, val in enumerate(alt_lengths)}
                gt_allele1 = 0 if abs_al1_len == ref_len else allele_to_index.get(abs_al1_len, 0)
                gt_allele2 = 0 if abs_al2_len == ref_len else allele_to_index.get(abs_al2_len, 0)
                gt = f"{gt_allele1}/{gt_allele2}"

                q_score = anchored + spanning
                if q_score > 99:
                    q_score = 99

                info = f"ID={locus_id};RU={repeatunit};START={pos};END={end_1based};PERIOD={period}"
                fmt = "GT:REPCN:Q:AD:DP"
                sample_data = f"{gt}:{abs_al1_rounded:.2f},{abs_al2_rounded:.2f}:{q_score}:{anchored},{spanning}:{depth}"

                # POS in VCF is 1-based. Catalog/STRling 'left' is 0-based start.
                ref_seq = build_repeat_sequence(repeatunit, ref_len)
                out.write(f"{chrom}\t{pos}\t{locus_id}\t{ref_seq}\t{alt_alleles}\t.\tPASS\t{info}\t{fmt}\t{sample_data}\n")

if __name__ == "__main__":
    main()
