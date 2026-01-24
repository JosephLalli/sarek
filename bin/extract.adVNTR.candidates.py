#!/usr/bin/env python3
import argparse, gzip, sys

def open_maybe_gz(path: str):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "rt")

def parse_tags(field4: str):
    """
    TRGT repeat-definition BED: 4th column is tag string with mandatory fields:
      ID=..., MOTIFS=comma,separated, STRUC=...
    (may contain additional tags too)
    """
    tags = {}
    for part in field4.strip().split(";"):
        part = part.strip()
        if not part:
            continue
        if "=" in part:
            k, v = part.split("=", 1)
            tags[k] = v
        else:
            tags[part] = True
    return tags

def main():
    ap = argparse.ArgumentParser(
        description="Filter TRGT repeat-definition BED for candidate VNTRs by motif length and >=2 copies (estimated from region length)."
    )
    ap.add_argument("trgt_bed", help="TRGT repeat definition BED (optionally .gz)")
    ap.add_argument("--min_motif", type=int, default=6, help="Minimum motif length (bp). Default: 6")
    ap.add_argument("--max_motif", type=int, default=100, help="Maximum motif length (bp) for 'core' set. Default: 100")
    ap.add_argument("--min_copies", type=float, default=2.0, help="Minimum estimated copies. Default: 2.0")
    ap.add_argument("--out_core_bed", default=None, help="Write passing (min_motif..max_motif) loci to BED (optional)")
    ap.add_argument("--out_relaxed_bed", default=None, help="Write passing (>=min_motif, no upper bound) loci to BED (optional)")
    args = ap.parse_args()

    n_total = 0
    n_core = 0
    n_relaxed = 0

    core_out = open(args.out_core_bed, "wt") if args.out_core_bed else None
    rel_out  = open(args.out_relaxed_bed, "wt") if args.out_relaxed_bed else None

    with open_maybe_gz(args.trgt_bed) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            chrom, start_s, end_s, field4, *rest = line.rstrip("\n").split("\t")
            start = int(start_s)
            end = int(end_s)
            region_len = end - start  # BED is 0-based, half-open; length=end-start

            tags = parse_tags(field4)
            motifs = tags.get("MOTIFS")
            if not motifs:
                continue

            motif_list = [m.strip() for m in motifs.split(",") if m.strip()]
            if not motif_list:
                continue

            # Decide locus passes if ANY motif meets criteria (works well for simple loci;
            # for complex STRUC, this is an approximation).
            passed_relaxed = False
            passed_core = False

            for m in motif_list:
                mlen = len(m)
                if mlen < args.min_motif:
                    continue
                copies_est = region_len / mlen if mlen > 0 else 0.0
                if copies_est < args.min_copies:
                    continue

                passed_relaxed = True
                if mlen <= args.max_motif:
                    passed_core = True

            n_total += 1
            if passed_relaxed:
                n_relaxed += 1
                if rel_out:
                    rel_out.write(line)
            if passed_core:
                n_core += 1
                if core_out:
                    core_out.write(line)

    if core_out: core_out.close()
    if rel_out: rel_out.close()

    print(f"Parsed loci (non-header BED lines with MOTIFS): {n_total}")
    print(f"Core candidates (copies>={args.min_copies} AND motif_len in [{args.min_motif},{args.max_motif}]): {n_core}")
    print(f"Relaxed candidates (copies>={args.min_copies} AND motif_len>={args.min_motif}, no upper bound): {n_relaxed}")
    print(f"Additional if allowing motif_len > {args.max_motif}: {n_relaxed - n_core}")

if __name__ == "__main__":
    main()
