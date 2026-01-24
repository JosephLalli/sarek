import os
import subprocess

outdir = "assets/test_data/pangenome/haplotype-sampling/"
os.makedirs(outdir, exist_ok=True)

# Reference snippets from grep output
# chr6: TTTCATTCAGTTGGCCACTGCTGAGCAGCTGAGAAGGTGGCGACGTAGGGGCCATGGGGCTGGGCCGGGTCCTGCTGTTT
# chr19: AATAACATCCTGTGCGCTGCTGAGCTGAGCTGGGGCGCGGCCGCCTGTCTGCACCGGCAGCACCATGTTGCTCATGGTCG

# VCF Creation
vcf_content = """##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##contig=<ID=chr6,length=13026>
##contig=<ID=chr19,length=14346>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1
chr6	1000	.	C	A	.	PASS	.	GT	0|1
chr19	1000	.	G	T	.	PASS	.	GT	1|0
"""

vcf_path = os.path.join(outdir, "test_pangenie_panel.vcf")
with open(vcf_path, "w") as f:
    f.write(vcf_content)

# Compress and Index VCF
subprocess.run(f"bgzip -f {vcf_path}", shell=True, check=True)
subprocess.run(f"tabix -f -p vcf {vcf_path}.gz", shell=True, check=True)

# EH Catalog (JSON)
# chr6:1000 is likely G (based on VCF above, assuming REF correct). 
# Actually, let's just use a region. 
# chr6:1020-1025 "TTGGC"
eh_json = """[
  {
    "LocusId": "TestSTR_chr6",
    "LocusStructure": "(T)*",
    "ReferenceRegion": "chr6:1020-1025",
    "VariantType": "Repeat"
  }
]"""
# chr6:1-5 is TTTC. Not perfect T repeat but close enough for testing pipeline mechanics.

with open(os.path.join(outdir, "test_str_catalog.json"), "w") as f:
    f.write(eh_json)

# STRling/GangSTR Catalog (BED)
# chr	start	end	name	period	motif
bed_content = "chr6\t1020\t1025\tTestSTR_chr6\t.\t.\tT\n"

with open(os.path.join(outdir, "test_str_catalog.bed"), "w") as f:
    f.write(bed_content)

print("Assets created successfully.")
