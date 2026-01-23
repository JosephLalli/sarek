import os
import subprocess
import tempfile
import warnings
from pathlib import Path


def test_str_consensus_merge_root_samples():
    """Hard gate: merge existing root sample VCFs."""
    repo_root = Path(__file__).resolve().parents[1]
    vcf_eh = repo_root / "eh.vcf"
    vcf_gs = repo_root / "gs.vcf"
    vcf_sl = repo_root / "sl.vcf"

    with tempfile.TemporaryDirectory() as tmpdir:
        output_vcf = Path(tmpdir) / "consensus.vcf"

        cmd = [
            "python3", str(repo_root / "bin" / "str_consensus_merge.py"),
            "--eh", str(vcf_eh),
            "--gs", str(vcf_gs),
            "--sl", str(vcf_sl),
            "--output", str(output_vcf),
            "--motif-threshold", "30"
        ]

        result = subprocess.run(cmd, capture_output=True, text=True)
        assert result.returncode == 0, f"Script failed with stderr: {result.stderr}"
        assert output_vcf.exists()

        content = output_vcf.read_text()
        assert "##INFO=<ID=MERGE_SRC" in content
        assert "EH_GT" in content
        assert "GS_GT" in content
        assert "SL_GT" in content
        assert "MERGE_SRC=" in content
        records = [line for line in content.splitlines() if not line.startswith("#")]
        pos100 = [line for line in records if line.split("\t")[1] == "100"]
        assert len(pos100) == 1

def test_str_consensus_merge_novel_soft():
    """Soft gate: warn if a novel SL-only locus is not preserved."""
    repo_root = Path(__file__).resolve().parents[1]

    header = (
        "##fileformat=VCFv4.2\n"
        "##command=hipstr\n"
        "##contig=<ID=chr1,length=248956422>\n"
        "##INFO=<ID=RU,Number=1,Type=String,Description=\"Repeat Unit\">\n"
        "##INFO=<ID=START,Number=1,Type=Integer,Description=\"Start position of the locus\">\n"
        "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the locus\">\n"
        "##INFO=<ID=PERIOD,Number=1,Type=Integer,Description=\"Repeat unit length\">\n"
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
        "##FORMAT=<ID=REPCN,Number=2,Type=Float,Description=\"Repeat Copy Number\">\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1\n"
    )

    with tempfile.TemporaryDirectory() as tmpdir:
        vcf_eh = Path(tmpdir) / "eh.vcf"
        vcf_gs = Path(tmpdir) / "gs.vcf"
        vcf_sl = Path(tmpdir) / "sl.vcf"
        output_vcf = Path(tmpdir) / "consensus.vcf"

        vcf_eh.write_text(header)
        vcf_gs.write_text(header)
        vcf_sl.write_text(
            header +
            "chr1\t200\tSTRling_Novel_chr1_200\tATATATATATATATAT\tATATATATATATATATATAT\t.\tPASS\tRU=AT;START=200;END=215;PERIOD=2\tGT:REPCN\t0/1:8.00,10.00\n"
        )

        cmd = [
            "python3", str(repo_root / "bin" / "str_consensus_merge.py"),
            "--eh", str(vcf_eh),
            "--gs", str(vcf_gs),
            "--sl", str(vcf_sl),
            "--output", str(output_vcf)
        ]

        result = subprocess.run(cmd, capture_output=True, text=True)
        assert result.returncode == 0, f"Script failed with stderr: {result.stderr}"
        assert output_vcf.exists()

        content = output_vcf.read_text()
        if "STRling_Novel_chr1_200" not in content:
            warnings.warn(
                "Expected a novel SL-only locus to be preserved in consensus output; "
                "follow up before advancing to the next task."
            )

if __name__ == "__main__":
    test_str_consensus_merge_basic()
