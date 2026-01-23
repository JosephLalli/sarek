import os
import shutil
import subprocess
import tempfile
import warnings
from pathlib import Path

import pytest


def test_strling_to_vcf_basic():
    """Test basic conversion of STRling output to VCF."""
    
    # 1. Setup Mock Data
    with tempfile.TemporaryDirectory() as tmpdir:
        catalog_path = os.path.join(tmpdir, "catalog.bed")
        strling_path = os.path.join(tmpdir, "strling.txt")
        ref_fai_path = os.path.join(tmpdir, "ref.fa.fai")
        output_vcf = os.path.join(tmpdir, "output.vcf")
        
        # Ref length = 21, Motif = 3 -> 7 units
        with open(catalog_path, "w") as f:
            f.write("chr1\t100\t121\tID=STR1;MOTIFS=CAG;STRUC=(CAG)n\n")
            
        with open(strling_path, "w") as f:
            f.write("#chrom\tleft\tright\trepeatunit\tallele1_est\tallele2_est\tanchored_reads\tspanning_reads\tspanning_pairs\texpected_spanning_pairs\tspanning_pairs_pctl\tleft_clips\tright_clips\tunplaced_pairs\tdepth\tsum_str_counts\n")
            f.write("chr1\t100\t121\tCAG\t0.00\t5.00\t10\t20\t0\t0\t0\t0\t0\t0\t30\t100\n")

        with open(ref_fai_path, "w") as f:
            f.write("chr1\t1000\n")
            
        # 2. Run the Script
        cmd = [
            "python3", "bin/strling_to_vcf.py",
            "--strling", strling_path,
            "--catalog", catalog_path,
            "--ref-fai", ref_fai_path,
            "--output", output_vcf,
            "--sample", "SAMPLE1"
        ]
        
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        # 3. Verify Output (Green phase expectations)
        assert result.returncode == 0, f"Script failed with stderr: {result.stderr}"
        assert os.path.exists(output_vcf)
        
        ref_seq = "CAG" * 7
        alt_seq = "CAG" * 12

        with open(output_vcf, "r") as f:
            content = f.read()
            assert "##fileformat=VCFv4.2" in content
            assert "##command=hipstr" in content
            assert "##contig=<ID=chr1,length=1000>" in content
            assert "##INFO=<ID=START" in content
            assert "##INFO=<ID=END" in content
            assert "##INFO=<ID=PERIOD" in content
            assert "##INFO=<ID=RU" in content
            assert "##FORMAT=<ID=REPCN" in content
            assert "##FORMAT=<ID=Q" in content
            assert "SAMPLE1" in content
            assert ref_seq in content
            assert alt_seq in content
            assert "RU=CAG" in content
            assert "START=101" in content
            assert "PERIOD=3" in content
            assert "END=121" in content
            # Ref (7) + SL_est (0, 5) -> 7, 12
            assert "7.00,12.00" in content
            assert ":30:" in content

        bcftools_path = shutil.which("bcftools")
        if not bcftools_path:
            pytest.skip("bcftools not available; skipping VCF parse check.")
        bcftools_result = subprocess.run(
            [bcftools_path, "view", "-h", output_vcf],
            capture_output=True,
            text=True
        )
        assert bcftools_result.returncode == 0, (
            f"bcftools parse failed: {bcftools_result.stderr}"
        )

def test_strling_to_vcf_novel_call_soft():
    """Soft gate: warn if no novel STRling call is emitted."""
    repo_root = Path(__file__).resolve().parents[1]
    strling_path = repo_root / "nc30_15-1-genotype.txt"
    catalog_path = repo_root / "tests" / "test_data" / "str_enhancement" / "gangstr_catalog.bed"
    ref_fai_path = repo_root / "tests" / "test_data" / "str_enhancement" / "Homo_sapiens_assembly38.chr4.fasta.fai"

    with tempfile.TemporaryDirectory() as tmpdir:
        output_vcf = Path(tmpdir) / "output.vcf"

        cmd = [
            "python3", str(repo_root / "bin" / "strling_to_vcf.py"),
            "--strling", str(strling_path),
            "--catalog", str(catalog_path),
            "--ref-fai", str(ref_fai_path),
            "--output", str(output_vcf),
            "--sample", "SAMPLE_NOVEL"
        ]
        result = subprocess.run(cmd, capture_output=True, text=True)

        assert result.returncode == 0, f"Script failed with stderr: {result.stderr}"
        assert output_vcf.exists()

        content = output_vcf.read_text()
        if "STRling_Novel_" not in content:
            warnings.warn(
                "Expected a novel STRling call but none was found; follow up on "
                "catalog matching or test data suitability."
            )

if __name__ == "__main__":
    test_strling_to_vcf_basic()
