import subprocess
import os
import pytest
import tempfile

def test_strling_to_vcf_basic():
    """Test basic conversion of STRling output to VCF."""
    
    # 1. Setup Mock Data
    with tempfile.TemporaryDirectory() as tmpdir:
        catalog_path = os.path.join(tmpdir, "catalog.bed")
        strling_path = os.path.join(tmpdir, "strling.txt")
        output_vcf = os.path.join(tmpdir, "output.vcf")
        
        # Ref length = 21, Motif = 3 -> 7 units
        with open(catalog_path, "w") as f:
            f.write("chr1\t100\t121\tID=STR1;MOTIFS=CAG;STRUC=(CAG)n\n")
            
        with open(strling_path, "w") as f:
            f.write("#chrom\tleft\tright\trepeatunit\tallele1_est\tallele2_est\tanchored_reads\tspanning_reads\tspanning_pairs\texpected_spanning_pairs\tspanning_pairs_pctl\tleft_clips\tright_clips\tunplaced_pairs\tdepth\tsum_str_counts\n")
            f.write("chr1\t100\t121\tCAG\t0.00\t5.00\t10\t20\t0\t0\t0\t0\t0\t0\t30\t100\n")
            
        # 2. Run the Script
        cmd = [
            "python3", "bin/strling_to_vcf.py",
            "--strling", strling_path,
            "--catalog", catalog_path,
            "--output", output_vcf,
            "--sample", "SAMPLE1"
        ]
        
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        # 3. Verify Output (Green phase expectations)
        assert result.returncode == 0, f"Script failed with stderr: {result.stderr}"
        assert os.path.exists(output_vcf)
        
        with open(output_vcf, "r") as f:
            content = f.read()
            assert "##fileformat=VCFv4.2" in content
            assert "SAMPLE1" in content
            assert "STR1" in content
            # Check for absolute copy number: Ref (7) + SL_est (0, 5) -> 7, 12
            # VCF usually represents this in REPCN or GT (as allele lengths)
            # For now let's just check if it contains some expected strings
            assert "GT" in content

if __name__ == "__main__":
    test_strling_to_vcf_basic()
