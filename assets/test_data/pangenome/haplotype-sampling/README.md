# Giraffe Integration Test Data

This directory contains test data for the vg giraffe pangenome aligner integration.

## Test Data Source

Data from vg repository haplotype-sampling test:
- Source: https://ucsc-ci.com/vgteam/vg/-/tree/e6f03265/test/haplotype-sampling
- Region: MICB-KIR3DL1 locus (chr6:31498140 and chr19:54816436)
- Sample: HG003
- Size: ~27 KB reference sequences, ~444 KB FASTQ

## Files

### Input Data
- `HG003_R1.fastq.gz` - Paired-end R1 reads (~781K, ~10k read pairs)
- `HG003_R2.fastq.gz` - Paired-end R2 reads (~816K)
- `HG003.kff` - K-mer file (375K)
- `micb-kir3dl1.gfa` - Pangenome graph source (1.5M)

Reads extracted from HG003 30X CRAM (CHM13-aligned) for MICB-KIR3DL1 regions:
- chr6:31340000-31375000 (MICB)
- chr19:18915000-18955000 (KIR3DL1)

### Giraffe Indices
- `micb-kir3dl1.giraffe.gbz` - GBZ graph index (79K)
- `micb-kir3dl1.dist` - Distance index (114K)
- `micb-kir3dl1.shortread.withzip.min` - Minimizer index (518K)

### Reference Sequences
- `micb-kir3dl1.GRCh38.fa/fai/dict` - GRCh38 reference paths
- `micb-kir3dl1.CHM13.fa/fai/dict` - CHM13 reference paths

### Test Configuration
- `samplesheet.csv` - Sarek input samplesheet for HG003
- `params_giraffe_test.json` - Parameters file for giraffe integration test

## Running the Integration Test

From the repository root directory:

```bash
nextflow run main.nf -params-file assets/test_data/pangenome/haplotype-sampling/params_giraffe_test.json
```

Or with explicit syntax mode:

```bash
NXF_SYNTAX_PARSER=v2 nextflow run main.nf -params-file assets/test_data/pangenome/haplotype-sampling/params_giraffe_test.json
```

## Expected Output

The test should:
1. Run fastq preprocessing on HG003.fq.gz
2. Align reads to pangenome using vg giraffe
3. Surject alignments to GRCh38 reference (BAM output)
4. Sort and index BAM files

Output will be in `test_output/giraffe/` directory.

## Index Generation

Indices were generated using vg v1.70.0:

```bash
# Generate giraffe indices from GFA
vg autoindex --workflow giraffe \
    --prefix micb-kir3dl1 \
    --gfa micb-kir3dl1.gfa

# Extract reference sequences
vg paths -x micb-kir3dl1.giraffe.gbz -F -Q GRCh38 > micb-kir3dl1.GRCh38.fa
vg paths -x micb-kir3dl1.giraffe.gbz -F -Q CHM13 > micb-kir3dl1.CHM13.fa

# Create auxiliary files
samtools faidx micb-kir3dl1.GRCh38.fa
samtools dict micb-kir3dl1.GRCh38.fa -o micb-kir3dl1.GRCh38.dict
samtools faidx micb-kir3dl1.CHM13.fa
samtools dict micb-kir3dl1.CHM13.fa -o micb-kir3dl1.CHM13.dict
```
