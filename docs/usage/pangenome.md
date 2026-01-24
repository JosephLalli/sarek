# nf-core/sarek: Pangenome Analysis

## Introduction

Sarek supports pangenome-aware analysis, leveraging graph-based alignment and variant calling to improve accuracy in complex genomic regions. This implementation integrates [VG Giraffe](https://github.com/vgteam/vg) for pangenome alignment and a specialized [DeepVariant](https://github.com/google/deepvariant) module for pangenome-aware variant calling.

## Pangenome Alignment (VG Giraffe)

VG Giraffe is a fast, pangenome-aware short-read aligner. In Sarek, Giraffe can be used as a drop-in replacement for standard linear aligners (like BWA-MEM).

### Key Parameters

* `--aligner giraffe`: Enables VG Giraffe alignment.
* `--pangenome_gbz`: Path to the pangenome graph in GBZ format.
* `--pangenome_dist`: Path to the distance index for the graph.
* `--pangenome_min`: Path to the minimizer index for the graph.
* `--pangenome_zipcodes`: (Optional) Path to the zipcodes file for the graph.
* `--pangenome_ref_paths`: (Optional) Text file containing names of paths in the graph to surject alignments onto. Defaults to the primary reference contigs.

### Reference Naming and Surjection

By default, Sarek surjects graph-based alignments to a linear reference (BAM/CRAM). 
Pangenome graphs often use complex path names (e.g., `CHM13#0#chr1`). 
Sarek includes logic to automatically strip common prefixes like `CHM13#0#` during surjection to ensure the resulting BAM files are compatible with standard downstream tools.

## Pangenome-Aware DeepVariant

Sarek includes a specialized subworkflow for Pangenome-Aware DeepVariant. This tool uses the pangenome graph topology during variant calling to better resolve reads in highly polymorphic regions.

### Usage

To enable pangenome-aware variant calling, add `deepvariant_pangenome` to your `--tools` parameter:

```bash
nextflow run nf-core/sarek \
    --input samplesheet.csv \
    --aligner giraffe \
    --tools deepvariant_pangenome \
    --pangenome_gbz path/to/graph.gbz \
    --pangenome_dist path/to/graph.dist \
    --pangenome_min path/to/graph.min \
    --genome CHM13 \
    -profile docker
```

### Important Notes

* **Mutual Exclusivity:** `deepvariant` and `deepvariant_pangenome` are mutually exclusive.
* **Reference Consistency:** The FASTA reference provided must match the coordinate system of the surjected BAM.
* **CHM13 Support:** When using CHM13-based pangenomes, set `--genome CHM13` to automatically apply optimized arguments to DeepVariant.

## Evaluation & Testing

For developers wishing to verify the pangenome implementation, Sarek includes specialized `nf-test` integration tests.

### Running Integration Tests

You can run the full CHM13 pangenome integration test using:

```bash
nf-test test tests/chm13_pangenome_full.nf.test --profile test,docker
```

This test verifies the entire flow from FASTQ to DeepVariant VCF using a rebased CHM13 subset.
