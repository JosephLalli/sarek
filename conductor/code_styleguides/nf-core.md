# nf-core Style Guide Summary

This document summarizes nf-core-specific conventions for Nextflow pipeline development. For core Nextflow syntax and coding conventions, see `nextflow.md`.

## 1. Module Organization

Each nf-core module requires these files:

```
modules/nf-core/<tool>/<subtool>/
├── main.nf              # Process definition
├── meta.yml             # Metadata and documentation (required)
├── environment.yml      # Conda environment specification
└── tests/
    ├── main.nf.test     # nf-test test file
    └── main.nf.test.snap # Snapshot file (auto-generated)
```

- **Module path mirrors process name:** `modules/nf-core/bwa/mem/main.nf` → `BWA_MEM`
- **One process per module:** Each module should contain a single atomic tool invocation.
- **Pre-commit hooks:** Install with `pre-commit install` to enforce formatting and linting.

## 2. The Meta Map

The meta map is an nf-core convention for carrying sample metadata through pipelines. It replaces creating separate channels for each metadata attribute.

- **Required field:** `meta.id` - unique sample identifier.
- **Common optional fields:** `single_end` (boolean), `strandedness`, custom fields as needed.
- **Standard structure:** `[ id: 'sample1', single_end: false ]`

**Passing through processes:**
```nextflow
input:
tuple val(meta), path(reads)

output:
tuple val(meta), path("*.bam"), emit: bam
```

**Numbered meta variables:** Use `meta`, `meta2`, `meta3` when inputs come from different sources (e.g., sample vs. reference).

**Modifying meta maps:**
```nextflow
meta + [ new_field: 'value' ]              // Add field
meta.subMap(['id', 'single_end'])          // Extract subset
meta.findAll { it.key != 'unwanted' }      // Remove field
```

## 3. ext Properties (modules.config)

nf-core modules use `task.ext` properties for configuration, keeping processes generic and reusable.

| Property | Purpose |
|----------|---------|
| `ext.args` | Primary tool arguments |
| `ext.args2`, `ext.args3` | Additional arguments for multi-tool modules |
| `ext.prefix` | Output filename prefix |
| `ext.when` | Conditional execution |

**In the module (main.nf):**
```nextflow
script:
def args   = task.ext.args   ?: ''
def prefix = task.ext.prefix ?: "${meta.id}"
```

**In conf/modules.config:**
```nextflow
withName: 'BWA_MEM' {
    ext.args   = { params.bwa_sort ? '-M' : '' }
    ext.prefix = { "${meta.id}.aligned" }
    publishDir = [
        path: { "${params.outdir}/aligned" },
        mode: params.publish_dir_mode,
        saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
    ]
}
```

- **Always use closures:** Wrap `ext.args` in `{ }` to access `params` and `meta` at runtime.
- **Never hardcode arguments:** Tool arguments belong in modules.config, not in the process itself.

## 4. Version Tracking

Every module must emit software versions for reproducibility.

```nextflow
output:
path "versions.yml", emit: versions

script:
"""
cat <<-END_VERSIONS > versions.yml
"${task.process}":
    tool: \$(tool --version | sed 's/tool //')
END_VERSIONS
"""
```

- **Quote `${task.process}`:** YAML requires quotes for strings containing colons.
- **Use heredoc with dash:** `<<-END_VERSIONS` allows indentation in the script.

## 5. Harshil Alignment

nf-core uses vertical alignment of punctuation for readability. Prettier handles basic formatting; Harshil alignment is manual.

**Curly bracket alignment (imports):**
```nextflow
include { SAMTOOLS_SORT      } from '../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_INDEX     } from '../modules/nf-core/samtools/index/main'
include { BAM_STATS_SAMTOOLS } from '../subworkflows/bam_stats_samtools/main'
```

**Equals sign alignment:**
```nextflow
stats    = TOOL.out.stats    // channel: [ val(meta), path(stats)    ]
flagstat = TOOL.out.flagstat // channel: [ val(meta), path(flagstat) ]
idxstats = TOOL.out.idxstats // channel: [ val(meta), path(idxstats) ]
```

**Comma alignment (outputs):**
```nextflow
tuple val(meta), path("*.bam"), emit: bam    , optional: true
tuple val(meta), path("*.log"), emit: log
path  "versions.yml"          , emit: versions
```

**Channel comment format:**
```nextflow
ch_input    // channel: [ val(meta), path(reads) ]
ch_fasta    // channel: [ val(meta), path(fasta) ]
```

## 6. Documentation (meta.yml)

Every module requires a `meta.yml` with structured metadata.

```yaml
name: "tool_subtool"
description: Brief description of what the module does
keywords:
  - keyword1
  - keyword2

tools:
  - toolname:
      description: Tool description
      homepage: https://tool-homepage.com
      documentation: https://tool-docs.com
      doi: "10.xxxx/xxxxx"
      licence: ["MIT"]

input:
  - - meta:
        type: map
        description: |
          Groovy Map containing sample information
          e.g. [ id:'test', single_end:false ]
    - reads:
        type: file
        description: Input FASTQ files
        pattern: "*.{fq,fastq,fq.gz,fastq.gz}"

output:
  - bam:
      - meta:
          type: map
          description: Sample information
      - "*.bam":
          type: file
          description: Aligned reads
          pattern: "*.bam"
  - versions:
      - "versions.yml":
          type: file
          description: Software versions
          pattern: "versions.yml"

authors:
  - "@github_username"
maintainers:
  - "@github_username"
```

**In-code documentation:** Use channel signature comments in workflows.
```nextflow
workflow ALIGN {
    take:
    ch_reads    // channel: [ val(meta), path(reads) ]
    ch_index    // channel: [ val(meta), path(index) ]

    emit:
    bam      = ALIGNER.out.bam      // channel: [ val(meta), path(bam) ]
    versions = ch_versions          // channel: [ path(versions.yml)  ]
}
```

## 7. Testing (nf-test)

nf-core uses nf-test for module and workflow testing.

**Test file location:** `modules/<tool>/tests/main.nf.test`

**Required tags:** `modules`, `modules_nfcore`, `<tool>`, `<tool>/<subtool>`

**Test naming:** `"organism - filetype - description"` format

**Basic test structure:**
```nextflow
nextflow_process {
    name "Test Process TOOL_SUBTOOL"
    script "../main.nf"
    process "TOOL_SUBTOOL"

    tag "modules"
    tag "modules_nfcore"
    tag "tool"
    tag "tool/subtool"

    test("sarscov2 - bam - basic") {
        when {
            process {
                """
                input[0] = [
                    [ id:'test', single_end:false ],
                    file(params.modules_testdata_base_path + 'path/to/test.bam', checkIfExists: true)
                ]
                """
            }
        }

        then {
            assertAll(
                { assert process.success },
                { assert snapshot(process.out).match() }
            )
        }
    }

    test("sarscov2 - bam - stub") {
        options "-stub"
        when {
            process {
                """
                input[0] = [
                    [ id:'test', single_end:false ],
                    file(params.modules_testdata_base_path + 'path/to/test.bam', checkIfExists: true)
                ]
                """
            }
        }

        then {
            assertAll(
                { assert process.success },
                { assert snapshot(process.out).match() }
            )
        }
    }
}
```

**Minimum assertions:**
```nextflow
assertAll(
    { assert process.success },
    { assert snapshot(process.out.versions).match("versions") }
)
```

**Non-deterministic outputs:** Use existence checks, size checks, or content matching instead of snapshots.
```nextflow
{ assert path(process.out.bam[0][1]).exists() }
{ assert path(process.out.bam[0][1]).size() > 1000 }
{ assert path(process.out.txt[0][1]).text.contains("Expected") }
```

**Running tests:**
```bash
nf-core modules test <tool>/<subtool>           # Run tests
nf-core modules test <tool>/<subtool> --update  # Update snapshots
```

## 8. Pipeline File Structure

Standard nf-core pipeline layout:

```
pipeline/
├── main.nf                 # Entry point
├── nextflow.config         # Main config
├── nextflow_schema.json    # Parameter schema (generated)
├── modules.json            # Installed module versions
├── conf/
│   ├── base.config         # Resource labels
│   ├── modules.config      # Module ext.args and publishDir
│   ├── igenomes.config     # Reference genome paths
│   └── test.config         # Test profile
├── modules/
│   ├── local/              # Pipeline-specific modules
│   └── nf-core/            # Installed nf-core modules
├── subworkflows/
│   ├── local/              # Pipeline-specific subworkflows
│   └── nf-core/            # Installed nf-core subworkflows
├── workflows/              # Main workflow definitions
├── lib/                    # Groovy helper classes
├── bin/                    # Scripts callable from processes
├── assets/                 # Additional resources
└── docs/                   # Documentation
```

## 9. Resource Labels

Use standard nf-core labels for process resources:

| Label | Typical Use |
|-------|-------------|
| `process_single` | Single-threaded, low memory |
| `process_low` | 2 CPUs, 12 GB memory |
| `process_medium` | 6 CPUs, 36 GB memory |
| `process_high` | 12 CPUs, 72 GB memory |
| `process_long` | Extended time limits |
| `process_high_memory` | Memory-intensive tasks |

Define in `conf/base.config`, reference in processes with `label 'process_medium'`.

**BE CONSISTENT.** Follow established nf-core patterns. When contributing modules, use the nf-core template exactly.

*Sources:*
- *[nf-core Contributing Guidelines](https://nf-co.re/docs/contributing)*
- *[nf-test Documentation](https://www.nf-test.com/docs/)*
