# Nextflow DSL2 Style Guide

This document defines code style rules and best practices for Nextflow DSL2 pipelines, emphasizing strict syntax compliance for maintainability and future compatibility.

## 1. Strict Syntax Mode

- **Enable strict syntax:** Set `NXF_SYNTAX_PARSER=v2` in your environment. This enables the strict parser that catches errors early and ensures forward compatibility.
  ```bash
  export NXF_SYNTAX_PARSER=v2
  ```
- **Why it matters:** Strict syntax restricts Groovy to specific patterns that Nextflow can analyze and optimize. Code that works in lenient mode may break in future versions.
- **Validation:** Run `nextflow lint` to check your code before execution.

## 2. Script Structure

Nextflow scripts must contain either declarations (processes, workflows, functions) OR executable statements—never both at the top level.

- **Separate declarations from execution:**
  ```groovy
  // WRONG - mixing declarations and statements
  process HELLO { /* ... */ }
  println 'Starting pipeline'  // Error: statement at top level

  // CORRECT - statements inside workflow block
  process HELLO { /* ... */ }
  workflow {
      println 'Starting pipeline'
      HELLO()
  }
  ```

- **File organization:**
  - `main.nf` - Entry point with main workflow
  - `modules/` - Process definitions (one process per file recommended)
  - `subworkflows/` - Reusable workflow compositions
  - `lib/` - Complex Groovy classes and helper code
  - `conf/` - Configuration profiles

- **No import statements:** Use fully qualified class names instead.
  ```groovy
  // WRONG
  import groovy.json.JsonSlurper
  def data = new JsonSlurper().parseText(text)

  // CORRECT
  def data = new groovy.json.JsonSlurper().parseText(text)
  ```

## 3. Variables and Types

### Variable Declarations

- **Always use `def`:** All variables must be declared with `def`. Type-specific declarations like `String x` are not allowed in strict mode.
  ```groovy
  // WRONG
  String name = 'sample'
  Integer count = 42

  // CORRECT
  def name = 'sample'
  def count = 42
  ```

- **Type annotations (Nextflow 25.10+):** Use colon syntax for optional type hints that aid documentation and IDE support.
  ```groovy
  def name: String = 'sample'
  def count: Integer = 42
  def (a: Integer, b: Integer) = [1, 2]
  ```

### Naming Conventions

- **Processes:** `UPPERCASE_WITH_UNDERSCORES` (e.g., `FASTQC`, `BWA_MEM`)
- **Workflows:** `UPPERCASE_WITH_UNDERSCORES` for subworkflows, `lowercase` acceptable for entry workflow
- **Channels:** `snake_case` with descriptive suffix (e.g., `reads_ch`, `bam_indexed_ch`)
- **Variables:** `snake_case` (e.g., `sample_id`, `output_dir`)
- **Parameters:** `snake_case` in params block (e.g., `params.input_dir`)
- **Functions:** `camelCase` or `snake_case` consistently

### String Handling

- **Supported string types:** Single-quoted, double-quoted, triple-quoted (multiline), and slashy strings.
  ```groovy
  def simple = 'no interpolation'
  def interpolated = "sample: ${sample_id}"
  def multiline = """
      Multiple
      lines
      """
  def regex = ~/pattern\d+/
  ```

- **Dollar-slashy strings are prohibited:** `$/.../$` syntax is not allowed.

- **Slashy strings cannot be interpolated:** Use double-quoted strings when you need variable substitution.

## 4. Process Definitions

Processes are the atomic units of computation in Nextflow. Each process should do one thing well.

### Basic Structure

```groovy
process PROCESS_NAME {
    tag "${meta.id}"
    label 'process_medium'

    container 'biocontainers/tool:version'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.bam"), emit: bam

    script:
    """
    tool command --input ${reads} --output ${meta.id}.bam
    """
}
```

### Input/Output Rules

- **Quote environment variable names:**
  ```groovy
  // WRONG
  input:
  env FOO

  // CORRECT
  input:
  env 'FOO'
  ```

- **Use explicit qualifiers:** Always specify `val()`, `path()`, `tuple()`, etc.
  ```groovy
  // WRONG - implicit typing
  input:
  sample_id
  reads

  // CORRECT - explicit qualifiers
  input:
  val(sample_id)
  path(reads)
  ```

- **Named outputs with `emit`:** Always name outputs for clarity in workflow wiring.
  ```groovy
  output:
  tuple val(meta), path("*.bam"), emit: aligned
  path("*.log"), emit: logs
  ```

- **Type annotations for inputs/outputs (25.10+):**
  ```groovy
  process FASTQC {
      input:
      (id, fastq_1, fastq_2): Tuple<String, Path, Path>

      output:
      file("fastqc_${id}_logs")
  }
  ```

### Script Section

- **The `script:` label can only be omitted** if the process has no other labeled sections (input, output, when, stub).
  ```groovy
  // OK - no other sections
  process SIMPLE {
      """
      echo "hello"
      """
  }

  // REQUIRED - has input section
  process WITH_INPUT {
      input:
      val(x)

      script:  // Required!
      """
      echo ${x}
      """
  }
  ```

- **Avoid `shell` section:** It is deprecated. Use `script` with proper escaping.

- **Use `stub` for testing:** Define stub blocks for dry-run testing.
  ```groovy
  process ALIGN {
      // ... input/output ...

      stub:
      """
      touch ${meta.id}.bam
      touch ${meta.id}.bam.bai
      """

      script:
      """
      bwa mem ... > ${meta.id}.bam
      samtools index ${meta.id}.bam
      """
  }
  ```

### Directives

- **Avoid `when` directive:** Implement conditional logic in the calling workflow instead.
  ```groovy
  // DISCOURAGED - logic inside process
  process OPTIONAL_STEP {
      when:
      params.run_optional
      // ...
  }

  // PREFERRED - logic in workflow
  workflow {
      if (params.run_optional) {
          OPTIONAL_STEP(input_ch)
      }
  }
  ```

- **Avoid `each` input qualifier:** Use `combine` operator instead for cross-product operations.
  ```groovy
  // DISCOURAGED
  process ALIGN {
      input:
      path(seq)
      each mode

      // ...
  }

  // PREFERRED
  process ALIGN {
      input:
      tuple path(seq), val(mode)
      // ...
  }

  workflow {
      sequences = channel.fromPath('*.fa')
      methods = channel.of('regular', 'espresso')
      ALIGN(sequences.combine(methods))
  }
  ```

## 5. Workflow Definitions

Workflows compose processes and subworkflows into pipelines.

### Structure

```groovy
workflow SUBWORKFLOW_NAME {
    take:
    input_ch      // Named inputs
    reference

    main:
    STEP_ONE(input_ch)
    STEP_TWO(STEP_ONE.out.results, reference)

    emit:
    final_output = STEP_TWO.out.data   // Named outputs
}
```

### Entry Workflow

- **Use `workflow { }` without a name** for the entry point.
- **Params should only be accessed** in the entry workflow and output blocks—pass values explicitly to subworkflows.
  ```groovy
  // WRONG - params deep in subworkflow
  workflow ALIGN {
      take: reads
      main:
      BWA_MEM(reads, params.reference)  // Avoid!
  }

  // CORRECT - pass params explicitly
  workflow ALIGN {
      take:
      reads
      reference   // Accept as input

      main:
      BWA_MEM(reads, reference)
  }

  workflow {
      ALIGN(reads_ch, params.reference)  // Pass from entry
  }
  ```

### Workflow Handlers (25.10+)

- **Define handlers inside workflows**, not as top-level config.
  ```groovy
  workflow {
      main:
      // pipeline logic

      onComplete {
          println "Pipeline completed at: ${workflow.complete}"
      }

      onError {
          println "Error: ${workflow.errorMessage}"
      }
  }
  ```

### Type Annotations for Workflows (25.10+)

```groovy
workflow RNASEQ {
    take:
    read_pairs_ch: Channel<Tuple<String, Path, Path>>
    transcriptome: Path

    main:
    // ...

    emit:
    fastqc: Channel<Tuple<String, Path>> = fastqc_ch
    quant: Channel<Tuple<String, Path>> = quant_ch
}
```

## 6. Channel Operations

Channels are the communication mechanism between processes. Use operators to transform data flowing through channels.

### Factory Methods

- **Use lowercase `channel`:** The capitalized `Channel` is deprecated.
  ```groovy
  // DEPRECATED
  Channel.of(1, 2, 3)
  Channel.fromPath('*.fq')

  // CORRECT
  channel.of(1, 2, 3)
  channel.fromPath('*.fq')
  ```

### Closure Parameters

- **Always declare closure parameters explicitly:** Never rely on implicit `it`.
  ```groovy
  // WRONG - implicit 'it'
  ch.map { it * 2 }
  ch.filter { it.size() > 100 }

  // CORRECT - explicit parameter
  ch.map { v -> v * 2 }
  ch.filter { file -> file.size() > 100 }

  // CORRECT - destructuring for tuples
  ch.map { meta, reads -> [meta.id, reads] }
  ```

### Avoid Deprecated Operators

- **No pipe operators:** Use explicit method calls and assignments.
  ```groovy
  // WRONG - pipe syntax
  channel.of('Hello', 'World')
      | greet
      | map { v -> v.toUpperCase() }

  // CORRECT - explicit calls
  ch_input = channel.of('Hello', 'World')
  ch_greet = greet(ch_input)
  ch_upper = ch_greet.map { v -> v.toUpperCase() }
  ```

- **No `set` or `tap` operators:** Use standard assignments.
  ```groovy
  // WRONG
  PROCESS(input).out.set { result_ch }

  // CORRECT
  result_ch = PROCESS(input).out
  ```

### No Spread Operator

- **Enumerate elements explicitly:** The `*` spread operator is not allowed.
  ```groovy
  // WRONG
  ch.map { meta, bambai -> [meta, *bambai] }

  // CORRECT
  ch.map { meta, bambai -> [meta, bambai[0], bambai[1]] }
  ```

### Higher-Order Functions Instead of Loops

- **Use `each` method instead of for/while loops:**
  ```groovy
  // WRONG - for loop
  for (module in ['fastqc', 'multiqc']) {
      // ...
  }

  // CORRECT - each method
  ['fastqc', 'multiqc'].each { module ->
      // ...
  }
  ```

## 7. Configuration

Configuration controls execution environment, resources, and pipeline behavior.

### File Structure

- `nextflow.config` - Main configuration file
- `conf/base.config` - Default process resources
- `conf/modules.config` - Per-process publish and container settings
- `conf/<profile>.config` - Environment-specific settings

### Syntax

```groovy
// Dot notation
params.input = 'data/'
process.cpus = 4

// Block notation (preferred for grouping)
params {
    input = 'data/'
    outdir = 'results/'
}

process {
    cpus = 4
    memory = '8.GB'
}
```

### Profiles

- **Use profiles for different environments:**
  ```groovy
  profiles {
      standard {
          process.executor = 'local'
      }

      slurm {
          process.executor = 'slurm'
          process.queue = 'batch'
      }

      docker {
          docker.enabled = true
      }

      singularity {
          singularity.enabled = true
          singularity.autoMounts = true
      }
  }
  ```

### Process Selectors

- **Use `withName` and `withLabel` for process-specific config:**
  ```groovy
  process {
      withLabel: 'process_high' {
          cpus = 16
          memory = '64.GB'
      }

      withName: 'FASTQC' {
          cpus = 2
          memory = '4.GB'
      }

      withName: '.*:SUBWORKFLOW:PROCESS' {
          // Selector with workflow context
      }
  }
  ```

### Environment Variables

- **Use `System.getenv()` or `env()` function:** Implicit env var access is not allowed.
  ```groovy
  // WRONG - implicit access
  println "Home: $HOME"

  // CORRECT
  println "Home: ${System.getenv('HOME')}"
  println "Home: ${env('HOME')}"
  ```

### Dynamic Configuration

- **Use dynamic includes for conditional config:**
  ```groovy
  includeConfig ({
      def hostname = 'hostname'.execute().text.trim()
      if (hostname.startsWith('hpc'))
          return 'conf/hpc.config'
      else
          return 'conf/local.config'
  }())
  ```

### Parameters Block (25.10+)

- **Declare typed parameters:**
  ```groovy
  params {
      input: Path
      outdir: Path = 'results'
      skip_qc: Boolean = false
      threads: Integer = 4
  }
  ```

## 8. Output Publishing

Nextflow 25.10 introduces workflow outputs as the preferred method for publishing results, replacing process-level `publishDir`.

### Why Workflow Outputs?

- **Centralized control:** All output rules in one place
- **Channel-aware:** Access metadata for dynamic paths
- **Index files:** Auto-generate samplesheets of outputs
- **No module pollution:** Processes remain reusable without hardcoded paths

### Basic Syntax

```groovy
workflow {
    main:
    FASTQC(reads_ch)
    ALIGN(reads_ch, reference)

    publish:
    fastqc_results = FASTQC.out.reports
    aligned_bams = ALIGN.out.bam
}

output {
    fastqc_results {
        path 'qc/fastqc'
    }

    aligned_bams {
        path 'aligned'
        mode 'copy'
    }
}
```

### Dynamic Paths

```groovy
output {
    samples {
        path { meta, files -> "results/${meta.id}" }
    }
}
```

### Index Files (Samplesheets)

```groovy
output {
    samples {
        path { sample -> "data/${sample.id}" }
        index {
            path 'samples.csv'
            header true
        }
    }
}
```

### Global Output Settings

```groovy
// In nextflow.config
params.outdir = 'results'

workflow.output {
    directory = params.outdir
    mode = 'copy'        // copy, symlink, link, move, rellink
    overwrite = true
}
```

### Legacy `publishDir` (Still Supported)

If using `publishDir` in processes, follow these patterns:

```groovy
process EXAMPLE {
    publishDir "${params.outdir}/subdir", mode: 'copy'

    // For multiple outputs to different locations:
    publishDir "${params.outdir}/logs", pattern: '*.log', mode: 'copy'
    publishDir "${params.outdir}/data", pattern: '*.bam', mode: 'copy'

    // ...
}
```

## 9. Anti-Patterns to Avoid

### Removed Syntax (Errors in Strict Mode)

| Pattern | Problem | Solution |
|---------|---------|----------|
| `import x.y.Z` | Import statements | Use fully qualified names |
| `class MyClass {}` | Class declarations | Move to `lib/` directory |
| `for (x in list)` | Loop statements | Use `.each { }` method |
| `switch (x) { }` | Switch statements | Use if-else chains |
| `[a, *list, b]` | Spread operator | Enumerate explicitly |
| `x = y = 1` | Chained assignment | Separate assignments |
| `foo(x = 1)` | Assignment in expression | Assign before call |

### Deprecated Syntax (Warnings Now, Errors Later)

| Pattern | Problem | Solution |
|---------|---------|----------|
| `Channel.of()` | Capitalized Channel | Use `channel.of()` |
| `{ it * 2 }` | Implicit `it` | Use `{ v -> v * 2 }` |
| `ch \| process` | Pipe operator | Use `process(ch)` |
| `ch.set { name }` | Set operator | Use `name = ch` |
| `shell` section | Deprecated | Use `script` section |
| `addParams` in include | Deprecated | Pass as workflow inputs |

### Structural Anti-Patterns

- **Params everywhere:** Only access `params` in entry workflow; pass values explicitly
- **Logic in processes:** Keep conditional logic in workflows, not process `when` blocks
- **Hardcoded paths:** Use `params` and `projectDir`/`launchDir` variables
- **Missing containers:** Always specify containers for reproducibility
- **Unnamed outputs:** Always use `emit:` for process outputs

## 10. Migration Checklist

When updating existing pipelines to strict syntax:

### Quick Checks

- [ ] Set `NXF_SYNTAX_PARSER=v2` and run `nextflow lint`
- [ ] Replace `Channel.` with `channel.`
- [ ] Add explicit closure parameters (no implicit `it`)
- [ ] Remove pipe operators (`|`, `&`)
- [ ] Remove `set` and `tap` operators
- [ ] Quote `env` input/output names

### Structural Changes

- [ ] Move class definitions to `lib/` directory
- [ ] Replace `import` with fully qualified names
- [ ] Convert loops to `.each { }` or channel operators
- [ ] Convert `switch` to if-else
- [ ] Remove spread operator usage
- [ ] Add `script:` label where required

### Modernization (25.10+)

- [ ] Add type annotations to params, inputs, outputs
- [ ] Convert `publishDir` to workflow `output` blocks
- [ ] Move workflow handlers inside workflow blocks
- [ ] Use typed `params { }` block

### Testing

- [ ] Run with `-stub` to test process stubs
- [ ] Verify outputs publish correctly
- [ ] Test all execution profiles

---

**BE CONSISTENT.** When editing existing code, match the established style. When starting new projects, follow this guide from the beginning.

*Sources:*
- [Nextflow Strict Syntax](https://nextflow.io/docs/latest/strict-syntax.html)
- [Static Types Tutorial](https://nextflow.io/docs/latest/tutorials/static-types.html)
- [Workflow Outputs Tutorial](https://nextflow.io/docs/latest/tutorials/workflow-outputs.html)
- [25.10 Migration Guide](https://nextflow.io/docs/latest/migrations/25-10.html)
- [Configuration Reference](https://nextflow.io/docs/latest/reference/config.html)
