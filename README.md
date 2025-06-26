# WDLab
WDLab is a practical, example-driven crash course to learn the Workflow Description Language (WDL) for bioinformatics workflows development.
## Table of Contents
1. [Introduction to WDL](#1.-Introduction-to-WDL)  
2. [Installation with Bioconda](#installation-with-bioconda)  
3. [Cromwell as an Engine to Run WDL](#cromwell-as-an-engine-to-run-wdl)  
4. [WDL Syntax Structure and Main Components](#wdl-syntax-structure-and-main-components)  
   - 4.1 [Tasks](#tasks)  
   - 4.2 [Workflows](#workflows)  
   - 4.3 [Calls](#calls)  
   - 4.4 [Inputs, Outputs, and Declarations](#inputs-outputs-and-declarations)  
5. [Data Types in WDL](#data-types-in-wdl)  
   - 5.1 [Primitive Types](#primitive-types)  
   - 5.2 [Array and Map Types](#array-and-map-types)  
   - 5.3 [Optional Types](#optional-types)  
6. [Control Structures in WDL](#control-structures-in-wdl)  
   - 6.1 [Scatter](#scatter)  
   - 6.2 [If (Conditional)](#if-conditional)  
7. [`runtime` Settings](#`runtime`-Settings)  
8. [Inputs and Outputs in WDL](#inputs-and-outputs-in-wdl)  
   - 8.1 [Task-Level Inputs/Outputs](#task-level-inputsoutputs)  
   - 8.2 [Workflow-Level Inputs/Outputs](#workflow-level-inputsoutputs)  
9. [Output Directory Structure](#output-directory-structure)  
   - 9.1 [Execution Directory Layout](#execution-directory-layout)  
   - 9.2 [Logging Directory and Files](#logging-directory-and-files)  
10. [JSON Input File Structure](#json-input-file-structure)  
    - 10.1 [Formatting a Valid Inputs JSON](#formatting-a-valid-inputs-json)  
    - 10.2 [Example: Bioinformatics Inputs JSON](#example-bioinformatics-inputs-json)   
11. [Parallelization in WDL and Multi-Task Mode on HPCs](#parallelization-in-wdl-and-multi-task-mode-on-hpcs)  
    - 12.1 [Parallelism via Scatter Blocks](#parallelism-via-scatter-blocks)  
    - 12.2 [Running Cromwell on HPC Backends (SLURM/SGE)](#running-cromwell-on-hpc-backends-slurm-sge)  
12. [Additional Topics & Tips](#additional-topics--tips)  

---

## 1. Introduction to WDL

**Definition.**  
WDL (Workflow Description Language) is a human-readable, declarative language designed to specify data-processing workflows—particularly for bioinformatics and scientific computing. A WDL file describes *tasks* (individual steps, e.g., running BWA alignment, sorting BAMs, calling variants) and a *workflow* (the orchestration of those tasks, their inputs, outputs, and the order in which they run).  

**History.**  
- Originally developed within the Broad Institute (around 2014–2015) to standardize pipeline descriptions for genomic data analysis.  
- Open-sourced via the [OpenWDL](https://github.com/openwdl/wdl) community, with active participation from Broad, DNAnexus, and other organizations.  
- The current specification (WDL 1.0.0) is maintained on GitHub. Over the years, several engines (Cromwell, miniWDL, dxWDL, toil‐wdl-runner) have implemented the spec.  

**Why WDL?**  
- **Clarity**: Clear separation of tasks vs. workflow logic.  
- **Portability**: A single WDL file can be executed with any compliant engine on your laptop, HPC, or cloud.  
- **Reusability**: Define tasks once (e.g., `bwa_mem`), then call them in multiple workflows.  
- **Extensibility**: Ability to specify Docker or Singularity images per task, or embed conditional logic (“if” blocks, scatter loops).  
- **Community**: Rich ecosystem of example pipelines (GATK, RNA-seq, ChIP-seq), WDL “best practices,” and auxiliary tools (linting, input generation, static analysis).

---

## 2. Installation with Bioconda

WDL itself is a plain text file (no installation required), but you need:

1. **Java 8+** (for Cromwell) or **Python 3.6+** (for miniWDL, wdl-runner).  
2. **Cromwell JAR** (for running WDL).  
3. **`wdltool`** (for validation and input JSON generation).  

Below are steps to install via [Bioconda](https://bioconda.github.io/):

1. **Ensure Conda and Bioconda Channels Are Configured**  
   ```
   # Create a new conda environment (recommended)
   conda create -n wdl_env python=3.8 -y
   conda activate wdl_env

   # Configure channels (if not already set)
   conda config --add channels defaults
   conda config --add channels bioconda
   conda config --add channels conda-forge

2. **Install `womtool`**
   ```bash
   conda install bioconda::womtool
   ```

   
4. **Verify Installation**

   ```bash
   # cromwell (should print usage)
   cromwell --help
   ```

---

## 3. Cromwell as an Engine to Run WDL

**What is Cromwell?**

* Cromwell is a workflow execution engine developed by the Broad Institute.
* Supports WDL (primary) and CWL (limited).
* Can run workflows locally (single-machine), on HPC clusters (SLURM, SGE, PBS), or in cloud environments (Google Cloud, AWS Batch, Azure).

### 3.1 Downloading and Running Cromwell
1. **Install Cromwell**
   Cromwell is a Java JAR you can download manually, but Bioconda also packages a wrapper:

   ```bash
   conda install -y cromwell
   ```

   * This installs a wrapper script named `cromwell` that invokes the latest stable JAR under the hood.

2. **Run a WDL with Cromwell**

   ```bash
   java -jar cromwell.jar run my_workflow.wdl --inputs my_inputs.json
   ```
   * This will execute locally, producing an “execution directory” in the same folder (`cromwell-executions/`).
 
## 4. WDL Syntax Structure and Main Components

A WDL script (.wdl) is composed of:

1. **Workflow Declaration**
2. **Task Declarations**
3. **Optional `import` Statements** (to reuse tasks from other WDL files)
4. **Metadata/Version Specification** (WDL 1.0 supports `version 1.0` at top)

### 4.1 Tasks

* A **task** defines a single command-line operation, including:

  * `command { ... }` block
  * `input` declarations
  * `output` declarations
  * `runtime` attributes (e.g., `docker`, `cpu`, `memory`)

#### Example: BWA-MEM Alignment Task

```wdl
task bwa_mem {
  input {
    File ref_fasta
    File read1_fastq
    File read2_fastq
    Int threads = 4
  }

  command {
    bwa mem -t ${threads} ${ref_fasta} ${read1_fastq} ${read2_fastq} > aligned.sam
  }

  output {
    File sam = "aligned.sam"
  }

  runtime {
    docker: "biocontainers/bwa:v0.7.17_cv1"
    cpu: threads
    memory: "8 GB"
  }
}
```

* **Explanation**:

  * `input` block: declares inputs to this task. Default `threads=4`.
  * `command`: runs BWA‐MEM on `ref_fasta` and FASTQ pairs, outputs `aligned.sam`.
  * `output`: uses a glob (literal in quotes) to capture the output file.
  * `runtime`: instructs Cromwell to pull the `bwa` Docker image, allocate `cpu` and `memory`.

### 4.2 Workflows

* A **workflow** ties together multiple tasks, specifying:

  * `input` declarations (global inputs).
  * Calls to tasks (with aliases, if desired).
  * `scatter`/`if` blocks for parallel/conditional execution.
  * Workflow‐level `output` declarations.

#### Example: Simple Alignment → Sorting Workflow

```wdl
workflow bwa_alignment_pipeline {
  input {
    File reference_fasta
    Array[File] fastq_pairs
    Int threads = 4
  }

  scatter (pair in fastq_pairs) {
    call bwa_mem {
      input:
        ref_fasta = reference_fasta,
        read1_fastq = pair[0],
        read2_fastq = pair[1],
        threads = threads
    }

    call samtools_sort {
      input:
        sam_file = bwa_mem.sam
    }
  }

  output {
    Array[File] sorted_bams = samtools_sort.sorted_bam
  }
}
```

* **Explanation**:

  * `scatter (pair in fastq_pairs)`: parallelizes alignment for each FASTQ pair.
  * Two calls inside scatter:

    1. `bwa_mem` (produces `aligned.sam`)
    2. `samtools_sort` (sorts `aligned.sam` to produce `aligned.sorted.bam`)
  * Workflow output aggregates an array of sorted BAMs.

### 4.3 Calls

* A **call** invokes a task or sub-workflow.
* You can rename outputs by specifying `output:` in the call block.
* Calls can be conditional (inside an `if` block) or parallel (inside a `scatter`).

#### Syntax:

```wdl
call task_name {
  input:
    var1 = expr1,
    var2 = expr2
}
```

### 4.4 Inputs, Outputs, and Declarations

* **`input { ... }`**:

  * In a `task`: declares parameters the task needs.
  * In a `workflow`: declares global inputs (to be filled by the JSON).

* **`output { ... }`**:

  * In a `task`: declares which files/values the task will emit.
  * In a `workflow`: declares the final outputs (can reference outputs of calls).

* **`runtime { ... }`**:

  * Declares resource requirements (CPU, memory, Docker image, time, etc.).
  * Interpreted by Cromwell or another engine to dispatch jobs accordingly.

---

## 5. Data Types in WDL

WDL supports a variety of data types that can be used in both tasks and workflows.

### 5.1 Primitive Types

* **`Int`**: 32-bit integer.
* **`Float`**: Double‐precision.
* **`String`**: Unicode text. Example: `"sample1"`.
* **`Boolean`**: `true` or `false`.
* **`File`**: Path to a file on disk (or on a cloud bucket).

#### Example:

```wdl
task example_primitives {
  input {
    Int count
    Float threshold
    String sample_name
    Boolean do_filter
    File data_table
  }
  command {
    echo "Sample: ${sample_name}"
    if [ ${do_filter} = true ]; then
      python filter_by_threshold.py --input ${data_table} --threshold ${threshold} --n ${count} > filtered.txt
    else
      cp ${data_table} filtered.txt
    fi
  }
  output {
    File out_table = "filtered.txt"
  }
  runtime {
    cpu: 1
    memory: "2 GB"
    docker: "python:3.8"
  }
}
```

### 5.2 Array and Map Types

* **`Array[T]`**: Ordered list of zero or more elements of type `T`.

  * Example: `Array[File] fastq_files = [ "R1.fastq", "R2.fastq", "R3.fastq" ]`.
* **`Map[K, V]`**: A key–value mapping.

  * Example: `Map[String, Int] sample_to_reps = { "SampleA": 3, "SampleB": 2 }`.

#### Usage Example (Array):

```wdl
workflow count_fastq_reads {
  input {
    Array[File] all_fastqs
  }

  scatter (fq in all_fastqs) {
    call count_reads {
      input:
        fastq = fq
    }
  }

  output {
    Array[Int] read_counts = count_reads.read_count
  }
}

task count_reads {
  input {
    File fastq
  }
  command {
    echo $(zcat ${fastq} | echo $((`wc -l`/4))) > ${basename(fastq)}.count.txt
  }
  output {
    Int read_count = read_int("${basename(fastq)}.count.txt")
  }
  runtime {
    cpu: 1
    memory: "1 GB"
    docker: "biocontainers/samtools:v1.9-4-deb_cv1"
  }
}
```

* Note: `read_int(file)` is a WDL standard library function that reads an integer from a file.

#### Usage Example (Map):

```wdl
task print_sample_info {
  input {
    Map[String, String] sample_to_fastq
  }
  command {
    printf "Sample\tFastq\n" > sample_info.tsv
    %{ for (sample, fq) in sample_to_fastq ~}
    printf "${sample}\t${fq}\n" >> sample_info.tsv
    %{ endfor ~}
  }
  output {
    File info = "sample_info.tsv"
  }
  runtime {
    cpu: 1
    memory: "1 GB"
    docker: "ubuntu:20.04"
  }
}
```

### 5.3 Optional Types

* Use `?` to denote an optional type. Example: `File? config_yaml`.
* If a user does not supply the input, the value will be `null`. Tasks can check for null with `if ... != null`.

#### Example:

```wdl
task optional_config {
  input {
    File? config_yaml
    File data
  }
  command {
    if [[ ${config_yaml} != "null" ]]; then
      python process_with_config.py --data ${data} --config ${config_yaml}
    else
      python process_without_config.py --data ${data}
    fi
  }
  output {
    File result = "output.txt"
  }
  runtime {
    cpu: 1
    memory: "2 GB"
    docker: "python:3.8"
  }
}
```

---

## 6. Control Structures in WDL

### 6.1 Scatter

* The `scatter` block parallelizes calls over an array.
* Within a `scatter`, each iteration creates separate task calls that can run concurrently (subject to resource limits).

#### Syntax:

```wdl
scatter (item in my_array) {
  call some_task {
    input:
      param = item
  }
}
```

#### Bioinformatics Example: Trim and Align Many Samples

```wdl
workflow trim_and_align {
  input {
    Array[File] fastq_R1_list
    Array[File] fastq_R2_list
    File reference_fasta
    Int threads = 4
  }

  scatter (idx in range(length(fastq_R1_list))) {
    call fastp_trim {
      input:
        read1 = fastq_R1_list[idx],
        read2 = fastq_R2_list[idx]
    }

    call bwa_mem {
      input:
        ref_fasta = reference_fasta,
        read1_fastq = fastp_trim.trimmed_R1,
        read2_fastq = fastp_trim.trimmed_R2,
        threads = threads
    }
  }

  output {
    Array[File] all_bams = bwa_mem.sam
  }
}

task fastp_trim {
  input {
    File read1
    File read2
  }
  command <<<
    fastp --in1 ${read1} --in2 ${read2} \
      --out1 trimmed_R1.fastq --out2 trimmed_R2.fastq \
      --thread 2
  >>>
  output {
    File trimmed_R1 = "trimmed_R1.fastq"
    File trimmed_R2 = "trimmed_R2.fastq"
  }
  runtime {
    docker: "quay.io/biocontainers/fastp:0.20.1--hdfd78af_0"
    cpu: 2
    memory: "4 GB"
  }
}

task bwa_mem {
  input {
    File ref_fasta
    File read1_fastq
    File read2_fastq
    Int threads = 4
  }
  command {
    bwa mem -t ${threads} ${ref_fasta} ${read1_fastq} ${read2_fastq} > aligned.sam
  }
  output {
    File sam = "aligned.sam"
  }
  runtime {
    docker: "biocontainers/bwa:v0.7.17_cv1"
    cpu: threads
    memory: "8 GB"
  }
}
```

### 6.2 If (Conditional)

* The `if` block conditionally executes calls or declarations.
* Condition must be a Boolean expression.

#### Syntax:

```wdl
if (some_condition) {
  call my_task {
    input:
      ...
  }
}
```

#### Example: Optional Variant Annotation

```wdl
workflow annotate_variants {
  input {
    File vcf_file
    Boolean do_annotate = true
  }

  call filter_variants {
    input:
      vcf = vcf_file
  }

  if (do_annotate) {
    call annotate_with_snpeff {
      input:
        vcf_in = filter_variants.filtered_vcf
    }
  }

  output {
    File final_vcf = if (do_annotate) then annotate_with_snpeff.annotated_vcf else filter_variants.filtered_vcf
  }
}

task filter_variants {
  input {
    File vcf
  }
  command {
    bcftools view -f PASS ${vcf} -Oz -o filtered.vcf.gz
  }
  output {
    File filtered_vcf = "filtered.vcf.gz"
  }
  runtime {
    docker: "quay.io/biocontainers/bcftools:1.9--h9ee0642_2"
    cpu: 1
    memory: "2 GB"
  }
}

task annotate_with_snpeff {
  input {
    File vcf_in
  }
  command {
    snpeff ann All ${vcf_in} > annotated.vcf
  }
  output {
    File annotated_vcf = "annotated.vcf"
  }
  runtime {
    docker: "biocontainers/snpeff:4.3.1t_cv2"
    cpu: 1
    memory: "4 GB"
  }
}
```

---

## 7. `runtime` Settings
* The `runtime` block in each task indicates to the execution engine how to allocate resources:

  * **`cpu`**: Number of CPU cores.
  * **`memory`**: Memory requirement (e.g., `"8 GB"`).
  * **`docker`** or **`singularity`**: Container image.
  * **`disks`** (optional): Disk allocation (e.g., `"local-disk 100 HDD"`).
  * **`time`** (optional, Cromwell‐only): Runtime limit (e.g., `"2 hours"`).

#### Example:

```wdl
task variant_call {
  input {
    File bam
    File ref_fasta
  }
  command {
    gatk --java-options "-Xmx${memory}" HaplotypeCaller \
      -R ${ref_fasta} \
      -I ${bam} \
      -O raw_variants.vcf
  }
  output {
    File raw_vcf = "raw_variants.vcf"
  }
  runtime {
    docker: "broadinstitute/gatk:4.1.4.1"
    cpu: 4
    memory: "16 GB"
    time: "2 hours"
  }
}
```

## 8. Inputs and Outputs in WDL

### 8.1 Task-Level Inputs/Outputs

* **`input { ... }`**: each declaration must have a type (and optionally a default).
* **`output { ... }`**: each line declares a name, an expression (often a file literal or a built-in function).

#### Example: Task Inputs/Outputs for Quality Control

```wdl
task fastqc {
  input {
    File fastq
  }

  command {
    fastqc ${fastq} --outdir .
  }

  output {
    File report_html = glob("*.html")
    File report_zip  = glob("*.zip")
  }

  runtime {
    docker: "biocontainers/fastqc:v0.11.9_cv8"
    cpu: 1
    memory: "2 GB"
  }
}
```

* Outputs use `glob("pattern")` to capture one or more files matching the pattern.

### 8.2 Workflow-Level Inputs/Outputs

* **Workflow `input { ... }`**: declares parameters that must be provided in the JSON.
* **Workflow `output { ... }`**: can reference any call’s outputs or intermediate expressions.

#### Example:

```wdl
workflow qc_and_align {
  input {
    Array[File] fastqs
    File reference_fasta
  }

  scatter (fq in fastqs) {
    call fastqc {
      input:
        fastq = fq
    }

    call bwa_mem {
      input:
        ref_fasta = reference_fasta
        read1_fastq = fq
        read2_fastq = fq  # if single-end, use same file twice
        threads = 2
    }
  }

  output {
    Array[File] html_reports = fastqc.report_html
    Array[File] aligned_sams   = bwa_mem.sam
  }
}
```

---

## 9. Output Directory Structure

When Cromwell (or another engine) runs a workflow, it organizes outputs and logs into a hierarchical directory:

### 9.1 Execution Directory Layout

By default, Cromwell creates a top‐level directory named `cromwell-executions/` alongside your WDL. Inside:

```
cromwell-executions/
└── <workflow_name>_<uuid>/          # Unique folder per workflow run
    ├── call-<task_name_1>/           # Each task call gets its own folder
    │   ├── script                 # Generated shell script
    │   ├── stdout                 # Captured standard output
    │   ├── stderr                 # Captured standard error
    │   ├── return_code            # Exit code of the command
    │   ├── <outputs>              # Any declared `output {}` files
    │   └── execution_metadata.json
    ├── call-<task_name_2>/           # Another task call
    │   └── ...
    ├── call-<subworkflow_call>/      # If sub‐workflow invoked
    │   └── ...
    └── workflow.log               # High‐level workflow log (optional)
```

* `<workflow_name>_<uuid>`: e.g., `qc_and_align_d71a0a9b-fefa-459a-b79f-0e09d422eb14`.
* Each `call-<task>` directory contains:

  * `script`: the `.sh` file that Cromwell generated to run your `command {}`.
  * `stdout`/`stderr`: logs from that particular invocation.
  * `return_code`: integer exit status.
  * Declared outputs, e.g., `aligned.sam`, `report.html`, etc.
  * `execution_metadata.json`: metadata about inputs, instance type, docker image, runtime values, timestamps, etc.

### 9.2 Logging Directory and Files

* **`workflow.log`** (if enabled via options): a consolidated log for the entire workflow.
* **Per‐task `stdout`/`stderr`** files: useful for debugging command failures.
* **`execution_metadata.json`**: includes full metadata for that call (start time, end time, inputs, etc.).

Example:

```
cromwell-executions/qc_and_align_d71a0a9b-fefa-459a-b79f-0e09d422eb14/
└── call-bwa_mem/
    ├── script
    ├── stdout
    ├── stderr
    ├── return_code
    ├── aligned.sam
    └── execution_metadata.json
```

You can archive or rename these directories post‐run to keep records.

---

## 10. JSON Input File Structure

To run a WDL, you must supply a JSON (or YAML) file specifying all required workflow inputs (as declared in `workflow input { ... }`).

### 10.1 Formatting a Valid Inputs JSON

* Every key in the JSON must follow the pattern:

  ```
  "<workflow_name>.<input_name>"
  ```
* Values must match the declared type (e.g., strings for File paths, arrays for `Array[File]`, etc.).
* Relative paths are usually interpreted relative to the directory from which you invoke Cromwell.

#### General Template:

```jsonc
{
  "my_workflow.reference_fasta": "/path/to/genome.fa",
  "my_workflow.fastq_R1_list": ["samples/sampleA_R1.fastq.gz", "samples/sampleB_R1.fastq.gz"],
  "my_workflow.fastq_R2_list": ["samples/sampleA_R2.fastq.gz", "samples/sampleB_R2.fastq.gz"],
  "my_workflow.threads": 4
}
```

### 10.2 Example: Bioinformatics Inputs JSON

Assume a workflow named `variant_calling_pipeline`:

```wdl
workflow variant_calling_pipeline {
  input {
    File reference_fasta
    Array[File] bam_files
    String output_prefix = "results"
    Int threads = 8
  }
  ...
}
```

A corresponding `inputs.json`:

```json
{
  "variant_calling_pipeline.reference_fasta": "/data/genomes/GRCh38.fa",
  "variant_calling_pipeline.bam_files": [
    "/data/bams/sample1.bam",
    "/data/bams/sample2.bam",
    "/data/bams/sample3.bam"
  ],
  "variant_calling_pipeline.output_prefix": "projectX_calls",
  "variant_calling_pipeline.threads": 8
}
```

* Note that `output_prefix` has a default, so you could omit it:

  ```json
  {
    "variant_calling_pipeline.reference_fasta": "/data/genomes/GRCh38.fa",
    "variant_calling_pipeline.bam_files": [
      "/data/bams/sample1.bam",
      "/data/bams/sample2.bam"
    ]
  }
  ```

---

## 11. Parallelization in WDL and Multi-Task Mode on HPCs

### 11.1 Parallelism via Scatter Blocks

* The simplest form of parallelism in WDL is a `scatter` over an `Array[T]`.
* Each iteration spins up an independent call to the task(s) inside the block.

#### Example: Depth-of-Coverage Calculation for Many BAMs

```wdl
workflow coverage_scatter {
  input {
    Array[File] bams
    File reference_fasta
  }

  scatter (bam_file in bams) {
    call compute_coverage {
      input:
        bam = bam_file,
        ref = reference_fasta
    }
  }

  output {
    Array[File] coverage_reports = compute_coverage.report
  }
}

task compute_coverage {
  input {
    File bam
    File ref
  }
  command {
    bedtools genomecov -ibam ${bam} -g ${ref}.fai > coverage.txt
  }
  output {
    File report = "coverage.txt"
  }
  runtime {
    docker: "quay.io/biocontainers/bedtools:2.29.2--he860b03_6"
    cpu: 1
    memory: "2 GB"
  }
}
```

* Cromwell will launch up to `workflow_max_concurrent` parallel jobs (default = number of CPU cores on local machine). Use an `options.json` (Cromwell) to adjust concurrency limits:

  ```json
  {
    "backend": {
      "default": {
        "max-concurrent-tasks": 100
      }
    }
  }
  ```

### 12.2 Running Cromwell on HPC Backends (SLURM/SGE)

To leverage an HPC scheduler:

1. **Prepare Backend Configuration**

   * Save a file named (e.g.) `slurm.conf`:

     ```hocon
     include required(classpath("application"))

     akka {
       loglevel = "INFO"
       stdout-loglevel = "INFO"
     }

     system {
       server {
         production = false
       }
     }

     backend {
       providers {
         SLURM {
           actor-factory = "cromwell.backend.impl.sfs.config.ConfigBackendLifecycleActorFactory"
           config {
             root = "/path/to/cromwell-executions"
             filesystems {
               local { }
             }
             run {
               cpu = default("1")
               memory = default("1 GB")
               time = default("01:00:00")
               queue = default("short")
             }
             submit = "sbatch --mem=${memory} --cpus-per-task=${cpu} --time=${time} --partition=${queue} --output=${cwd}/slurm-%j.out"
             kill = "scancel ${job_id}"
             check-alive = "squeue -j ${job_id}"
             job-id-regex = "Submitted batch job (\\d+)"
             queue-commands {
               query = "sinfo --format=\"%.100P %.100c %.100D\" -h"
               filter = "grep \"${queue}\""
             }
           }
         }
       }
       default = "SLURM"
     }
     ```
   * Adjust `root`, `queue`, default `cpu`, `memory`, etc., for your HPC.

2. **Launch Cromwell with SLURM Backend**

   ```bash
   java -Dconfig.file=slurm.conf -jar cromwell.jar run my_workflow.wdl --inputs inputs.json
   ```

   * On each `call`, Cromwell will submit a job via `sbatch`.
   * You will see HPC job IDs printed in Cromwell logs, and each Cromwell `script` file will exist in the respective `call-<task>` directory.

3. **Multi‐Task Mode**

   * Cromwell will queue tasks to SLURM as long as you haven’t exceeded quotas.
   * To run multiple workflows concurrently, launch multiple Cromwell processes (or use a single Cromwell server with `server.play.allow-workflow-launches = true` and send HTTP API calls).
   * **Note**: If your HPC restricts number of concurrent submissions, coordinate with cluster admins or set `"max-concurrent-job-submissions": <N>` in Cromwell’s `options.json`.

---

## 13. Additional Topics & Tips

### 13.1 Importing and Reusing WDL Files

* You can `import "common_tasks.wdl"` at the top of your WDL to reuse tasks defined elsewhere:

  ```wdl
  version 1.0

  import "common_tasks.wdl" as Common

  workflow my_workflow {
    input { ... }

    call Common.bwa_mem {
      input:
        ref_fasta = reference_fasta,
        read1_fastq = fastq_pairs[0],
        read2_fastq = fastq_pairs[1]
    }

    ...
  }
  ```

### 13.2 Aliasing Calls

* You can rename a call to avoid name collisions or to indicate purpose:

  ```wdl
  call bwa_mem as align_sampleA {
    input: ...
  }
  call bwa_mem as align_sampleB {
    input: ...
  }
  ```

### 13.3 Task Placement and Default Values

* Variables declared in one scope are visible only in that scope or downstream.
* You may specify default values in `input` declarations. If omitted in the JSON, WDL uses the default.

### 13.4 Common WDL Best Practices

* **Keep tasks small and modular**: One well-defined command per task.
* **Pin Docker versions**: Use a fully qualified container tag (e.g., `biocontainers/bwa:v0.7.17_cv1`).
* **Include `metadata` and logging**: Always look at `execution_metadata.json` to debug.
* **Use `wdltool` linting**: Catch common mistakes early.
* **Document inputs/outputs**: Add comments above task inputs and workflow inputs, explaining purpose.

### 13.5 Debugging Tips

* **Local Dry Run** (miniWDL):

  ```bash
  miniwdl run --sleep 0 my_workflow.wdl my_inputs.json --backend local
  ```

  Prints planned calls without executing.
* **Verbose Cromwell Logs**:

  ```bash
  java -Dconfig.file=slurm.conf -Dlogback.resource=cromwell-logback-debug.xml \
       -jar cromwell.jar run my_workflow.wdl --inputs inputs.json
  ```

  Use a debug `logback` configuration to capture detailed logs for troubleshooting.

### 13.6 WDL Style Guidelines for Bioinformatics

1. **Naming Conventions**:

   * Use lowercase and underscores for tasks/workflows (e.g., `bwa_mem`, `variant_calling`).
   * Prefix workflow-level inputs with `wf_` if desired, to distinguish from task-level inputs.
2. **Resource Estimation**:

   * Run small test data to measure runtime, CPU, memory, then adjust `runtime` fields.
   * Over‐requesting resources may waste cluster quotas; under‐requesting may cause job failures.
3. **Containerization**:

   * Always reference well‐maintained Docker images (e.g., Biocontainers).
   * For proprietary software, build a custom Docker image and host it on Docker Hub or a private registry.
4. **Data Localization**:

   * If using cloud storage (e.g., `gs://` or `s3://` paths), ensure Cromwell’s backend is configured to localize and delocalize automatically.
   * For HPC, provide a shared filesystem path (NFS, Lustre) that all nodes can access.

---

### Appendix: Complete Example WDL Pipeline for RNA-Seq

Below is a consolidated WDL pipeline to:

1. Trim adapters with `fastp`.
2. Align trimmed reads to a reference with `STAR`.
3. Sort and index BAM with `samtools`.
4. Quantify expression (e.g., `featureCounts`).
5. Merge results.

```wdl
version 1.0
# ========================
# Workflow: rna_seq_pipeline
# ========================
workflow rna_seq_pipeline {
  input {
    Array[File] fastq_R1_list
    Array[File] fastq_R2_list
    File star_index_dir
    File gtf
    Int threads_trim = 4
    Int threads_align = 8
    Int threads_fc = 4
  }

  scatter (i in range(length(fastq_R1_list))) {
    call fastp_trim {
      input:
        read1 = fastq_R1_list[i],
        read2 = fastq_R2_list[i]
    }

    call star_align {
      input:
        trimmed_R1 = fastp_trim.trimmed_R1,
        trimmed_R2 = fastp_trim.trimed_R2,
        star_index_dir = star_index_dir,
        threads = threads_align
    }

    call samtools_index {
      input:
        bam = star_align.aligned_bam
    }

    call featurecounts {
      input:
        bam = star_align.aligned_bam,
        gtf = gtf,
        threads = threads_fc
    }
  }

  output {
    Array[File] bam_files    = star_align.aligned_bam
    Array[File] bam_indexes = samtools_index.bai
    Array[File] count_tables = featurecounts.counts
  }
}

# ========================
# Task: fastp_trim
# ========================
task fastp_trim {
  input {
    File read1
    File read2
  }

  command <<<
    fastp \
      --in1 ${read1} --in2 ${read2} \
      --out1 trimmed_R1.fastq.gz --out2 trimmed_R2.fastq.gz \
      --thread 4 \
      --html fastp_report.html \
      --json fastp_report.json
  >>>

  output {
    File trimmed_R1 = "trimmed_R1.fastq.gz"
    File trimmed_R2 = "trimmed_R2.fastq.gz"
    File html_report = "fastp_report.html"
    File json_report = "fastp_report.json"
  }

  runtime {
    docker: "quay.io/biocontainers/fastp:0.23.1--hf7d754f_0"
    cpu: 4
    memory: "8 GB"
  }
}

# ========================
# Task: star_align
# ========================
task star_align {
  input {
    File trimmed_R1
    File trimmed_R2
    File star_index_dir  # Path to STAR index directory
    Int threads = 8
  }

  command <<<
    STAR \
      --runThreadN ${threads} \
      --genomeDir ${star_index_dir} \
      --readFilesIn ${trimmed_R1} ${trimmed_R2} \
      --readFilesCommand zcat \
      --outSAMtype BAM SortedByCoordinate \
      --outFileNamePrefix sample_
  >>>

  output {
    File aligned_bam = "sample_Aligned.sortedByCoord.out.bam"
    File log_final = "sample_Log.final.out"
  }

  runtime {
    docker: "quay.io/biocontainers/star:2.7.9a--0"
    cpu: threads
    memory: "32 GB"
  }
}

# ========================
# Task: samtools_index
# ========================
task samtools_index {
  input {
    File bam
  }

  command {
    samtools index ${bam}
  }

  output {
    File bai = "${bam}.bai"
  }

  runtime {
    docker: "quay.io/biocontainers/samtools:1.14--h2e538c0_0"
    cpu: 1
    memory: "2 GB"
  }
}

# ========================
# Task: featurecounts
# ========================
task featurecounts {
  input {
    File bam
    File gtf
    Int threads = 4
  }

  command <<<
    featureCounts \
      -T ${threads} \
      -a ${gtf} \
      -o counts.txt \
      ${bam}
  >>>

  output {
    File counts = "counts.txt"
  }

  runtime {
    docker: "quay.io/biocontainers/subread:2.0.1--h9ee0642_0"
    cpu: threads
    memory: "8 GB"
  }
}

```

#### Example `inputs.json` for RNA-Seq Pipeline

```json
{
  "rna_seq_pipeline.fastq_R1_list": [
    "/data/reads/sample1_R1.fastq.gz",
    "/data/reads/sample2_R1.fastq.gz"
  ],
  "rna_seq_pipeline.fastq_R2_list": [
    "/data/reads/sample1_R2.fastq.gz",
    "/data/reads/sample2_R2.fastq.gz"
  ],
  "rna_seq_pipeline.star_index_dir": "/data/indices/STAR_hg38/",
  "rna_seq_pipeline.gtf": "/data/annotations/hg38.gtf",
  "rna_seq_pipeline.threads_trim": 4,
  "rna_seq_pipeline.threads_align": 8,
  "rna_seq_pipeline.threads_fc": 4
}
```

---
