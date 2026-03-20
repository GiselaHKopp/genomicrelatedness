# nf-core/genomicrelatedness: Output

## Introduction

This document describes the output produced by the nfcore/genomicrelatedness pipeline. The output primarily consists of the following main components: output files (e.g. CRAM, BAM or VCF files), and summary statistics of the whole run presented in a [`MultiQC`](https://multiqc.info) report. Intermediate files and module-specific statistics files are also retained. 
The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory. The results directory of the pipeline needs to be specified with the `--outdir`flag when running the pipeline. During the run, intermediate files will be written to the work directory, which can be specified with the `-work-dir`parameter.

```text
{outdir}
├── bootstrapping
│   ├── round_1
│   │   ├── bqsr
│   │   |   ├── cram
|   |   |   |   └── merged
│   │   |   └── qc
|   |   ├── stats
│   │   └── variants
│   │       ├── db
│   │       ├── filtered
│   │       └── merged
|   |
│   ├── round_2
│   │   ├── …
|   ⋮    ⋮
|   ⋮
│   └── round_3
│       ├── …
|       ⋮
|
├── intervals
|
├── multiqc
|
├── pipeline_info
|
├── preprocessing
│   ├── alignment
|   |    ├── bam
|   |    |   └── bwamem2
│   |    └── cram
|   |
│   ├── coverage
│   ├── fastp
│   ├── preseq
│   └── stats
|
├── variant_calling
│   ├── bcftools
|   |    ├── merged
│   |    └── bam
|   |
│   └── gatk
|        ├── I<interval>_joint
|        ⋮
|        ├── merged
│        └── stats
|
└── relatedness_estimation
    ├── exclude
    |
    ├── intersection
    |
    ├── thinned
    |
    └── angsd_ngsrelate
work/
.nextflow.log
```

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

- [Preprocessing](#preprocessing)
  - [Prepare Reference Genome](#prepare-reference-genome)
  - [Prepare Intervals](#prepare-intervals)
  - [Prepare Input Files](#prepare-input-files)
  - [Map to Reference](#map-to-reference)
  - [Mark Duplicates](#mark-duplicates)
- [Bootstrapping](#bootstrapping)
  - [Call Variants](#call-variants)
  - [Hard Filter Variants](#hard-filter-variants)
  - [Base Quality Score Recalibration](#base-quality-score-recalibration)
- [Variant Calling](#variant-calling)
- [Relatedness Estimation](#relatedness-estimation)

## Preprocessing

### Prepare Reference Genome

### Prepare Intervals

### Prepare Input Files

### Map to Reference

### Mark Duplicates

## Bootstrapping

### Call Variants

### Hard Filter Variants

### Base Quality Score Recalibration

## Variant Calling
In the third section, **genotyping**, the recalibrated CRAMs are processed to generate genotype likelihoods and multi-sample VCFs. Two variant calling approaches, bcftools and GATK4, are chosen to mitigate caller-specific biases. The subworkflow *call_variants_gatk* is executed like in the previous section,  the subworkflow *call_variants_bcftools* converts the recalibrated CRAM with samtools `convert`, employs bcftools `mpileup` (--output-type z -d 100) and `call` (--output-type z -m -v –write-index=tbi) per scaffold, and concatenates results into a full cohort VCF with `concat` (--output-type z –write-index=tbi). The callsets from both subworkflows are accompanied by variant-level quality summaries from bcftools `stats`. This stage outputs two harmonized, multi-sample VCFs — one from GATK and one from bcftools. Both subworkflows run in parallel. Additionally, summary statistics are produced such as transition/transversion ratios that can be used to judge the quality of/improvement in the called variants and whether subsequent rounds of BQSR could be beneficial. Finally, in the subworkflow `intersect_variants`, the two callsets are intersected using bcftools `isec` to retain only sites called by both subworkflows and filtered using bcftools `exclude` (using parameter -include_scaffolds or -exclude_scaffolds), e.g. to exclude mitochondrial or gonosomal scaffolds or only include autosomal scaffolds, and `thin` (--remove-filtered-all –remove-indels –maf 0.025 –recode –recode-INFO-all –max-missing 0.75) producing the final variant set. First,  Optionally, mitochondrial and gonosomal scaffolds can be excluded with vcftools. The shared variant set is then filtered and thinned with VCFtools (Danecek et al., 2011)  (removing indels, sites with minor allele frequency below 0.025, and sites missing in more than 25% of samples; default values) to reduce linkage disequilibrium and mitigate ascertainment bias.

<details markdown="1">
<summary>Output files</summary>

- `variant_calling/`
  - `bcftools/`: directory containing the called variants as vcf and tbi files for each interval, the individual bam files per sample (in the subdirectory bam) and the merged vcf and tbi file (in the subdirectory merged).
  - `gatk/`: directory containing the called variants for each individuals as vcf and tbi files for each interval, the merged vcf and tbi file (in the subdirectory merged) and a text file with variant statistics (in the subfolder stats).
  
</details>

## Relatedness Estimation
The fourth section, **relatedness estimation**, uses the final variant set to produce robust estimates of pairwise relatedness suitable for low-coverage whole-genome sequencing data. It infers relatedness and inbreeding from genotype likelihood data.using NgsRelate v2 as implemented in ANGSD  (Korneliussen and Moltke, 2015; Hanghøj et al., 2019). The final output is a pairwise relatedness matrix.

<details markdown="1">
<summary>Output files</summary>

- `relatedness_estimation/`
  - `angsd_ngsrelate/`: directory containing the tab-delimited file with the pairwise relatedness matrix of the analyzed samples.
  - `excluded/`: directory containing the vcf file that was filtered to only included variants of the scaffolds defined to be included.
  - `include_contigs/`: directory containing bed file with the scaffolds to define which variants to use for relatedness estimation.
  - `intersection/`: directory containing the vcf files produced by bcftools `isec`, a README.txt to explain the content of the four different vcf files, and a text file with the sites.
  - `thinned/`: a directory containing the vcf file with the filtered variant data used as input for the relatedness estimation.

</details>

### MultiQC

<details markdown="1">
<summary>Output files</summary>

- `multiqc/`
  - `multiqc_report.html`: a standalone HTML file, which includes the most important overview statistics and can be viewed in your web browser.
  - `multiqc_data/`: directory containing parsed statistics from the different tools used in the pipeline.
  - `multiqc_plots/`: directory containing static images from the report in various formats.

</details>

[MultiQC](http://multiqc.info) is a visualization tool that generates a single HTML report summarising all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in the report data directory.

Results generated by MultiQC collate pipeline QC from supported tools e.g. FastP. The pipeline has special steps which also allow the software versions to be reported in the MultiQC output for future traceability. For more information about how to use MultiQC reports, see <http://multiqc.info>.

### Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  - Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.yml`. The `pipeline_report*` files will only be present if the `--email` / `--email_on_fail` parameter's are used when running the pipeline.
  - Reformatted samplesheet files used as input to the pipeline: `samplesheet.valid.csv`.
  - Parameters used by the pipeline run: `params.json`.

</details>

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.
