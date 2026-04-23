# Metagenome_Analysis_SOP

This repository contains the Beman Lab metagenome processing pipeline from raw paired-end reads through assembly, binning, DAS Tool refinement, and downstream QC summaries.

## Default workflow

The default entrypoint is [`run_pipeline.sh`](./run_pipeline.sh). It wraps the existing stage scripts, centralizes configuration, runs preflight validation, captures run metadata, and generates a final aggregate report.

1. Edit [`config/pipeline.env`](./config/pipeline.env) for your run defaults.
2. Run a preflight check:

```bash
./run_pipeline.sh --preflight-only
```

3. Launch a run:

```bash
./run_pipeline.sh --run-dir /scratch/$USER/mag_runs/project_001
```

## Config

The shared config file lives at [`config/pipeline.env`](./config/pipeline.env).

Common settings:

- `RAW_READS_DIR`: input directory of paired-end raw reads
- `RQC_LABELING_PATTERN_FWD` / `RQC_LABELING_PATTERN_REV`: read suffix patterns
- `RUN_ROOT`: default output root
- `RUN_RQCFILTER`, `RUN_ASSEMBLY`, `RUN_MEGAHIT`, `RUN_BINNING`, `RUN_DAS_TOOL`, `RUN_BIN_ANALYSIS`, `RUN_DAS_BINS_QC`, `RUN_GTDBTK`: stage toggles
- `BBTOOLS_IMAGE`, `SPADES_IMAGE`, `MEGAHIT_IMAGE`, `METABAT2_IMAGE`: container paths
- `MAXBIN2_ENV`, `CONCOCT_ENV`, `DAS_TOOL_ENV`, `BUSCO_ENV`, `GTDBTK_ENV`: conda environments

For development or partial reruns, you can override the config per run:

```bash
./run_pipeline.sh --config /path/to/custom.env --from-step binning --to-step report
```

## Output layout

Each pipeline run writes a structured output tree under `RUN_ROOT`:

- `01_rqcfilter/`: raw RQCFilter outputs
- `02_filtered_reads/`: flattened `_interleaved_filtered.fastq.gz` links used by downstream steps
- `03_assemblies/`: SPAdes and MEGAHIT outputs plus flat FASTA links for binning
- `04_binning/`: binning outputs, DAS Tool outputs, and QC CSVs
- `05_gtdbtk/`: optional GTDB-Tk results
- `06_reports/`: aggregate Markdown and JSON run summaries
- `metadata/`: per-stage logs, resolved config, and run manifest

## Reports

At the end of a run, the wrapper generates:

- `06_reports/run_report.md`: human-readable run summary
- `06_reports/run_summary.json`: machine-readable summary of key artifacts and counts

You can regenerate just the final report from an existing run directory:

```bash
./run_pipeline.sh --run-dir /scratch/$USER/mag_runs/project_001 --report-only
```

## Testing

Run the repository test suite with:

```bash
./tests/run_tests.sh
```
