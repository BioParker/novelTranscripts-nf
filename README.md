# Pipeline_novelTranscripts

Updated, cleaner version of Pipeline_novelTranscripts, using containers rather than a conda environment for ease of deployment on cloud services. Nextflow pipeline for identifying transcripts that contain novel cassette exons from IsoQuant outputs

## Parameters

All essential parameters can be set in config file.

### Key files

  * models: Path to transcript models gtf output from IsoQuant
  * model_tpms: Path to model transcript tpms output from IsoQuant
  * refgtf: Path to GENCODE reference gtf

### Paths

  * outdir: Directory in which to place output files

## Dockerfiles

Two named Dockerfiles are provided.

  * bior_docker: Dockerfile for generating container image for running the R scripts in workflow
  * suppatools: Dockerfile for generating container image required for running suppa2 and bedtools

