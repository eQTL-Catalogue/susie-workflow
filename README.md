# susie-workflow
Nextflow workflow for applying [Sum of Single Effects (SuSiE)](https://stephenslab.github.io/susieR/) finemapping model to a large number of molecular QTLs.

# Usage

```bash
nextflow run main.nf \
 --qtl_results Fairfax_2014.tsv\
 --cisdistance 500000\
 --n_batches 200\
 --eqtlutils '/gpfs/hpc/home/a72094/projects/eQTLUtils'\
 -resume\
 -profile finemapping\
 -executor.queueSize 1
```

