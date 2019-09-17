
nextflow run main.nf --qtl_results /gpfs/hpc/home/a72094/datasets/controlled_access/SampleArcheology/finemapping/Alasoo_2018.tsv\
 --cisdistance 1000000 \
 --n_batches 200\
 --eqtlutils '/gpfs/hpc/home/a72094/projects/eQTLUtils'\
 -resume\
 -profile finemapping\
 -executor.queueSize 1

