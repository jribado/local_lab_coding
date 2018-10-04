~/miniconda3/bin/snakemake -p \
--jobs 100 \
--configfile config.yaml \
--cluster-config clusterconfig.yaml \
--restart-times 0 \
-s Snakefile \
--profile scg
