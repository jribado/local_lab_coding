snakemake \
  --jobs 100 \
  --configfile /srv/gsfs0/projects/bhatt/ribado/tools/qc_snakemake/config.json \
  --cluster "sbatch --account=asbhatt \
                    	  --job-name {params.job_name} \
                        -t {params.max_time} \
                        --mem={params.max_mem} \
                        -n {threads} \
                        -o {params.log_files}.log \
                        -e {params.log_files}.error"
