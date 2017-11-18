# dry run
snakemake -np

# Run
snakemake -p --latency-wait 999999

# Submit jobs
snakemake --cluster-config cluster.yml --cluster "sbatch --time {cluster.time} --mem {cluster.mem} --cpus-per-task {cluster.cpus}  -o logs/slurm-%A_%a.out --mail-user=arushiv@umich.edu --mail-type=FAIL --parsable" -j 60 -p --latency-wait 400

# print workflow
snakemake --dag | dot -Tsvg > workflow.svg
