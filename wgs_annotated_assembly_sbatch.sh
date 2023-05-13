#!/bin/bash

#SBATCH --partition=synergy,cpu2019,cpu2021,cpu2022
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=7-00:00:00
#SBATCH --mem=5G
#SBATCH --error=wgs_annotated_assembly_run.%J.err
#SBATCH --output=wgs_annotated_assembly_run.%J.out

log_dir="$(pwd)"
log_file="logs/wgs_annotated_assembly.log.txt"

# The number of jobs for the snakemake command.
num_jobs=60

# The number of seconds to wait before checking if the output file of a snakemake rule is created.
latency_wait=15

# The number of times to restart a job if it fails.
restart_times=10

echo "started at: `date`"

# Load the ~/.bashrc file as source.
source ~/.bashrc

# Activate the snakemake conda environment.
conda activate snakemake

snakemake --cluster-config cluster.json --cluster 'sbatch --partition={cluster.partition} --cpus-per-task={cluster.cpus-per-task} --nodes={cluster.nodes} --ntasks={cluster.ntasks} --time={cluster.time} --mem={cluster.mem} --output={cluster.output} --error={cluster.error}' --jobs $num_jobs --latency-wait $latency_wait --restart-times $restart_times --rerun-incomplete --use-conda &> $log_dir/$log_file

output_dir=$(grep "output_dir" < config.yaml | grep -v "#" | cut -d ' ' -f2 | sed 's/"//g')
list_files=$(grep "list_files" < config.yaml | grep -v "#" | cut -d ' ' -f2 | sed 's/"//g')

snakemake_file_dir="${output_dir}/snakemake_files"
mkdir -p $snakemake_file_dir

cp $list_files $snakemake_file_dir

cp Snakefile $snakemake_file_dir
cp config.yaml $snakemake_file_dir
cp cluster.json $snakemake_file_dir
cp wgs_annotated_assembly_sbatch.sh $snakemake_file_dir 
cp wgs_annotated_assembly_run* $snakemake_file_dir

cp -rf logs $snakemake_file_dir
cp -rf utils $snakemake_file_dir

#python utils/scripts/parse_snakemake_command_logs.py --log_infile $log_dir/$log_file --output_dir $output_dir

echo "finished with exit code $? at: `date`"

