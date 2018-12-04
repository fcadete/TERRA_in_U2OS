#!/bin/bash
#SBATCH --job-name=examine_counts
#SBATCH --ntasks=24
#SBATCH --mem-per-cpu=8G
#SBATCH --cpus-per-task=1
#SBATCH --image=docker:rocker/tidyverse:latest

echo "#!/bin/bash" > scripts/extract_sequences_to_map_subtelomeres.sbatch
echo "#SBATCH --job-name=Extract_repeats" >> scripts/extract_sequences_to_map_subtelomeres.sbatch
echo "#SBATCH --ntasks=12" >> scripts/extract_sequences_to_map_subtelomeres.sbatch
echo "#SBATCH --mem-per-cpu=1G" >> scripts/extract_sequences_to_map_subtelomeres.sbatch
echo "#SBATCH --cpus-per-task=1" >> scripts/extract_sequences_to_map_subtelomeres.sbatch
echo "mkdir extracted" >> scripts/extract_sequences_to_map_subtelomeres.sbatch
echo "mkdir extracted_separated" >> scripts/extract_sequences_to_map_subtelomeres.sbatch


for sample_ID in $(ls raw_data/*1.fastq | sed -e 's/_1.fastq//' | sed -e 's/raw_data\///'); do

   srun -n1 -N1 --exclusive shifter Rscript scripts/analyse_repeat_counts_parallel.R $sample_ID

done
wait

echo "wait" >> scripts/extract_sequences_to_map_subtelomeres.sbatch

echo "Statistics for job $SLURM_JOB_ID:"
sacct --format="JOBID,NodeList,NNodes,Start,End,Elapsed,AllocCPUs,CPUTime,AveDiskRead,AveDiskWrite,MaxRSS,MaxVMSize,exitcode,derivedexitcode" -j $SLURM_JOB_ID


