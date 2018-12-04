#!/bin/bash
#SBATCH --job-name=download_datasets
#SBATCH --time=6:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --ntasks=24
#SBATCH --cpus-per-task=1
#SBATCH --image=docker:argrosso/htspreprocessing:0.1.2
#SBATCH --workdir=/mnt/beegfs/scratch/SALMEIDA/fcadete/TERRA_in_U2OS/

mkdir raw_data

while read line; do
  array=($line)
  SRA_ID=${array[3]}
  srun --exclusive -n1 -N1 shifter fastq-dump --split-files -O raw_data $SRA_ID &
done < info/datasets.txt
wait

#There is a mixup of files in the databases, correct them here
mv raw_data/SRR7880262_2.fastq raw_data/temp_SRR7880261_2.fastq
mv raw_data/SRR7880261_2.fastq raw_data/SRR7880262_2.fastq
mv raw_data/temp_SRR7880261_2.fastq raw_data/SRR7880261_2.fastq

echo "Statistics for job $SLURM_JOB_ID:"
sacct --format="JOBID,Start,End,Elapsed,AllocCPUs,CPUTime,AveDiskRead,AveDiskWrite,MaxRSS,MaxVMSize,exitcode,derivedexitcode" -j $SLURM_JOB_ID

