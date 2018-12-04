#!/bin/bash
#SBATCH --job-name=Rep_Count
#SBATCH --ntasks=24
#SBATCH --mem-per-cpu=1G
#SBATCH --cpus-per-task=1

mkdir repeat_counts

for file in raw_data/*fastq
do
   filename=$(basename "${file%.*}")
   srun -n1 -N1 --exclusive bash -c "grep -o -n 'TTAGGG\|CCCTAA' $file | cut -d : -f 1 | uniq -c | sed -e 's/  \+//' > repeat_counts/${filename}.repeat_count" &
done
wait

echo "Statistics for job $SLURM_JOB_ID:"
sacct --format="JOBID,NodeList,NNodes,Start,End,Elapsed,AllocCPUs,CPUTime,AveDiskRead,AveDiskWrite,MaxRSS,MaxVMSize,exitcode,derivedexitcode" -j $SLURM_JOB_ID

