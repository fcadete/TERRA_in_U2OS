#!/bin/bash
#SBATCH --job-name=count_chromosome_hits
#SBATCH --time=4:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --image=docker:ummidock/bowtie2_samtools:latest
#SBATCH --workdir=/mnt/beegfs/scratch/SALMEIDA/fcadete/TERRA_in_U2OS/

rm chromosome_hits.tsv
touch chromosome_hits.tsv

for f in alignments/*bam; do
	filename=$(basename "${f%.*}")
	srun --exclusive -n1 -N1 shifter bash -c "samtools idxstats $f | sed -e 's/$/\t$filename/' >> chromosome_hits.tsv" 
done

echo "Statistics for job $SLURM_JOB_ID:"
sacct --format="JOBID,NodeList,NNodes,Start,End,Elapsed,AllocCPUs,CPUTime,AveDiskRead,AveDiskWrite,MaxRSS,MaxVMSize,exitcode,derivedexitcode" -j $SLURM_JOB_ID

