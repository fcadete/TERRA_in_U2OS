#!/bin/bash
#SBATCH --job-name=TTAGGG_reads-bowtie2
#SBATCH --time=4:00:00
#SBATCH --mem-per-cpu=2G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --image=docker:ummidock/bowtie2_samtools:latest
#SBATCH --workdir=/mnt/nfs/lobo/SALMEIDA-NFS/fcadete/scratch/TERRA_in_U2OS/

mkdir alignments

for f in extracted/*1_extracted.fastq; do
	filename=$(basename "${f%.*}")
	srun --exclusive -n1 -N1 shifter bash -c "echo 'Running $filename'; bowtie2 --very-sensitive -p32 -x reference/genome -1 $f -2 ${f//_1_/_2_}| samtools view -bSF4 -@ 32 - | samtools sort -@ 32 - > alignments/${filename//_1_/_paired_}_TTAGG_reads_hg38.bam" &
done
wait

for f in extracted/*1_extracted.fastq; do
	filename=$(basename "${f%.*}")
	srun --exclusive -n1 -N1 --cpus-per-task=1 --mem-per-cpu=4G shifter samtools index alignments/${filename//_1_/_paired_}_TTAGG_reads_hg38.bam
        srun --exclusive -n1 -N1 --cpus-per-task=1 --mem-per-cpu=4Gshifter --image=docker:argrosso/htstools:0.2.1 bedtools genomecov -g reference/genome.fa -ibam alignments/${filename//_1_/_paired_}_TTAGG_reads_hg38.bam -bg -split > alignments/${filename//_1_/_paired_}_TTAGG_reads_hg38_pileup.bedGraph
done

echo "Statistics for job $SLURM_JOB_ID:"
sacct --format="JOBID,NodeList,NNodes,Start,End,Elapsed,AllocCPUs,CPUTime,AveDiskRead,AveDiskWrite,MaxRSS,MaxVMSize,exitcode,derivedexitcode" -j $SLURM_JOB_ID

