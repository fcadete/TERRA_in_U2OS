#!/bin/bash
#SBATCH --job-name=TTAGGG_reads-bowtie2
#SBATCH --time=4:00:00
#SBATCH --mem-per-cpu=2G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --image=docker:ummidock/bowtie2_samtools:latest
#SBATCH --workdir=/mnt/nfs/lobo/SALMEIDA-NFS/fcadete/scratch/TERRA_in_U2OS/

mkdir alignments

for f in extracted_separated/*fastq; do
        filename=$(basename "${f%.*}")
        srun --exclusive -n1 -N1 shifter bash -c "echo 'Running $filename'; bowtie2 --very-sensitive -p32 -x reference/ConcatenatedFASTAAassemblies_hTel.txt -U $f | samtools view -bSF4 -@ 32 - | samtools sort -@ 32 - > alignments/${filename}_TTAGG_reads_telomeres.bam" &
done

wait

for f in extracted_separated/*fastq; do
        filename=$(basename "${f%.*}")
        srun --exclusive -n1 -N1 --cpus-per-task=1 --mem-per-cpu=4G shifter samtools index alignments/${filename}_TTAGG_reads_telomeres.bam &
        srun --exclusive -n1 -N1 --cpus-per-task=1 --mem-per-cpu=4G shifter --image=docker:argrosso/htstools:0.2.1 bedtools genomecov -g reference/ConcatenatedFASTAAassemblies_hTel.txt.fa -ibam alignments/${filename}_TTAGG_reads_telomeres.bam -bg > alignments/${filename}_TTAGG_reads_telomeres_pileup.bedGraph &
done

wait

echo "Statistics for job $SLURM_JOB_ID:"
sacct --format="JOBID,NodeList,NNodes,Start,End,Elapsed,AllocCPUs,CPUTime,AveDiskRead,AveDiskWrite,MaxRSS,MaxVMSize,exitcode,derivedexitcode" -j $SLURM_JOB_ID

