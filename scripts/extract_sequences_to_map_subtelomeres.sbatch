#!/bin/bash
#SBATCH --job-name=Extract_repeats
#SBATCH --ntasks=48
#SBATCH --mem-per-cpu=1G
#SBATCH --cpus-per-task=1
#SBATCH --time=36:00:00
mkdir extracted
mkdir extracted_separated
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_1_paired_extract.sed raw_data/SRR1246789_1.fastq > extracted/SRR1246789_1_1_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_1_paired_extract.sed raw_data/SRR1246789_2.fastq > extracted/SRR1246789_2_1_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_2_paired_extract.sed raw_data/SRR1246789_1.fastq > extracted/SRR1246789_1_2_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_2_paired_extract.sed raw_data/SRR1246789_2.fastq > extracted/SRR1246789_2_2_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_3_paired_extract.sed raw_data/SRR1246789_1.fastq > extracted/SRR1246789_1_3_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_3_paired_extract.sed raw_data/SRR1246789_2.fastq > extracted/SRR1246789_2_3_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_4_paired_extract.sed raw_data/SRR1246789_1.fastq > extracted/SRR1246789_1_4_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_4_paired_extract.sed raw_data/SRR1246789_2.fastq > extracted/SRR1246789_2_4_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_5_paired_extract.sed raw_data/SRR1246789_1.fastq > extracted/SRR1246789_1_5_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_5_paired_extract.sed raw_data/SRR1246789_2.fastq > extracted/SRR1246789_2_5_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_6_paired_extract.sed raw_data/SRR1246789_1.fastq > extracted/SRR1246789_1_6_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_6_paired_extract.sed raw_data/SRR1246789_2.fastq > extracted/SRR1246789_2_6_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_7_paired_extract.sed raw_data/SRR1246789_1.fastq > extracted/SRR1246789_1_7_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_7_paired_extract.sed raw_data/SRR1246789_2.fastq > extracted/SRR1246789_2_7_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_8_paired_extract.sed raw_data/SRR1246789_1.fastq > extracted/SRR1246789_1_8_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_8_paired_extract.sed raw_data/SRR1246789_2.fastq > extracted/SRR1246789_2_8_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_9_paired_extract.sed raw_data/SRR1246789_1.fastq > extracted/SRR1246789_1_9_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_9_paired_extract.sed raw_data/SRR1246789_2.fastq > extracted/SRR1246789_2_9_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_10_paired_extract.sed raw_data/SRR1246789_1.fastq > extracted/SRR1246789_1_10_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_10_paired_extract.sed raw_data/SRR1246789_2.fastq > extracted/SRR1246789_2_10_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_11_paired_extract.sed raw_data/SRR1246789_1.fastq > extracted/SRR1246789_1_11_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_11_paired_extract.sed raw_data/SRR1246789_2.fastq > extracted/SRR1246789_2_11_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_12_paired_extract.sed raw_data/SRR1246789_1.fastq > extracted/SRR1246789_1_12_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_12_paired_extract.sed raw_data/SRR1246789_2.fastq > extracted/SRR1246789_2_12_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_13_paired_extract.sed raw_data/SRR1246789_1.fastq > extracted/SRR1246789_1_13_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_13_paired_extract.sed raw_data/SRR1246789_2.fastq > extracted/SRR1246789_2_13_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_14_paired_extract.sed raw_data/SRR1246789_1.fastq > extracted/SRR1246789_1_14_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_14_paired_extract.sed raw_data/SRR1246789_2.fastq > extracted/SRR1246789_2_14_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_15_paired_extract.sed raw_data/SRR1246789_1.fastq > extracted/SRR1246789_1_15_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_15_paired_extract.sed raw_data/SRR1246789_2.fastq > extracted/SRR1246789_2_15_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_16_paired_extract.sed raw_data/SRR1246789_1.fastq > extracted/SRR1246789_1_16_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_16_paired_extract.sed raw_data/SRR1246789_2.fastq > extracted/SRR1246789_2_16_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_17_paired_extract.sed raw_data/SRR1246789_1.fastq > extracted/SRR1246789_1_17_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_17_paired_extract.sed raw_data/SRR1246789_2.fastq > extracted/SRR1246789_2_17_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_18_paired_extract.sed raw_data/SRR1246789_1.fastq > extracted/SRR1246789_1_18_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_18_paired_extract.sed raw_data/SRR1246789_2.fastq > extracted/SRR1246789_2_18_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_19_paired_extract.sed raw_data/SRR1246789_1.fastq > extracted/SRR1246789_1_19_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_19_paired_extract.sed raw_data/SRR1246789_2.fastq > extracted/SRR1246789_2_19_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_20_paired_extract.sed raw_data/SRR1246789_1.fastq > extracted/SRR1246789_1_20_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_20_paired_extract.sed raw_data/SRR1246789_2.fastq > extracted/SRR1246789_2_20_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_21_paired_extract.sed raw_data/SRR1246789_1.fastq > extracted/SRR1246789_1_21_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_21_paired_extract.sed raw_data/SRR1246789_2.fastq > extracted/SRR1246789_2_21_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_22_paired_extract.sed raw_data/SRR1246789_1.fastq > extracted/SRR1246789_1_22_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_22_paired_extract.sed raw_data/SRR1246789_2.fastq > extracted/SRR1246789_2_22_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_23_paired_extract.sed raw_data/SRR1246789_1.fastq > extracted/SRR1246789_1_23_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_23_paired_extract.sed raw_data/SRR1246789_2.fastq > extracted/SRR1246789_2_23_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_24_paired_extract.sed raw_data/SRR1246789_1.fastq > extracted/SRR1246789_1_24_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_24_paired_extract.sed raw_data/SRR1246789_2.fastq > extracted/SRR1246789_2_24_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_25_paired_extract.sed raw_data/SRR1246789_1.fastq > extracted/SRR1246789_1_25_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_25_paired_extract.sed raw_data/SRR1246789_2.fastq > extracted/SRR1246789_2_25_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_26_paired_extract.sed raw_data/SRR1246789_1.fastq > extracted/SRR1246789_1_26_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_26_paired_extract.sed raw_data/SRR1246789_2.fastq > extracted/SRR1246789_2_26_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_27_paired_extract.sed raw_data/SRR1246789_1.fastq > extracted/SRR1246789_1_27_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_27_paired_extract.sed raw_data/SRR1246789_2.fastq > extracted/SRR1246789_2_27_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_28_paired_extract.sed raw_data/SRR1246789_1.fastq > extracted/SRR1246789_1_28_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_28_paired_extract.sed raw_data/SRR1246789_2.fastq > extracted/SRR1246789_2_28_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_29_paired_extract.sed raw_data/SRR1246789_1.fastq > extracted/SRR1246789_1_29_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_29_paired_extract.sed raw_data/SRR1246789_2.fastq > extracted/SRR1246789_2_29_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_30_paired_extract.sed raw_data/SRR1246789_1.fastq > extracted/SRR1246789_1_30_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_30_paired_extract.sed raw_data/SRR1246789_2.fastq > extracted/SRR1246789_2_30_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_31_paired_extract.sed raw_data/SRR1246789_1.fastq > extracted/SRR1246789_1_31_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_31_paired_extract.sed raw_data/SRR1246789_2.fastq > extracted/SRR1246789_2_31_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_32_paired_extract.sed raw_data/SRR1246789_1.fastq > extracted/SRR1246789_1_32_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_32_paired_extract.sed raw_data/SRR1246789_2.fastq > extracted/SRR1246789_2_32_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_33_paired_extract.sed raw_data/SRR1246789_1.fastq > extracted/SRR1246789_1_33_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_33_paired_extract.sed raw_data/SRR1246789_2.fastq > extracted/SRR1246789_2_33_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_34_paired_extract.sed raw_data/SRR1246789_1.fastq > extracted/SRR1246789_1_34_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_34_paired_extract.sed raw_data/SRR1246789_2.fastq > extracted/SRR1246789_2_34_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_35_paired_extract.sed raw_data/SRR1246789_1.fastq > extracted/SRR1246789_1_35_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_35_paired_extract.sed raw_data/SRR1246789_2.fastq > extracted/SRR1246789_2_35_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_36_paired_extract.sed raw_data/SRR1246789_1.fastq > extracted/SRR1246789_1_36_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_36_paired_extract.sed raw_data/SRR1246789_2.fastq > extracted/SRR1246789_2_36_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_37_paired_extract.sed raw_data/SRR1246789_1.fastq > extracted/SRR1246789_1_37_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_37_paired_extract.sed raw_data/SRR1246789_2.fastq > extracted/SRR1246789_2_37_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_38_paired_extract.sed raw_data/SRR1246789_1.fastq > extracted/SRR1246789_1_38_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_38_paired_extract.sed raw_data/SRR1246789_2.fastq > extracted/SRR1246789_2_38_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_39_paired_extract.sed raw_data/SRR1246789_1.fastq > extracted/SRR1246789_1_39_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_39_paired_extract.sed raw_data/SRR1246789_2.fastq > extracted/SRR1246789_2_39_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_40_paired_extract.sed raw_data/SRR1246789_1.fastq > extracted/SRR1246789_1_40_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_40_paired_extract.sed raw_data/SRR1246789_2.fastq > extracted/SRR1246789_2_40_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_41_paired_extract.sed raw_data/SRR1246789_1.fastq > extracted/SRR1246789_1_41_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_41_paired_extract.sed raw_data/SRR1246789_2.fastq > extracted/SRR1246789_2_41_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_42_paired_extract.sed raw_data/SRR1246789_1.fastq > extracted/SRR1246789_1_42_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_42_paired_extract.sed raw_data/SRR1246789_2.fastq > extracted/SRR1246789_2_42_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_43_paired_extract.sed raw_data/SRR1246789_1.fastq > extracted/SRR1246789_1_43_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_43_paired_extract.sed raw_data/SRR1246789_2.fastq > extracted/SRR1246789_2_43_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_44_paired_extract.sed raw_data/SRR1246789_1.fastq > extracted/SRR1246789_1_44_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_44_paired_extract.sed raw_data/SRR1246789_2.fastq > extracted/SRR1246789_2_44_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_45_paired_extract.sed raw_data/SRR1246789_1.fastq > extracted/SRR1246789_1_45_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_45_paired_extract.sed raw_data/SRR1246789_2.fastq > extracted/SRR1246789_2_45_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_46_paired_extract.sed raw_data/SRR1246789_1.fastq > extracted/SRR1246789_1_46_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_46_paired_extract.sed raw_data/SRR1246789_2.fastq > extracted/SRR1246789_2_46_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_47_paired_extract.sed raw_data/SRR1246789_1.fastq > extracted/SRR1246789_1_47_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_47_paired_extract.sed raw_data/SRR1246789_2.fastq > extracted/SRR1246789_2_47_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_48_paired_extract.sed raw_data/SRR1246789_1.fastq > extracted/SRR1246789_1_48_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_48_paired_extract.sed raw_data/SRR1246789_2.fastq > extracted/SRR1246789_2_48_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_49_paired_extract.sed raw_data/SRR1246789_1.fastq > extracted/SRR1246789_1_49_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_49_paired_extract.sed raw_data/SRR1246789_2.fastq > extracted/SRR1246789_2_49_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_50_paired_extract.sed raw_data/SRR1246789_1.fastq > extracted/SRR1246789_1_50_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_50_paired_extract.sed raw_data/SRR1246789_2.fastq > extracted/SRR1246789_2_50_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_51_paired_extract.sed raw_data/SRR1246789_1.fastq > extracted/SRR1246789_1_51_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_51_paired_extract.sed raw_data/SRR1246789_2.fastq > extracted/SRR1246789_2_51_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_52_paired_extract.sed raw_data/SRR1246789_1.fastq > extracted/SRR1246789_1_52_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_52_paired_extract.sed raw_data/SRR1246789_2.fastq > extracted/SRR1246789_2_52_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_53_paired_extract.sed raw_data/SRR1246789_1.fastq > extracted/SRR1246789_1_53_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_53_paired_extract.sed raw_data/SRR1246789_2.fastq > extracted/SRR1246789_2_53_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_54_paired_extract.sed raw_data/SRR1246789_1.fastq > extracted/SRR1246789_1_54_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_54_paired_extract.sed raw_data/SRR1246789_2.fastq > extracted/SRR1246789_2_54_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_55_paired_extract.sed raw_data/SRR1246789_1.fastq > extracted/SRR1246789_1_55_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_55_paired_extract.sed raw_data/SRR1246789_2.fastq > extracted/SRR1246789_2_55_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_56_paired_extract.sed raw_data/SRR1246789_1.fastq > extracted/SRR1246789_1_56_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_56_paired_extract.sed raw_data/SRR1246789_2.fastq > extracted/SRR1246789_2_56_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_57_paired_extract.sed raw_data/SRR1246789_1.fastq > extracted/SRR1246789_1_57_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_57_paired_extract.sed raw_data/SRR1246789_2.fastq > extracted/SRR1246789_2_57_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_58_paired_extract.sed raw_data/SRR1246789_1.fastq > extracted/SRR1246789_1_58_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_58_paired_extract.sed raw_data/SRR1246789_2.fastq > extracted/SRR1246789_2_58_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_59_paired_extract.sed raw_data/SRR1246789_1.fastq > extracted/SRR1246789_1_59_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_59_paired_extract.sed raw_data/SRR1246789_2.fastq > extracted/SRR1246789_2_59_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_60_paired_extract.sed raw_data/SRR1246789_1.fastq > extracted/SRR1246789_1_60_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_60_paired_extract.sed raw_data/SRR1246789_2.fastq > extracted/SRR1246789_2_60_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_61_paired_extract.sed raw_data/SRR1246789_1.fastq > extracted/SRR1246789_1_61_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_61_paired_extract.sed raw_data/SRR1246789_2.fastq > extracted/SRR1246789_2_61_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_62_paired_extract.sed raw_data/SRR1246789_1.fastq > extracted/SRR1246789_1_62_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_62_paired_extract.sed raw_data/SRR1246789_2.fastq > extracted/SRR1246789_2_62_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_63_paired_extract.sed raw_data/SRR1246789_1.fastq > extracted/SRR1246789_1_63_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_63_paired_extract.sed raw_data/SRR1246789_2.fastq > extracted/SRR1246789_2_63_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_64_paired_extract.sed raw_data/SRR1246789_1.fastq > extracted/SRR1246789_1_64_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_64_paired_extract.sed raw_data/SRR1246789_2.fastq > extracted/SRR1246789_2_64_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_65_paired_extract.sed raw_data/SRR1246789_1.fastq > extracted/SRR1246789_1_65_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_65_paired_extract.sed raw_data/SRR1246789_2.fastq > extracted/SRR1246789_2_65_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_66_paired_extract.sed raw_data/SRR1246789_1.fastq > extracted/SRR1246789_1_66_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_66_paired_extract.sed raw_data/SRR1246789_2.fastq > extracted/SRR1246789_2_66_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_67_paired_extract.sed raw_data/SRR1246789_1.fastq > extracted/SRR1246789_1_67_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_67_paired_extract.sed raw_data/SRR1246789_2.fastq > extracted/SRR1246789_2_67_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_68_paired_extract.sed raw_data/SRR1246789_1.fastq > extracted/SRR1246789_1_68_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_68_paired_extract.sed raw_data/SRR1246789_2.fastq > extracted/SRR1246789_2_68_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_69_paired_extract.sed raw_data/SRR1246789_1.fastq > extracted/SRR1246789_1_69_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_69_paired_extract.sed raw_data/SRR1246789_2.fastq > extracted/SRR1246789_2_69_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_70_paired_extract.sed raw_data/SRR1246789_1.fastq > extracted/SRR1246789_1_70_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_70_paired_extract.sed raw_data/SRR1246789_2.fastq > extracted/SRR1246789_2_70_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_71_paired_extract.sed raw_data/SRR1246789_1.fastq > extracted/SRR1246789_1_71_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_71_paired_extract.sed raw_data/SRR1246789_2.fastq > extracted/SRR1246789_2_71_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_72_paired_extract.sed raw_data/SRR1246789_1.fastq > extracted/SRR1246789_1_72_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_72_paired_extract.sed raw_data/SRR1246789_2.fastq > extracted/SRR1246789_2_72_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_73_paired_extract.sed raw_data/SRR1246789_1.fastq > extracted/SRR1246789_1_73_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_73_paired_extract.sed raw_data/SRR1246789_2.fastq > extracted/SRR1246789_2_73_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_74_paired_extract.sed raw_data/SRR1246789_1.fastq > extracted/SRR1246789_1_74_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_74_paired_extract.sed raw_data/SRR1246789_2.fastq > extracted/SRR1246789_2_74_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_75_paired_extract.sed raw_data/SRR1246789_1.fastq > extracted/SRR1246789_1_75_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_75_paired_extract.sed raw_data/SRR1246789_2.fastq > extracted/SRR1246789_2_75_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_76_paired_extract.sed raw_data/SRR1246789_1.fastq > extracted/SRR1246789_1_76_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_76_paired_extract.sed raw_data/SRR1246789_2.fastq > extracted/SRR1246789_2_76_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_77_paired_extract.sed raw_data/SRR1246789_1.fastq > extracted/SRR1246789_1_77_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_77_paired_extract.sed raw_data/SRR1246789_2.fastq > extracted/SRR1246789_2_77_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_78_paired_extract.sed raw_data/SRR1246789_1.fastq > extracted/SRR1246789_1_78_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_78_paired_extract.sed raw_data/SRR1246789_2.fastq > extracted/SRR1246789_2_78_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_79_paired_extract.sed raw_data/SRR1246789_1.fastq > extracted/SRR1246789_1_79_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_79_paired_extract.sed raw_data/SRR1246789_2.fastq > extracted/SRR1246789_2_79_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_80_paired_extract.sed raw_data/SRR1246789_1.fastq > extracted/SRR1246789_1_80_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_80_paired_extract.sed raw_data/SRR1246789_2.fastq > extracted/SRR1246789_2_80_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_81_paired_extract.sed raw_data/SRR1246789_1.fastq > extracted/SRR1246789_1_81_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_81_paired_extract.sed raw_data/SRR1246789_2.fastq > extracted/SRR1246789_2_81_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_82_paired_extract.sed raw_data/SRR1246789_1.fastq > extracted/SRR1246789_1_82_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_82_paired_extract.sed raw_data/SRR1246789_2.fastq > extracted/SRR1246789_2_82_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_83_paired_extract.sed raw_data/SRR1246789_1.fastq > extracted/SRR1246789_1_83_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_83_paired_extract.sed raw_data/SRR1246789_2.fastq > extracted/SRR1246789_2_83_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_84_paired_extract.sed raw_data/SRR1246789_1.fastq > extracted/SRR1246789_1_84_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246789_84_paired_extract.sed raw_data/SRR1246789_2.fastq > extracted/SRR1246789_2_84_extracted.fastq &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246789_with_3_repeats_1_1_extract.sed raw_data/SRR1246789_1.fastq > extracted_separated/SRR1246789_with_3_repeats_1.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246789_with_3_repeats_1_2_extract.sed raw_data/SRR1246789_1.fastq > extracted_separated/SRR1246789_with_3_repeats_2.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246789_with_3_repeats_1_3_extract.sed raw_data/SRR1246789_1.fastq > extracted_separated/SRR1246789_with_3_repeats_3.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246789_with_3_repeats_1_4_extract.sed raw_data/SRR1246789_1.fastq > extracted_separated/SRR1246789_with_3_repeats_4.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246789_with_3_repeats_1_5_extract.sed raw_data/SRR1246789_1.fastq > extracted_separated/SRR1246789_with_3_repeats_5.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246789_with_3_repeats_1_6_extract.sed raw_data/SRR1246789_1.fastq > extracted_separated/SRR1246789_with_3_repeats_6.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246789_with_3_repeats_1_7_extract.sed raw_data/SRR1246789_1.fastq > extracted_separated/SRR1246789_with_3_repeats_7.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246789_with_3_repeats_1_8_extract.sed raw_data/SRR1246789_1.fastq > extracted_separated/SRR1246789_with_3_repeats_8.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246789_with_3_repeats_1_9_extract.sed raw_data/SRR1246789_1.fastq > extracted_separated/SRR1246789_with_3_repeats_9.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246789_with_3_repeats_1_10_extract.sed raw_data/SRR1246789_1.fastq > extracted_separated/SRR1246789_with_3_repeats_10.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246789_with_3_repeats_1_11_extract.sed raw_data/SRR1246789_1.fastq > extracted_separated/SRR1246789_with_3_repeats_11.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246789_with_3_repeats_1_12_extract.sed raw_data/SRR1246789_1.fastq > extracted_separated/SRR1246789_with_3_repeats_12.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246789_with_3_repeats_1_13_extract.sed raw_data/SRR1246789_1.fastq > extracted_separated/SRR1246789_with_3_repeats_13.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246789_with_3_repeats_1_14_extract.sed raw_data/SRR1246789_1.fastq > extracted_separated/SRR1246789_with_3_repeats_14.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246789_with_3_repeats_1_15_extract.sed raw_data/SRR1246789_1.fastq > extracted_separated/SRR1246789_with_3_repeats_15.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246789_with_3_repeats_1_16_extract.sed raw_data/SRR1246789_1.fastq > extracted_separated/SRR1246789_with_3_repeats_16.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246789_with_3_repeats_1_17_extract.sed raw_data/SRR1246789_1.fastq > extracted_separated/SRR1246789_with_3_repeats_17.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246789_with_3_repeats_1_18_extract.sed raw_data/SRR1246789_1.fastq > extracted_separated/SRR1246789_with_3_repeats_18.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246789_with_3_repeats_1_19_extract.sed raw_data/SRR1246789_1.fastq > extracted_separated/SRR1246789_with_3_repeats_19.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246789_with_3_repeats_1_20_extract.sed raw_data/SRR1246789_1.fastq > extracted_separated/SRR1246789_with_3_repeats_20.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246789_with_3_repeats_1_21_extract.sed raw_data/SRR1246789_1.fastq > extracted_separated/SRR1246789_with_3_repeats_21.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246789_with_3_repeats_1_22_extract.sed raw_data/SRR1246789_1.fastq > extracted_separated/SRR1246789_with_3_repeats_22.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246789_with_3_repeats_1_23_extract.sed raw_data/SRR1246789_1.fastq > extracted_separated/SRR1246789_with_3_repeats_23.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246789_with_3_repeats_1_24_extract.sed raw_data/SRR1246789_1.fastq > extracted_separated/SRR1246789_with_3_repeats_24.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246789_with_3_repeats_1_25_extract.sed raw_data/SRR1246789_1.fastq > extracted_separated/SRR1246789_with_3_repeats_25.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246789_with_3_repeats_1_26_extract.sed raw_data/SRR1246789_1.fastq > extracted_separated/SRR1246789_with_3_repeats_26.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246789_with_3_repeats_1_27_extract.sed raw_data/SRR1246789_1.fastq > extracted_separated/SRR1246789_with_3_repeats_27.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246789_with_3_repeats_1_28_extract.sed raw_data/SRR1246789_1.fastq > extracted_separated/SRR1246789_with_3_repeats_28.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246789_with_3_repeats_1_29_extract.sed raw_data/SRR1246789_1.fastq > extracted_separated/SRR1246789_with_3_repeats_29.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246789_with_3_repeats_1_30_extract.sed raw_data/SRR1246789_1.fastq > extracted_separated/SRR1246789_with_3_repeats_30.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246789_with_3_repeats_1_31_extract.sed raw_data/SRR1246789_1.fastq > extracted_separated/SRR1246789_with_3_repeats_31.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246789_with_3_repeats_1_32_extract.sed raw_data/SRR1246789_1.fastq > extracted_separated/SRR1246789_with_3_repeats_32.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246789_with_3_repeats_1_33_extract.sed raw_data/SRR1246789_1.fastq > extracted_separated/SRR1246789_with_3_repeats_33.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246789_with_3_repeats_1_34_extract.sed raw_data/SRR1246789_1.fastq > extracted_separated/SRR1246789_with_3_repeats_34.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246789_with_3_repeats_1_35_extract.sed raw_data/SRR1246789_1.fastq > extracted_separated/SRR1246789_with_3_repeats_35.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246789_with_3_repeats_1_36_extract.sed raw_data/SRR1246789_1.fastq > extracted_separated/SRR1246789_with_3_repeats_36.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246789_with_3_repeats_1_37_extract.sed raw_data/SRR1246789_1.fastq > extracted_separated/SRR1246789_with_3_repeats_37.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246789_with_3_repeats_1_38_extract.sed raw_data/SRR1246789_1.fastq > extracted_separated/SRR1246789_with_3_repeats_38.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246789_with_3_repeats_1_39_extract.sed raw_data/SRR1246789_1.fastq > extracted_separated/SRR1246789_with_3_repeats_39.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246789_with_3_repeats_1_40_extract.sed raw_data/SRR1246789_1.fastq > extracted_separated/SRR1246789_with_3_repeats_40.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246789_with_3_repeats_1_41_extract.sed raw_data/SRR1246789_1.fastq > extracted_separated/SRR1246789_with_3_repeats_41.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246789_with_3_repeats_1_42_extract.sed raw_data/SRR1246789_1.fastq > extracted_separated/SRR1246789_with_3_repeats_42.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246789_with_3_repeats_1_43_extract.sed raw_data/SRR1246789_1.fastq > extracted_separated/SRR1246789_with_3_repeats_43.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246789_with_3_repeats_1_44_extract.sed raw_data/SRR1246789_1.fastq > extracted_separated/SRR1246789_with_3_repeats_44.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246789_with_3_repeats_1_45_extract.sed raw_data/SRR1246789_1.fastq > extracted_separated/SRR1246789_with_3_repeats_45.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246789_with_3_repeats_1_46_extract.sed raw_data/SRR1246789_1.fastq > extracted_separated/SRR1246789_with_3_repeats_46.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246789_with_3_repeats_1_47_extract.sed raw_data/SRR1246789_1.fastq > extracted_separated/SRR1246789_with_3_repeats_47.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246789_with_3_repeats_1_48_extract.sed raw_data/SRR1246789_1.fastq > extracted_separated/SRR1246789_with_3_repeats_48.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246789_with_3_repeats_1_49_extract.sed raw_data/SRR1246789_1.fastq > extracted_separated/SRR1246789_with_3_repeats_49.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246789_with_3_repeats_1_50_extract.sed raw_data/SRR1246789_1.fastq > extracted_separated/SRR1246789_with_3_repeats_50.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246789_with_3_repeats_1_51_extract.sed raw_data/SRR1246789_1.fastq > extracted_separated/SRR1246789_with_3_repeats_51.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246789_with_3_repeats_2_1_extract.sed raw_data/SRR1246789_2.fastq > extracted_separated/SRR1246789_with_3_repeats_1.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246789_with_3_repeats_2_2_extract.sed raw_data/SRR1246789_2.fastq > extracted_separated/SRR1246789_with_3_repeats_2.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246789_with_3_repeats_2_3_extract.sed raw_data/SRR1246789_2.fastq > extracted_separated/SRR1246789_with_3_repeats_3.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246789_with_3_repeats_2_4_extract.sed raw_data/SRR1246789_2.fastq > extracted_separated/SRR1246789_with_3_repeats_4.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246789_with_3_repeats_2_5_extract.sed raw_data/SRR1246789_2.fastq > extracted_separated/SRR1246789_with_3_repeats_5.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246789_with_3_repeats_2_6_extract.sed raw_data/SRR1246789_2.fastq > extracted_separated/SRR1246789_with_3_repeats_6.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246789_with_3_repeats_2_7_extract.sed raw_data/SRR1246789_2.fastq > extracted_separated/SRR1246789_with_3_repeats_7.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246789_with_3_repeats_2_8_extract.sed raw_data/SRR1246789_2.fastq > extracted_separated/SRR1246789_with_3_repeats_8.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246789_with_3_repeats_2_9_extract.sed raw_data/SRR1246789_2.fastq > extracted_separated/SRR1246789_with_3_repeats_9.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246789_with_3_repeats_2_10_extract.sed raw_data/SRR1246789_2.fastq > extracted_separated/SRR1246789_with_3_repeats_10.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246789_with_3_repeats_2_11_extract.sed raw_data/SRR1246789_2.fastq > extracted_separated/SRR1246789_with_3_repeats_11.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246789_with_3_repeats_2_12_extract.sed raw_data/SRR1246789_2.fastq > extracted_separated/SRR1246789_with_3_repeats_12.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246789_with_3_repeats_2_13_extract.sed raw_data/SRR1246789_2.fastq > extracted_separated/SRR1246789_with_3_repeats_13.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246789_with_3_repeats_2_14_extract.sed raw_data/SRR1246789_2.fastq > extracted_separated/SRR1246789_with_3_repeats_14.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246789_with_3_repeats_2_15_extract.sed raw_data/SRR1246789_2.fastq > extracted_separated/SRR1246789_with_3_repeats_15.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246789_with_3_repeats_2_16_extract.sed raw_data/SRR1246789_2.fastq > extracted_separated/SRR1246789_with_3_repeats_16.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246789_with_3_repeats_2_17_extract.sed raw_data/SRR1246789_2.fastq > extracted_separated/SRR1246789_with_3_repeats_17.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246789_with_3_repeats_2_18_extract.sed raw_data/SRR1246789_2.fastq > extracted_separated/SRR1246789_with_3_repeats_18.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246789_with_3_repeats_2_19_extract.sed raw_data/SRR1246789_2.fastq > extracted_separated/SRR1246789_with_3_repeats_19.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246789_with_3_repeats_2_20_extract.sed raw_data/SRR1246789_2.fastq > extracted_separated/SRR1246789_with_3_repeats_20.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246789_with_3_repeats_2_21_extract.sed raw_data/SRR1246789_2.fastq > extracted_separated/SRR1246789_with_3_repeats_21.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246789_with_3_repeats_2_22_extract.sed raw_data/SRR1246789_2.fastq > extracted_separated/SRR1246789_with_3_repeats_22.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246789_with_3_repeats_2_23_extract.sed raw_data/SRR1246789_2.fastq > extracted_separated/SRR1246789_with_3_repeats_23.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246789_with_3_repeats_2_24_extract.sed raw_data/SRR1246789_2.fastq > extracted_separated/SRR1246789_with_3_repeats_24.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246789_with_3_repeats_2_25_extract.sed raw_data/SRR1246789_2.fastq > extracted_separated/SRR1246789_with_3_repeats_25.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246789_with_3_repeats_2_26_extract.sed raw_data/SRR1246789_2.fastq > extracted_separated/SRR1246789_with_3_repeats_26.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246789_with_3_repeats_2_27_extract.sed raw_data/SRR1246789_2.fastq > extracted_separated/SRR1246789_with_3_repeats_27.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246789_with_3_repeats_2_28_extract.sed raw_data/SRR1246789_2.fastq > extracted_separated/SRR1246789_with_3_repeats_28.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246789_with_3_repeats_2_29_extract.sed raw_data/SRR1246789_2.fastq > extracted_separated/SRR1246789_with_3_repeats_29.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246789_with_3_repeats_2_30_extract.sed raw_data/SRR1246789_2.fastq > extracted_separated/SRR1246789_with_3_repeats_30.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246789_with_3_repeats_2_31_extract.sed raw_data/SRR1246789_2.fastq > extracted_separated/SRR1246789_with_3_repeats_31.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246789_with_3_repeats_2_32_extract.sed raw_data/SRR1246789_2.fastq > extracted_separated/SRR1246789_with_3_repeats_32.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246789_with_3_repeats_2_33_extract.sed raw_data/SRR1246789_2.fastq > extracted_separated/SRR1246789_with_3_repeats_33.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246789_with_3_repeats_2_34_extract.sed raw_data/SRR1246789_2.fastq > extracted_separated/SRR1246789_with_3_repeats_34.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246789_with_3_repeats_2_35_extract.sed raw_data/SRR1246789_2.fastq > extracted_separated/SRR1246789_with_3_repeats_35.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246789_with_3_repeats_2_36_extract.sed raw_data/SRR1246789_2.fastq > extracted_separated/SRR1246789_with_3_repeats_36.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246789_with_3_repeats_2_37_extract.sed raw_data/SRR1246789_2.fastq > extracted_separated/SRR1246789_with_3_repeats_37.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246789_with_3_repeats_2_38_extract.sed raw_data/SRR1246789_2.fastq > extracted_separated/SRR1246789_with_3_repeats_38.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246789_with_3_repeats_2_39_extract.sed raw_data/SRR1246789_2.fastq > extracted_separated/SRR1246789_with_3_repeats_39.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246789_with_3_repeats_2_40_extract.sed raw_data/SRR1246789_2.fastq > extracted_separated/SRR1246789_with_3_repeats_40.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246789_with_3_repeats_2_41_extract.sed raw_data/SRR1246789_2.fastq > extracted_separated/SRR1246789_with_3_repeats_41.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246789_with_3_repeats_2_42_extract.sed raw_data/SRR1246789_2.fastq > extracted_separated/SRR1246789_with_3_repeats_42.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246789_with_3_repeats_2_43_extract.sed raw_data/SRR1246789_2.fastq > extracted_separated/SRR1246789_with_3_repeats_43.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246789_with_3_repeats_2_44_extract.sed raw_data/SRR1246789_2.fastq > extracted_separated/SRR1246789_with_3_repeats_44.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246789_without_3_repeats_1_extract.sed raw_data/SRR1246789_1.fastq > extracted_separated/SRR1246789_without_3_repeats.fastq ;
                                            sed -n -f scripts/sed_scripts/SRR1246789_without_3_repeats_2_extract.sed raw_data/SRR1246789_2.fastq >> extracted_separated/SRR1246789_without_3_repeats.fastq" &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246790_1_paired_extract.sed raw_data/SRR1246790_1.fastq > extracted/SRR1246790_1_1_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246790_1_paired_extract.sed raw_data/SRR1246790_2.fastq > extracted/SRR1246790_2_1_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246790_2_paired_extract.sed raw_data/SRR1246790_1.fastq > extracted/SRR1246790_1_2_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246790_2_paired_extract.sed raw_data/SRR1246790_2.fastq > extracted/SRR1246790_2_2_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246790_3_paired_extract.sed raw_data/SRR1246790_1.fastq > extracted/SRR1246790_1_3_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246790_3_paired_extract.sed raw_data/SRR1246790_2.fastq > extracted/SRR1246790_2_3_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246790_4_paired_extract.sed raw_data/SRR1246790_1.fastq > extracted/SRR1246790_1_4_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246790_4_paired_extract.sed raw_data/SRR1246790_2.fastq > extracted/SRR1246790_2_4_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246790_5_paired_extract.sed raw_data/SRR1246790_1.fastq > extracted/SRR1246790_1_5_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246790_5_paired_extract.sed raw_data/SRR1246790_2.fastq > extracted/SRR1246790_2_5_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246790_6_paired_extract.sed raw_data/SRR1246790_1.fastq > extracted/SRR1246790_1_6_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246790_6_paired_extract.sed raw_data/SRR1246790_2.fastq > extracted/SRR1246790_2_6_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246790_7_paired_extract.sed raw_data/SRR1246790_1.fastq > extracted/SRR1246790_1_7_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246790_7_paired_extract.sed raw_data/SRR1246790_2.fastq > extracted/SRR1246790_2_7_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246790_8_paired_extract.sed raw_data/SRR1246790_1.fastq > extracted/SRR1246790_1_8_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246790_8_paired_extract.sed raw_data/SRR1246790_2.fastq > extracted/SRR1246790_2_8_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246790_9_paired_extract.sed raw_data/SRR1246790_1.fastq > extracted/SRR1246790_1_9_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246790_9_paired_extract.sed raw_data/SRR1246790_2.fastq > extracted/SRR1246790_2_9_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246790_10_paired_extract.sed raw_data/SRR1246790_1.fastq > extracted/SRR1246790_1_10_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246790_10_paired_extract.sed raw_data/SRR1246790_2.fastq > extracted/SRR1246790_2_10_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246790_11_paired_extract.sed raw_data/SRR1246790_1.fastq > extracted/SRR1246790_1_11_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246790_11_paired_extract.sed raw_data/SRR1246790_2.fastq > extracted/SRR1246790_2_11_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246790_12_paired_extract.sed raw_data/SRR1246790_1.fastq > extracted/SRR1246790_1_12_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246790_12_paired_extract.sed raw_data/SRR1246790_2.fastq > extracted/SRR1246790_2_12_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246790_13_paired_extract.sed raw_data/SRR1246790_1.fastq > extracted/SRR1246790_1_13_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246790_13_paired_extract.sed raw_data/SRR1246790_2.fastq > extracted/SRR1246790_2_13_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246790_14_paired_extract.sed raw_data/SRR1246790_1.fastq > extracted/SRR1246790_1_14_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246790_14_paired_extract.sed raw_data/SRR1246790_2.fastq > extracted/SRR1246790_2_14_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246790_15_paired_extract.sed raw_data/SRR1246790_1.fastq > extracted/SRR1246790_1_15_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246790_15_paired_extract.sed raw_data/SRR1246790_2.fastq > extracted/SRR1246790_2_15_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246790_16_paired_extract.sed raw_data/SRR1246790_1.fastq > extracted/SRR1246790_1_16_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246790_16_paired_extract.sed raw_data/SRR1246790_2.fastq > extracted/SRR1246790_2_16_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246790_17_paired_extract.sed raw_data/SRR1246790_1.fastq > extracted/SRR1246790_1_17_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246790_17_paired_extract.sed raw_data/SRR1246790_2.fastq > extracted/SRR1246790_2_17_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246790_18_paired_extract.sed raw_data/SRR1246790_1.fastq > extracted/SRR1246790_1_18_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246790_18_paired_extract.sed raw_data/SRR1246790_2.fastq > extracted/SRR1246790_2_18_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246790_19_paired_extract.sed raw_data/SRR1246790_1.fastq > extracted/SRR1246790_1_19_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246790_19_paired_extract.sed raw_data/SRR1246790_2.fastq > extracted/SRR1246790_2_19_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246790_20_paired_extract.sed raw_data/SRR1246790_1.fastq > extracted/SRR1246790_1_20_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246790_20_paired_extract.sed raw_data/SRR1246790_2.fastq > extracted/SRR1246790_2_20_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246790_21_paired_extract.sed raw_data/SRR1246790_1.fastq > extracted/SRR1246790_1_21_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246790_21_paired_extract.sed raw_data/SRR1246790_2.fastq > extracted/SRR1246790_2_21_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246790_22_paired_extract.sed raw_data/SRR1246790_1.fastq > extracted/SRR1246790_1_22_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246790_22_paired_extract.sed raw_data/SRR1246790_2.fastq > extracted/SRR1246790_2_22_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246790_23_paired_extract.sed raw_data/SRR1246790_1.fastq > extracted/SRR1246790_1_23_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246790_23_paired_extract.sed raw_data/SRR1246790_2.fastq > extracted/SRR1246790_2_23_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246790_24_paired_extract.sed raw_data/SRR1246790_1.fastq > extracted/SRR1246790_1_24_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246790_24_paired_extract.sed raw_data/SRR1246790_2.fastq > extracted/SRR1246790_2_24_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246790_25_paired_extract.sed raw_data/SRR1246790_1.fastq > extracted/SRR1246790_1_25_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246790_25_paired_extract.sed raw_data/SRR1246790_2.fastq > extracted/SRR1246790_2_25_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246790_26_paired_extract.sed raw_data/SRR1246790_1.fastq > extracted/SRR1246790_1_26_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246790_26_paired_extract.sed raw_data/SRR1246790_2.fastq > extracted/SRR1246790_2_26_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246790_27_paired_extract.sed raw_data/SRR1246790_1.fastq > extracted/SRR1246790_1_27_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246790_27_paired_extract.sed raw_data/SRR1246790_2.fastq > extracted/SRR1246790_2_27_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246790_28_paired_extract.sed raw_data/SRR1246790_1.fastq > extracted/SRR1246790_1_28_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246790_28_paired_extract.sed raw_data/SRR1246790_2.fastq > extracted/SRR1246790_2_28_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246790_29_paired_extract.sed raw_data/SRR1246790_1.fastq > extracted/SRR1246790_1_29_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246790_29_paired_extract.sed raw_data/SRR1246790_2.fastq > extracted/SRR1246790_2_29_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246790_30_paired_extract.sed raw_data/SRR1246790_1.fastq > extracted/SRR1246790_1_30_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246790_30_paired_extract.sed raw_data/SRR1246790_2.fastq > extracted/SRR1246790_2_30_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246790_31_paired_extract.sed raw_data/SRR1246790_1.fastq > extracted/SRR1246790_1_31_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246790_31_paired_extract.sed raw_data/SRR1246790_2.fastq > extracted/SRR1246790_2_31_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246790_32_paired_extract.sed raw_data/SRR1246790_1.fastq > extracted/SRR1246790_1_32_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246790_32_paired_extract.sed raw_data/SRR1246790_2.fastq > extracted/SRR1246790_2_32_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246790_33_paired_extract.sed raw_data/SRR1246790_1.fastq > extracted/SRR1246790_1_33_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246790_33_paired_extract.sed raw_data/SRR1246790_2.fastq > extracted/SRR1246790_2_33_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246790_34_paired_extract.sed raw_data/SRR1246790_1.fastq > extracted/SRR1246790_1_34_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246790_34_paired_extract.sed raw_data/SRR1246790_2.fastq > extracted/SRR1246790_2_34_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246790_35_paired_extract.sed raw_data/SRR1246790_1.fastq > extracted/SRR1246790_1_35_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246790_35_paired_extract.sed raw_data/SRR1246790_2.fastq > extracted/SRR1246790_2_35_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246790_36_paired_extract.sed raw_data/SRR1246790_1.fastq > extracted/SRR1246790_1_36_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246790_36_paired_extract.sed raw_data/SRR1246790_2.fastq > extracted/SRR1246790_2_36_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246790_37_paired_extract.sed raw_data/SRR1246790_1.fastq > extracted/SRR1246790_1_37_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246790_37_paired_extract.sed raw_data/SRR1246790_2.fastq > extracted/SRR1246790_2_37_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246790_38_paired_extract.sed raw_data/SRR1246790_1.fastq > extracted/SRR1246790_1_38_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246790_38_paired_extract.sed raw_data/SRR1246790_2.fastq > extracted/SRR1246790_2_38_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246790_39_paired_extract.sed raw_data/SRR1246790_1.fastq > extracted/SRR1246790_1_39_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246790_39_paired_extract.sed raw_data/SRR1246790_2.fastq > extracted/SRR1246790_2_39_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246790_40_paired_extract.sed raw_data/SRR1246790_1.fastq > extracted/SRR1246790_1_40_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246790_40_paired_extract.sed raw_data/SRR1246790_2.fastq > extracted/SRR1246790_2_40_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246790_41_paired_extract.sed raw_data/SRR1246790_1.fastq > extracted/SRR1246790_1_41_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246790_41_paired_extract.sed raw_data/SRR1246790_2.fastq > extracted/SRR1246790_2_41_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246790_42_paired_extract.sed raw_data/SRR1246790_1.fastq > extracted/SRR1246790_1_42_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246790_42_paired_extract.sed raw_data/SRR1246790_2.fastq > extracted/SRR1246790_2_42_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246790_43_paired_extract.sed raw_data/SRR1246790_1.fastq > extracted/SRR1246790_1_43_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246790_43_paired_extract.sed raw_data/SRR1246790_2.fastq > extracted/SRR1246790_2_43_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246790_44_paired_extract.sed raw_data/SRR1246790_1.fastq > extracted/SRR1246790_1_44_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246790_44_paired_extract.sed raw_data/SRR1246790_2.fastq > extracted/SRR1246790_2_44_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246790_45_paired_extract.sed raw_data/SRR1246790_1.fastq > extracted/SRR1246790_1_45_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246790_45_paired_extract.sed raw_data/SRR1246790_2.fastq > extracted/SRR1246790_2_45_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246790_46_paired_extract.sed raw_data/SRR1246790_1.fastq > extracted/SRR1246790_1_46_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246790_46_paired_extract.sed raw_data/SRR1246790_2.fastq > extracted/SRR1246790_2_46_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246790_47_paired_extract.sed raw_data/SRR1246790_1.fastq > extracted/SRR1246790_1_47_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246790_47_paired_extract.sed raw_data/SRR1246790_2.fastq > extracted/SRR1246790_2_47_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246790_48_paired_extract.sed raw_data/SRR1246790_1.fastq > extracted/SRR1246790_1_48_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246790_48_paired_extract.sed raw_data/SRR1246790_2.fastq > extracted/SRR1246790_2_48_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246790_49_paired_extract.sed raw_data/SRR1246790_1.fastq > extracted/SRR1246790_1_49_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246790_49_paired_extract.sed raw_data/SRR1246790_2.fastq > extracted/SRR1246790_2_49_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246790_50_paired_extract.sed raw_data/SRR1246790_1.fastq > extracted/SRR1246790_1_50_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246790_50_paired_extract.sed raw_data/SRR1246790_2.fastq > extracted/SRR1246790_2_50_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246790_51_paired_extract.sed raw_data/SRR1246790_1.fastq > extracted/SRR1246790_1_51_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246790_51_paired_extract.sed raw_data/SRR1246790_2.fastq > extracted/SRR1246790_2_51_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246790_52_paired_extract.sed raw_data/SRR1246790_1.fastq > extracted/SRR1246790_1_52_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246790_52_paired_extract.sed raw_data/SRR1246790_2.fastq > extracted/SRR1246790_2_52_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246790_53_paired_extract.sed raw_data/SRR1246790_1.fastq > extracted/SRR1246790_1_53_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246790_53_paired_extract.sed raw_data/SRR1246790_2.fastq > extracted/SRR1246790_2_53_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246790_54_paired_extract.sed raw_data/SRR1246790_1.fastq > extracted/SRR1246790_1_54_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246790_54_paired_extract.sed raw_data/SRR1246790_2.fastq > extracted/SRR1246790_2_54_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246790_55_paired_extract.sed raw_data/SRR1246790_1.fastq > extracted/SRR1246790_1_55_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246790_55_paired_extract.sed raw_data/SRR1246790_2.fastq > extracted/SRR1246790_2_55_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246790_56_paired_extract.sed raw_data/SRR1246790_1.fastq > extracted/SRR1246790_1_56_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246790_56_paired_extract.sed raw_data/SRR1246790_2.fastq > extracted/SRR1246790_2_56_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246790_57_paired_extract.sed raw_data/SRR1246790_1.fastq > extracted/SRR1246790_1_57_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246790_57_paired_extract.sed raw_data/SRR1246790_2.fastq > extracted/SRR1246790_2_57_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246790_58_paired_extract.sed raw_data/SRR1246790_1.fastq > extracted/SRR1246790_1_58_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246790_58_paired_extract.sed raw_data/SRR1246790_2.fastq > extracted/SRR1246790_2_58_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246790_59_paired_extract.sed raw_data/SRR1246790_1.fastq > extracted/SRR1246790_1_59_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246790_59_paired_extract.sed raw_data/SRR1246790_2.fastq > extracted/SRR1246790_2_59_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246790_60_paired_extract.sed raw_data/SRR1246790_1.fastq > extracted/SRR1246790_1_60_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246790_60_paired_extract.sed raw_data/SRR1246790_2.fastq > extracted/SRR1246790_2_60_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246790_61_paired_extract.sed raw_data/SRR1246790_1.fastq > extracted/SRR1246790_1_61_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246790_61_paired_extract.sed raw_data/SRR1246790_2.fastq > extracted/SRR1246790_2_61_extracted.fastq &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246790_with_3_repeats_1_1_extract.sed raw_data/SRR1246790_1.fastq > extracted_separated/SRR1246790_with_3_repeats_1.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246790_with_3_repeats_1_2_extract.sed raw_data/SRR1246790_1.fastq > extracted_separated/SRR1246790_with_3_repeats_2.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246790_with_3_repeats_1_3_extract.sed raw_data/SRR1246790_1.fastq > extracted_separated/SRR1246790_with_3_repeats_3.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246790_with_3_repeats_1_4_extract.sed raw_data/SRR1246790_1.fastq > extracted_separated/SRR1246790_with_3_repeats_4.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246790_with_3_repeats_1_5_extract.sed raw_data/SRR1246790_1.fastq > extracted_separated/SRR1246790_with_3_repeats_5.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246790_with_3_repeats_1_6_extract.sed raw_data/SRR1246790_1.fastq > extracted_separated/SRR1246790_with_3_repeats_6.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246790_with_3_repeats_1_7_extract.sed raw_data/SRR1246790_1.fastq > extracted_separated/SRR1246790_with_3_repeats_7.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246790_with_3_repeats_1_8_extract.sed raw_data/SRR1246790_1.fastq > extracted_separated/SRR1246790_with_3_repeats_8.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246790_with_3_repeats_1_9_extract.sed raw_data/SRR1246790_1.fastq > extracted_separated/SRR1246790_with_3_repeats_9.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246790_with_3_repeats_1_10_extract.sed raw_data/SRR1246790_1.fastq > extracted_separated/SRR1246790_with_3_repeats_10.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246790_with_3_repeats_1_11_extract.sed raw_data/SRR1246790_1.fastq > extracted_separated/SRR1246790_with_3_repeats_11.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246790_with_3_repeats_1_12_extract.sed raw_data/SRR1246790_1.fastq > extracted_separated/SRR1246790_with_3_repeats_12.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246790_with_3_repeats_1_13_extract.sed raw_data/SRR1246790_1.fastq > extracted_separated/SRR1246790_with_3_repeats_13.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246790_with_3_repeats_1_14_extract.sed raw_data/SRR1246790_1.fastq > extracted_separated/SRR1246790_with_3_repeats_14.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246790_with_3_repeats_1_15_extract.sed raw_data/SRR1246790_1.fastq > extracted_separated/SRR1246790_with_3_repeats_15.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246790_with_3_repeats_1_16_extract.sed raw_data/SRR1246790_1.fastq > extracted_separated/SRR1246790_with_3_repeats_16.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246790_with_3_repeats_1_17_extract.sed raw_data/SRR1246790_1.fastq > extracted_separated/SRR1246790_with_3_repeats_17.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246790_with_3_repeats_1_18_extract.sed raw_data/SRR1246790_1.fastq > extracted_separated/SRR1246790_with_3_repeats_18.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246790_with_3_repeats_1_19_extract.sed raw_data/SRR1246790_1.fastq > extracted_separated/SRR1246790_with_3_repeats_19.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246790_with_3_repeats_1_20_extract.sed raw_data/SRR1246790_1.fastq > extracted_separated/SRR1246790_with_3_repeats_20.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246790_with_3_repeats_1_21_extract.sed raw_data/SRR1246790_1.fastq > extracted_separated/SRR1246790_with_3_repeats_21.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246790_with_3_repeats_1_22_extract.sed raw_data/SRR1246790_1.fastq > extracted_separated/SRR1246790_with_3_repeats_22.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246790_with_3_repeats_1_23_extract.sed raw_data/SRR1246790_1.fastq > extracted_separated/SRR1246790_with_3_repeats_23.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246790_with_3_repeats_1_24_extract.sed raw_data/SRR1246790_1.fastq > extracted_separated/SRR1246790_with_3_repeats_24.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246790_with_3_repeats_1_25_extract.sed raw_data/SRR1246790_1.fastq > extracted_separated/SRR1246790_with_3_repeats_25.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246790_with_3_repeats_1_26_extract.sed raw_data/SRR1246790_1.fastq > extracted_separated/SRR1246790_with_3_repeats_26.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246790_with_3_repeats_1_27_extract.sed raw_data/SRR1246790_1.fastq > extracted_separated/SRR1246790_with_3_repeats_27.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246790_with_3_repeats_1_28_extract.sed raw_data/SRR1246790_1.fastq > extracted_separated/SRR1246790_with_3_repeats_28.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246790_with_3_repeats_1_29_extract.sed raw_data/SRR1246790_1.fastq > extracted_separated/SRR1246790_with_3_repeats_29.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246790_with_3_repeats_1_30_extract.sed raw_data/SRR1246790_1.fastq > extracted_separated/SRR1246790_with_3_repeats_30.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246790_with_3_repeats_1_31_extract.sed raw_data/SRR1246790_1.fastq > extracted_separated/SRR1246790_with_3_repeats_31.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246790_with_3_repeats_1_32_extract.sed raw_data/SRR1246790_1.fastq > extracted_separated/SRR1246790_with_3_repeats_32.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246790_with_3_repeats_1_33_extract.sed raw_data/SRR1246790_1.fastq > extracted_separated/SRR1246790_with_3_repeats_33.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246790_with_3_repeats_1_34_extract.sed raw_data/SRR1246790_1.fastq > extracted_separated/SRR1246790_with_3_repeats_34.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246790_with_3_repeats_1_35_extract.sed raw_data/SRR1246790_1.fastq > extracted_separated/SRR1246790_with_3_repeats_35.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246790_with_3_repeats_1_36_extract.sed raw_data/SRR1246790_1.fastq > extracted_separated/SRR1246790_with_3_repeats_36.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246790_with_3_repeats_1_37_extract.sed raw_data/SRR1246790_1.fastq > extracted_separated/SRR1246790_with_3_repeats_37.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246790_with_3_repeats_1_38_extract.sed raw_data/SRR1246790_1.fastq > extracted_separated/SRR1246790_with_3_repeats_38.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246790_with_3_repeats_1_39_extract.sed raw_data/SRR1246790_1.fastq > extracted_separated/SRR1246790_with_3_repeats_39.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246790_with_3_repeats_2_1_extract.sed raw_data/SRR1246790_2.fastq > extracted_separated/SRR1246790_with_3_repeats_1.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246790_with_3_repeats_2_2_extract.sed raw_data/SRR1246790_2.fastq > extracted_separated/SRR1246790_with_3_repeats_2.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246790_with_3_repeats_2_3_extract.sed raw_data/SRR1246790_2.fastq > extracted_separated/SRR1246790_with_3_repeats_3.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246790_with_3_repeats_2_4_extract.sed raw_data/SRR1246790_2.fastq > extracted_separated/SRR1246790_with_3_repeats_4.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246790_with_3_repeats_2_5_extract.sed raw_data/SRR1246790_2.fastq > extracted_separated/SRR1246790_with_3_repeats_5.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246790_with_3_repeats_2_6_extract.sed raw_data/SRR1246790_2.fastq > extracted_separated/SRR1246790_with_3_repeats_6.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246790_with_3_repeats_2_7_extract.sed raw_data/SRR1246790_2.fastq > extracted_separated/SRR1246790_with_3_repeats_7.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246790_with_3_repeats_2_8_extract.sed raw_data/SRR1246790_2.fastq > extracted_separated/SRR1246790_with_3_repeats_8.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246790_with_3_repeats_2_9_extract.sed raw_data/SRR1246790_2.fastq > extracted_separated/SRR1246790_with_3_repeats_9.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246790_with_3_repeats_2_10_extract.sed raw_data/SRR1246790_2.fastq > extracted_separated/SRR1246790_with_3_repeats_10.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246790_with_3_repeats_2_11_extract.sed raw_data/SRR1246790_2.fastq > extracted_separated/SRR1246790_with_3_repeats_11.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246790_with_3_repeats_2_12_extract.sed raw_data/SRR1246790_2.fastq > extracted_separated/SRR1246790_with_3_repeats_12.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246790_with_3_repeats_2_13_extract.sed raw_data/SRR1246790_2.fastq > extracted_separated/SRR1246790_with_3_repeats_13.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246790_with_3_repeats_2_14_extract.sed raw_data/SRR1246790_2.fastq > extracted_separated/SRR1246790_with_3_repeats_14.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246790_with_3_repeats_2_15_extract.sed raw_data/SRR1246790_2.fastq > extracted_separated/SRR1246790_with_3_repeats_15.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246790_with_3_repeats_2_16_extract.sed raw_data/SRR1246790_2.fastq > extracted_separated/SRR1246790_with_3_repeats_16.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246790_with_3_repeats_2_17_extract.sed raw_data/SRR1246790_2.fastq > extracted_separated/SRR1246790_with_3_repeats_17.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246790_with_3_repeats_2_18_extract.sed raw_data/SRR1246790_2.fastq > extracted_separated/SRR1246790_with_3_repeats_18.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246790_with_3_repeats_2_19_extract.sed raw_data/SRR1246790_2.fastq > extracted_separated/SRR1246790_with_3_repeats_19.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246790_with_3_repeats_2_20_extract.sed raw_data/SRR1246790_2.fastq > extracted_separated/SRR1246790_with_3_repeats_20.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246790_with_3_repeats_2_21_extract.sed raw_data/SRR1246790_2.fastq > extracted_separated/SRR1246790_with_3_repeats_21.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246790_with_3_repeats_2_22_extract.sed raw_data/SRR1246790_2.fastq > extracted_separated/SRR1246790_with_3_repeats_22.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246790_with_3_repeats_2_23_extract.sed raw_data/SRR1246790_2.fastq > extracted_separated/SRR1246790_with_3_repeats_23.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246790_with_3_repeats_2_24_extract.sed raw_data/SRR1246790_2.fastq > extracted_separated/SRR1246790_with_3_repeats_24.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246790_with_3_repeats_2_25_extract.sed raw_data/SRR1246790_2.fastq > extracted_separated/SRR1246790_with_3_repeats_25.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246790_with_3_repeats_2_26_extract.sed raw_data/SRR1246790_2.fastq > extracted_separated/SRR1246790_with_3_repeats_26.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246790_with_3_repeats_2_27_extract.sed raw_data/SRR1246790_2.fastq > extracted_separated/SRR1246790_with_3_repeats_27.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246790_with_3_repeats_2_28_extract.sed raw_data/SRR1246790_2.fastq > extracted_separated/SRR1246790_with_3_repeats_28.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246790_with_3_repeats_2_29_extract.sed raw_data/SRR1246790_2.fastq > extracted_separated/SRR1246790_with_3_repeats_29.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246790_with_3_repeats_2_30_extract.sed raw_data/SRR1246790_2.fastq > extracted_separated/SRR1246790_with_3_repeats_30.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246790_with_3_repeats_2_31_extract.sed raw_data/SRR1246790_2.fastq > extracted_separated/SRR1246790_with_3_repeats_31.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246790_without_3_repeats_1_extract.sed raw_data/SRR1246790_1.fastq > extracted_separated/SRR1246790_without_3_repeats.fastq ;
                                            sed -n -f scripts/sed_scripts/SRR1246790_without_3_repeats_2_extract.sed raw_data/SRR1246790_2.fastq >> extracted_separated/SRR1246790_without_3_repeats.fastq" &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_1_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_1_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_1_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_1_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_2_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_2_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_2_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_2_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_3_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_3_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_3_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_3_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_4_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_4_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_4_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_4_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_5_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_5_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_5_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_5_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_6_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_6_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_6_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_6_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_7_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_7_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_7_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_7_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_8_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_8_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_8_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_8_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_9_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_9_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_9_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_9_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_10_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_10_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_10_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_10_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_11_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_11_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_11_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_11_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_12_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_12_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_12_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_12_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_13_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_13_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_13_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_13_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_14_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_14_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_14_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_14_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_15_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_15_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_15_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_15_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_16_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_16_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_16_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_16_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_17_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_17_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_17_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_17_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_18_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_18_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_18_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_18_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_19_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_19_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_19_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_19_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_20_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_20_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_20_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_20_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_21_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_21_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_21_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_21_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_22_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_22_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_22_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_22_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_23_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_23_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_23_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_23_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_24_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_24_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_24_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_24_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_25_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_25_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_25_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_25_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_26_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_26_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_26_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_26_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_27_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_27_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_27_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_27_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_28_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_28_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_28_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_28_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_29_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_29_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_29_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_29_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_30_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_30_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_30_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_30_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_31_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_31_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_31_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_31_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_32_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_32_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_32_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_32_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_33_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_33_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_33_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_33_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_34_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_34_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_34_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_34_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_35_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_35_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_35_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_35_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_36_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_36_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_36_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_36_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_37_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_37_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_37_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_37_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_38_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_38_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_38_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_38_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_39_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_39_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_39_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_39_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_40_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_40_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_40_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_40_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_41_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_41_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_41_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_41_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_42_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_42_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_42_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_42_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_43_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_43_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_43_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_43_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_44_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_44_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_44_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_44_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_45_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_45_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_45_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_45_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_46_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_46_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_46_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_46_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_47_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_47_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_47_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_47_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_48_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_48_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_48_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_48_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_49_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_49_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_49_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_49_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_50_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_50_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_50_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_50_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_51_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_51_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_51_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_51_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_52_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_52_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_52_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_52_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_53_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_53_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_53_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_53_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_54_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_54_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_54_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_54_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_55_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_55_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_55_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_55_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_56_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_56_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_56_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_56_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_57_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_57_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_57_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_57_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_58_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_58_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_58_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_58_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_59_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_59_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_59_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_59_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_60_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_60_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_60_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_60_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_61_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_61_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_61_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_61_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_62_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_62_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_62_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_62_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_63_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_63_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_63_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_63_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_64_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_64_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_64_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_64_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_65_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_65_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_65_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_65_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_66_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_66_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_66_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_66_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_67_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_67_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_67_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_67_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_68_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_68_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_68_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_68_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_69_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_69_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_69_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_69_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_70_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_70_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_70_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_70_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_71_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_71_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_71_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_71_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_72_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_72_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_72_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_72_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_73_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_73_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_73_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_73_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_74_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_74_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_74_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_74_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_75_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_75_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_75_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_75_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_76_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_76_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_76_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_76_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_77_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_77_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_77_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_77_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_78_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_78_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_78_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_78_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_79_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_79_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_79_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_79_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_80_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_80_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_80_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_80_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_81_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_81_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_81_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_81_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_82_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_82_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_82_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_82_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_83_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_83_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_83_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_83_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_84_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_84_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_84_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_84_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_85_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_85_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_85_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_85_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_86_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_86_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_86_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_86_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_87_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_87_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_87_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_87_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_88_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_88_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_88_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_88_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_89_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_89_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_89_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_89_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_90_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_90_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_90_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_90_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_91_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_91_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_91_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_91_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_92_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_92_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_92_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_92_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_93_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_93_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_93_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_93_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_94_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_94_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_94_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_94_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_95_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_95_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_95_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_95_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_96_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_96_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_96_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_96_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_97_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_97_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_97_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_97_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_98_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_98_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_98_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_98_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_99_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_99_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_99_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_99_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_100_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_100_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_100_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_100_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_101_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_101_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_101_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_101_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_102_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_102_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_102_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_102_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_103_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_103_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_103_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_103_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_104_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_104_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_104_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_104_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_105_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_105_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_105_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_105_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_106_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_106_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_106_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_106_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_107_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_107_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_107_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_107_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_108_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_108_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_108_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_108_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_109_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_109_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_109_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_109_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_110_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_110_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_110_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_110_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_111_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_111_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_111_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_111_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_112_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_112_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_112_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_112_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_113_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_113_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_113_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_113_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_114_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_114_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_114_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_114_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_115_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_115_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_115_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_115_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_116_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_116_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_116_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_116_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_117_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_117_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_117_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_117_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_118_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_118_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_118_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_118_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_119_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_119_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_119_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_119_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_120_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_120_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_120_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_120_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_121_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_121_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_121_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_121_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_122_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_122_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_122_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_122_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_123_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_123_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_123_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_123_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_124_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_124_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_124_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_124_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_125_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_125_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_125_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_125_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_126_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_126_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_126_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_126_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_127_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_127_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_127_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_127_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_128_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_128_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_128_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_128_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_129_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_129_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_129_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_129_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_130_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_130_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_130_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_130_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_131_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_131_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_131_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_131_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_132_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_132_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_132_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_132_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_133_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_133_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_133_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_133_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_134_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_134_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_134_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_134_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_135_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_135_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_135_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_135_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_136_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_136_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_136_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_136_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_137_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_137_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_137_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_137_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_138_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_138_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_138_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_138_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_139_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_139_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_139_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_139_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_140_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_140_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_140_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_140_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_141_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_141_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_141_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_141_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_142_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_142_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_142_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_142_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_143_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_143_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_143_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_143_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_144_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_144_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_144_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_144_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_145_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_145_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_145_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_145_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_146_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_146_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_146_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_146_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_147_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_147_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_147_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_147_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_148_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_148_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_148_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_148_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_149_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_149_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_149_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_149_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_150_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_150_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_150_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_150_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_151_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_151_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_151_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_151_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_152_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_152_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_152_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_152_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_153_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_153_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_153_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_153_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_154_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_154_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_154_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_154_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_155_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_155_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_155_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_155_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_156_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_156_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_156_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_156_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_157_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_157_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_157_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_157_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_158_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_158_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_158_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_158_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_159_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_159_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_159_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_159_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_160_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_160_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_160_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_160_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_161_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_161_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_161_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_161_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_162_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_162_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_162_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_162_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_163_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_163_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_163_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_163_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_164_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_164_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_164_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_164_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_165_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_165_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_165_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_165_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_166_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_166_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_166_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_166_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_167_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_167_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_167_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_167_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_168_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_168_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_168_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_168_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_169_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_169_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_169_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_169_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_170_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_170_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_170_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_170_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_171_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_171_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_171_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_171_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_172_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_172_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_172_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_172_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_173_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_173_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_173_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_173_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_174_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_174_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_174_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_174_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_175_paired_extract.sed raw_data/SRR1246791_1.fastq > extracted/SRR1246791_1_175_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246791_175_paired_extract.sed raw_data/SRR1246791_2.fastq > extracted/SRR1246791_2_175_extracted.fastq &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_1_1_extract.sed raw_data/SRR1246791_1.fastq > extracted_separated/SRR1246791_with_3_repeats_1.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_1_2_extract.sed raw_data/SRR1246791_1.fastq > extracted_separated/SRR1246791_with_3_repeats_2.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_1_3_extract.sed raw_data/SRR1246791_1.fastq > extracted_separated/SRR1246791_with_3_repeats_3.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_1_4_extract.sed raw_data/SRR1246791_1.fastq > extracted_separated/SRR1246791_with_3_repeats_4.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_1_5_extract.sed raw_data/SRR1246791_1.fastq > extracted_separated/SRR1246791_with_3_repeats_5.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_1_6_extract.sed raw_data/SRR1246791_1.fastq > extracted_separated/SRR1246791_with_3_repeats_6.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_1_7_extract.sed raw_data/SRR1246791_1.fastq > extracted_separated/SRR1246791_with_3_repeats_7.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_1_8_extract.sed raw_data/SRR1246791_1.fastq > extracted_separated/SRR1246791_with_3_repeats_8.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_1_9_extract.sed raw_data/SRR1246791_1.fastq > extracted_separated/SRR1246791_with_3_repeats_9.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_1_10_extract.sed raw_data/SRR1246791_1.fastq > extracted_separated/SRR1246791_with_3_repeats_10.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_1_11_extract.sed raw_data/SRR1246791_1.fastq > extracted_separated/SRR1246791_with_3_repeats_11.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_1_12_extract.sed raw_data/SRR1246791_1.fastq > extracted_separated/SRR1246791_with_3_repeats_12.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_1_13_extract.sed raw_data/SRR1246791_1.fastq > extracted_separated/SRR1246791_with_3_repeats_13.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_1_14_extract.sed raw_data/SRR1246791_1.fastq > extracted_separated/SRR1246791_with_3_repeats_14.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_1_15_extract.sed raw_data/SRR1246791_1.fastq > extracted_separated/SRR1246791_with_3_repeats_15.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_1_16_extract.sed raw_data/SRR1246791_1.fastq > extracted_separated/SRR1246791_with_3_repeats_16.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_1_17_extract.sed raw_data/SRR1246791_1.fastq > extracted_separated/SRR1246791_with_3_repeats_17.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_1_18_extract.sed raw_data/SRR1246791_1.fastq > extracted_separated/SRR1246791_with_3_repeats_18.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_1_19_extract.sed raw_data/SRR1246791_1.fastq > extracted_separated/SRR1246791_with_3_repeats_19.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_1_20_extract.sed raw_data/SRR1246791_1.fastq > extracted_separated/SRR1246791_with_3_repeats_20.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_1_21_extract.sed raw_data/SRR1246791_1.fastq > extracted_separated/SRR1246791_with_3_repeats_21.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_1_22_extract.sed raw_data/SRR1246791_1.fastq > extracted_separated/SRR1246791_with_3_repeats_22.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_1_23_extract.sed raw_data/SRR1246791_1.fastq > extracted_separated/SRR1246791_with_3_repeats_23.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_1_24_extract.sed raw_data/SRR1246791_1.fastq > extracted_separated/SRR1246791_with_3_repeats_24.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_1_25_extract.sed raw_data/SRR1246791_1.fastq > extracted_separated/SRR1246791_with_3_repeats_25.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_1_26_extract.sed raw_data/SRR1246791_1.fastq > extracted_separated/SRR1246791_with_3_repeats_26.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_1_27_extract.sed raw_data/SRR1246791_1.fastq > extracted_separated/SRR1246791_with_3_repeats_27.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_1_28_extract.sed raw_data/SRR1246791_1.fastq > extracted_separated/SRR1246791_with_3_repeats_28.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_1_29_extract.sed raw_data/SRR1246791_1.fastq > extracted_separated/SRR1246791_with_3_repeats_29.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_1_30_extract.sed raw_data/SRR1246791_1.fastq > extracted_separated/SRR1246791_with_3_repeats_30.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_1_31_extract.sed raw_data/SRR1246791_1.fastq > extracted_separated/SRR1246791_with_3_repeats_31.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_1_32_extract.sed raw_data/SRR1246791_1.fastq > extracted_separated/SRR1246791_with_3_repeats_32.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_1_33_extract.sed raw_data/SRR1246791_1.fastq > extracted_separated/SRR1246791_with_3_repeats_33.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_1_34_extract.sed raw_data/SRR1246791_1.fastq > extracted_separated/SRR1246791_with_3_repeats_34.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_1_35_extract.sed raw_data/SRR1246791_1.fastq > extracted_separated/SRR1246791_with_3_repeats_35.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_1_36_extract.sed raw_data/SRR1246791_1.fastq > extracted_separated/SRR1246791_with_3_repeats_36.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_1_37_extract.sed raw_data/SRR1246791_1.fastq > extracted_separated/SRR1246791_with_3_repeats_37.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_1_38_extract.sed raw_data/SRR1246791_1.fastq > extracted_separated/SRR1246791_with_3_repeats_38.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_1_39_extract.sed raw_data/SRR1246791_1.fastq > extracted_separated/SRR1246791_with_3_repeats_39.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_1_40_extract.sed raw_data/SRR1246791_1.fastq > extracted_separated/SRR1246791_with_3_repeats_40.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_1_41_extract.sed raw_data/SRR1246791_1.fastq > extracted_separated/SRR1246791_with_3_repeats_41.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_1_42_extract.sed raw_data/SRR1246791_1.fastq > extracted_separated/SRR1246791_with_3_repeats_42.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_1_43_extract.sed raw_data/SRR1246791_1.fastq > extracted_separated/SRR1246791_with_3_repeats_43.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_1_44_extract.sed raw_data/SRR1246791_1.fastq > extracted_separated/SRR1246791_with_3_repeats_44.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_1_45_extract.sed raw_data/SRR1246791_1.fastq > extracted_separated/SRR1246791_with_3_repeats_45.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_1_46_extract.sed raw_data/SRR1246791_1.fastq > extracted_separated/SRR1246791_with_3_repeats_46.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_1_47_extract.sed raw_data/SRR1246791_1.fastq > extracted_separated/SRR1246791_with_3_repeats_47.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_1_48_extract.sed raw_data/SRR1246791_1.fastq > extracted_separated/SRR1246791_with_3_repeats_48.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_1_49_extract.sed raw_data/SRR1246791_1.fastq > extracted_separated/SRR1246791_with_3_repeats_49.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_1_50_extract.sed raw_data/SRR1246791_1.fastq > extracted_separated/SRR1246791_with_3_repeats_50.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_1_51_extract.sed raw_data/SRR1246791_1.fastq > extracted_separated/SRR1246791_with_3_repeats_51.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_1_52_extract.sed raw_data/SRR1246791_1.fastq > extracted_separated/SRR1246791_with_3_repeats_52.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_1_53_extract.sed raw_data/SRR1246791_1.fastq > extracted_separated/SRR1246791_with_3_repeats_53.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_1_54_extract.sed raw_data/SRR1246791_1.fastq > extracted_separated/SRR1246791_with_3_repeats_54.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_1_55_extract.sed raw_data/SRR1246791_1.fastq > extracted_separated/SRR1246791_with_3_repeats_55.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_1_56_extract.sed raw_data/SRR1246791_1.fastq > extracted_separated/SRR1246791_with_3_repeats_56.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_1_57_extract.sed raw_data/SRR1246791_1.fastq > extracted_separated/SRR1246791_with_3_repeats_57.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_1_58_extract.sed raw_data/SRR1246791_1.fastq > extracted_separated/SRR1246791_with_3_repeats_58.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_1_59_extract.sed raw_data/SRR1246791_1.fastq > extracted_separated/SRR1246791_with_3_repeats_59.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_1_60_extract.sed raw_data/SRR1246791_1.fastq > extracted_separated/SRR1246791_with_3_repeats_60.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_1_61_extract.sed raw_data/SRR1246791_1.fastq > extracted_separated/SRR1246791_with_3_repeats_61.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_1_62_extract.sed raw_data/SRR1246791_1.fastq > extracted_separated/SRR1246791_with_3_repeats_62.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_1_63_extract.sed raw_data/SRR1246791_1.fastq > extracted_separated/SRR1246791_with_3_repeats_63.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_1_64_extract.sed raw_data/SRR1246791_1.fastq > extracted_separated/SRR1246791_with_3_repeats_64.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_1_65_extract.sed raw_data/SRR1246791_1.fastq > extracted_separated/SRR1246791_with_3_repeats_65.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_1_66_extract.sed raw_data/SRR1246791_1.fastq > extracted_separated/SRR1246791_with_3_repeats_66.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_1_67_extract.sed raw_data/SRR1246791_1.fastq > extracted_separated/SRR1246791_with_3_repeats_67.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_1_68_extract.sed raw_data/SRR1246791_1.fastq > extracted_separated/SRR1246791_with_3_repeats_68.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_1_69_extract.sed raw_data/SRR1246791_1.fastq > extracted_separated/SRR1246791_with_3_repeats_69.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_1_70_extract.sed raw_data/SRR1246791_1.fastq > extracted_separated/SRR1246791_with_3_repeats_70.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_1_71_extract.sed raw_data/SRR1246791_1.fastq > extracted_separated/SRR1246791_with_3_repeats_71.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_1_72_extract.sed raw_data/SRR1246791_1.fastq > extracted_separated/SRR1246791_with_3_repeats_72.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_1_73_extract.sed raw_data/SRR1246791_1.fastq > extracted_separated/SRR1246791_with_3_repeats_73.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_1_74_extract.sed raw_data/SRR1246791_1.fastq > extracted_separated/SRR1246791_with_3_repeats_74.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_1_75_extract.sed raw_data/SRR1246791_1.fastq > extracted_separated/SRR1246791_with_3_repeats_75.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_1_76_extract.sed raw_data/SRR1246791_1.fastq > extracted_separated/SRR1246791_with_3_repeats_76.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_1_77_extract.sed raw_data/SRR1246791_1.fastq > extracted_separated/SRR1246791_with_3_repeats_77.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_1_78_extract.sed raw_data/SRR1246791_1.fastq > extracted_separated/SRR1246791_with_3_repeats_78.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_1_79_extract.sed raw_data/SRR1246791_1.fastq > extracted_separated/SRR1246791_with_3_repeats_79.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_1_80_extract.sed raw_data/SRR1246791_1.fastq > extracted_separated/SRR1246791_with_3_repeats_80.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_1_81_extract.sed raw_data/SRR1246791_1.fastq > extracted_separated/SRR1246791_with_3_repeats_81.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_1_82_extract.sed raw_data/SRR1246791_1.fastq > extracted_separated/SRR1246791_with_3_repeats_82.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_1_83_extract.sed raw_data/SRR1246791_1.fastq > extracted_separated/SRR1246791_with_3_repeats_83.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_1_84_extract.sed raw_data/SRR1246791_1.fastq > extracted_separated/SRR1246791_with_3_repeats_84.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_1_85_extract.sed raw_data/SRR1246791_1.fastq > extracted_separated/SRR1246791_with_3_repeats_85.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_1_86_extract.sed raw_data/SRR1246791_1.fastq > extracted_separated/SRR1246791_with_3_repeats_86.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_1_87_extract.sed raw_data/SRR1246791_1.fastq > extracted_separated/SRR1246791_with_3_repeats_87.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_1_88_extract.sed raw_data/SRR1246791_1.fastq > extracted_separated/SRR1246791_with_3_repeats_88.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_1_89_extract.sed raw_data/SRR1246791_1.fastq > extracted_separated/SRR1246791_with_3_repeats_89.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_1_90_extract.sed raw_data/SRR1246791_1.fastq > extracted_separated/SRR1246791_with_3_repeats_90.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_1_91_extract.sed raw_data/SRR1246791_1.fastq > extracted_separated/SRR1246791_with_3_repeats_91.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_1_92_extract.sed raw_data/SRR1246791_1.fastq > extracted_separated/SRR1246791_with_3_repeats_92.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_1_93_extract.sed raw_data/SRR1246791_1.fastq > extracted_separated/SRR1246791_with_3_repeats_93.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_1_94_extract.sed raw_data/SRR1246791_1.fastq > extracted_separated/SRR1246791_with_3_repeats_94.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_1_95_extract.sed raw_data/SRR1246791_1.fastq > extracted_separated/SRR1246791_with_3_repeats_95.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_1_96_extract.sed raw_data/SRR1246791_1.fastq > extracted_separated/SRR1246791_with_3_repeats_96.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_1_97_extract.sed raw_data/SRR1246791_1.fastq > extracted_separated/SRR1246791_with_3_repeats_97.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_1_98_extract.sed raw_data/SRR1246791_1.fastq > extracted_separated/SRR1246791_with_3_repeats_98.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_1_99_extract.sed raw_data/SRR1246791_1.fastq > extracted_separated/SRR1246791_with_3_repeats_99.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_1_100_extract.sed raw_data/SRR1246791_1.fastq > extracted_separated/SRR1246791_with_3_repeats_100.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_1_101_extract.sed raw_data/SRR1246791_1.fastq > extracted_separated/SRR1246791_with_3_repeats_101.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_1_102_extract.sed raw_data/SRR1246791_1.fastq > extracted_separated/SRR1246791_with_3_repeats_102.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_1_103_extract.sed raw_data/SRR1246791_1.fastq > extracted_separated/SRR1246791_with_3_repeats_103.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_1_104_extract.sed raw_data/SRR1246791_1.fastq > extracted_separated/SRR1246791_with_3_repeats_104.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_1_105_extract.sed raw_data/SRR1246791_1.fastq > extracted_separated/SRR1246791_with_3_repeats_105.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_1_106_extract.sed raw_data/SRR1246791_1.fastq > extracted_separated/SRR1246791_with_3_repeats_106.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_1_107_extract.sed raw_data/SRR1246791_1.fastq > extracted_separated/SRR1246791_with_3_repeats_107.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_1_108_extract.sed raw_data/SRR1246791_1.fastq > extracted_separated/SRR1246791_with_3_repeats_108.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_1_109_extract.sed raw_data/SRR1246791_1.fastq > extracted_separated/SRR1246791_with_3_repeats_109.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_1_110_extract.sed raw_data/SRR1246791_1.fastq > extracted_separated/SRR1246791_with_3_repeats_110.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_1_111_extract.sed raw_data/SRR1246791_1.fastq > extracted_separated/SRR1246791_with_3_repeats_111.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_2_1_extract.sed raw_data/SRR1246791_2.fastq > extracted_separated/SRR1246791_with_3_repeats_1.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_2_2_extract.sed raw_data/SRR1246791_2.fastq > extracted_separated/SRR1246791_with_3_repeats_2.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_2_3_extract.sed raw_data/SRR1246791_2.fastq > extracted_separated/SRR1246791_with_3_repeats_3.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_2_4_extract.sed raw_data/SRR1246791_2.fastq > extracted_separated/SRR1246791_with_3_repeats_4.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_2_5_extract.sed raw_data/SRR1246791_2.fastq > extracted_separated/SRR1246791_with_3_repeats_5.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_2_6_extract.sed raw_data/SRR1246791_2.fastq > extracted_separated/SRR1246791_with_3_repeats_6.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_2_7_extract.sed raw_data/SRR1246791_2.fastq > extracted_separated/SRR1246791_with_3_repeats_7.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_2_8_extract.sed raw_data/SRR1246791_2.fastq > extracted_separated/SRR1246791_with_3_repeats_8.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_2_9_extract.sed raw_data/SRR1246791_2.fastq > extracted_separated/SRR1246791_with_3_repeats_9.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_2_10_extract.sed raw_data/SRR1246791_2.fastq > extracted_separated/SRR1246791_with_3_repeats_10.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_2_11_extract.sed raw_data/SRR1246791_2.fastq > extracted_separated/SRR1246791_with_3_repeats_11.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_2_12_extract.sed raw_data/SRR1246791_2.fastq > extracted_separated/SRR1246791_with_3_repeats_12.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_2_13_extract.sed raw_data/SRR1246791_2.fastq > extracted_separated/SRR1246791_with_3_repeats_13.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_2_14_extract.sed raw_data/SRR1246791_2.fastq > extracted_separated/SRR1246791_with_3_repeats_14.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_2_15_extract.sed raw_data/SRR1246791_2.fastq > extracted_separated/SRR1246791_with_3_repeats_15.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_2_16_extract.sed raw_data/SRR1246791_2.fastq > extracted_separated/SRR1246791_with_3_repeats_16.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_2_17_extract.sed raw_data/SRR1246791_2.fastq > extracted_separated/SRR1246791_with_3_repeats_17.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_2_18_extract.sed raw_data/SRR1246791_2.fastq > extracted_separated/SRR1246791_with_3_repeats_18.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_2_19_extract.sed raw_data/SRR1246791_2.fastq > extracted_separated/SRR1246791_with_3_repeats_19.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_2_20_extract.sed raw_data/SRR1246791_2.fastq > extracted_separated/SRR1246791_with_3_repeats_20.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_2_21_extract.sed raw_data/SRR1246791_2.fastq > extracted_separated/SRR1246791_with_3_repeats_21.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_2_22_extract.sed raw_data/SRR1246791_2.fastq > extracted_separated/SRR1246791_with_3_repeats_22.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_2_23_extract.sed raw_data/SRR1246791_2.fastq > extracted_separated/SRR1246791_with_3_repeats_23.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_2_24_extract.sed raw_data/SRR1246791_2.fastq > extracted_separated/SRR1246791_with_3_repeats_24.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_2_25_extract.sed raw_data/SRR1246791_2.fastq > extracted_separated/SRR1246791_with_3_repeats_25.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_2_26_extract.sed raw_data/SRR1246791_2.fastq > extracted_separated/SRR1246791_with_3_repeats_26.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_2_27_extract.sed raw_data/SRR1246791_2.fastq > extracted_separated/SRR1246791_with_3_repeats_27.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_2_28_extract.sed raw_data/SRR1246791_2.fastq > extracted_separated/SRR1246791_with_3_repeats_28.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_2_29_extract.sed raw_data/SRR1246791_2.fastq > extracted_separated/SRR1246791_with_3_repeats_29.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_2_30_extract.sed raw_data/SRR1246791_2.fastq > extracted_separated/SRR1246791_with_3_repeats_30.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_2_31_extract.sed raw_data/SRR1246791_2.fastq > extracted_separated/SRR1246791_with_3_repeats_31.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_2_32_extract.sed raw_data/SRR1246791_2.fastq > extracted_separated/SRR1246791_with_3_repeats_32.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_2_33_extract.sed raw_data/SRR1246791_2.fastq > extracted_separated/SRR1246791_with_3_repeats_33.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_2_34_extract.sed raw_data/SRR1246791_2.fastq > extracted_separated/SRR1246791_with_3_repeats_34.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_2_35_extract.sed raw_data/SRR1246791_2.fastq > extracted_separated/SRR1246791_with_3_repeats_35.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_2_36_extract.sed raw_data/SRR1246791_2.fastq > extracted_separated/SRR1246791_with_3_repeats_36.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_2_37_extract.sed raw_data/SRR1246791_2.fastq > extracted_separated/SRR1246791_with_3_repeats_37.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_2_38_extract.sed raw_data/SRR1246791_2.fastq > extracted_separated/SRR1246791_with_3_repeats_38.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_2_39_extract.sed raw_data/SRR1246791_2.fastq > extracted_separated/SRR1246791_with_3_repeats_39.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_2_40_extract.sed raw_data/SRR1246791_2.fastq > extracted_separated/SRR1246791_with_3_repeats_40.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_2_41_extract.sed raw_data/SRR1246791_2.fastq > extracted_separated/SRR1246791_with_3_repeats_41.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_2_42_extract.sed raw_data/SRR1246791_2.fastq > extracted_separated/SRR1246791_with_3_repeats_42.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_2_43_extract.sed raw_data/SRR1246791_2.fastq > extracted_separated/SRR1246791_with_3_repeats_43.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_2_44_extract.sed raw_data/SRR1246791_2.fastq > extracted_separated/SRR1246791_with_3_repeats_44.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_2_45_extract.sed raw_data/SRR1246791_2.fastq > extracted_separated/SRR1246791_with_3_repeats_45.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_2_46_extract.sed raw_data/SRR1246791_2.fastq > extracted_separated/SRR1246791_with_3_repeats_46.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_2_47_extract.sed raw_data/SRR1246791_2.fastq > extracted_separated/SRR1246791_with_3_repeats_47.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_2_48_extract.sed raw_data/SRR1246791_2.fastq > extracted_separated/SRR1246791_with_3_repeats_48.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_2_49_extract.sed raw_data/SRR1246791_2.fastq > extracted_separated/SRR1246791_with_3_repeats_49.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_2_50_extract.sed raw_data/SRR1246791_2.fastq > extracted_separated/SRR1246791_with_3_repeats_50.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_2_51_extract.sed raw_data/SRR1246791_2.fastq > extracted_separated/SRR1246791_with_3_repeats_51.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_2_52_extract.sed raw_data/SRR1246791_2.fastq > extracted_separated/SRR1246791_with_3_repeats_52.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_2_53_extract.sed raw_data/SRR1246791_2.fastq > extracted_separated/SRR1246791_with_3_repeats_53.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_2_54_extract.sed raw_data/SRR1246791_2.fastq > extracted_separated/SRR1246791_with_3_repeats_54.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_2_55_extract.sed raw_data/SRR1246791_2.fastq > extracted_separated/SRR1246791_with_3_repeats_55.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_2_56_extract.sed raw_data/SRR1246791_2.fastq > extracted_separated/SRR1246791_with_3_repeats_56.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_2_57_extract.sed raw_data/SRR1246791_2.fastq > extracted_separated/SRR1246791_with_3_repeats_57.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_2_58_extract.sed raw_data/SRR1246791_2.fastq > extracted_separated/SRR1246791_with_3_repeats_58.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_2_59_extract.sed raw_data/SRR1246791_2.fastq > extracted_separated/SRR1246791_with_3_repeats_59.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_2_60_extract.sed raw_data/SRR1246791_2.fastq > extracted_separated/SRR1246791_with_3_repeats_60.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_2_61_extract.sed raw_data/SRR1246791_2.fastq > extracted_separated/SRR1246791_with_3_repeats_61.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_2_62_extract.sed raw_data/SRR1246791_2.fastq > extracted_separated/SRR1246791_with_3_repeats_62.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_2_63_extract.sed raw_data/SRR1246791_2.fastq > extracted_separated/SRR1246791_with_3_repeats_63.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_2_64_extract.sed raw_data/SRR1246791_2.fastq > extracted_separated/SRR1246791_with_3_repeats_64.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_2_65_extract.sed raw_data/SRR1246791_2.fastq > extracted_separated/SRR1246791_with_3_repeats_65.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_2_66_extract.sed raw_data/SRR1246791_2.fastq > extracted_separated/SRR1246791_with_3_repeats_66.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_2_67_extract.sed raw_data/SRR1246791_2.fastq > extracted_separated/SRR1246791_with_3_repeats_67.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_2_68_extract.sed raw_data/SRR1246791_2.fastq > extracted_separated/SRR1246791_with_3_repeats_68.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_2_69_extract.sed raw_data/SRR1246791_2.fastq > extracted_separated/SRR1246791_with_3_repeats_69.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_2_70_extract.sed raw_data/SRR1246791_2.fastq > extracted_separated/SRR1246791_with_3_repeats_70.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_2_71_extract.sed raw_data/SRR1246791_2.fastq > extracted_separated/SRR1246791_with_3_repeats_71.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_2_72_extract.sed raw_data/SRR1246791_2.fastq > extracted_separated/SRR1246791_with_3_repeats_72.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_2_73_extract.sed raw_data/SRR1246791_2.fastq > extracted_separated/SRR1246791_with_3_repeats_73.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_2_74_extract.sed raw_data/SRR1246791_2.fastq > extracted_separated/SRR1246791_with_3_repeats_74.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_2_75_extract.sed raw_data/SRR1246791_2.fastq > extracted_separated/SRR1246791_with_3_repeats_75.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_2_76_extract.sed raw_data/SRR1246791_2.fastq > extracted_separated/SRR1246791_with_3_repeats_76.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_2_77_extract.sed raw_data/SRR1246791_2.fastq > extracted_separated/SRR1246791_with_3_repeats_77.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_2_78_extract.sed raw_data/SRR1246791_2.fastq > extracted_separated/SRR1246791_with_3_repeats_78.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_2_79_extract.sed raw_data/SRR1246791_2.fastq > extracted_separated/SRR1246791_with_3_repeats_79.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_2_80_extract.sed raw_data/SRR1246791_2.fastq > extracted_separated/SRR1246791_with_3_repeats_80.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_2_81_extract.sed raw_data/SRR1246791_2.fastq > extracted_separated/SRR1246791_with_3_repeats_81.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_2_82_extract.sed raw_data/SRR1246791_2.fastq > extracted_separated/SRR1246791_with_3_repeats_82.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_2_83_extract.sed raw_data/SRR1246791_2.fastq > extracted_separated/SRR1246791_with_3_repeats_83.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_2_84_extract.sed raw_data/SRR1246791_2.fastq > extracted_separated/SRR1246791_with_3_repeats_84.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_2_85_extract.sed raw_data/SRR1246791_2.fastq > extracted_separated/SRR1246791_with_3_repeats_85.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_2_86_extract.sed raw_data/SRR1246791_2.fastq > extracted_separated/SRR1246791_with_3_repeats_86.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_2_87_extract.sed raw_data/SRR1246791_2.fastq > extracted_separated/SRR1246791_with_3_repeats_87.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_2_88_extract.sed raw_data/SRR1246791_2.fastq > extracted_separated/SRR1246791_with_3_repeats_88.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_2_89_extract.sed raw_data/SRR1246791_2.fastq > extracted_separated/SRR1246791_with_3_repeats_89.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_2_90_extract.sed raw_data/SRR1246791_2.fastq > extracted_separated/SRR1246791_with_3_repeats_90.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_2_91_extract.sed raw_data/SRR1246791_2.fastq > extracted_separated/SRR1246791_with_3_repeats_91.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_2_92_extract.sed raw_data/SRR1246791_2.fastq > extracted_separated/SRR1246791_with_3_repeats_92.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_2_93_extract.sed raw_data/SRR1246791_2.fastq > extracted_separated/SRR1246791_with_3_repeats_93.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_2_94_extract.sed raw_data/SRR1246791_2.fastq > extracted_separated/SRR1246791_with_3_repeats_94.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_2_95_extract.sed raw_data/SRR1246791_2.fastq > extracted_separated/SRR1246791_with_3_repeats_95.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_2_96_extract.sed raw_data/SRR1246791_2.fastq > extracted_separated/SRR1246791_with_3_repeats_96.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_2_97_extract.sed raw_data/SRR1246791_2.fastq > extracted_separated/SRR1246791_with_3_repeats_97.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_2_98_extract.sed raw_data/SRR1246791_2.fastq > extracted_separated/SRR1246791_with_3_repeats_98.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_2_99_extract.sed raw_data/SRR1246791_2.fastq > extracted_separated/SRR1246791_with_3_repeats_99.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_2_100_extract.sed raw_data/SRR1246791_2.fastq > extracted_separated/SRR1246791_with_3_repeats_100.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_2_101_extract.sed raw_data/SRR1246791_2.fastq > extracted_separated/SRR1246791_with_3_repeats_101.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_2_102_extract.sed raw_data/SRR1246791_2.fastq > extracted_separated/SRR1246791_with_3_repeats_102.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_2_103_extract.sed raw_data/SRR1246791_2.fastq > extracted_separated/SRR1246791_with_3_repeats_103.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_2_104_extract.sed raw_data/SRR1246791_2.fastq > extracted_separated/SRR1246791_with_3_repeats_104.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_2_105_extract.sed raw_data/SRR1246791_2.fastq > extracted_separated/SRR1246791_with_3_repeats_105.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_2_106_extract.sed raw_data/SRR1246791_2.fastq > extracted_separated/SRR1246791_with_3_repeats_106.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_2_107_extract.sed raw_data/SRR1246791_2.fastq > extracted_separated/SRR1246791_with_3_repeats_107.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246791_with_3_repeats_2_108_extract.sed raw_data/SRR1246791_2.fastq > extracted_separated/SRR1246791_with_3_repeats_108.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246791_without_3_repeats_1_extract.sed raw_data/SRR1246791_1.fastq > extracted_separated/SRR1246791_without_3_repeats.fastq ;
                                            sed -n -f scripts/sed_scripts/SRR1246791_without_3_repeats_2_extract.sed raw_data/SRR1246791_2.fastq >> extracted_separated/SRR1246791_without_3_repeats.fastq" &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_1_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_1_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_1_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_1_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_2_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_2_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_2_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_2_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_3_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_3_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_3_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_3_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_4_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_4_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_4_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_4_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_5_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_5_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_5_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_5_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_6_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_6_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_6_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_6_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_7_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_7_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_7_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_7_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_8_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_8_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_8_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_8_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_9_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_9_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_9_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_9_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_10_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_10_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_10_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_10_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_11_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_11_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_11_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_11_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_12_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_12_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_12_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_12_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_13_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_13_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_13_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_13_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_14_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_14_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_14_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_14_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_15_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_15_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_15_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_15_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_16_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_16_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_16_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_16_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_17_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_17_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_17_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_17_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_18_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_18_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_18_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_18_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_19_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_19_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_19_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_19_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_20_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_20_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_20_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_20_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_21_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_21_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_21_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_21_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_22_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_22_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_22_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_22_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_23_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_23_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_23_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_23_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_24_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_24_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_24_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_24_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_25_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_25_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_25_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_25_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_26_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_26_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_26_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_26_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_27_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_27_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_27_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_27_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_28_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_28_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_28_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_28_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_29_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_29_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_29_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_29_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_30_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_30_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_30_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_30_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_31_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_31_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_31_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_31_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_32_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_32_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_32_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_32_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_33_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_33_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_33_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_33_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_34_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_34_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_34_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_34_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_35_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_35_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_35_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_35_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_36_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_36_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_36_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_36_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_37_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_37_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_37_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_37_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_38_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_38_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_38_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_38_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_39_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_39_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_39_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_39_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_40_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_40_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_40_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_40_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_41_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_41_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_41_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_41_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_42_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_42_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_42_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_42_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_43_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_43_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_43_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_43_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_44_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_44_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_44_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_44_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_45_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_45_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_45_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_45_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_46_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_46_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_46_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_46_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_47_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_47_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_47_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_47_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_48_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_48_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_48_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_48_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_49_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_49_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_49_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_49_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_50_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_50_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_50_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_50_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_51_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_51_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_51_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_51_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_52_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_52_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_52_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_52_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_53_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_53_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_53_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_53_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_54_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_54_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_54_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_54_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_55_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_55_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_55_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_55_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_56_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_56_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_56_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_56_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_57_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_57_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_57_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_57_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_58_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_58_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_58_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_58_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_59_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_59_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_59_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_59_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_60_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_60_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_60_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_60_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_61_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_61_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_61_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_61_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_62_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_62_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_62_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_62_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_63_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_63_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_63_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_63_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_64_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_64_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_64_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_64_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_65_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_65_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_65_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_65_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_66_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_66_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_66_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_66_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_67_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_67_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_67_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_67_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_68_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_68_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_68_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_68_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_69_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_69_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_69_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_69_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_70_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_70_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_70_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_70_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_71_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_71_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_71_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_71_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_72_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_72_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_72_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_72_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_73_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_73_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_73_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_73_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_74_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_74_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_74_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_74_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_75_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_75_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_75_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_75_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_76_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_76_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_76_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_76_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_77_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_77_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_77_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_77_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_78_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_78_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_78_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_78_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_79_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_79_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_79_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_79_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_80_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_80_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_80_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_80_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_81_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_81_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_81_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_81_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_82_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_82_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_82_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_82_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_83_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_83_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_83_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_83_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_84_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_84_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_84_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_84_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_85_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_85_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_85_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_85_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_86_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_86_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_86_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_86_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_87_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_87_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_87_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_87_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_88_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_88_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_88_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_88_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_89_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_89_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_89_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_89_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_90_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_90_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_90_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_90_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_91_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_91_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_91_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_91_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_92_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_92_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_92_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_92_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_93_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_93_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_93_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_93_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_94_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_94_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_94_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_94_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_95_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_95_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_95_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_95_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_96_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_96_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_96_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_96_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_97_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_97_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_97_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_97_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_98_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_98_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_98_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_98_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_99_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_99_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_99_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_99_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_100_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_100_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_100_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_100_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_101_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_101_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_101_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_101_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_102_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_102_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_102_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_102_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_103_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_103_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_103_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_103_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_104_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_104_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_104_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_104_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_105_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_105_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_105_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_105_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_106_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_106_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_106_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_106_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_107_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_107_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_107_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_107_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_108_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_108_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_108_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_108_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_109_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_109_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_109_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_109_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_110_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_110_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_110_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_110_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_111_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_111_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_111_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_111_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_112_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_112_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_112_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_112_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_113_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_113_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_113_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_113_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_114_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_114_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_114_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_114_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_115_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_115_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_115_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_115_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_116_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_116_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_116_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_116_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_117_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_117_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_117_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_117_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_118_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_118_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_118_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_118_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_119_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_119_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_119_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_119_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_120_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_120_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_120_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_120_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_121_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_121_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_121_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_121_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_122_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_122_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_122_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_122_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_123_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_123_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_123_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_123_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_124_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_124_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_124_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_124_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_125_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_125_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_125_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_125_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_126_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_126_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_126_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_126_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_127_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_127_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_127_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_127_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_128_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_128_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_128_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_128_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_129_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_129_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_129_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_129_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_130_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_130_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_130_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_130_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_131_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_131_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_131_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_131_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_132_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_132_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_132_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_132_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_133_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_133_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_133_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_133_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_134_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_134_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_134_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_134_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_135_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_135_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_135_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_135_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_136_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_136_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_136_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_136_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_137_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_137_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_137_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_137_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_138_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_138_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_138_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_138_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_139_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_139_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_139_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_139_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_140_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_140_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_140_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_140_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_141_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_141_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_141_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_141_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_142_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_142_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_142_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_142_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_143_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_143_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_143_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_143_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_144_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_144_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_144_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_144_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_145_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_145_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_145_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_145_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_146_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_146_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_146_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_146_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_147_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_147_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_147_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_147_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_148_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_148_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_148_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_148_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_149_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_149_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_149_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_149_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_150_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_150_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_150_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_150_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_151_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_151_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_151_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_151_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_152_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_152_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_152_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_152_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_153_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_153_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_153_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_153_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_154_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_154_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_154_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_154_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_155_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_155_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_155_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_155_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_156_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_156_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_156_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_156_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_157_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_157_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_157_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_157_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_158_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_158_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_158_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_158_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_159_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_159_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_159_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_159_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_160_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_160_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_160_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_160_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_161_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_161_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_161_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_161_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_162_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_162_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_162_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_162_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_163_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_163_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_163_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_163_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_164_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_164_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_164_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_164_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_165_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_165_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_165_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_165_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_166_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_166_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_166_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_166_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_167_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_167_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_167_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_167_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_168_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_168_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_168_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_168_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_169_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_169_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_169_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_169_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_170_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_170_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_170_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_170_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_171_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_171_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_171_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_171_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_172_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_172_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_172_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_172_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_173_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_173_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_173_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_173_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_174_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_174_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_174_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_174_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_175_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_175_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_175_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_175_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_176_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_176_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_176_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_176_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_177_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_177_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_177_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_177_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_178_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_178_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_178_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_178_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_179_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_179_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_179_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_179_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_180_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_180_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_180_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_180_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_181_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_181_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_181_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_181_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_182_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_182_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_182_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_182_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_183_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_183_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_183_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_183_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_184_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_184_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_184_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_184_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_185_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_185_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_185_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_185_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_186_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_186_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_186_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_186_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_187_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_187_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_187_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_187_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_188_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_188_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_188_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_188_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_189_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_189_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_189_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_189_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_190_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_190_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_190_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_190_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_191_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_191_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_191_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_191_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_192_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_192_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_192_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_192_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_193_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_193_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_193_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_193_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_194_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_194_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_194_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_194_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_195_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_195_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_195_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_195_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_196_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_196_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_196_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_196_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_197_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_197_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_197_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_197_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_198_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_198_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_198_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_198_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_199_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_199_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_199_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_199_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_200_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_200_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_200_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_200_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_201_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_201_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_201_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_201_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_202_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_202_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_202_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_202_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_203_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_203_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_203_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_203_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_204_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_204_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_204_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_204_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_205_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_205_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_205_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_205_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_206_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_206_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_206_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_206_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_207_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_207_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_207_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_207_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_208_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_208_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_208_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_208_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_209_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_209_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_209_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_209_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_210_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_210_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_210_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_210_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_211_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_211_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_211_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_211_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_212_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_212_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_212_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_212_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_213_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_213_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_213_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_213_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_214_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_214_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_214_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_214_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_215_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_215_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_215_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_215_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_216_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_216_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_216_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_216_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_217_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_217_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_217_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_217_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_218_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_218_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_218_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_218_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_219_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_219_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_219_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_219_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_220_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_220_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_220_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_220_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_221_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_221_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_221_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_221_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_222_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_222_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_222_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_222_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_223_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_223_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_223_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_223_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_224_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_224_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_224_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_224_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_225_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_225_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_225_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_225_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_226_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_226_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_226_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_226_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_227_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_227_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_227_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_227_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_228_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_228_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_228_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_228_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_229_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_229_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_229_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_229_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_230_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_230_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_230_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_230_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_231_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_231_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_231_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_231_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_232_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_232_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_232_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_232_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_233_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_233_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_233_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_233_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_234_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_234_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_234_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_234_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_235_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_235_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_235_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_235_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_236_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_236_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_236_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_236_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_237_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_237_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_237_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_237_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_238_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_238_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_238_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_238_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_239_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_239_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_239_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_239_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_240_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_240_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_240_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_240_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_241_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_241_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_241_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_241_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_242_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_242_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_242_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_242_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_243_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_243_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_243_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_243_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_244_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_244_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_244_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_244_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_245_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_245_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_245_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_245_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_246_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_246_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_246_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_246_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_247_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_247_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_247_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_247_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_248_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_248_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_248_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_248_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_249_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_249_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_249_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_249_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_250_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_250_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_250_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_250_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_251_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_251_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_251_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_251_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_252_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_252_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_252_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_252_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_253_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_253_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_253_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_253_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_254_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_254_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_254_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_254_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_255_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_255_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_255_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_255_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_256_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_256_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_256_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_256_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_257_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_257_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_257_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_257_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_258_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_258_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_258_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_258_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_259_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_259_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_259_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_259_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_260_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_260_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_260_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_260_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_261_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_261_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_261_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_261_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_262_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_262_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_262_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_262_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_263_paired_extract.sed raw_data/SRR1246792_1.fastq > extracted/SRR1246792_1_263_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR1246792_263_paired_extract.sed raw_data/SRR1246792_2.fastq > extracted/SRR1246792_2_263_extracted.fastq &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_1_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_1.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_2_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_2.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_3_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_3.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_4_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_4.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_5_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_5.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_6_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_6.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_7_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_7.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_8_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_8.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_9_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_9.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_10_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_10.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_11_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_11.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_12_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_12.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_13_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_13.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_14_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_14.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_15_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_15.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_16_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_16.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_17_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_17.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_18_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_18.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_19_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_19.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_20_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_20.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_21_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_21.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_22_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_22.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_23_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_23.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_24_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_24.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_25_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_25.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_26_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_26.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_27_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_27.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_28_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_28.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_29_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_29.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_30_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_30.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_31_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_31.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_32_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_32.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_33_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_33.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_34_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_34.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_35_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_35.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_36_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_36.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_37_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_37.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_38_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_38.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_39_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_39.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_40_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_40.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_41_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_41.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_42_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_42.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_43_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_43.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_44_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_44.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_45_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_45.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_46_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_46.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_47_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_47.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_48_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_48.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_49_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_49.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_50_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_50.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_51_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_51.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_52_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_52.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_53_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_53.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_54_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_54.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_55_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_55.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_56_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_56.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_57_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_57.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_58_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_58.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_59_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_59.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_60_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_60.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_61_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_61.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_62_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_62.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_63_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_63.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_64_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_64.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_65_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_65.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_66_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_66.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_67_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_67.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_68_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_68.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_69_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_69.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_70_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_70.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_71_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_71.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_72_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_72.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_73_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_73.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_74_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_74.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_75_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_75.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_76_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_76.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_77_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_77.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_78_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_78.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_79_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_79.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_80_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_80.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_81_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_81.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_82_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_82.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_83_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_83.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_84_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_84.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_85_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_85.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_86_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_86.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_87_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_87.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_88_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_88.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_89_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_89.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_90_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_90.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_91_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_91.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_92_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_92.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_93_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_93.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_94_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_94.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_95_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_95.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_96_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_96.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_97_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_97.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_98_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_98.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_99_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_99.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_100_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_100.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_101_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_101.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_102_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_102.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_103_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_103.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_104_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_104.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_105_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_105.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_106_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_106.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_107_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_107.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_108_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_108.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_109_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_109.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_110_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_110.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_111_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_111.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_112_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_112.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_113_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_113.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_114_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_114.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_115_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_115.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_116_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_116.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_117_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_117.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_118_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_118.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_119_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_119.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_120_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_120.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_121_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_121.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_122_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_122.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_123_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_123.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_124_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_124.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_125_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_125.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_126_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_126.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_127_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_127.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_128_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_128.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_129_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_129.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_130_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_130.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_131_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_131.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_132_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_132.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_133_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_133.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_134_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_134.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_135_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_135.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_136_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_136.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_137_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_137.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_138_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_138.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_139_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_139.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_140_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_140.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_141_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_141.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_142_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_142.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_143_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_143.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_144_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_144.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_145_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_145.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_146_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_146.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_147_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_147.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_148_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_148.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_149_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_149.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_150_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_150.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_151_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_151.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_152_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_152.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_153_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_153.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_154_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_154.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_155_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_155.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_156_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_156.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_157_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_157.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_158_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_158.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_159_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_159.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_160_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_160.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_161_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_161.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_162_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_162.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_163_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_163.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_164_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_164.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_165_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_165.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_166_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_166.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_1_167_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_with_3_repeats_167.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_1_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_1.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_2_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_2.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_3_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_3.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_4_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_4.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_5_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_5.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_6_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_6.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_7_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_7.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_8_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_8.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_9_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_9.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_10_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_10.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_11_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_11.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_12_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_12.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_13_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_13.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_14_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_14.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_15_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_15.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_16_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_16.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_17_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_17.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_18_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_18.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_19_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_19.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_20_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_20.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_21_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_21.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_22_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_22.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_23_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_23.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_24_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_24.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_25_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_25.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_26_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_26.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_27_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_27.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_28_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_28.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_29_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_29.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_30_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_30.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_31_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_31.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_32_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_32.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_33_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_33.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_34_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_34.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_35_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_35.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_36_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_36.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_37_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_37.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_38_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_38.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_39_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_39.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_40_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_40.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_41_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_41.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_42_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_42.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_43_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_43.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_44_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_44.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_45_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_45.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_46_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_46.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_47_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_47.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_48_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_48.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_49_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_49.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_50_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_50.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_51_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_51.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_52_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_52.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_53_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_53.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_54_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_54.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_55_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_55.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_56_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_56.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_57_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_57.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_58_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_58.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_59_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_59.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_60_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_60.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_61_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_61.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_62_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_62.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_63_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_63.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_64_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_64.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_65_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_65.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_66_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_66.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_67_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_67.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_68_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_68.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_69_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_69.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_70_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_70.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_71_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_71.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_72_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_72.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_73_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_73.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_74_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_74.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_75_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_75.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_76_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_76.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_77_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_77.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_78_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_78.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_79_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_79.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_80_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_80.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_81_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_81.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_82_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_82.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_83_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_83.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_84_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_84.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_85_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_85.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_86_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_86.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_87_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_87.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_88_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_88.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_89_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_89.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_90_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_90.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_91_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_91.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_92_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_92.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_93_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_93.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_94_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_94.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_95_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_95.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_96_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_96.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_97_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_97.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_98_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_98.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_99_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_99.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_100_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_100.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_101_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_101.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_102_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_102.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_103_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_103.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_104_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_104.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_105_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_105.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_106_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_106.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_107_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_107.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_108_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_108.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_109_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_109.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_110_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_110.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_111_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_111.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_112_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_112.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_113_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_113.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_114_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_114.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_115_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_115.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_116_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_116.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_117_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_117.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_118_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_118.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_119_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_119.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_120_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_120.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_121_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_121.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_122_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_122.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_123_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_123.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_124_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_124.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_125_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_125.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_126_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_126.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_127_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_127.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_128_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_128.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_129_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_129.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_130_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_130.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_131_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_131.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_132_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_132.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_133_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_133.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_134_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_134.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_135_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_135.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_136_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_136.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_137_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_137.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_138_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_138.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_139_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_139.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_140_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_140.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_141_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_141.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_142_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_142.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_143_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_143.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_144_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_144.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_145_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_145.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_146_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_146.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_147_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_147.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_148_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_148.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_149_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_149.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_150_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_150.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_151_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_151.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_152_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_152.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_153_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_153.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_154_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_154.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_155_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_155.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR1246792_with_3_repeats_2_156_extract.sed raw_data/SRR1246792_2.fastq > extracted_separated/SRR1246792_with_3_repeats_156.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR1246792_without_3_repeats_1_extract.sed raw_data/SRR1246792_1.fastq > extracted_separated/SRR1246792_without_3_repeats.fastq ;
                                            sed -n -f scripts/sed_scripts/SRR1246792_without_3_repeats_2_extract.sed raw_data/SRR1246792_2.fastq >> extracted_separated/SRR1246792_without_3_repeats.fastq" &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR3637136_1_paired_extract.sed raw_data/SRR3637136_1.fastq > extracted/SRR3637136_1_1_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR3637136_1_paired_extract.sed raw_data/SRR3637136_2.fastq > extracted/SRR3637136_2_1_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR3637136_2_paired_extract.sed raw_data/SRR3637136_1.fastq > extracted/SRR3637136_1_2_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR3637136_2_paired_extract.sed raw_data/SRR3637136_2.fastq > extracted/SRR3637136_2_2_extracted.fastq &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR3637136_with_3_repeats_1_1_extract.sed raw_data/SRR3637136_1.fastq > extracted_separated/SRR3637136_with_3_repeats_1.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR3637136_with_3_repeats_2_1_extract.sed raw_data/SRR3637136_2.fastq > extracted_separated/SRR3637136_with_3_repeats_1.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR3637136_without_3_repeats_1_extract.sed raw_data/SRR3637136_1.fastq > extracted_separated/SRR3637136_without_3_repeats.fastq ;
                                            sed -n -f scripts/sed_scripts/SRR3637136_without_3_repeats_2_extract.sed raw_data/SRR3637136_2.fastq >> extracted_separated/SRR3637136_without_3_repeats.fastq" &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR3637137_1_paired_extract.sed raw_data/SRR3637137_1.fastq > extracted/SRR3637137_1_1_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR3637137_1_paired_extract.sed raw_data/SRR3637137_2.fastq > extracted/SRR3637137_2_1_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR3637137_2_paired_extract.sed raw_data/SRR3637137_1.fastq > extracted/SRR3637137_1_2_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR3637137_2_paired_extract.sed raw_data/SRR3637137_2.fastq > extracted/SRR3637137_2_2_extracted.fastq &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR3637137_with_3_repeats_1_1_extract.sed raw_data/SRR3637137_1.fastq > extracted_separated/SRR3637137_with_3_repeats_1.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR3637137_with_3_repeats_1_2_extract.sed raw_data/SRR3637137_1.fastq > extracted_separated/SRR3637137_with_3_repeats_2.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR3637137_with_3_repeats_2_1_extract.sed raw_data/SRR3637137_2.fastq > extracted_separated/SRR3637137_with_3_repeats_1.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR3637137_with_3_repeats_2_2_extract.sed raw_data/SRR3637137_2.fastq > extracted_separated/SRR3637137_with_3_repeats_2.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR3637137_without_3_repeats_1_extract.sed raw_data/SRR3637137_1.fastq > extracted_separated/SRR3637137_without_3_repeats.fastq ;
                                            sed -n -f scripts/sed_scripts/SRR3637137_without_3_repeats_2_extract.sed raw_data/SRR3637137_2.fastq >> extracted_separated/SRR3637137_without_3_repeats.fastq" &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR3745014_1_paired_extract.sed raw_data/SRR3745014_1.fastq > extracted/SRR3745014_1_1_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR3745014_1_paired_extract.sed raw_data/SRR3745014_2.fastq > extracted/SRR3745014_2_1_extracted.fastq &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR3745014_with_3_repeats_1_1_extract.sed raw_data/SRR3745014_1.fastq > extracted_separated/SRR3745014_with_3_repeats_1.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR3745014_with_3_repeats_2_1_extract.sed raw_data/SRR3745014_2.fastq > extracted_separated/SRR3745014_with_3_repeats_1.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR3745014_without_3_repeats_1_extract.sed raw_data/SRR3745014_1.fastq > extracted_separated/SRR3745014_without_3_repeats.fastq ;
                                            sed -n -f scripts/sed_scripts/SRR3745014_without_3_repeats_2_extract.sed raw_data/SRR3745014_2.fastq >> extracted_separated/SRR3745014_without_3_repeats.fastq" &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR3745015_1_paired_extract.sed raw_data/SRR3745015_1.fastq > extracted/SRR3745015_1_1_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR3745015_1_paired_extract.sed raw_data/SRR3745015_2.fastq > extracted/SRR3745015_2_1_extracted.fastq &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR3745015_with_3_repeats_1_1_extract.sed raw_data/SRR3745015_1.fastq > extracted_separated/SRR3745015_with_3_repeats_1.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR3745015_with_3_repeats_2_1_extract.sed raw_data/SRR3745015_2.fastq > extracted_separated/SRR3745015_with_3_repeats_1.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR3745015_without_3_repeats_1_extract.sed raw_data/SRR3745015_1.fastq > extracted_separated/SRR3745015_without_3_repeats.fastq ;
                                            sed -n -f scripts/sed_scripts/SRR3745015_without_3_repeats_2_extract.sed raw_data/SRR3745015_2.fastq >> extracted_separated/SRR3745015_without_3_repeats.fastq" &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR3745016_1_paired_extract.sed raw_data/SRR3745016_1.fastq > extracted/SRR3745016_1_1_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR3745016_1_paired_extract.sed raw_data/SRR3745016_2.fastq > extracted/SRR3745016_2_1_extracted.fastq &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR3745016_with_3_repeats_1_1_extract.sed raw_data/SRR3745016_1.fastq > extracted_separated/SRR3745016_with_3_repeats_1.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR3745016_with_3_repeats_2_1_extract.sed raw_data/SRR3745016_2.fastq > extracted_separated/SRR3745016_with_3_repeats_1.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR3745016_without_3_repeats_1_extract.sed raw_data/SRR3745016_1.fastq > extracted_separated/SRR3745016_without_3_repeats.fastq ;
                                            sed -n -f scripts/sed_scripts/SRR3745016_without_3_repeats_2_extract.sed raw_data/SRR3745016_2.fastq >> extracted_separated/SRR3745016_without_3_repeats.fastq" &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR5015865_1_paired_extract.sed raw_data/SRR5015865_1.fastq > extracted/SRR5015865_1_1_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR5015865_1_paired_extract.sed raw_data/SRR5015865_2.fastq > extracted/SRR5015865_2_1_extracted.fastq &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR5015865_with_3_repeats_1_1_extract.sed raw_data/SRR5015865_1.fastq > extracted_separated/SRR5015865_with_3_repeats_1.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR5015865_with_3_repeats_2_1_extract.sed raw_data/SRR5015865_2.fastq > extracted_separated/SRR5015865_with_3_repeats_1.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR5015865_without_3_repeats_1_extract.sed raw_data/SRR5015865_1.fastq > extracted_separated/SRR5015865_without_3_repeats.fastq ;
                                            sed -n -f scripts/sed_scripts/SRR5015865_without_3_repeats_2_extract.sed raw_data/SRR5015865_2.fastq >> extracted_separated/SRR5015865_without_3_repeats.fastq" &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR5015866_1_paired_extract.sed raw_data/SRR5015866_1.fastq > extracted/SRR5015866_1_1_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR5015866_1_paired_extract.sed raw_data/SRR5015866_2.fastq > extracted/SRR5015866_2_1_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR5015866_2_paired_extract.sed raw_data/SRR5015866_1.fastq > extracted/SRR5015866_1_2_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR5015866_2_paired_extract.sed raw_data/SRR5015866_2.fastq > extracted/SRR5015866_2_2_extracted.fastq &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR5015866_with_3_repeats_1_1_extract.sed raw_data/SRR5015866_1.fastq > extracted_separated/SRR5015866_with_3_repeats_1.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR5015866_with_3_repeats_1_2_extract.sed raw_data/SRR5015866_1.fastq > extracted_separated/SRR5015866_with_3_repeats_2.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR5015866_with_3_repeats_2_1_extract.sed raw_data/SRR5015866_2.fastq > extracted_separated/SRR5015866_with_3_repeats_1.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR5015866_with_3_repeats_2_2_extract.sed raw_data/SRR5015866_2.fastq > extracted_separated/SRR5015866_with_3_repeats_2.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR5015866_without_3_repeats_1_extract.sed raw_data/SRR5015866_1.fastq > extracted_separated/SRR5015866_without_3_repeats.fastq ;
                                            sed -n -f scripts/sed_scripts/SRR5015866_without_3_repeats_2_extract.sed raw_data/SRR5015866_2.fastq >> extracted_separated/SRR5015866_without_3_repeats.fastq" &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR5015867_1_paired_extract.sed raw_data/SRR5015867_1.fastq > extracted/SRR5015867_1_1_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR5015867_1_paired_extract.sed raw_data/SRR5015867_2.fastq > extracted/SRR5015867_2_1_extracted.fastq &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR5015867_with_3_repeats_1_1_extract.sed raw_data/SRR5015867_1.fastq > extracted_separated/SRR5015867_with_3_repeats_1.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR5015867_with_3_repeats_2_1_extract.sed raw_data/SRR5015867_2.fastq > extracted_separated/SRR5015867_with_3_repeats_1.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR5015867_without_3_repeats_1_extract.sed raw_data/SRR5015867_1.fastq > extracted_separated/SRR5015867_without_3_repeats.fastq ;
                                            sed -n -f scripts/sed_scripts/SRR5015867_without_3_repeats_2_extract.sed raw_data/SRR5015867_2.fastq >> extracted_separated/SRR5015867_without_3_repeats.fastq" &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR5015868_1_paired_extract.sed raw_data/SRR5015868_1.fastq > extracted/SRR5015868_1_1_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR5015868_1_paired_extract.sed raw_data/SRR5015868_2.fastq > extracted/SRR5015868_2_1_extracted.fastq &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR5015868_with_3_repeats_1_1_extract.sed raw_data/SRR5015868_1.fastq > extracted_separated/SRR5015868_with_3_repeats_1.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR5015868_with_3_repeats_2_1_extract.sed raw_data/SRR5015868_2.fastq > extracted_separated/SRR5015868_with_3_repeats_1.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR5015868_without_3_repeats_1_extract.sed raw_data/SRR5015868_1.fastq > extracted_separated/SRR5015868_without_3_repeats.fastq ;
                                            sed -n -f scripts/sed_scripts/SRR5015868_without_3_repeats_2_extract.sed raw_data/SRR5015868_2.fastq >> extracted_separated/SRR5015868_without_3_repeats.fastq" &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR5015870_1_paired_extract.sed raw_data/SRR5015870_1.fastq > extracted/SRR5015870_1_1_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR5015870_1_paired_extract.sed raw_data/SRR5015870_2.fastq > extracted/SRR5015870_2_1_extracted.fastq &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR5015870_with_3_repeats_1_1_extract.sed raw_data/SRR5015870_1.fastq > extracted_separated/SRR5015870_with_3_repeats_1.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR5015870_with_3_repeats_2_1_extract.sed raw_data/SRR5015870_2.fastq > extracted_separated/SRR5015870_with_3_repeats_1.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR5015870_without_3_repeats_1_extract.sed raw_data/SRR5015870_1.fastq > extracted_separated/SRR5015870_without_3_repeats.fastq ;
                                            sed -n -f scripts/sed_scripts/SRR5015870_without_3_repeats_2_extract.sed raw_data/SRR5015870_2.fastq >> extracted_separated/SRR5015870_without_3_repeats.fastq" &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR6122935_1_paired_extract.sed raw_data/SRR6122935_1.fastq > extracted/SRR6122935_1_1_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR6122935_1_paired_extract.sed raw_data/SRR6122935_2.fastq > extracted/SRR6122935_2_1_extracted.fastq &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR6122935_with_3_repeats_1_1_extract.sed raw_data/SRR6122935_1.fastq > extracted_separated/SRR6122935_with_3_repeats_1.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR6122935_with_3_repeats_2_1_extract.sed raw_data/SRR6122935_2.fastq > extracted_separated/SRR6122935_with_3_repeats_1.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR6122935_without_3_repeats_1_extract.sed raw_data/SRR6122935_1.fastq > extracted_separated/SRR6122935_without_3_repeats.fastq ;
                                            sed -n -f scripts/sed_scripts/SRR6122935_without_3_repeats_2_extract.sed raw_data/SRR6122935_2.fastq >> extracted_separated/SRR6122935_without_3_repeats.fastq" &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR7880259_1_paired_extract.sed raw_data/SRR7880259_1.fastq > extracted/SRR7880259_1_1_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR7880259_1_paired_extract.sed raw_data/SRR7880259_2.fastq > extracted/SRR7880259_2_1_extracted.fastq &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR7880259_with_3_repeats_2_extract.sed raw_data/SRR7880259_2.fastq > extracted_separated/SRR7880259_with_3_repeats.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR7880259_without_3_repeats_1_extract.sed raw_data/SRR7880259_1.fastq > extracted_separated/SRR7880259_without_3_repeats.fastq" &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR7880260_1_paired_extract.sed raw_data/SRR7880260_1.fastq > extracted/SRR7880260_1_1_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR7880260_1_paired_extract.sed raw_data/SRR7880260_2.fastq > extracted/SRR7880260_2_1_extracted.fastq &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR7880260_with_3_repeats_1_1_extract.sed raw_data/SRR7880260_1.fastq > extracted_separated/SRR7880260_with_3_repeats_1.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR7880260_with_3_repeats_2_1_extract.sed raw_data/SRR7880260_2.fastq > extracted_separated/SRR7880260_with_3_repeats_1.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR7880260_without_3_repeats_1_extract.sed raw_data/SRR7880260_1.fastq > extracted_separated/SRR7880260_without_3_repeats.fastq ;
                                            sed -n -f scripts/sed_scripts/SRR7880260_without_3_repeats_2_extract.sed raw_data/SRR7880260_2.fastq >> extracted_separated/SRR7880260_without_3_repeats.fastq" &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR7880261_1_paired_extract.sed raw_data/SRR7880261_1.fastq > extracted/SRR7880261_1_1_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR7880261_1_paired_extract.sed raw_data/SRR7880261_2.fastq > extracted/SRR7880261_2_1_extracted.fastq &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR7880261_with_3_repeats_1_1_extract.sed raw_data/SRR7880261_1.fastq > extracted_separated/SRR7880261_with_3_repeats_1.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR7880261_with_3_repeats_2_1_extract.sed raw_data/SRR7880261_2.fastq > extracted_separated/SRR7880261_with_3_repeats_1.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR7880261_without_3_repeats_1_extract.sed raw_data/SRR7880261_1.fastq > extracted_separated/SRR7880261_without_3_repeats.fastq ;
                                            sed -n -f scripts/sed_scripts/SRR7880261_without_3_repeats_2_extract.sed raw_data/SRR7880261_2.fastq >> extracted_separated/SRR7880261_without_3_repeats.fastq" &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR7880262_1_paired_extract.sed raw_data/SRR7880262_1.fastq > extracted/SRR7880262_1_1_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR7880262_1_paired_extract.sed raw_data/SRR7880262_2.fastq > extracted/SRR7880262_2_1_extracted.fastq &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR7880262_with_3_repeats_1_1_extract.sed raw_data/SRR7880262_1.fastq > extracted_separated/SRR7880262_with_3_repeats_1.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR7880262_with_3_repeats_2_1_extract.sed raw_data/SRR7880262_2.fastq > extracted_separated/SRR7880262_with_3_repeats_1.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR7880262_without_3_repeats_1_extract.sed raw_data/SRR7880262_1.fastq > extracted_separated/SRR7880262_without_3_repeats.fastq ;
                                            sed -n -f scripts/sed_scripts/SRR7880262_without_3_repeats_2_extract.sed raw_data/SRR7880262_2.fastq >> extracted_separated/SRR7880262_without_3_repeats.fastq" &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR7880263_1_paired_extract.sed raw_data/SRR7880263_1.fastq > extracted/SRR7880263_1_1_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR7880263_1_paired_extract.sed raw_data/SRR7880263_2.fastq > extracted/SRR7880263_2_1_extracted.fastq &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR7880263_with_3_repeats_1_1_extract.sed raw_data/SRR7880263_1.fastq > extracted_separated/SRR7880263_with_3_repeats_1.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR7880263_with_3_repeats_2_1_extract.sed raw_data/SRR7880263_2.fastq > extracted_separated/SRR7880263_with_3_repeats_1.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR7880263_without_3_repeats_1_extract.sed raw_data/SRR7880263_1.fastq > extracted_separated/SRR7880263_without_3_repeats.fastq ;
                                            sed -n -f scripts/sed_scripts/SRR7880263_without_3_repeats_2_extract.sed raw_data/SRR7880263_2.fastq >> extracted_separated/SRR7880263_without_3_repeats.fastq" &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR7880264_1_paired_extract.sed raw_data/SRR7880264_1.fastq > extracted/SRR7880264_1_1_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR7880264_1_paired_extract.sed raw_data/SRR7880264_2.fastq > extracted/SRR7880264_2_1_extracted.fastq &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR7880264_with_3_repeats_1_1_extract.sed raw_data/SRR7880264_1.fastq > extracted_separated/SRR7880264_with_3_repeats_1.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR7880264_with_3_repeats_2_1_extract.sed raw_data/SRR7880264_2.fastq > extracted_separated/SRR7880264_with_3_repeats_1.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR7880264_without_3_repeats_1_extract.sed raw_data/SRR7880264_1.fastq > extracted_separated/SRR7880264_without_3_repeats.fastq ;
                                            sed -n -f scripts/sed_scripts/SRR7880264_without_3_repeats_2_extract.sed raw_data/SRR7880264_2.fastq >> extracted_separated/SRR7880264_without_3_repeats.fastq" &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR7880265_1_paired_extract.sed raw_data/SRR7880265_1.fastq > extracted/SRR7880265_1_1_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR7880265_1_paired_extract.sed raw_data/SRR7880265_2.fastq > extracted/SRR7880265_2_1_extracted.fastq &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR7880265_with_3_repeats_1_1_extract.sed raw_data/SRR7880265_1.fastq > extracted_separated/SRR7880265_with_3_repeats_1.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR7880265_with_3_repeats_2_1_extract.sed raw_data/SRR7880265_2.fastq > extracted_separated/SRR7880265_with_3_repeats_1.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR7880265_without_3_repeats_1_extract.sed raw_data/SRR7880265_1.fastq > extracted_separated/SRR7880265_without_3_repeats.fastq ;
                                            sed -n -f scripts/sed_scripts/SRR7880265_without_3_repeats_2_extract.sed raw_data/SRR7880265_2.fastq >> extracted_separated/SRR7880265_without_3_repeats.fastq" &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR7880266_1_paired_extract.sed raw_data/SRR7880266_1.fastq > extracted/SRR7880266_1_1_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR7880266_1_paired_extract.sed raw_data/SRR7880266_2.fastq > extracted/SRR7880266_2_1_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR7880266_2_paired_extract.sed raw_data/SRR7880266_1.fastq > extracted/SRR7880266_1_2_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR7880266_2_paired_extract.sed raw_data/SRR7880266_2.fastq > extracted/SRR7880266_2_2_extracted.fastq &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR7880266_with_3_repeats_1_1_extract.sed raw_data/SRR7880266_1.fastq > extracted_separated/SRR7880266_with_3_repeats_1.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR7880266_with_3_repeats_2_1_extract.sed raw_data/SRR7880266_2.fastq > extracted_separated/SRR7880266_with_3_repeats_1.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR7880266_with_3_repeats_2_2_extract.sed raw_data/SRR7880266_2.fastq > extracted_separated/SRR7880266_with_3_repeats_2.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR7880266_without_3_repeats_1_extract.sed raw_data/SRR7880266_1.fastq > extracted_separated/SRR7880266_without_3_repeats.fastq ;
                                            sed -n -f scripts/sed_scripts/SRR7880266_without_3_repeats_2_extract.sed raw_data/SRR7880266_2.fastq >> extracted_separated/SRR7880266_without_3_repeats.fastq" &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR7880267_1_paired_extract.sed raw_data/SRR7880267_1.fastq > extracted/SRR7880267_1_1_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR7880267_1_paired_extract.sed raw_data/SRR7880267_2.fastq > extracted/SRR7880267_2_1_extracted.fastq &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR7880267_with_3_repeats_1_1_extract.sed raw_data/SRR7880267_1.fastq > extracted_separated/SRR7880267_with_3_repeats_1.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR7880267_with_3_repeats_2_1_extract.sed raw_data/SRR7880267_2.fastq > extracted_separated/SRR7880267_with_3_repeats_1.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR7880267_without_3_repeats_1_extract.sed raw_data/SRR7880267_1.fastq > extracted_separated/SRR7880267_without_3_repeats.fastq ;
                                            sed -n -f scripts/sed_scripts/SRR7880267_without_3_repeats_2_extract.sed raw_data/SRR7880267_2.fastq >> extracted_separated/SRR7880267_without_3_repeats.fastq" &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR7880268_1_paired_extract.sed raw_data/SRR7880268_1.fastq > extracted/SRR7880268_1_1_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR7880268_1_paired_extract.sed raw_data/SRR7880268_2.fastq > extracted/SRR7880268_2_1_extracted.fastq &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR7880268_with_3_repeats_1_1_extract.sed raw_data/SRR7880268_1.fastq > extracted_separated/SRR7880268_with_3_repeats_1.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR7880268_with_3_repeats_2_1_extract.sed raw_data/SRR7880268_2.fastq > extracted_separated/SRR7880268_with_3_repeats_1.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR7880268_without_3_repeats_1_extract.sed raw_data/SRR7880268_1.fastq > extracted_separated/SRR7880268_without_3_repeats.fastq ;
                                            sed -n -f scripts/sed_scripts/SRR7880268_without_3_repeats_2_extract.sed raw_data/SRR7880268_2.fastq >> extracted_separated/SRR7880268_without_3_repeats.fastq" &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR7880269_1_paired_extract.sed raw_data/SRR7880269_1.fastq > extracted/SRR7880269_1_1_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR7880269_1_paired_extract.sed raw_data/SRR7880269_2.fastq > extracted/SRR7880269_2_1_extracted.fastq &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR7880269_with_3_repeats_1_1_extract.sed raw_data/SRR7880269_1.fastq > extracted_separated/SRR7880269_with_3_repeats_1.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR7880269_with_3_repeats_2_1_extract.sed raw_data/SRR7880269_2.fastq > extracted_separated/SRR7880269_with_3_repeats_1.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR7880269_without_3_repeats_1_extract.sed raw_data/SRR7880269_1.fastq > extracted_separated/SRR7880269_without_3_repeats.fastq ;
                                            sed -n -f scripts/sed_scripts/SRR7880269_without_3_repeats_2_extract.sed raw_data/SRR7880269_2.fastq >> extracted_separated/SRR7880269_without_3_repeats.fastq" &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR7880270_1_paired_extract.sed raw_data/SRR7880270_1.fastq > extracted/SRR7880270_1_1_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR7880270_1_paired_extract.sed raw_data/SRR7880270_2.fastq > extracted/SRR7880270_2_1_extracted.fastq &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR7880270_with_3_repeats_1_1_extract.sed raw_data/SRR7880270_1.fastq > extracted_separated/SRR7880270_with_3_repeats_1.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR7880270_with_3_repeats_2_1_extract.sed raw_data/SRR7880270_2.fastq > extracted_separated/SRR7880270_with_3_repeats_1.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR7880270_without_3_repeats_1_extract.sed raw_data/SRR7880270_1.fastq > extracted_separated/SRR7880270_without_3_repeats.fastq ;
                                            sed -n -f scripts/sed_scripts/SRR7880270_without_3_repeats_2_extract.sed raw_data/SRR7880270_2.fastq >> extracted_separated/SRR7880270_without_3_repeats.fastq" &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR7880271_1_paired_extract.sed raw_data/SRR7880271_1.fastq > extracted/SRR7880271_1_1_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR7880271_1_paired_extract.sed raw_data/SRR7880271_2.fastq > extracted/SRR7880271_2_1_extracted.fastq &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR7880271_with_3_repeats_1_1_extract.sed raw_data/SRR7880271_1.fastq > extracted_separated/SRR7880271_with_3_repeats_1.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR7880271_with_3_repeats_2_1_extract.sed raw_data/SRR7880271_2.fastq > extracted_separated/SRR7880271_with_3_repeats_1.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR7880271_without_3_repeats_1_extract.sed raw_data/SRR7880271_1.fastq > extracted_separated/SRR7880271_without_3_repeats.fastq ;
                                            sed -n -f scripts/sed_scripts/SRR7880271_without_3_repeats_2_extract.sed raw_data/SRR7880271_2.fastq >> extracted_separated/SRR7880271_without_3_repeats.fastq" &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR7880272_1_paired_extract.sed raw_data/SRR7880272_1.fastq > extracted/SRR7880272_1_1_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR7880272_1_paired_extract.sed raw_data/SRR7880272_2.fastq > extracted/SRR7880272_2_1_extracted.fastq &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR7880272_with_3_repeats_1_1_extract.sed raw_data/SRR7880272_1.fastq > extracted_separated/SRR7880272_with_3_repeats_1.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR7880272_with_3_repeats_2_1_extract.sed raw_data/SRR7880272_2.fastq > extracted_separated/SRR7880272_with_3_repeats_1.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR7880272_without_3_repeats_1_extract.sed raw_data/SRR7880272_1.fastq > extracted_separated/SRR7880272_without_3_repeats.fastq ;
                                            sed -n -f scripts/sed_scripts/SRR7880272_without_3_repeats_2_extract.sed raw_data/SRR7880272_2.fastq >> extracted_separated/SRR7880272_without_3_repeats.fastq" &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR7880273_1_paired_extract.sed raw_data/SRR7880273_1.fastq > extracted/SRR7880273_1_1_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR7880273_1_paired_extract.sed raw_data/SRR7880273_2.fastq > extracted/SRR7880273_2_1_extracted.fastq &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR7880273_with_3_repeats_1_1_extract.sed raw_data/SRR7880273_1.fastq > extracted_separated/SRR7880273_with_3_repeats_1.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR7880273_with_3_repeats_2_1_extract.sed raw_data/SRR7880273_2.fastq > extracted_separated/SRR7880273_with_3_repeats_1.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR7880273_without_3_repeats_1_extract.sed raw_data/SRR7880273_1.fastq > extracted_separated/SRR7880273_without_3_repeats.fastq ;
                                            sed -n -f scripts/sed_scripts/SRR7880273_without_3_repeats_2_extract.sed raw_data/SRR7880273_2.fastq >> extracted_separated/SRR7880273_without_3_repeats.fastq" &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR7880274_1_paired_extract.sed raw_data/SRR7880274_1.fastq > extracted/SRR7880274_1_1_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR7880274_1_paired_extract.sed raw_data/SRR7880274_2.fastq > extracted/SRR7880274_2_1_extracted.fastq &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR7880274_with_3_repeats_1_1_extract.sed raw_data/SRR7880274_1.fastq > extracted_separated/SRR7880274_with_3_repeats_1.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR7880274_with_3_repeats_2_1_extract.sed raw_data/SRR7880274_2.fastq > extracted_separated/SRR7880274_with_3_repeats_1.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR7880274_without_3_repeats_1_extract.sed raw_data/SRR7880274_1.fastq > extracted_separated/SRR7880274_without_3_repeats.fastq ;
                                            sed -n -f scripts/sed_scripts/SRR7880274_without_3_repeats_2_extract.sed raw_data/SRR7880274_2.fastq >> extracted_separated/SRR7880274_without_3_repeats.fastq" &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR7880275_1_paired_extract.sed raw_data/SRR7880275_1.fastq > extracted/SRR7880275_1_1_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR7880275_1_paired_extract.sed raw_data/SRR7880275_2.fastq > extracted/SRR7880275_2_1_extracted.fastq &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR7880275_with_3_repeats_1_1_extract.sed raw_data/SRR7880275_1.fastq > extracted_separated/SRR7880275_with_3_repeats_1.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR7880275_with_3_repeats_2_1_extract.sed raw_data/SRR7880275_2.fastq > extracted_separated/SRR7880275_with_3_repeats_1.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR7880275_without_3_repeats_1_extract.sed raw_data/SRR7880275_1.fastq > extracted_separated/SRR7880275_without_3_repeats.fastq ;
                                            sed -n -f scripts/sed_scripts/SRR7880275_without_3_repeats_2_extract.sed raw_data/SRR7880275_2.fastq >> extracted_separated/SRR7880275_without_3_repeats.fastq" &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR7880276_1_paired_extract.sed raw_data/SRR7880276_1.fastq > extracted/SRR7880276_1_1_extracted.fastq &
srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/SRR7880276_1_paired_extract.sed raw_data/SRR7880276_2.fastq > extracted/SRR7880276_2_1_extracted.fastq &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR7880276_with_3_repeats_1_1_extract.sed raw_data/SRR7880276_1.fastq > extracted_separated/SRR7880276_with_3_repeats_1.fastq" &
srun -n1 -N1 --exclusive bash -c "
                                            sed -n -f scripts/sed_scripts/SRR7880276_with_3_repeats_2_1_extract.sed raw_data/SRR7880276_2.fastq > extracted_separated/SRR7880276_with_3_repeats_1.fastq" &
srun -n1 -N1 --exclusive bash -c "sed -n -f scripts/sed_scripts/SRR7880276_without_3_repeats_1_extract.sed raw_data/SRR7880276_1.fastq > extracted_separated/SRR7880276_without_3_repeats.fastq ;
                                            sed -n -f scripts/sed_scripts/SRR7880276_without_3_repeats_2_extract.sed raw_data/SRR7880276_2.fastq >> extracted_separated/SRR7880276_without_3_repeats.fastq" &
wait
