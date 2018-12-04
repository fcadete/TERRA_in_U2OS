
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)

sample_ID <- args[1]

write("", file=paste0("repeat_count_summaries/", sample_ID, "_repeat_count_summary.txt"))


sample_1 <- read_delim(paste0("repeat_counts/", sample_ID, "_1.repeat_count"),
                           delim = " ",
                           col_names = c("sample_1_counts", "sample_ID"))

sample_1$sample <- sample_ID

sample_2 <- read_delim(paste0("repeat_counts/", sample_ID, "_2.repeat_count"),
                           delim = " ",
                           col_names = c("sample_2_counts", "sample_ID"))

sample_2$sample <- sample_ID

sample_counts <- full_join(sample_1, sample_2)

sample_counts[is.na(sample_counts$sample_1_counts), "sample_1_counts"] <- 0
sample_counts[is.na(sample_counts$sample_2_counts), "sample_2_counts"] <- 0

write(sample_ID, file=paste0("repeat_count_summaries/", sample_ID, "_repeat_count_summary.txt"))
write(with(sample_counts, table(sample_1_counts, sample_2_counts)), file =paste0("repeat_count_summaries/", sample_ID, "_repeat_count_summary.txt"), append = TRUE)




sample_1 <- read_delim(paste0("repeat_counts/", sample_ID, "_1.repeat_count"),
                           delim = " ",
                           col_names = c("sample_1_counts", "sample_ID"))

sample_1$sample <- sample_ID

sample_2 <- read_delim(paste0("repeat_counts/", sample_ID, "_2.repeat_count"),
                           delim = " ",
                           col_names = c("sample_2_counts", "sample_ID"))

sample_2$sample <- sample_ID

sample_counts <- full_join(sample_1, sample_2)

sample_counts[is.na(sample_counts$sample_1_counts), "sample_1_counts"] <- 0
sample_counts[is.na(sample_counts$sample_2_counts), "sample_2_counts"] <- 0

sed_commands <- paste0("srun -n1 -N1 --exclusive sed -n '",
       with(sample_counts %>% filter(sample_1_counts >= 3 | sample_2_counts >= 3),
            paste0(paste(paste0(sample_ID - 1, ",", sample_ID + 2, "p"),
                         collapse = ";"), ";", max(sample_ID) + 3, "q")),
       "' raw_data/", sample_ID, "_1.fastq > extracted/", sample_ID, "_1_extracted.fastq &")

sed_commands <- append(sed_commands, paste0("srun -n1 -N1 --exclusive sed -n '",
       with(sample_counts %>% filter(sample_1_counts >= 3 | sample_2_counts >= 3),
            paste0(paste(paste0(sample_ID - 1, ",", sample_ID + 2, "p"),
                         collapse = ";"), ";", max(sample_ID) + 3, "q")),
       "' raw_data/", sample_ID, "_2.fastq > extracted/", sample_ID, "_2_extracted.fastq &"))

sed_commands <- append(sed_commands, paste0("srun -n1 -N1 --exclusive bash -c \"sed -n '",
                                            with(sample_counts %>% filter(sample_1_counts >= 3 | sample_2_counts >= 3) %>%
                                                   filter(sample_1_counts >= 3),
            paste0(paste(paste0(sample_ID - 1, ",", sample_ID + 2, "p"),
                         collapse = ";"), ";", max(sample_ID) + 3, "q")),
       "' raw_data/", sample_ID, "_1.fastq > extracted_separated/", sample_ID, "_with_3_repeats.fastq | sed -n '",
                                            with(sample_counts %>% filter(sample_1_counts >= 3 | sample_2_counts >= 3) %>%
                                                   filter(sample_2_counts >= 3),
            paste0(paste(paste0(sample_ID - 1, ",", sample_ID + 2, "p"),
                         collapse = ";"), ";", max(sample_ID) + 3, "q")),
       "' raw_data/", sample_ID, "_2.fastq >> extracted_separated/", sample_ID, "_with_3_repeats.fastq\" &"))


sed_commands <- append(sed_commands, paste0("srun -n1 -N1 --exclusive bash -c \"sed -n '",
                                            with(sample_counts %>% filter(sample_1_counts >= 3 | sample_2_counts >= 3) %>%
                                                   filter(sample_1_counts < 3),
            paste0(paste(paste0(sample_ID - 1, ",", sample_ID + 2, "p"),
                         collapse = ";"), ";", max(sample_ID) + 3, "q")),
       "' raw_data/", sample_ID, "_1.fastq > extracted_separated/", sample_ID, "_without_3_repeats.fastq | sed -n '",
                                            with(sample_counts %>% filter(sample_1_counts >= 3 | sample_2_counts >= 3) %>%
                                                   filter(sample_2_counts < 3),
            paste0(paste(paste0(sample_ID - 1, ",", sample_ID + 2, "p"),
                         collapse = ";"), ";", max(sample_ID) + 3, "q")),
       "' raw_data/", sample_ID, "_2.fastq >> extracted_separated/", sample_ID, "_without_3_repeats.fastq\" &"))



write(sed_commands, file = "scripts/extract_sequences_to_map_subtelomeres.sbatch", append = TRUE)


