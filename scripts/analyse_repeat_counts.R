
library(tidyverse)

write("", file="repeat_count_summary.txt")

for(sample_ID in gsub("_1.repeat_count", "", list.files("repeat_counts", pattern = "_1.repeat_count"))) {

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

write(sample_ID, file = "repeat_count_summary.txt", append = TRUE)
write.table(with(sample_counts, table(sample_1_counts, sample_2_counts)), file = "repeat_count_summary.txt", append = TRUE)

}


sed_commands <- c("#!/bin/bash", "#SBATCH --job-name=Extract_repeats",
                  "#SBATCH --ntasks=48", "#SBATCH --mem-per-cpu=1G",
                  "#SBATCH --cpus-per-task=1",
                  "mkdir extracted", "mkdir extracted_separated")

for(sample_ID in gsub("_1.repeat_count", "", list.files("repeat_counts", pattern = "_1.repeat_count"))) {

print(sample_ID)
print(length(sed_commands))

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


more_than_3_repeats_in_either <- sample_counts %>% filter(sample_1_counts >= 3 | sample_2_counts >= 3)

if (is.infinite(max(more_than_3_repeats_in_either$sample_ID)) == FALSE) {

split_sample_IDs <- split(more_than_3_repeats_in_either$sample_ID, ceiling(seq_along(more_than_3_repeats_in_either$sample_ID)/1000))

for (i in 1:length(split_sample_IDs)) {

sed_commands <- append(sed_commands, paste0("srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/",
                                            sample_ID, "_", i, "_paired_extract.sed raw_data/",
                                            sample_ID, "_1.fastq > extracted/",
                                            sample_ID, "_1_", i, "_extracted.fastq &"))

sed_commands <- append(sed_commands, paste0("srun -n1 -N1 --exclusive sed -n -f scripts/sed_scripts/",
                                            sample_ID, "_", i, "_paired_extract.sed raw_data/",
                                            sample_ID, "_2.fastq > extracted/",
                                            sample_ID, "_2_", i, "_extracted.fastq &"))

this_sample_IDs <- split_sample_IDs[[i]]

write(paste0(paste0(paste(paste0(format(this_sample_IDs - 1, scientific = FALSE, trim = TRUE),
                                      ",",
                                      format(this_sample_IDs + 2, scientific = FALSE, trim = TRUE),
                                      "p"),
                         collapse = "\n"), "\n", format(max(this_sample_IDs) + 3, scientific = FALSE, trim = TRUE), "q")),
      file = paste0("scripts/sed_scripts/", sample_ID, "_", i, "_paired_extract.sed"))

}

}

more_than_3_repeats_in_first_member <- sample_counts %>% filter(sample_1_counts >= 3 | sample_2_counts >= 3) %>%
                                                   filter(sample_1_counts >= 3)
more_than_3_repeats_in_second_member <- sample_counts %>% filter(sample_1_counts >= 3 | sample_2_counts >= 3) %>%
                                                   filter(sample_2_counts >= 3)

if ((is.infinite(max(more_than_3_repeats_in_first_member$sample_ID)) == FALSE) & (is.infinite(max(more_than_3_repeats_in_second_member$sample_ID)) == FALSE)) {

split_sample_IDs <- split(more_than_3_repeats_in_first_member$sample_ID, ceiling(seq_along(more_than_3_repeats_in_first_member$sample_ID)/1000))

for (i in 1:length(split_sample_IDs)) {

sed_commands <- append(sed_commands, paste0("srun -n1 -N1 --exclusive bash -c \"sed -n -f scripts/sed_scripts/",
                                            sample_ID, "_with_3_repeats_1_", i, "_extract.sed raw_data/",
                                            sample_ID, "_1.fastq > extracted_separated/",
                                            sample_ID, "_with_3_repeats_", i, ".fastq\" &"))

this_sample_IDs <- split_sample_IDs[[i]]

write(paste0(paste0(paste(paste0(format(this_sample_IDs - 1, scientific = FALSE, trim = TRUE),
                                ",",
                                format(this_sample_IDs + 2, scientific = FALSE, trim = TRUE),
                                "p"),
                         collapse = "\n"), "\n", format(max(this_sample_IDs) + 3, scientific = FALSE, trim = TRUE), "q")),
      file = paste0("scripts/sed_scripts/", sample_ID, "_with_3_repeats_1_", i, "_extract.sed"))

}

split_sample_IDs <- split(more_than_3_repeats_in_second_member$sample_ID, ceiling(seq_along(more_than_3_repeats_in_second_member$sample_ID)/1000))

for (i in 1:length(split_sample_IDs)) {


sed_commands <- append(sed_commands, paste0("srun -n1 -N1 --exclusive bash -c \"
                                            sed -n -f scripts/sed_scripts/",
                                            sample_ID, "_with_3_repeats_2_", i, "_extract.sed raw_data/",
                                            sample_ID, "_2.fastq > extracted_separated/",
                                            sample_ID, "_with_3_repeats_", i, ".fastq\" &"))

this_sample_IDs <- split_sample_IDs[[i]]

write(paste0(paste0(paste(paste0(format(this_sample_IDs - 1, scientific = FALSE, trim = TRUE),
                                ",",
                                format(this_sample_IDs + 2, scientific = FALSE, trim = TRUE),
                                "p"),
                         collapse = "\n"), "\n", format(max(this_sample_IDs) + 3, scientific = FALSE, trim = TRUE), "q")),
      file = paste0("scripts/sed_scripts/", sample_ID, "_with_3_repeats_2_", i, "_extract.sed"))
}

} else if (is.infinite(max(more_than_3_repeats_in_first_member$sample_ID))) {

sed_commands <- append(sed_commands, paste0("srun -n1 -N1 --exclusive bash -c \"sed -n -f scripts/sed_scripts/",
                                            sample_ID, "_with_3_repeats_2_extract.sed raw_data/",
                                            sample_ID, "_2.fastq > extracted_separated/",
                                            sample_ID, "_with_3_repeats.fastq\" &"))

write(paste0(with(sample_counts %>% filter(sample_1_counts >= 3 | sample_2_counts >= 3) %>%
                                                   filter(sample_2_counts >= 3),
            paste0(paste(paste0(format(sample_ID - 1, scientific = FALSE, trim = TRUE),
                                ",",
                                format(sample_ID + 2, scientific = FALSE, trim = TRUE),
                                "p"),
                         collapse = "\n"), "\n", format(max(sample_ID) + 3, scientific = FALSE, trim = TRUE), "q"))),
      file = paste0("scripts/sed_scripts/", sample_ID, "_with_3_repeats_2_extract.sed"))


} else if (is.infinite(max(more_than_3_repeats_in_second_member$sample_ID))) {

sed_commands <- append(sed_commands, paste0("srun -n1 -N1 --exclusive bash -c \"sed -n -f scripts/sed_scripts/",
                                            sample_ID, "_with_3_repeats_1_extract.sed raw_data/",
                                            sample_ID, "_1.fastq > extracted_separated/",
                                            sample_ID, "_with_3_repeats.fastq\" &"))

write(paste0(with(sample_counts %>% filter(sample_1_counts >= 3 | sample_2_counts >= 3) %>%
                                                   filter(sample_1_counts >= 3),
            paste0(paste(paste0(format(sample_ID - 1, scientific = FALSE, trim = TRUE),
                                ",",
                                format(sample_ID + 2, scientific = FALSE, trim = TRUE),
                                "p"),
                         collapse = "\n"), "\n", format(max(sample_ID) + 3, scientific = FALSE, trim = TRUE), "q"))),
      file = paste0("scripts/sed_scripts/", sample_ID, "_with_3_repeats_1_extract.sed"))


}

less_than_3_repeats_in_first_member <- sample_counts %>% filter(sample_1_counts >= 3 | sample_2_counts >= 3) %>%
                                                   filter(sample_1_counts < 3)
less_than_3_repeats_in_second_member <- sample_counts %>% filter(sample_1_counts >= 3 | sample_2_counts >= 3) %>%
                                                   filter(sample_2_counts < 3)

if ((is.infinite(max(less_than_3_repeats_in_first_member$sample_ID)) == FALSE) & (is.infinite(max(less_than_3_repeats_in_second_member$sample_ID)) == FALSE)) {

sed_commands <- append(sed_commands, paste0("srun -n1 -N1 --exclusive bash -c \"sed -n -f scripts/sed_scripts/",
                                            sample_ID, "_without_3_repeats_1_extract.sed raw_data/",
                                            sample_ID, "_1.fastq > extracted_separated/",
                                            sample_ID, "_without_3_repeats.fastq ;
                                            sed -n -f scripts/sed_scripts/",
                                            sample_ID, "_without_3_repeats_2_extract.sed raw_data/",
                                            sample_ID, "_2.fastq >> extracted_separated/",
                                            sample_ID, "_without_3_repeats.fastq\" &"))


write(paste0(with(sample_counts %>% filter(sample_1_counts >= 3 | sample_2_counts >= 3) %>%
                                                   filter(sample_1_counts < 3),
            paste0(paste(paste0(format(sample_ID - 1, scientific = FALSE, trim = TRUE),
                                ",",
                                format(sample_ID + 2, scientific = FALSE, trim = TRUE),
                                "p"),
                         collapse = "\n"), "\n", format(max(sample_ID) + 3, scientific = FALSE, trim = TRUE), "q"))),
       file = paste0("scripts/sed_scripts/", sample_ID, "_without_3_repeats_1_extract.sed"))
                                            
write(paste0(with(sample_counts %>% filter(sample_1_counts >= 3 | sample_2_counts >= 3) %>%
                                                   filter(sample_2_counts < 3),
            paste0(paste(paste0(format(sample_ID - 1, scientific =  FALSE, trim = TRUE),
                                ",",
                                format(sample_ID + 2, scientific = FALSE, trim = TRUE),
                                "p"),
                         collapse = "\n"), "\n", format(max(sample_ID) + 3, scientific = FALSE, trim = TRUE), "q"))),
       file = paste0("scripts/sed_scripts/", sample_ID, "_without_3_repeats_2_extract.sed"))

} else if (is.infinite(max(less_than_3_repeats_in_first_member$sample_ID))) {

sed_commands <- append(sed_commands, paste0("srun -n1 -N1 --exclusive bash -c \"sed -n -f scripts/sed_scripts/",
                                            sample_ID, "_without_3_repeats_2_extract.sed raw_data/",
                                            sample_ID, "_2.fastq > extracted_separated/",
                                            sample_ID, "_without_3_repeats.fastq\" &"))

write(paste0(with(sample_counts %>% filter(sample_1_counts >= 3 | sample_2_counts >= 3) %>%
                                                   filter(sample_2_counts < 3),
            paste0(paste(paste0(format(sample_ID - 1, scientific = FALSE, trim = TRUE),
                                ",",
                                format(sample_ID + 2, scientific = FALSE, trim = TRUE),
                                "p"),
                         collapse = "\n"), "\n", format(max(sample_ID) + 3, scientific = FALSE, trim = TRUE), "q"))),
       file = paste0("scripts/sed_scripts/", sample_ID, "_without_3_repeats_2_extract.sed"))

} else if (is.infinite(max(less_than_3_repeats_in_second_member$sample_ID))) {

sed_commands <- append(sed_commands, paste0("srun -n1 -N1 --exclusive bash -c \"sed -n -f scripts/sed_scripts/",
                                            sample_ID, "_without_3_repeats_1_extract.sed raw_data/",
                                            sample_ID, "_1.fastq > extracted_separated/",
                                            sample_ID, "_without_3_repeats.fastq\" &"))

write(paste0(with(sample_counts %>% filter(sample_1_counts >= 3 | sample_2_counts >= 3) %>%
                                                   filter(sample_1_counts < 3),
            paste0(paste(paste0(format(sample_ID - 1, scientific = FALSE, trim = TRUE),
                                ",",
                                format(sample_ID + 2, scientific = FALSE, trim = TRUE),
                                "p"),
                         collapse = "\n"), "\n", format(max(sample_ID) + 3, scientific = FALSE, trim = TRUE), "q"))),
       file = paste0("scripts/sed_scripts/", sample_ID, "_without_3_repeats_1_extract.sed"))
            
}


}

sed_commands <- append(sed_commands, c("wait"))

write(sed_commands, file = "scripts/extract_sequences_to_map_subtelomeres.sbatch")

