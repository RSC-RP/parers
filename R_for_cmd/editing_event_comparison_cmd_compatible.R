### Script for making bar graphs with contextualized types of editing occurring at individual editing sites




# command line arguments
args <- commandArgs(trailingOnly = TRUE)
primer_f <- args[1]
primer_r <- args[2]
data_list <- strsplit(args[3], ",")[[1]]
sample_names <- strsplit(args[4], ",")[[1]]
wrdir <- args[5]
today <- args[6]
gene <- args[7]





# What size would you like your output images to be? h = height, w = width
hi <- 10
wi <- 5 * length(sample_names)



#### I advise against changing any of the code below ####

#install.packages("tidyverse")
library("tidyverse")

file_type <- ".pdf"

if (wi < 15 ){
  wi == 25
}

remove_cell_name <- function(name) {
  new_name <- gsub(gene, "", name) %>%
    trimws()
  return(new_name)
}
sample2 <- lapply(sample_names, remove_cell_name)
sample <- data.frame(sample = unlist(sample2))

#remove Ts from primer sequences and extract primer lengths
pre_primer_length_f <- gsub("T", "", primer_f) 
primer_length_f <- nchar(pre_primer_length_f)

pre_primer_length_r <- gsub("T", "", primer_r) 
primer_length_r <- nchar(pre_primer_length_r)

#create a list of the paths to call the ordered counts csv data frames
data_paths2 <- as.list(data_list)

#create a list of the data frames
genes_bar <- lapply(data_paths2, read.csv)

# determine max value in NE column
max_ne <- max(genes_bar[[1]]$NE)

# take out any forward primer bases that experienced editing
process_dataframe <- function(df) {
  pre_for <- df %>% 
    slice((1:primer_length_f))
  for_slice2 <- filter(pre_for, NE != max(NE))
  return(for_slice2)
}
for_slice = lapply(genes_bar, process_dataframe)

# take out any reverse primer bases that experienced editing
nrows <- nrow(genes_bar[[1]]) # determine length of transcript
rev_primer_index <- nrows - primer_length_r + 1 # find where reverse primer starts

process_rev <- function(df, rev_primer_index) {
  nrows <- nrow(df)
  pre_rev <- df %>% 
    slice(rev_primer_index:nrows)
  rev_slice2 <- filter(pre_rev, NE != max(NE))
  return(rev_slice2)
}
rev_slice <- lapply(genes_bar, process_rev, rev_primer_index)

max_rows_f <- max(sapply(for_slice, nrow))
max_rows_r <- max(sapply(rev_slice, nrow))

forward_trim <- primer_length_f - max_rows_f
reverse_trim <- primer_length_r - max_rows_r

#remove forward primers from sequences
slice_dataframe_f <- function(df, primer_length_f) {
  df %>% slice(-(1:primer_length_f))
}
genes_subset_pre <- lapply(genes_bar, slice_dataframe_f, forward_trim)

#remove reverse primers from sequences
slice_dataframe_r <- function(df, primer_length_r) {
  n <- nrow(df)
  df %>% slice(1:(n - primer_length_r))
}
genes_subset <- lapply(genes_subset_pre, slice_dataframe_r, reverse_trim)

#determine the length of each transcript, needed to maintain sequence order later
lengt <- map(genes_subset, ~ .x %>%
               select(Nucleotide) %>%
               nrow())
lengt <- lengt[[1]]


# pull out the nucleotide and order columns
consensus_columns <- list()
for (df in genes_subset) {
  selected_columns_df <- df[, c("Nucleotide", "Order")]
  consensus_columns <- c(consensus_columns, list(selected_columns_df))
}

# create new data frame for the fifth level of categories
level1 <- list()
for (df in genes_subset) {
  level1_df <- df[, c("NCI", "NCD", "NE", "CI", "CII", "CEI", "DCI", "CINE", "CD", "CID", "CED", "ICD", "CDNE")]
  level1 <- c(level1, list(level1_df))
}
level1 <- map2(consensus_columns, level1, ~ bind_cols(.x, .y))

# create new data frame for the fifth level of categories
level2 <- list()
for (df in genes_subset) {
  level2_df <- df[, c("N.inc.ins", "N.inc.del", "No.edit1", "Ins.ed", "I.inc", "Del.ed", "D.inc")]
  level2 <- c(level2, list(level2_df))
}
level2 <- map2(consensus_columns, level2, ~ bind_cols(.x, .y))

# create new data frame for the fifth level of categories
level3 <- list()
for (df in genes_subset) {
  level3_df <- df[, c("Insertion", "Deletion", "No.edit2")]
  level3 <- c(level3, list(level3_df))
}
level3 <- map2(consensus_columns, level3, ~ bind_cols(.x, .y))

# create new data frame for the fifth level of categories
level4 <- list()
for (df in genes_subset) {
  level4_df <- df[, c("NC.edit", "PE", "FE")]
  level4 <- c(level4, list(level4_df))
}
level4 <- map2(consensus_columns, level4, ~ bind_cols(.x, .y))

# create new data frame for the fifth level of categories
level5 <- list()
for (df in genes_subset) {
  level5_df <- df[, c("No.edit3", "Insertion1", "Deletion1", "NC.Insertion", "NC.Deletion")]
  level5 <- c(level5, list(level5_df))
}
level5 <- map2(consensus_columns, level5, ~ bind_cols(.x, .y))

# extract the editing events of interest and sum up their total counts
sum_col <- function(df) {
  df_trunc <- df[, -c(1, 2)]
  df_summed <- colSums(df_trunc)
  return(df_summed)
}
summed_cols1 <- lapply(level1, sum_col)
summed_cols2 <- lapply(level2, sum_col)
summed_cols3 <- lapply(level3, sum_col)
summed_cols4 <- lapply(level4, sum_col)
summed_cols5 <- lapply(level5, sum_col)


# find the percent occurrence of each editing event in each cell line
quick_maths <- function(dat_list) {
  perc <- lapply(dat_list, function(x) x / (1000000 * lengt))
  perc_sum <- sum(unlist(perc))
  perc_norm <- lapply(perc, function(x) x / perc_sum)
  perc_100 <- lapply(perc_norm, function(x) x * 100)
  return(perc_100)
}
perc_count1 <- lapply(summed_cols1, quick_maths)
perc_count2 <- lapply(summed_cols2, quick_maths)
perc_count3 <- lapply(summed_cols3, quick_maths)
perc_count4 <- lapply(summed_cols4, quick_maths)
perc_count5 <- lapply(summed_cols5, quick_maths)

# convert list of lists to a data frame
perc_df1 <- do.call(rbind, lapply(perc_count1, as.data.frame))
perc_df2 <- do.call(rbind, lapply(perc_count2, as.data.frame))
perc_df3 <- do.call(rbind, lapply(perc_count3, as.data.frame))
perc_df4 <- do.call(rbind, lapply(perc_count4, as.data.frame))
perc_df5 <- do.call(rbind, lapply(perc_count5, as.data.frame))

# add cell line names to the dataframe
perc_df1 <- cbind(sample, perc_df1)
perc_df2 <- cbind(sample, perc_df2)
perc_df3 <- cbind(sample, perc_df3)
perc_df4 <- cbind(sample, perc_df4)
perc_df5 <- cbind(sample, perc_df5)

# Convert to long format
perc_df_long1 <- perc_df1 %>%
  pivot_longer(
    cols = -sample,          
    names_to = "Group",     
    values_to = "Value"     
  )
perc_df_long2 <- perc_df2 %>%
  pivot_longer(
    cols = -sample,          
    names_to = "Group",     
    values_to = "Value"     
  )

perc_df_long3 <- perc_df3 %>%
  pivot_longer(
    cols = -sample,          
    names_to = "Group",     
    values_to = "Value"     
  )

perc_df_long4 <- perc_df4 %>%
  pivot_longer(
    cols = -sample,          
    names_to = "Group",     
    values_to = "Value"     
  )

perc_df_long5 <- perc_df5 %>%
  pivot_longer(
    cols = -sample,          
    names_to = "Group",     
    values_to = "Value"     
  )


### start of level 1 plot ###
ggplot(data = perc_df_long1, aes(x = sample, y = Value, fill = Group)) +
  geom_bar(stat = "identity", position = "stack", width = .5) +
  labs(title = paste0("Comparison of Editing Events\nBetween ", gene, " Cell Lines"), x = "", y = "% of Editing Events") +
  scale_fill_manual(values = c(NE = "#e6e7e8", NCI = "#000000", NCD = "#888787", CII = "#FFFF00", CID = "#33FFFF",
                               CI = "#ff0000", CD = "#0066ff", CEI = "#ffc000", DCI = "#0fb2f6",
                               CINE = "#654321", CED = "#9932cc", ICD = "#6ff59b", CDNE = "#ffe4b5"),
                    labels = c(NCI = "NCI", NCD = "NCD", NE = "NE", CI = "CI", CII = "CII", CEI = "CEI", DCI = "DCI", CINE = "CINE", 
                               CD = "CD", CID = "CID", CED = "CED", ICD = "ICD", CDNE = "CDNE")) +
  scale_y_continuous(limits = c(0, 101),
                     breaks = seq(0, 100, by = 10), expand = c(0,0)) +
  theme_minimal() + 
  theme(
    plot.title.position = "plot",
    plot.title = element_text(hjust = .5, size = 45, face = "bold"),
    legend.position = "bottom",
    legend.title = element_blank(),
    axis.title.x = element_text(size = 50, face = "bold"),
    axis.title.y = element_text(size = 50, face = "bold"),
    axis.text.x = element_text(size = 35, face = "bold", color = "#000000"),
    axis.text.y = element_text(size = 35, face = "bold", color = "#000000"),
    legend.text = element_text(size = 25, face = "bold")
  )


# save the plot
ggsave(paste0(wrdir, "/", gene, "_editing_event_comparison_bar_lvl1_", today, file_type), dpi = 600, width = wi, height = hi, limitsize = FALSE)
### end of level 1 plot ###

### start of level 2 plot ###
ggplot(data = perc_df_long2, aes(x = sample, y = Value, fill = Group)) +
  geom_bar(stat = "identity", position = "stack", width = .5) +
  labs(title = paste0("Comparison of Editing Events\nBetween ", gene, " Cell Lines"), x = "", y = "% of Editing Events") +
  scale_fill_manual(values = c(No.edit1 = "#e6e7e8", N.inc.ins = "#000000", N.inc.del = "#888787", I.inc = "#ffc000", D.inc = "#6ff59b",
                               Ins.ed = "#ff0000", Del.ed = "#0066ff"),
                    labels = c(No.edit1 = "No Edit", Del.ed = "Deletion", Ins.ed = "Insertion", D.inc = "Incomplete Deletion", 
                               I.inc = "Incomplete Insertion", N.inc.ins = "NC Incomplete Insertion", N.inc.del = "NC Incomplete Deletion")) +
  scale_y_continuous(limits = c(0, 101),
                     breaks = seq(0, 100, by = 10), expand = c(0,0)) +
  theme_minimal() + 
  theme(
    plot.title.position = "plot",
    plot.title = element_text(hjust = .5, size = 45, face = "bold"),
    legend.position = "bottom",
    legend.title = element_blank(),
    axis.title.x = element_text(size = 50, face = "bold"),
    axis.title.y = element_text(size = 50, face = "bold"),
    axis.text.x = element_text(size = 35, face = "bold", color = "#000000"),
    axis.text.y = element_text(size = 35, face = "bold", color = "#000000"),
    legend.text = element_text(size = 25, face = "bold")
  )


# save the plot
ggsave(paste0(wrdir, "/", gene, "_editing_event_comparison_bar_lvl2_", today, file_type), dpi = 600, width = wi, height = hi, limitsize = FALSE)
### end of level 2 plot ###

### start of level 3 plot ###
ggplot(data = perc_df_long3, aes(x = sample, y = Value, fill = Group)) +
  geom_bar(stat = "identity", position = "stack", width = .5) +
  labs(title = paste0("Comparison of Editing Events\nBetween ", gene, " Cell Lines"), x = "", y = "% of Editing Events") +
  scale_fill_manual(values = c(No.edit2 = "#e6e7e8", "Insertion" = "#ff0000", "Deletion" = "#0066ff"),
                    labels = c(No.edit2 = "No Edit", "Insertion", "Deletion")) +
  scale_y_continuous(limits = c(0, 101),
                     breaks = seq(0, 100, by = 10), expand = c(0,0)) +
  theme_minimal() + 
  theme(
    plot.title.position = "plot",
    plot.title = element_text(hjust = .5, size = 45, face = "bold"),
    legend.position = "bottom",
    legend.title = element_blank(),
    axis.title.x = element_text(size = 50, face = "bold"),
    axis.title.y = element_text(size = 50, face = "bold"),
    axis.text.x = element_text(size = 35, face = "bold", color = "#000000"),
    axis.text.y = element_text(size = 35, face = "bold", color = "#000000"),
    legend.text = element_text(size = 25, face = "bold")
  )


# save the plot
ggsave(paste0(wrdir, "/", gene, "_editing_event_comparison_bar_lvl3_", today, file_type), dpi = 600, width = wi, height = hi, limitsize = FALSE)
### end of level 3 plot ###

### start of level 4 plot ###
ggplot(data = perc_df_long4, aes(x = sample, y = Value, fill = Group)) +
  geom_bar(stat = "identity", position = "stack", width = .5) +
  labs(title = paste0("Comparison of Editing Events\nBetween ", gene, " Cell Lines"), x = "", y = "% of Editing Events") +
  scale_fill_manual(values = c("FE" = "#000000", "NC.edit" = "#E07EFD", "PE" = "#c0c0c0"),
                    labels = c("PE" = "Pre-edited", "FE" = "Fully Edited", "NC.edit" = "Non-canonical Edit")) +
  scale_y_continuous(limits = c(0, 101),
                     breaks = seq(0, 100, by = 10), expand = c(0,0)) +
  theme_minimal() + 
  theme(
    plot.title.position = "plot",
    plot.title = element_text(hjust = .5, size = 45, face = "bold"),
    legend.position = "bottom",
    legend.title = element_blank(),
    axis.title.x = element_text(size = 50, face = "bold"),
    axis.title.y = element_text(size = 50, face = "bold"),
    axis.text.x = element_text(size = 35, face = "bold", color = "#000000"),
    axis.text.y = element_text(size = 35, face = "bold", color = "#000000"),
    legend.text = element_text(size = 25, face = "bold")
  )


# save the plot
ggsave(paste0(wrdir, "/", gene, "_editing_event_comparison_bar_lvl4_", today, file_type), dpi = 600, width = wi, height = hi, limitsize = FALSE)
### end of level 4 plot ###

### start of level 5 plot ###
ggplot(data = perc_df_long5, aes(x = sample, y = Value, fill = Group)) +
  geom_bar(stat = "identity", position = "stack", width = .5) +
  labs(title = paste0("Comparison of Editing Events\nBetween ", gene, " Cell Lines"), x = "", y = "% of Editing Events") +
  scale_fill_manual(values = c(NC.Insertion = "#ffc000", NC.Deletion = "#6ff59b",
                               Insertion1 = "#ff0000", Deletion1 = "#0066ff", No.edit3 = "#e6e7e8"),
                    labels = c(Deletion1 = "Deletion", Insertion1 = "Insertion", NC.Deletion = "Incomplete Deletion",
                               NC.Insertion = "Incomplete Insertion", No.edit3 = "No edit")) +
  scale_y_continuous(limits = c(0, 101),
                     breaks = seq(0, 100, by = 10), expand = c(0,0)) +
  theme_minimal() + 
  theme(
    plot.title.position = "plot",
    plot.title = element_text(hjust = .5, size = 45, face = "bold"),
    legend.position = "bottom",
    legend.title = element_blank(),
    axis.title.x = element_text(size = 50, face = "bold"),
    axis.title.y = element_text(size = 50, face = "bold"),
    axis.text.x = element_text(size = 35, face = "bold", color = "#000000"),
    axis.text.y = element_text(size = 35, face = "bold", color = "#000000"),
    legend.text = element_text(size = 25, face = "bold")
  )


# save the plot
ggsave(paste0(wrdir, "/", gene, "_editing_event_comparison_bar_lvl5_", today, file_type), dpi = 600, width = wi, height = hi, limitsize = FALSE)
### end of level 5 plot ###

### start of making csv summary table for level 1 data ###
# calculate the total number of editing events in each sample
sum_edit_events <- function(df) {
  total_sum <- sum(df[, c(3:4, 6:9, 11:14)])
  return(total_sum)
}

sum_events21 <- lapply(level1, sum_edit_events)
sum_events1 <- data.frame(total_events = unlist(sum_events21))

avg_events1 <- lapply(sum_events1, function(x) x / 1000000)

table_df1 <- cbind(perc_df1, sum_events1, avg_events1 = avg_events1)
column_headers1 <- c("Sample", "% NCI", "% NCD", "% NE", "% CI", "% CII", "% CEI", "% DCI", "% CINE", "% CD", "% CID", "% CED", "% ICD", "% CDNE",
                     "Total number of events", "Average events per normalized read")
table_df1 <- setNames(table_df1, column_headers1)


# write out data frame to a csv file
write.csv(table_df1, paste0(wrdir, "/", gene, "_editing_event_table_lvl1_", today, ".csv"), row.names = FALSE)
### end of making csv summary table for level 1 data ###

### start of making csv summary table for level 2 data ###
# calculate the total number of editing events in each sample
sum_edit_events <- function(df) {
  total_sum <- sum(df[, c(3:4, 6:ncol(df))])
  return(total_sum)
}

sum_events22 <- lapply(level2, sum_edit_events)
sum_events2 <- data.frame(total_events = unlist(sum_events22))

avg_events2 <- lapply(sum_events2, function(x) x / 1000000)

table_df2 <- cbind(perc_df2, sum_events2, avg_events2 = avg_events2)
column_headers2 <- c("Sample", "% NC Incomplete Insertion", "% NC Incomplete Deletion", "% No Edit", "% Insertion", "% Incomplete Insertion",
                     "% Deletion", "% Incomplete Deletion", "Total number of events", "Average events per normalized read")
table_df2 <- setNames(table_df2, column_headers2)

# write out data frame to a csv file
write.csv(table_df2, paste0(wrdir, "/", gene, "_editing_event_table_lvl2_", today, ".csv"), row.names = FALSE)
### end of making csv summary table for level 2 data ###

### start of making csv summary table for level 3 data ###
# calculate the total number of editing events in each sample
sum_edit_events <- function(df) {
  total_sum <- sum(df[,3:4])
  return(total_sum)
}

sum_events23 <- lapply(level3, sum_edit_events)
sum_events3 <- data.frame(total_events = unlist(sum_events23))

avg_events3 <- lapply(sum_events3, function(x) x / 1000000)

table_df3 <- cbind(perc_df3, sum_events3, avg_events3 = avg_events3)
column_headers3 <- c("Sample", "% Insertion", "% Deletion", "% No Edit", "% Total number of events", "% Average events per normalized read")
table_df3 <- setNames(table_df3, column_headers3)

# write out data frame to a csv file
write.csv(table_df3, paste0(wrdir, "/", gene, "_editing_event_table_lvl3_", today, ".csv"), row.names = FALSE)
### end of making csv summary table for level 3 data ###

### start of making csv summary table for level 4 data ###
#since this is a measure of C vs NC we do not sum the number of events
table_df4 <- cbind(perc_df4)
column_headers4 <- c("Sample", "% Canonical", "% Non-canonical")
table_df4 <- setNames(table_df4, column_headers4)

# write out data frame to a csv file
write.csv(table_df4, paste0(wrdir, "/", gene, "_editing_event_table_lvl4_", today, ".csv"), row.names = FALSE)
### end of making csv summary table for level 4 data ###

### start of making csv summary table for level 5 data ###
# calculate the total number of editing events in each sample
sum_edit_events <- function(df) {
  total_sum <- sum(df[,4:7])
  return(total_sum)
}

sum_events25 <- lapply(level5, sum_edit_events)
sum_events5 <- data.frame(total_events = unlist(sum_events25))

avg_events5 <- lapply(sum_events5, function(x) x / 1000000)

table_df5 <- cbind(perc_df5, sum_events5, avg_events5 = avg_events5)
column_headers5 <- c("Sample", "% No Edit", "% Insertion", "% Deletion", "% Incomplete Insertion", "% Incomplete Deletion",
                     "Total number of events", "Average events per normalized read")
table_df5 <- setNames(table_df5, column_headers5)


# write out data frame to a csv file
write.csv(table_df5, paste0(wrdir, "/", gene, "_editing_event_table_lvl5_", today, ".csv"), row.names = FALSE)
### end of making csv summary table for level 5 data ###






