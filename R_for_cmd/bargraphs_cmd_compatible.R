### Script for making bar graphs with contextualized types of editing occurring at individual editing sites


args <- commandArgs(trailingOnly = TRUE)
sample <- strsplit(args[1], ",")[[1]]
primer_f <- args[2]
primer_r <- args[3]
data_list <- strsplit(args[4], ",")[[1]]
wrdir <- strsplit(args[5], ",")[[1]]
today <- args[6]
gene_name <- args[7]


gene_name_title <- paste("Position in", gene_name, "mRNA", sep = " ")


# set height of figure
hi <- 10

#### I advise against changing any of the code below ####

#install.packages("tidyverse")
library("tidyverse")

file_type <- ".pdf"

storage_paths <- as.list(wrdir)

#remove Ts from primer sequences and extract primer lengths
pre_primer_length_f <- gsub("T", "", primer_f) 
primer_length_f <- nchar(pre_primer_length_f)

pre_primer_length_r <- gsub("T", "", primer_r) 
primer_length_r <- nchar(pre_primer_length_r)

#create a list of the paths to call the ordered counts csv data frames
data_paths2 <- as.list(data_list)

#create a list of the data frames
genes_bar <- lapply(data_paths2, read.csv)

#extract cell line information from the input file, raw values used for chart titles
bar_titles <- as.list(sample, "Level")

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

#add some extra characters to cell line information to generate usable file names
bar_file_titles <- bar_titles %>%
  gsub(" ", "_", .) %>%
  trimws() %>%
  paste0("_bargraphs")

#remove forward primers from sequences
slice_dataframe_f <- function(df, primer_length_f) {
  df %>% slice(-(1:primer_length_f))
}
genes_subset_pre <- lapply(genes_bar, slice_dataframe_f, forward_trim)

#remove reverse primers from sequences
#slice_dataframe_r <- function(df, primer_length_r) {
#  n <- nrow(df)
#  df %>% slice(1:(n - primer_length_r))
#}
#genes_subset <- lapply(genes_subset_pre, slice_dataframe_r, reverse_trim)
slice_dataframe_r <- function(df, primer_length_r) {
  n <- nrow(df)
  df %>% slice(1:(n - primer_length_r))
}
genes_subset_pre2 <- lapply(genes_subset_pre, slice_dataframe_r, reverse_trim)

divide_by_thousand <- function(df) {
  df <- df %>%
    mutate(across(-c(1, 2), ~ . / 10000))
}
genes_subset <- lapply(genes_subset_pre2, divide_by_thousand)


#determine the length of each transcript, needed to maintain sequence order later
lengt <- map(genes_subset, ~ .x %>%
               select(Nucleotide) %>%
               nrow())

# pull out the nucleotide and order columns
consensus_columns <- list()
for (df in genes_subset) {
  selected_columns_df <- df[, c("Nucleotide", "Order")]
  consensus_columns <- c(consensus_columns, list(selected_columns_df))
}

# create new data frame for the first level of categories
level1 <- list()
for (df in genes_subset) {
  level1_df <- df[, c("NCI", "NCD", "NE", "CI", "CII", "CEI", "DCI", "CINE", "CD", "CID", "CED", "ICD", "CDNE")]
  level1 <- c(level1, list(level1_df))
}
level1 <- map2(consensus_columns, level1, ~ bind_cols(.x, .y))

# create new data frame for the second level of categories
level2 <- list()
for (df in genes_subset) {
  level2_df <- df[, c("N.inc.ins", "N.inc.del", "No.edit1", "Ins.ed", "I.inc", "Del.ed", "D.inc")]
  level2 <- c(level2, list(level2_df))
}
level2 <- map2(consensus_columns, level2, ~ bind_cols(.x, .y))

# create new data frame for the third level of categories
level3 <- list()
for (df in genes_subset) {
  level3_df <- df[, c("Insertion", "Deletion", "No.edit2")]
  level3 <- c(level3, list(level3_df))
}
level3 <- map2(consensus_columns, level3, ~ bind_cols(.x, .y))

# create new data frame for the fourth level of categories
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

create_file_paths <- function(dirs, titles, date, level, file_type) {
  file_paths <- vector("list", length(dirs))
  
  for (i in seq_along(dirs)) {
    name <- paste0(dirs[[i]], "/", titles[[i]], "_", date, level, file_type)
    file_paths[[i]] <- name
  }
  
  return(file_paths)
}

wi <- 0.75 * lengt[[1]]

######################## Creating bar plots for level 1 data ########################
bar_file_titles1 <- create_file_paths(wrdir, bar_file_titles, today, "_level1", file_type)

#create a function that transforms data into a form that is usable by geom_bar
genes_bar_clean_function_level1 <- function(genes_data, lengt) {
  genes_data %>%
    mutate(lengt = 1:lengt) %>%
    rename("lengt" = "lengt") %>%
    unite(x, Nucleotide, lengt) %>%
    pivot_longer(cols = c(NCI, NCD, NE, CI, CII, CEI, DCI, CINE, CD, CID, CED, ICD, CDNE),
                 names_to = "d",
                 values_to = "c")
}

# Store function results in a list
results_list1 <- list()

# Iterate through the data frames in genes_bar
for (i in seq_along(level1)) {
  result <- genes_bar_clean_function_level1(level1[[i]], lengt[[i]])
  results_list1[[i]] <- result
}

#get the sequence for each gene
ordered_counts_unique_b1 <- map(results_list1, ~ .x %>%
                                  pivot_wider(names_from = d, values_from = c) %>%
                                  separate(x, c("nucleotide", "number")) %>%
                                  select(nucleotide))

#pivot data back to a wide form, useful for ordering the edit types on the bar charts
ordered_counts_unique_c1 <- map(results_list1, ~ .x %>%
                                  pivot_wider(names_from = d, values_from = c)) 

#create bar charts for each gene
create_bar_plot1 <- function(bar, i) {
  ggplot(bar, aes(x = x, y = y)) +
    geom_bar(data = results_list1[[i]], aes(x = factor(x, levels = ordered_counts_unique_c1[[i]]$x), 
                                            fill = factor(d, levels = c("NE", "CINE", "CDNE", "NCI", "CEI", "ICD", "NCD", "DCI", "CED", "CI", "CII", "CD", "CID")), y = c),
             position = "stack", stat = "identity", width = .65) +
    labs(x = gene_name_title, y = "% Reads", title = bar_titles[[i]]) +
    scale_x_discrete(labels = ordered_counts_unique_b1[[i]]$nucleotide) +
    scale_fill_manual(name = "Event",
                      values = c(NE = "#e6e7e8", NCI = "#000000", NCD = "#6ff59b",
                                 CII = "#FFFF00", CID = "#33FFFF",
                                 CI = "#ff0000", CD = "#0066ff", CEI = "#FF7A0D", DCI = "#0fb2f6",
                                 CINE = "#654321", CED = "#9932cc", ICD = "#32cd32", CDNE = "#ffe4b5")) +
    scale_y_continuous(expand = c(0,0)) +
    theme_classic() +
    theme(plot.title.position = "plot",
          plot.title = element_text(hjust = .5, size = 45, face = "bold"),
          axis.title.x = element_text(size = 40, face = "bold"),
          axis.title.y = element_text(size = 40, face = "bold"),
          axis.text.x = element_text(size = 22, face = "bold", color = "#000000", hjust = 1.5),
          axis.text.y = element_text(size = 22, face = "bold", color = "#000000"),
          legend.text = element_text(size = 22),
          legend.title = element_text(size = 30, face = "bold"))
}

#store the bar charts in a list
plot_list1 <- list()
for (i in seq_along(results_list1)) {
  plot_list1[[i]] <- create_bar_plot1(results_list1[[i]], i)
  ggsave(paste(bar_file_titles1[[i]], sep = ""), dpi = 600, width = wi, height = hi, limitsize = FALSE)
}
######################## End for level 1 data ########################

######################## Creating bar plots for level 2 data ########################
bar_file_titles2 <- create_file_paths(wrdir, bar_file_titles, today, "_level2", file_type)

#create a function that transforms data into a form that is usable by geom_bar
genes_bar_clean_function_level2 <- function(genes_data, lengt) {
  genes_data %>%
    mutate(lengt = 1:lengt) %>%
    rename("lengt" = "lengt") %>%
    unite(x, Nucleotide, lengt) %>%
    pivot_longer(cols = c(D.inc, Del.ed, I.inc, Ins.ed, N.inc.ins, No.edit1, N.inc.del),
                 names_to = "d",
                 values_to = "c")
}

# Store function results in a list
results_list2 <- list()

# Iterate through the data frames in genes_bar
for (i in seq_along(level2)) {
  result <- genes_bar_clean_function_level2(level2[[i]], lengt[[i]])
  results_list2[[i]] <- result
}

#get the sequence for each gene
ordered_counts_unique_b2 <- map(results_list2, ~ .x %>%
                                  pivot_wider(names_from = d, values_from = c) %>%
                                  separate(x, c("nucleotide", "number")) %>%
                                  select(nucleotide))

#pivot data back to a wide form, useful for ordering the edit types on the bar charts
ordered_counts_unique_c2 <- map(results_list2, ~ .x %>%
                                  pivot_wider(names_from = d, values_from = c)) 

#create bar charts for each gene
create_bar_plot2 <- function(bar, i) {
  ggplot(bar, aes(x = x, y = y)) +
    geom_bar(data = results_list2[[i]], aes(x = factor(x, levels = ordered_counts_unique_c2[[i]]$x), 
                                            fill = factor(d, levels = c("No.edit1", "Del.ed", "Ins.ed", "D.inc", "I.inc", "N.inc.ins", "N.inc.del")), y = c),
             position = "stack", stat = "identity", width = .65) +
    labs(x = gene_name_title, y = "% Reads", title = bar_titles[[i]]) +
    scale_x_discrete(labels = ordered_counts_unique_b2[[i]]$nucleotide) +
    scale_fill_manual(name = "Event",
                      values = c(No.edit1 = "#e6e7e8", N.inc.ins = "#FF7A0D", N.inc.del = "#888787",
                                 I.inc = "#ffc000", D.inc = "#6ff59b",
                                 Ins.ed = "#ff0000", Del.ed = "#0066ff"),
                      labels = c("No.edit1" = "No Edit", "Del.ed" = "Deletion", "Ins.ed" = "Insertion", "D.inc" = "Incomplete Deletion", 
                                 "I.inc" = "Incomplete Insertion", "N.inc.ins" = "NC Incomplete Insertion", "N.inc.del" = "NC Incomplete Deletion")) +
    scale_y_continuous(expand = c(0,0)) +
    theme_classic() +
    theme(plot.title.position = "plot",
          plot.title = element_text(hjust = .5, size = 45, face = "bold"),
          axis.title.x = element_text(size = 40, face = "bold"),
          axis.title.y = element_text(size = 40, face = "bold"),
          axis.text.x = element_text(size = 22, face = "bold", color = "#000000", hjust = 1.5),
          axis.text.y = element_text(size = 22, face = "bold", color = "#000000"),
          legend.text = element_text(size = 22),
          legend.title = element_text(size = 30, face = "bold"))
}

#store the bar charts in a list
plot_list2 <- list()
for (i in seq_along(results_list2)) {
  plot_list2[[i]] <- create_bar_plot2(results_list2[[i]], i)
  ggsave(paste(bar_file_titles2[[i]], sep = ""), dpi = 600, width = wi, height = hi, limitsize = FALSE)
}
######################## End for level 2 data ########################

######################## Creating bar plots for level 3 data ########################
bar_file_titles3 <- create_file_paths(wrdir, bar_file_titles, today, "_level3", file_type)

#create a function that transforms data into a form that is usable by geom_bar
genes_bar_clean_function_level3 <- function(genes_data, lengt) {
  genes_data %>%
    mutate(lengt = 1:lengt) %>%
    rename("lengt" = "lengt") %>%
    unite(x, Nucleotide, lengt) %>%
    rename("Insertion" = Insertion, "Deletion" = Deletion, "No Edit" = No.edit2) %>%
    pivot_longer(cols = c("Insertion", "Deletion", "No Edit"),
                 names_to = "d",
                 values_to = "c")
}

# Store function results in a list
results_list3 <- list()

# Iterate through the data frames in genes_bar
for (i in seq_along(level3)) {
  result <- genes_bar_clean_function_level3(level3[[i]], lengt[[i]])
  results_list3[[i]] <- result
}

#get the sequence for each gene
ordered_counts_unique_b3 <- map(results_list3, ~ .x %>%
                                  pivot_wider(names_from = d, values_from = c) %>%
                                  separate(x, c("nucleotide", "number")) %>%
                                  select(nucleotide))

#pivot data back to a wide form, useful for ordering the edit types on the bar charts
ordered_counts_unique_c3 <- map(results_list3, ~ .x %>%
                                  pivot_wider(names_from = d, values_from = c)) 

#create bar charts for each gene
create_bar_plot3 <- function(bar, i) {
  ggplot(bar, aes(x = x, y = y)) +
    geom_bar(data = results_list3[[i]], aes(x = factor(x, levels = ordered_counts_unique_c3[[i]]$x), 
                                            fill = factor(d, levels = c("No Edit", "Insertion", "Deletion")), y = c),
             position = "stack", stat = "identity", width = .65) +
    labs(x = gene_name_title, y = "% Reads", title = bar_titles[[i]]) +
    scale_x_discrete(labels = ordered_counts_unique_b3[[i]]$nucleotide) +
    scale_fill_manual(name = "Event",
                      values = c("No Edit" = "#e6e7e8", "Insertion" = "#ff0000", "Deletion" = "#0066ff"),
                      labels = c("No Edit" = "No Edit", "Insertion" = "Insertion", "Deletion" = "Deletion")) +
    scale_y_continuous(expand = c(0,0)) +
    theme_classic() +
    theme(plot.title.position = "plot",
          plot.title = element_text(hjust = .5, size = 45, face = "bold"),
          axis.title.x = element_text(size = 40, face = "bold"),
          axis.title.y = element_text(size = 40, face = "bold"),
          axis.text.x = element_text(size = 22, face = "bold", color = "#000000", hjust = 1.5),
          axis.text.y = element_text(size = 22, face = "bold", color = "#000000"),
          legend.text = element_text(size = 22),
          legend.title = element_text(size = 30, face = "bold"))
}

#store the bar charts in a list
plot_list3 <- list()
for (i in seq_along(results_list3)) {
  plot_list3[[i]] <- create_bar_plot3(results_list3[[i]], i)
  ggsave(paste(bar_file_titles3[[i]], sep = ""), dpi = 600, width = wi, height = hi, limitsize = FALSE)
}
######################## End for level 3 data ########################

######################## Creating bar plots for level 4 data ########################
bar_file_titles4 <- create_file_paths(wrdir, bar_file_titles, today, "_level4", file_type)

#create a function that transforms data into a form that is usable by geom_bar
genes_bar_clean_function_level4 <- function(genes_data, lengt) {
  genes_data %>%
    mutate(lengt = 1:lengt) %>%
    rename("lengt" = "lengt") %>%
    unite(x, Nucleotide, lengt) %>%
    pivot_longer(cols = c("NC.edit", "PE", "FE"),
                 names_to = "d",
                 values_to = "c")
}

# Store function results in a list
results_list4 <- list()

# Iterate through the data frames in genes_bar
for (i in seq_along(level4)) {
  result <- genes_bar_clean_function_level4(level4[[i]], lengt[[i]])
  results_list4[[i]] <- result
}

#get the sequence for each gene
ordered_counts_unique_b4 <- map(results_list4, ~ .x %>%
                                  pivot_wider(names_from = d, values_from = c) %>%
                                  separate(x, c("nucleotide", "number")) %>%
                                  select(nucleotide))

#pivot data back to a wide form, useful for ordering the edit types on the bar charts
ordered_counts_unique_c4 <- map(results_list4, ~ .x %>%
                                  pivot_wider(names_from = d, values_from = c)) 

#create bar charts for each gene
create_bar_plot4 <- function(bar, i) {
  ggplot(bar, aes(x = x, y = y)) +
    geom_bar(data = results_list4[[i]], aes(x = factor(x, levels = ordered_counts_unique_c4[[i]]$x), 
                                            fill = factor(d, levels = c("PE", "NC.edit", "FE")), y = c),
             position = "stack", stat = "identity", width = .65) +
    labs(x = gene_name_title, y = "% Reads", title = bar_titles[[i]]) +
    scale_x_discrete(labels = ordered_counts_unique_b4[[i]]$nucleotide) +
    scale_fill_manual(name = "Event",
                      values = c("FE" = "#000000", "NC.edit" = "#E07EFD", "PE" = "#c0c0c0"),
                      labels = c("PE" = "Pre-edited", "FE" = "Fully Edited", "NC.edit" = "Non-canonical Edit")) +
    scale_y_continuous(expand = c(0,0)) +
    theme_classic() +
    theme(plot.title.position = "plot",
          plot.title = element_text(hjust = .5, size = 45, face = "bold"),
          axis.title.x = element_text(size = 40, face = "bold"),
          axis.title.y = element_text(size = 40, face = "bold"),
          axis.text.x = element_text(size = 22, face = "bold", color = "#000000", hjust = 1.5),
          axis.text.y = element_text(size = 22, face = "bold", color = "#000000"),
          legend.text = element_text(size = 22),
          legend.title = element_text(size = 30, face = "bold"))
}

#store the bar charts in a list
plot_list4 <- list()
for (i in seq_along(results_list4)) {
  plot_list4[[i]] <- create_bar_plot4(results_list4[[i]], i)
  ggsave(paste(bar_file_titles4[[i]], sep = ""), dpi = 600, width = wi, height = hi, limitsize = FALSE)
}
######################## End for level 4 data ########################

######################## Creating bar plots for level 5 data ########################
bar_file_titles5 <- create_file_paths(wrdir, bar_file_titles, today, "_level5", file_type)

#create a function that transforms data into a form that is usable by geom_bar
genes_bar_clean_function_level5 <- function(genes_data, lengt) {
  genes_data %>%
    mutate(lengt = 1:lengt) %>%
    rename("lengt" = "lengt") %>%
    unite(x, Nucleotide, lengt) %>%
    #rename(No.edit3 = "No Edit", Insertion1 = "Insertion", Deletion1 = "Deletion", NC.Insertion = "Incomplete Insertion", NC.Deletion = "Incomplete Deletion") %>%
    pivot_longer(cols = c(No.edit3, Insertion1, Deletion1, NC.Insertion, NC.Deletion),
                 names_to = "d",
                 values_to = "c")
}

# Store function results in a list
results_list5 <- list()

# Iterate through the data frames in genes_bar
for (i in seq_along(level5)) {
  result <- genes_bar_clean_function_level5(level5[[i]], lengt[[i]])
  results_list5[[i]] <- result
}

#get the sequence for each gene
ordered_counts_unique_b5 <- map(results_list5, ~ .x %>%
                                  pivot_wider(names_from = d, values_from = c) %>%
                                  separate(x, c("nucleotide", "number")) %>%
                                  select(nucleotide))

#pivot data back to a wide form, useful for ordering the edit types on the bar charts
ordered_counts_unique_c5 <- map(results_list5, ~ .x %>%
                                  pivot_wider(names_from = d, values_from = c)) 

#create bar charts for each gene
create_bar_plot5 <- function(bar, i) {
  ggplot(bar, aes(x = x, y = y)) +
    geom_bar(data = results_list5[[i]], aes(x = factor(x, levels = ordered_counts_unique_c5[[i]]$x), 
                                            fill = factor(d, levels = c("No.edit3", "Insertion1", "Deletion1", "NC.Insertion", "NC.Deletion")), y = c),
             position = "stack", stat = "identity", width = .65) +
    labs(x = gene_name_title, y = "% Reads", title = bar_titles[[i]]) +
    scale_x_discrete(labels = ordered_counts_unique_b5[[i]]$nucleotide) +
    scale_fill_manual(name = "Event",
                      values = c(No.edit3 = "#e6e7e8", NC.Insertion = "#ffc000", NC.Deletion = "#6ff59b",
                                 Insertion1 = "#ff0000", Deletion1 = "#0066ff"),
                      labels = c(No.edit3 = "No Edit", Insertion1 = "Insertion", Deletion1 = "Deletion",
                                 NC.Insertion = "Incomplete Insertion", NC.Deletion = "Incomplete Deletion")) +
    scale_y_continuous(expand = c(0,0)) +
    theme_minimal() +
    theme(plot.title.position = "plot",
          plot.title = element_text(hjust = .5, size = 45, face = "bold"),
          axis.title.x = element_text(size = 40, face = "bold"),
          axis.title.y = element_text(size = 40, face = "bold"),
          axis.text.x = element_text(size = 22, face = "bold", color = "#000000", hjust = 1.5),
          axis.text.y = element_text(size = 22, face = "bold", color = "#000000"),
          legend.text = element_text(size = 22),
          legend.title = element_text(size = 30, face = "bold"),
          panel.grid.major.x = element_blank(),  
          panel.grid.minor.x = element_blank())
}

#store the bar charts in a list
plot_list5 <- list()
for (i in seq_along(results_list5)) {
  plot_list5[[i]] <- create_bar_plot5(results_list5[[i]], i)
  ggsave(paste(bar_file_titles5[[i]], sep = ""), dpi = 600, width = wi, height = hi, limitsize = FALSE)
}

######################## End for level 5 data ########################