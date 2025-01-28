### Script for making difference plots with contextualized types of editing occurring at individual editing sites

####10/30/24  I need to make the following changes: add the "Position in [????] mRNA" x-axis



# command line arguments
args <- commandArgs(trailingOnly = TRUE)
sample <- strsplit(args[1], ",")[[1]]
primer_f <- args[2]
primer_r <- args[3]
data_list <- strsplit(args[4], ",")[[1]]
wrdir <- args[5]
today <- args[6]
control <- args[7]
y_labels <- strsplit(args[8], ",")[[1]] 
axis_increment <- as.numeric(args[9])
min_bound <- as.numeric(args[10]) 
max_bound <- as.numeric(args[11])
gene_name <- args[12]


gene_name_title <- paste("Position in", gene_name, "mRNA", sep = " ")

# specify the height of the plot
hi <- 15


#### I advise against changing any of the code below ####


control_index <- as.integer(control) + 1


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
genes_bar <- lapply(data_paths2, read.csv)

#extract cell line information from the input file, raw values used for chart titles
bar_titles <- as.list(sample, "Level")

# extract the names
control_name <- bar_titles[[control_index]]
exp_names <- bar_titles[-control_index]

# create new list of plot names
new_names_func <- function(cn, en) {
  new_name <- paste(en, " - ", cn)
  return(new_name)
}
updated_names <- lapply(exp_names, function(en) {
  new_names_func(control_name, en)
})


# determine max value in NE column
max_ne <- max(genes_bar[[1]]$NE)

# take out any forward primer bases that experienced editing
process_dataframe <- function(df) {
  pre_for <- df %>% 
    slice((1:primer_length_f))
  for_slice2 <- filter(pre_for, NE != max(NE))
  return(for_slice2)
}
for_slice <- lapply(genes_bar, process_dataframe)

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
bar_file_titles <- updated_names %>%
  gsub(" ", "_", .) %>%
  trimws() %>%
  paste0("_difference_plot")

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
genes_subset2 <- lapply(genes_subset_pre, slice_dataframe_r, reverse_trim)

#create a list of the data frames
genes_subset_control <- genes_subset2[[control_index]]
genes_subset_experiment <- genes_subset2[-control_index]

subtract_dataframe <- function(df_control, df_experiment) {
  df_control_subset <- df_control[, -c(1, 2)]
  df_experiment_subset <- df_experiment[, -c(1, 2)]
  
  # Perform the subtraction 
  df_result_subset <- (df_experiment_subset - df_control_subset)/10000
  
  # Add the first two columns back to the result from the experiment dataframe
  df_result <- cbind(df_experiment[, 1:2], df_result_subset)
  
  return(df_result)
}
genes_subset <- lapply(genes_subset_experiment, function(df_experiment) {
  subtract_dataframe(genes_subset_control, df_experiment)
})

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


# create new data frame for the fifth level of categories
level5 <- list()
for (df in genes_subset) {
  level5_df <- df[, c("No.edit3", "Insertion1", "Deletion1", "NC.Insertion", "NC.Deletion")]
  level5 <- c(level5, list(level5_df))
}
level5 <- map2(consensus_columns, level5, ~ bind_cols(.x, .y))

create_file_paths <- function(dirs, titles, date, level, file_type) {
  file_paths <- vector("list", length(dirs))
  
  for (i in seq_along(titles)) {
    name <- paste0(dirs, "/", titles[[i]], "_", date, level, file_type)
    file_paths[[i]] <- name
  }
  
  return(file_paths)
}

# determine necessary width for the plot
wi <- 0.75 * lengt[[1]]


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
    labs(x = gene_name_title, y = "% Reads", title = updated_names[[i]]) +
    scale_x_discrete(labels = ordered_counts_unique_b5[[i]]$nucleotide) +
    scale_fill_manual(values = c(No.edit3 = "#e6e7e8", NC.Insertion = "#ffc000", NC.Deletion = "#6ff59b",
                                 Insertion1 = "#ff0000", Deletion1 = "#0066ff"),
                      labels = c("No Edit", "Insertion", "Deletion", "Incomplete Insertion", "Incomplete Deletion")) +
    scale_y_continuous(limits = c(min_bound, max_bound),
                       breaks = seq(min_bound , max_bound, by = axis_increment), labels = y_labels) +
    theme_minimal() +
    theme(plot.title.position = "plot",
          plot.title = element_text(hjust = .5, size = 45, face = "bold"),
          axis.title.x = element_text(size = 40, face = "bold"),
          axis.title.y = element_text(size = 40, face = "bold"),
          axis.text.x = element_text(size = 22, face = "bold", color = "#000000", hjust = 1.5),
          axis.text.y = element_text(size = 22, face = "bold", color = "#000000"),
          legend.text = element_text(size = 22),
          legend.position = "bottom",
          legend.title = element_blank(),
          panel.grid.major.x = element_blank(),  
          panel.grid.minor.x = element_blank()
    )
}

#store the bar charts in a list
plot_list5 <- list()
for (i in seq_along(results_list5)) {
  plot_list5[[i]] <- create_bar_plot5(results_list5[[i]], i)
  ggsave(paste(bar_file_titles5[[i]], sep = ""), dpi = 600, width = wi, height = hi, limitsize = FALSE)
}

######################## End for level 5 data ########################