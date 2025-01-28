### Script for generating plots that display type of editing by percent




# command line arguments
args <- commandArgs(trailingOnly = TRUE)
pe_seq <- args[1]
fe_seq <- args[2]
data_list <- strsplit(args[3], ",")[[1]]
wrdir <- args[4]
today <- args[5]
gene <- args[6]



# specify height and width of plot
hi <- 15
wi <- 5 * length(data_list)


#### I advise against changing any of the code below ####

if (wi < 15 ){
  wi == 25
}

#install.packages("tidyverse")
library("tidyverse")

file_type <- ".pdf"

#create a list of the paths to call the ordered counts csv data frames
data_paths2 <- as.list(data_list)
genes_bar <- lapply(data_paths2, read.csv)

# extract the sample names from each cir file
extract_name <- function(df) {
  col_and_row <- df[1,1]
  sliced_string <- sub(":.*", "", col_and_row)
  return(sliced_string)
}
sample_names2 <- lapply(genes_bar, extract_name)

remove_cell_name <- function(name) {
  new_name <- gsub(gene, "", name) %>%
    trimws()
  return(new_name)
}
sample_names <- lapply(sample_names2, remove_cell_name)

# extract the frequency relating to the pe sequence and remove that row, forming a new dataframe
find_pe <- function(df) {
  matched_rows <- df[, 2] == pe_seq
  
  if (any(matched_rows)) {
    result <- df[matched_rows, "Frequency"]
    new_df <- df[!matched_rows, ]  
  } else {
    result <- 0
    new_df <- df  
  }
  
  return(list(result = result, new_df = new_df))
}
pe_data_lists <- lapply(genes_bar, find_pe)
# extract the new dataframe and frequencies
pe_freq_list <- lapply(pe_data_lists, function(x) x$result)
non_pe_dfs <- lapply(pe_data_lists, function(x) x$new_df)


# repeat the above for the fe seq
find_fe <- function(df) {
  matched_rows <- df[, 2] == fe_seq
  
  if (any(matched_rows)) {
    result <- df[matched_rows, "Frequency"]
    new_df <- df[!matched_rows, ]  
  } else {
    result <- 0
    new_df <- df  
  }
  
  return(list(result = result, new_df = new_df))
}
fe_data_lists <- lapply(non_pe_dfs, find_fe)
# extract the new dataframe and frequencies
fe_freq_list <- lapply(fe_data_lists, function(x) x$result)
no_pe_no_fe_df <- lapply(fe_data_lists, function(x) x$new_df)

# sum up the frequencies of all the partially edited sequences
sum_others <- function(df) {
  summed <- sum(df[, 5])
  return(summed)
}
summed_others <- lapply(no_pe_no_fe_df, sum_others)

# convert the lists to a dataframe
bar_start2 <- data.frame(unlist(sample_names), unlist(pe_freq_list), unlist(fe_freq_list), unlist(summed_others))
names(bar_start2) <- c("Samples", "PE", "FE", "Others")

# Convert to long format
bar_start_long <- bar_start2 %>%
  as.data.frame() %>%
  pivot_longer(
    cols = -Samples,          # Pivot all columns except the 'Sample' column
    names_to = "Group",     # This will be the original column names (Pre-Edited, Fully Edited, etc.)
    values_to = "Value"     # Values will be the actual numeric values
  )

# make the plot
ggplot(data = bar_start_long, aes(x = Samples, y = Value, fill = Group)) +
  geom_bar(stat = "identity", position = "stack", width = .5) +
  labs(title = paste0("Comparison of Editing\nBetween ", gene, " Cell Lines"), x = "", y = "% of Sequences") +
  scale_fill_manual(values = c(PE = "#BFBFBF", FE = "#FF0000", Others = "#8FAADC"),
                    labels = c("% Fully\nEdited   ", "% Partially\nEdited   ", "% Pre-edited   ")) +
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
    legend.text = element_text(size = 35)
  )


# save the plot
ggsave(paste0(wrdir, "/", gene, "_editing_extent_plot_", today, file_type), dpi = 600, width = wi, height = hi, limitsize = FALSE)

