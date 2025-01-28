### script for generating bubble plots from Ordered Counts csv output from PARERS2



#install.packages("tidyverse")
library("tidyverse")

#arguments for the command line
args <- commandArgs(trailingOnly = TRUE)
sample <- strsplit(args[1], ",")[[1]]
primer_f <- args[2]
primer_r <- args[3]
data_list <- strsplit(args[4], ",")[[1]]
wrdir <- strsplit(args[5], ",")[[1]]
today <- args[6]
gene_name <- args[7]

gene_name_title <- paste("Position in", gene_name, "mRNA", sep = " ")

# set height of plot. 10 is sufficient for the vast majority of cases
hi <- 10







##### I advise against changing any of the code below #####

file_type <- ".pdf"

#remove Ts from primer sequences and acquire the lengths, used for remove rows containing primer sequences
pre_primer_length_f <- gsub("T", "", primer_f) 
primer_length_f <- nchar(pre_primer_length_f)

pre_primer_length_r <- gsub("T", "", primer_r) 
primer_length_r <- nchar(pre_primer_length_r)

#create a list of the paths to call the ordered counts csv data frames
data_paths <- as.list(data_list) 

#create a list of the data frames
ordered_counts_b <- lapply(data_paths, read.csv)

#extract cell line information from the input file, raw values used for chart titles
chart_titles <- as.list(sample)

bubble_file_titles <- chart_titles %>%
  gsub(" ", "_", .) %>%
  trimws() %>%
  paste0("_bubble-plots")

#add some extra characters to cell line information to generate usable file names
create_file_paths <- function(dirs, titles, date, level, file_type) {
  file_paths <- vector("list", length(dirs))
  
  for (i in seq_along(dirs)) {
    name <- paste0(dirs[[i]], "/", titles[[i]], "_", date, level, file_type)
    file_paths[[i]] <- name
  }
  
  return(file_paths)
}

max_orders <- sapply(ordered_counts_b, function(df) max(df$Order, na.rm = TRUE))
overall_max <- max(max_orders, na.rm = TRUE)
rev_prim_index <- overall_max - primer_length_r


# take out those rows that are supposedly reverse primers only
oc_unique_rev <- function(df, rev_prim_index) {
  df %>%
    subset(Order >= rev_prim_index) %>%
    group_by(Order) %>%
    filter(n() > 1)
}
rev_not_uni_pre <- lapply(ordered_counts_b, oc_unique_rev, rev_prim_index = rev_prim_index)
rev_max_order <- max(sapply(rev_not_uni_pre, function(df) max(df[["Order"]], na.rm = TRUE)), na.rm = TRUE) # find rev primer cutoff

# make new rev_not_uni that are normalized across samples
oc_unique_rev2 <- function(df, rev_prim_index) {
  df %>%
    subset(Order >= rev_prim_index) %>%
    group_by(Order) %>%
    filter(Order <= rev_max_order)
}
rev_not_uni <- lapply(ordered_counts_b, oc_unique_rev2, rev_prim_index = rev_prim_index)

# take out those rows that are supposedly forward primers only
oc_subset_for <- function(df, primer_length_f) {
  df %>%
    subset(Order <= primer_length_f) %>%
    group_by(Order) %>%
    filter(n() > 1)
}
for_not_uni_pre <- lapply(ordered_counts_b, oc_subset_for, primer_length_f = primer_length_f)
for_min_order <- min(sapply(for_not_uni_pre, function(df) min(df[["Order"]], na.rm = TRUE)), na.rm = TRUE) # find forward primer cutoff

# make new for_not_uni that are normalized across samples
oc_subset_for2 <- function(df, primer_length_f) {
  df %>%
    subset(Order <= primer_length_f) %>%
    group_by(Order) %>%
    filter(Order >= for_min_order)
}
for_not_uni <- lapply(ordered_counts_b, oc_subset_for2, primer_length_f = primer_length_f)

#remove rows that belong to the primers while retaining all the rest of the data
filter_dataframe <- function(df, primer_length_f, rev_prim_index) {
  df %>%
    subset(Order > primer_length_f & Order < rev_prim_index)
}
filtered_df_list_pre <- lapply(ordered_counts_b, filter_dataframe, primer_length_f = primer_length_f, rev_prim_index = rev_prim_index)

# bind rows back on
filtered_df_list <- mapply(function(df1, df2, df3) {
  rbind(df1, df2, df3)
}, for_not_uni, filtered_df_list_pre, rev_not_uni, SIMPLIFY = FALSE)

#determine the non-T nucleotide sequence for each gene, this includes the primers
oc_unique <- filtered_df_list[[1]] %>%
  distinct(Order, .keep_all = TRUE) %>%
  select(Nucleotide)

# number of non-t nucleotides
nont_count <- nrow(oc_unique)

max_nums<- sapply(ordered_counts_b, function(df) max(df$Number, na.rm = TRUE))
y_max <- max(max_nums, na.rm = TRUE) 

size_breaks <- c(750000, 500000, 250000, 125000, 62500, 30000, 10000, 5000)
size_labels <- c("75%", "50%", "25%", "12.5%", "6.25%", "3%", "1%", "0.5%")

wi <- .75 * nrow(oc_unique)

############## Start of processing for level 1 ##############

bubble_file_titles1 <- create_file_paths(wrdir, bubble_file_titles, today, "_level1", file_type)

#create dot plots for each gene
create_dot_plot1 <- function(df) {
  ggplot(df, aes(x = x, y = y)) +
    geom_point(data=filtered_df_list[[i]], aes(x = factor(Order), y = Number, size = NewCount, color = factor(Lvl1))) + 
    scale_size_area(breaks = size_breaks, labels = size_labels, max_size = 13) + 
    scale_x_discrete(labels = oc_unique$Nucleotide) +
    scale_y_continuous(breaks=c(0:max(filtered_df_list[[i]]$Number)), limits=c(0, y_max)) +
    labs(x = gene_name_title, y = "# of U", color = "Observation Frequency", 
         size = "Observation Frequency", title = chart_titles[[i]]) +
    scale_color_manual(name = "Event",
                       values = c(NE = "#e6e7e8", NCI = "#000000", NCD = "#6ff59b",
                                  CII = "#FFFF00", CID = "#33FFFF",
                                  CI = "#ff0000", CD = "#0066ff", CEI = "#FF7A0D", DCI = "#0fb2f6",
                                  CINE = "#654321", CED = "#9932cc", ICD = "#32cd32", CDNE = "#ffe4b5")) +
    theme_classic() +
    theme(plot.title.position = "plot",
          plot.title = element_text(hjust = .5, size = 45, face = "bold"),
          axis.title.x = element_text(size = 40, face = "bold"),
          axis.title.y = element_text(size = 40, face = "bold"),
          axis.text.x = element_text(size = 22, face = "bold", color = "#000000", hjust = 1.5),
          axis.text.y = element_text(size = 22, face = "bold", color = "#000000"),
          axis.ticks.x = element_blank(),
          legend.text = element_text(size = 22),
          legend.title = element_text(size = 30, face = "bold")) +
    guides(color = guide_legend(title = "Event", override.aes = list(size = 6)),
           size = guide_legend(title = "Observation\nFrequency", ncol = 2))
}


#store the dot plots in a list
plot_list1 <- list()

#iterate through the plot list and automatically save to the file location designated in line 27
for (i in 1:length(ordered_counts_b)) {
  plot_list1[[i]] <- create_dot_plot1(ordered_counts_b[[i]])
  
  ggsave(paste(bubble_file_titles1[[i]], sep = ""), dpi = 600, width = wi, height = hi, limitsize = FALSE)
}
############## End of processing for level 1 ##############

############## Start of processing for level 2 ##############
bubble_file_titles2 <- create_file_paths(wrdir, bubble_file_titles, today, "_level2", file_type)

create_dot_plot2 <- function(df) {
  ggplot(df, aes(x = x, y = y)) +
    geom_point(data=filtered_df_list[[i]], aes(x = factor(Order), y = Number, size = NewCount, color = factor(Lvl2))) + 
    scale_size_area(breaks = size_breaks, labels = size_labels, max_size = 13) + 
    scale_x_discrete(labels = oc_unique$Nucleotide) +
    scale_y_continuous(breaks=c(0:max(filtered_df_list[[i]]$Number)), limits=c(0, y_max)) +
    labs(x = gene_name_title, y = "# of U", color = "Observation Frequency", 
         size = "Observation Frequency", title = chart_titles[[i]]) +
    scale_color_manual(name = "Event",
                       values = c("No edit1" = "#e6e7e8", "N-inc-ins" = "#FF7A0D", "N-inc-del" = "#888787",
                                  "I-inc" = "#ffc000", "D-inc" = "#6ff59b",
                                  "Ins ed" = "#ff0000", "Del ed" = "#0066ff"),
                       labels = c("D-inc" = "Incomplete Deletion", "Del ed" = "Deletion", "I-inc" = "Incomplete Insertion", "Ins ed" = "Insertion", 
                                  "N-inc-del" = "NC Incomplete Deletion", "N-inc-ins" = "NC Incomplete Insertion", "No edit1" = "No Edit")) +
    theme_classic() +
    theme(plot.title.position = "plot",
          plot.title = element_text(hjust = .5, size = 45, face = "bold"),
          axis.title.x = element_text(size = 40, face = "bold"),
          axis.title.y = element_text(size = 40, face = "bold"),
          axis.text.x = element_text(size = 22, face = "bold", color = "#000000"),
          axis.text.y = element_text(size = 22, face = "bold", color = "#000000"),
          legend.text = element_text(size = 22),
          legend.title = element_text(size = 30, face = "bold")) +
    guides(color = guide_legend(title = "Event", override.aes = list(size = 6)),
           size = guide_legend(title = "Observation\nFrequency", ncol = 2))
}


#store the dot plots in a list
plot_list2 <- list()

#iterate through the plot list and automatically save to the file location designated in line 27
for (i in 1:length(ordered_counts_b)) {
  plot_list2[[i]] <- create_dot_plot2(ordered_counts_b[[i]])
  
  ggsave(paste(bubble_file_titles2[[i]], sep = ""), dpi = 600, width = wi, height = hi, limitsize = FALSE)
}
############## End of processing for level 2 ##############

############## Start of processing for level 3 ##############
bubble_file_titles3 <- create_file_paths(wrdir, bubble_file_titles, today, "_level3", file_type)

create_dot_plot3 <- function(df) {
  ggplot(df, aes(x = x, y = y)) +
    geom_point(data=filtered_df_list[[i]], aes(x = factor(Order), y = Number, size = NewCount, color = factor(Lvl3))) + 
    scale_size_area(breaks = size_breaks, labels = size_labels, max_size = 13) + 
    scale_x_discrete(labels = oc_unique$Nucleotide) +
    scale_y_continuous(breaks=c(0:max(filtered_df_list[[i]]$Number)), limits=c(0, y_max)) +
    labs(x = gene_name_title, y = "# of U", color = "Observation Frequency", 
         size = "Observation Frequency", title = chart_titles[[i]]) +
    scale_color_manual(name = "Event",
                       values = c("No edit2" = "#e6e7e8", "Insertion" = "#ff0000", "Deletion" = "#0066ff"),
                       labels = c("Deletion", "Insertion", "No Edit")) +
    theme_classic() +
    theme(plot.title.position = "plot",
          plot.title = element_text(hjust = .5, size = 45, face = "bold"),
          axis.title.x = element_text(size = 40, face = "bold"),
          axis.title.y = element_text(size = 40, face = "bold"),
          axis.text.x = element_text(size = 22, face = "bold", color = "#000000"),
          axis.text.y = element_text(size = 22, face = "bold", color = "#000000"),
          legend.text = element_text(size = 22),
          legend.title = element_text(size = 30, face = "bold")) +
    guides(color = guide_legend(title = "Event", override.aes = list(size = 6)),
           size = guide_legend(title = "Observation\nFrequency", ncol = 2))
}


#store the dot plots in a list
plot_list3 <- list()

#iterate through the plot list and automatically save to the file location designated in line 27
for (i in 1:length(ordered_counts_b)) {
  plot_list3[[i]] <- create_dot_plot3(ordered_counts_b[[i]])
  
  ggsave(paste(bubble_file_titles3[[i]], sep = ""), dpi = 600, width = wi, height = hi, limitsize = FALSE)
}
############## End of processing for level 3 ##############

############## Start of processing for level 4 ##############
bubble_file_titles4 <- create_file_paths(wrdir, bubble_file_titles, today, "_level4", file_type)

create_dot_plot4 <- function(df) {
  ggplot(df, aes(x = x, y = y)) +
    geom_point(data=filtered_df_list[[i]], aes(x = factor(Order), y = Number, size = NewCount, color = factor(Lvl4))) + 
    scale_size_area(breaks = size_breaks, labels = size_labels, max_size = 13) + 
    scale_x_discrete(labels = oc_unique$Nucleotide) +
    scale_y_continuous(breaks=c(0:max(filtered_df_list[[i]]$Number)), limits=c(0, y_max)) +
    labs(x = gene_name_title, y = "# of U", color = "Observation Frequency", 
         size = "Observation Frequency", title = chart_titles[[i]]) +
    scale_color_manual(name = "Event",
                       values = c("FE" = "#000000", "PE" = "#c0c0c0", "NC edit" = "#E07EFD"),
                       labels = c("NC edit" = "Non-canonical Edit", "PE" = "Pre-edited", "FE" = "Fully Edited")) +
    theme_classic() +
    theme(plot.title.position = "plot",
          plot.title = element_text(hjust = .5, size = 45, face = "bold"),
          axis.title.x = element_text(size = 40, face = "bold"),
          axis.title.y = element_text(size = 40, face = "bold"),
          axis.text.x = element_text(size = 22, face = "bold", color = "#000000"),
          axis.text.y = element_text(size = 22, face = "bold", color = "#000000"),
          legend.text = element_text(size = 22),
          legend.title = element_text(size = 30, face = "bold")) +
    guides(color = guide_legend(title = "Event", override.aes = list(size = 6)),
           size = guide_legend(title = "Observation\nFrequency", ncol = 2))
}


#store the dot plots in a list
plot_list4 <- list()

#iterate through the plot list and automatically save to the file location designated in line 27
for (i in 1:length(ordered_counts_b)) {
  plot_list4[[i]] <- create_dot_plot4(ordered_counts_b[[i]])
  
  ggsave(paste(bubble_file_titles4[[i]], sep = ""), dpi = 600, width = wi, height = hi, limitsize = FALSE)
}
############## End of processing for level 4 ##############

############## Start of processing for level 5 ##############
bubble_file_titles5 <- create_file_paths(wrdir, bubble_file_titles, today, "_level5", file_type)

create_dot_plot5 <- function(df) {
  ggplot(df, aes(x = x, y = y)) +
    geom_point(data=filtered_df_list[[i]], aes(x = factor(Order), y = Number, size = NewCount, color = factor(Lvl5))) + 
    scale_size_area(breaks = size_breaks, labels = size_labels, max_size = 13) + 
    scale_x_discrete(labels = oc_unique$Nucleotide) +
    scale_y_continuous(breaks=c(0:y_max), limits=c(0, y_max)) +
    labs(x = gene_name_title, y = "# of U", size = "Frequency", color = "Frequency", 
         title = chart_titles[[i]]) +
    scale_color_manual(name = "Event",
                       values = c("No edit3" = "#e6e7e8", "NC Insertion" = "#000000", "NC Deletion" = "#000000",
                                  "Insertion1" = "#ff0000", "Deletion1" = "#0066ff"),
                       labels = c("Deletion1" = "Deletion", "Insertion1" = "Insertion", "NC Insertion" = "Incomplete Insertion",
                                  "NC Deletion" = "Incomplete Deletion", "No edit3" = "No Edit")) +
    theme_classic() +
    theme(plot.title.position = "plot",
          plot.title = element_text(hjust = .5, size = 45, face = "bold"),
          axis.title.x = element_text(size = 40, face = "bold"),
          axis.title.y = element_text(size = 40, face = "bold"),
          axis.text.x = element_text(size = 22, face = "bold", color = "#000000", hjust = 1.5),
          axis.text.y = element_text(size = 22, face = "bold", color = "#000000"),
          axis.ticks.x = element_blank(),
          legend.text = element_text(size = 22),
          legend.title = element_text(size = 30, face = "bold")) +
    guides(color = guide_legend(title = "Event", override.aes = list(size = 6)),
           size = guide_legend(title = "Observation\nFrequency", ncol = 2))
}


#store the dot plots in a list
plot_list5 <- list()

#iterate through the plot list and automatically save to the file location designated in line 27
for (i in 1:length(ordered_counts_b)) {
  plot_list5[[i]] <- create_dot_plot5(ordered_counts_b[[i]])
  
  ggsave(paste(bubble_file_titles5[[i]], sep = ""), dpi = 600, width = wi, height = hi, limitsize = FALSE)
}
############## End of processing for level 5 ##############
