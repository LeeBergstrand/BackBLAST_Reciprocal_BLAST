#!/usr/bin/env Rscript
# CombineBlastTables.R
# Copyright Jackson M. Tsuji, 2018
# Neufeld Research Group
# See required R packages below.

#####################################################
## User variables: #################################
RUN_COMMAND_LINE <- TRUE # If selected, all user input here is ignored, and terminal-based input is expected instead.

# Set other user variables here
if (RUN_COMMAND_LINE == FALSE) {
  setwd("/Users/JTsuji/Downloads") # your working directory where files are stored
  input_filenames <- c("subject1.csv", "subject2.csv", "subject3.csv") # NOTE: filename minux extesion will become the subject_name!
  output_filename <- "combined_blast_tables.csv"
}
#####################################################

#####################################################
## Load required packages: ##########################
library(getopt)
library(plyr)
suppressMessages(library(dplyr))
#####################################################

# Hard-coded variables
header_names <- c("qseqid", "sseqid", "pident", "evalue", "qcovhsp", "bitscore")


# Function to assign command line input to variables, or throw help message
parse_command_line_input <- function() {
  
  inputArgs <- commandArgs(trailingOnly = TRUE)
  if (length(inputArgs) < 2) {
    cat("CombineBlastTables.R: Combines BLAST tables and makes new column for sample ID based on the filenames.\n")
    cat("Copyright Neufeld Lab, 2018\n")
    cat("Contact Jackson M. Tsuji (jackson.@uwaterloo.ca) for bug reports or feature requests.\n\n")
    
    cat("Usage: CombineBlastTables.R subject1.csv subject2.csv ... subjectN.csv combined_blast_tables.csv\n\n")
    
    cat("Details:\n",
        "subject1.csv, etc.: CSV BLAST tables for all individual samples (subjects for BLAST). Name should be the subject name [Required]\n",
        "combined_blast_tables.csv: the output CSV file. [Required]\n\n")
        
    quit(status = 1)
  }
  
  num_args <- length(inputArgs)
  
  # Make variables from provided input and save as global variables (<<-)
  input_filenames <<- inputArgs[1:(num_args-1)]
  output_filename <<- inputArgs[num_args]
  
}

# Helper function for timestamps
ts <- function() {
  
  datetime_utc <- format(as.POSIXlt(Sys.time(), tz = "UTC"), "%a %d %b %Y %H:%M:%S %Z")
  
  date_message <- paste("[ ", datetime_utc, " ]: ", sep = "")
  
  return(date_message)
}


# Function to load table of data, add header IDs, and add query/subject name by parsing filename
load_individual_table <- function(table_filename, header_names) {
  cat(paste(ts(), "Reading table '", table_filename, "'\n", sep = ""))
  
  # Read the table and add columns
  data_table <- read.table(table_filename, header = FALSE, sep = ",", stringsAsFactors = FALSE, comment.char = "")
  colnames(data_table) <- header_names
  
  # Remove CSV ending
  table_filename_base <- gsub(pattern = ".csv$", replacement = "", x = table_filename)
  
  # Add new columns to table
  # TODO - consider also adding query_name
  # data_table$query_name <- query_name
  cat(paste(ts(), "Adding subject_name of '", subject_name, "'\n", sep = ""))
  data_table$subject_name <- subject_name
  
  # Re-order columns to be more user friendly
  data_table <- cbind(data_table[,ncol(data_table)], data_table[,1:(ncol(data_table)-1)])
  
  return(data_table)
  
}

main <- function() {
  # Run command line version if requested
  if (RUN_COMMAND_LINE == TRUE) {
    parse_command_line_input()
  }
  
  # Startup messages
  cat(paste(ts(), "Running CombineBlastTables.R\n", sep = ""))
  cat(paste(ts(), "Input tables: ", length(input_filenames), " total\n", sep = ""))
  cat(paste(ts(), "Output table: ", output_filename, "\n", sep = ""))
  
  # Load all tables
  cat(paste(ts(), "Loading BLAST tables\n", sep = ""))
  blast_tables <- lapply(input_tables, function(x) { load_individual_table(x, header_names) })
  
  cat(paste(ts(), "Combining BLAST tables\n", sep = ""))
  output_table <- dplyr::bin_rows(blast_tables)
  
  cat(paste(ts(), "Writing combining BLAST table to file (no headers, for convention)\n", sep = ""))
  write.table(output_table, file = output_filename, sep = ",", row.names = FALSE, col.names = FALSE)
  
  cat(paste(ts(), "Done\n", sep = ""))
}

main()

  