#!/usr/bin/env Rscript

# BackBLAST_combine_BLAST_tables.R
# Copyright Lee H. Bergstrand and Jackson M. Tsuji, 2019
# Combines BLAST output (CSV) tables into a single output table with headers
# Part of the BackBLAST pipeline

# Load libraries
library(getopt)
library(futile.logger)
library(tools)
library(plyr)
library(dplyr, warn.conflicts = FALSE)

# Hard-coded variables
HEADER_NAMES <- c("qseqid", "sseqid", "pident", "evalue", "qcovhsp", "bitscore")

# Function to assign command line input to variables, or throw help message
parse_command_line_input <- function(commandArgs) {
  
  if (length(commandArgs) < 2) {
    cat("BackBLAST_combine_BLAST_tables.R: Combines BLAST tables and makes new column for sample ID based on the filenames.\n")
    cat("Copyright Lee Bergstrand and Jackson M. Tsuji, 2019\n")
    cat("Part of the BackBLAST pipeline.\n\n")
    
    cat("Usage: BackBLAST_combine_BLAST_tables.R subject1.csv subject2.csv ... subjectN.csv combined_blast_tables.csv\n\n")
    
    cat("Details:\n",
        "subject1.csv, etc.: CSV BLAST tables for all individual samples (subjects for BLAST). Name should be the subject name [Required]\n",
        "combined_blast_tables.csv: the output CSV file. [Required]\n\n")
        
    quit(status = 1)
  }
  
  num_args <- length(commandArgs)
  
  # Return output params
  params <- list()
  params$input_filenames <- commandArgs[1:(num_args-1)]
  params$output_filename <- commandArgs[num_args]
  
  return(params)
}

# Function to load table of data, add header IDs, and add query/subject name by parsing filename
load_individual_table <- function(table_filename, header_names) {
  flog.info(paste("Reading table '", table_filename, "'", sep = ""))
  
  # Read the table and add columns
  data_table <- read.table(table_filename, header = FALSE, sep = ",", stringsAsFactors = FALSE, comment.char = "")
  colnames(data_table) <- header_names
  
  # Remove CSV ending and folder path. Set to subject_name.
  subject_name <- tools::file_path_sans_ext(basename(table_filename))
  
  # Add new columns to table
  # TODO - consider also adding query_name
  # data_table$query_name <- query_name
  flog.info(paste("Adding subject_name of '", subject_name, "'", sep = ""))
  data_table$subject_name <- subject_name
  
  # Re-order columns to be more user friendly
  data_table <- dplyr::select(data_table, subject_name, everything())
  
  return(data_table)
  
}

main <- function() {
  # Parse command line input
  params <- parse_command_line_input(commandArgs(trailingOnly = TRUE))
  
  # Startup messages
  flog.info("Running BackBLAST_combine_BLAST_tables.R")
  flog.info(paste("Input tables: ", length(params$input_filenames), " total", sep = ""))
  flog.info(paste("Output table: '", params$output_filename, "'", sep = ""))
  
  # Load all tables
  flog.info("Loading BLAST tables")
  blast_tables <- lapply(params$input_filenames, function(x) { load_individual_table(x, HEADER_NAMES) })
  
  flog.info("Combining BLAST tables")
  output_table <- dplyr::bind_rows(blast_tables)
  
  flog.info("Writing combining BLAST table to file (**with headers**)")
  write.table(output_table, file = params$output_filename, sep = ",", row.names = FALSE, col.names = TRUE)
  
  flog.info("BackBLAST_combine_BLAST_tables.R: Done.")
}

main()
