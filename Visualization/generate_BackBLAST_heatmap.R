#!/usr/bin/env Rscript
# generate_BackBLAST_heatmap.R
# Copyright Lee Bergstrand and Jackson M. Tsuji, 2018
# Plots a newick treefile and BLAST table together as a phylogenetic tree and heatmap
# See required R packages below.

#####################################################
## User variables: #################################
RUN_COMMAND_LINE <- TRUE # If selected, all user input here is ignored, and terminal-based input is expected instead.

# Set other user variables here
if (RUN_COMMAND_LINE == FALSE) {
  # See proper descriptions of variables under parse_command_line_input
  
  # Required inputs
  setwd("/home/jmtsuji/Research_General/Bioinformatics/02_git/BackBLAST_Reciprocal_BLAST/test/")
  input_phylogenetic_tree_filepath <- "GSB_riboproteins_tree_vs3.treefile"
  input_blast_table_filename <- "GSB_pathway_vs7_e40_best_hits_MOD3.csv"
  output_pdf_name <- "test.pdf"
  
  # Optional inputs (set to 'NA' to ignore)
  tree_metadata_filename <- NA
  tree_decorator_colname <- NA
  map_genome_names <- NA
  gene_naming_table_filename <- NA
  bootstrap_cutoff <- NA
  root_name <- NA #"Ignavibacterium_album_JCM_16511_NC_017464.1" # Optional; set to NA if you want to use the tree as-is.
    
  # To add...
  # tree_naming_mapping <- ""
  
}
#####################################################

#####################################################
## Load required packages: ##########################
library(getopt)
library(plyr)
suppressMessages(library(tidyverse))
suppressMessages(library(ggtree))
suppressMessages(library(glue))
suppressMessages(library(egg))
# library(RColorBrewer) # Loaded as part of tidyverse?
# library(phytools) # for re-rooting; optional
#####################################################

# Function: assign command line input to variables, or throw help message
parse_command_line_input <- function() {
  
  # Define flags
  params <- matrix(c('tree_filepath', 'i', 1, "character",
                     'blast_table_filepath', 'j', 1, "character",
                     'output_filepath', 'o', 1, "character",
                     'tree_metadata_filename', 'm', 2, "character",
                     'tree_decorator_colname', 'd', 2, "character",
                     'genome_plotting_names', 'n', 2, "logical",
                     'gene_naming_table_filename', 'g', 2, "character",
                     'bootstrap_cutoff', 'b', 2, "character",
                     'root_name', 'r', 2, "character",
                     'help', 'h', 2, "character"), byrow=TRUE, ncol=4)
  
  opt <- getopt(params)
  
  # If help was called, print help message and exit
  if ( !is.null(opt$help) ) {
    
    cat("generate_BackBLAST_heatmap.R: Binds a phylogenetic tree to a BLAST table heatmap.\n")
    cat("Copyright Lee Bergstrand and Jackson M. Tsuji, 2018\n\n")
    
    cat(getopt(params, usage = TRUE))
    
    cat("\n")
    
    message(glue::glue("
              Required inputs:
                --tree_filepath               Filepath for newick-format phylogenetic tree of the BLAST subject organisms
                --blast_table_filepath        Filepath for CSV-format BLAST hit table from CombineBlastTables.R
                --output_filepath             Output filepath for the PDF

              Optional inputs:
                --tree_metadata_filename      Filepath for TSV-format metadata file for the phylogenetic tree. Details below.
                --tree_decorator_colname      Column name from the metadata to map onto the tree as fill colours.
                                                Requires that --tree_metadata_filename is set. Details below.
                --genome_plotting_names       Set this flag to include a vector of genome names to plot in the tree_metadata_filename
                                                called 'plotting_name'. Details below.
                --gene_naming_table_filename  Filepath for TSV-format table linking query gene IDs in the BLAST table to proper
                                                gene names. Details below.
                --bootstrap_cutoff            A percentage at or above which to display the bootstrap values on the tree (e.g., 80)
                --root_name                   Exact name of the tip you wish to use as the root of the tree, if your tree is not 
                                                already rooted

              Tree vs. the BLAST table
                  Note that the subject organism names MUST be EXACTLY the same between the tree tips and the BLAST table 
                  'subject_name' column. Otherwise, the script will fail.

              Tree metadata:
                  You can optionally provide a tree_metadata table to overlay additional information onto the phylogenetic tree.
                  
                  The first column of this TSV (tab-separated) table MUST be called 'subject_name' and include the EXACT names
                  of all genomes in the phylogenetic tree. You can then optionally include the following:

                  1. genome_plotting_names: add a column called 'plotting_name' to include the names of the organisms that you 
                      want to appear on the final plot. You can use spaces, most special characters, and so on.

                  2. tree_decorator_colname: you can overlay characteristics of the organisms/genomes as fill colours on the tips
                      of the tree. Add a column to the metadata table with any name you'd like, e.g., 'GC_content', 
                      'predicted_metabolism', 'pathogenicity', and so on. Fill with meaningful data to you. Then, specify the
                      'tree_decorator_colname' to EXACTLY match the name of ONE of the additional columns. That info will then
                      be plotted on the tree. Enjoy! (We might expand this in the future to allow for font colours and so on to 
                      be varied.)

              Gene naming for BLAST table:
                  You can provide a gene_naming_table to provide custom names and ordering for query genes in your BLAST search, 
                  in place of the qseqid from NCBI, which may not be very human-readable.

                  The gene naming table must meet the following criteria:
                    - First column: 'qseqid' - the EXACT ID of ALL of the unique query proteins in the BLAST table must be included.
                          You can include additional qseqid's here if you'd like, they just won't be used.
                    - Second column: 'gene_name' - a corresponding name of your choice (e.g., rpoB, dsrA, and so on)
                  
                  The order of the rows in this table will dictate the order of the genes in the heatmap.

              "))
    
    quit(status = 1)
  }
  
  # Exit if required inputs are not provided
  if ( is.null(opt$tree_filepath) ) {
    stop("Filepath to phylogenetic tree required (-i). Try -h for help message.")
  } else if ( is.null(opt$blast_table_filepath) ) {
    stop("Filepath to BLAST hit table required (-j). Try -h for help message.")
  } else if ( is.null(opt$output_filepath) ) {
    stop("Output filepath required (-o). Try -h for help message.")
  }
  
  # Check on the tree_metadata_filename and tree_decorator_colname were provided, as needed
  if ( is.null(opt$tree_metadata_filename) ) {
    # Annul everything
    # TODO - throw a warning if some of the other flags were set
    opt$tree_metadata_filename <- NA
    opt$tree_decorator_colname <- NA
    opt$genome_plotting_names <- NA
  } else if ( is.null(opt$tree_metadata_filename) == FALSE && is.null(opt$genome_plotting_names) == FALSE
              && is.null(opt$tree_decorator_colname) == FALSE ) {
    # All are present
    opt$genome_plotting_names <- TRUE
  } else if ( is.null(opt$tree_metadata_filename) == FALSE && is.null(opt$genome_plotting_names) == FALSE ) {
    # tree_decorator_colname must not be set
    opt$tree_decorator_colname <- NA
    opt$genome_plotting_names <- TRUE
  } else if ( is.null(opt$tree_metadata_filename) == FALSE && is.null(opt$tree_decorator_colname) == FALSE ) {
    # genome_plotting_names must not be set
    opt$genome_plotting_names <- NA
  }
  
  # Set defaults
  if ( is.null(opt$gene_naming_table_filename) == TRUE ) {
    opt$gene_naming_table_filename <- NA
  }
  if ( is.null(opt$bootstrap_cutoff) == TRUE ) {
    opt$bootstrap_cutoff <- NA
  }
  if ( is.null(opt$root_name) == TRUE ) {
    opt$root_name <- NA
  }
  
  # Make variables from provided input and save as global variables (<<-)
  input_phylogenetic_tree_filepath <<- opt$tree_filepath
  input_blast_table_filename <<- opt$blast_table_filepath
  output_pdf_name <<- opt$output_filepath
  tree_metadata_filename <<- opt$tree_metadata_filename
  tree_decorator_colname <<- opt$tree_decorator_colname
  map_genome_names <<- opt$genome_plotting_names
  gene_naming_table_filename <<- opt$gene_naming_table_filename
  bootstrap_cutoff <<- opt$bootstrap_cutoff
  root_name <<- opt$root_name
  
}


# Function: helper function for timestamps. Returns a length 1 character vector.
ts <- function() {
  
  datetime_utc <- format(as.POSIXlt(Sys.time(), tz = "UTC"), "%a %d %b %Y %H:%M:%S %Z")
  
  date_message <- paste("[ ", datetime_utc, " ]: ", sep = "")
  
  return(date_message)
}


# Function: re-roots tree to terminal node name matching 'root_name'
# Input: 'phylo_tree' - ggtree-format tree; 'root_name' - character vector (length 1) with exact name of the to-be root
# Return: ggtree-format re-rooted tree
reroot_tree <- function(phylo_tree, root_name) {
  
  # Extract the data from the tree in tabular format
  tree_data <- ggtree(phylo_tree)$data
  
  # Check that the root_name exists
  if ( (root_name %in% tree_data$label) == FALSE ) {
    stop(ts(), "Could not find the root_name '", root_name, "' in the provided tree. Cannot re-root. Exiting...")
  }
  
  # Get the ggtree ID # of the parent node of the to-be root
  root_node_num <- as.numeric(dplyr::filter(tree_data, label == root_name)$parent)
  
  # Check that only one matching parent node exists
  if (length(root_node_num) != 1) {
    stop(ts(), "The provided root_name '", root_name, "' appears to have more than one parent node. Cannot re-root. Exiting...")
  }
  
  # Re-root
  tree_rooted <- suppressMessages(reroot(phylo_tree, root_node_num))
  
  return(tree_rooted)
}


# Function: extracts tree label (bootstrap) data and applies bootstrap cutoff if desired
# Input: 'phylo_tree' - ggtree-format tree; 'bootstrap_cutoff' - numeric vector (length 1) with minimum bootstrap cutoff in %
# Return: extracted tree data as a data frame with labels adjusted to only contain bootstraps (numeric) (optionally) above bootstrap cutoff for branches.
        # To be mapped onto the main tree during plotting.
generate_bootstrap_labels <- function(phylo_tree, bootstrap_cutoff) {
  # For adding bootstrap labels like this, see https://guangchuangyu.github.io/software/ggtree/faq/ (accessed Sept 14, 2018)
  # Note that 'labels' encompasses tip labels and node labels. You can to extract the tip labels alone.
  
  # Extract info from the tree into data frame format
  bootstrap_label_data <- ggtree(phylo_tree)$data
  
  # Set the labels to numeric (instead of character; will make non character go 'NA', like perhaps tip labels)
  bootstrap_label_data$label <- as.numeric(bootstrap_label_data$label)
  # TODO - make this more elegant to avoid warning
  
  # Filter out tip labels and 'NA' labels
  bootstrap_label_data <- dplyr::filter(bootstrap_label_data, is.na(label) == FALSE & isTip == FALSE)
  
  # Set bootstrap cutoff if desired
  if ( is.na(bootstrap_cutoff) == FALSE ) {
    message(ts(), "Applying bootstrap display cutoff of ", bootstrap_cutoff)
    bootstrap_label_data <- dplyr::filter(bootstrap_label_data, label >= bootstrap_cutoff)
  }
  
  # Return whole tree data frame (helps for mapping the labels onto the existing tree later)
  return(bootstrap_label_data)
  
}


# Function: chooses which discrete colour scale would work best for a plot based on the number of entries to be plotted
choose_discrete_colour_scale <- function(number_of_entries, data_type) {
  
  # Choose the best colour scale based on the number of entries to be plotted
  if ( number_of_entries == 2 ) {
    colour_palette <- RColorBrewer::brewer.pal(n = 3, name = "Dark2")[c(1,2)]
  } else if ( number_of_entries <= 8 ) {
    colour_palette <- RColorBrewer::brewer.pal(n = number_of_entries, name = "Dark2")
  } else if ( number_of_entries <= 12 ) {
    colour_palette <- RColorBrewer::brewer.pal(n = number_of_entries, name = "Set3")
  } else if ( number_of_entries > 12 ) {
    colour_palette <- hue_pal(h = c(20,290))(number_of_entries)
  } else {
    stop("Something is wrong with the number_of_entries ('", number_of_entries, "'). Is it non-numeric? Exiting...")
  }
  
  return(colour_palette)
  
}


# Function: makes an initial plot of the tree with overlaid bootstrap labels
# Inputs: 'phylo_tree' - ggtree-format tree; 'bootstrap_label_data' - data frame extracted from tree and possibly modified, generated in 'generate_bootstrap_labels'
# Return: ggtree plot
make_tree_plot <- function(phylo_tree, bootstrap_label_data) {
  
  tree_plot <- ggtree(phylo_tree, size = 1.5, colour = "black", ladderize = TRUE) +
    
    # Add the bootstrap labels from the external file
    # TODO - fine-tune the 'nudge' value so that it works more generically
    geom_text(data = bootstrap_label_data, aes(label = label), nudge_x = -0.035, nudge_y = 0.35, size = 3) +
    
    # Add scale bar
    # TODO - set y based on tree topography?
    geom_treescale(x = 0, y = 8, linesize = 1, fontsize = 3, offset = 0.5) +
    
    # Manually tune the y-axis boundaries to match the heatmap: https://stackoverflow.com/a/18718196 (accessed Sept. 15, 2018)
    # TODO - change this to be a function of the number of entries
    scale_y_discrete(expand = c(0,0.6)) +
  
    ## Manually set the righthand cutoff for the tree
    # TODO - change this to be a function of the length of text
    xlim_tree(1.8)
    # xlim etc.: https://guangchuangyu.github.io/2016/10/xlim_tree-set-x-axis-limits-for-only-tree-panel/ (accessed Sept 15, 2018)
  
  return(tree_plot)
  
}


# Function: loads tree metadata file and checks for errors
# Inputs: 'tree_metadata_filename'
          # 'phylo_tree' - loaded (not plotted) ggtree treefile. Just for checking columns match.
# Return: tree metadata (data frame)
load_tree_metadata <- function(tree_metadata_filename, phylo_tree) {
  
  # Load metadata to decorate the tree
  message(ts(), "Loading tree metadata")
  metadata <- read.table(tree_metadata_filename, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  
  # Check first column of metadata matches expected
  if ( colnames(metadata)[1] != "subject_name" ) {
    stop(ts(), "First column of the tree metadata file must be 'subject_name'; yours is '", colnames(metadata)[1], "'. Exiting...")
  }
  if ( identical(sort(unique(metadata$subject_name)), 
                 sort(dplyr::filter(ggtree(phylo_tree)$data, isTip == TRUE)$label)) == FALSE ) {
    stop(ts(), "Entries in 'subject_name' of the provided metadata file do not exactly match the tip labels on the tree. Exiting...")
  }
  
  return(metadata)
}


# Function: adds tip labels to the tree plot, either the standard ones or as specified in metadata
# Inputs: 'tree_plot' - ggtree plot; '
# 'metadata' - OPTIONAL data frame of user-supplied information for 'plotting_name' (required column)
            # If not provided, will just use standard tip labels
# Return: ggtree plot with labels
add_tip_labels_to_tree_plot <- function(tree_plot, metadata_with_plotting_name = NULL) {
  
  if ( is.null(metadata_with_plotting_name) ) {
    message(ts(), "Adding standard tip labels to tree")
    
    tree_plot <- tree_plot +
      # Align tips to far right
      # TODO - make the size a function of total plot size? Or is fixed okay? But need to get it to match the heatmap font on the x-axis
      # TODO - make the offset a function of total plot size
      # TODO - ideally, make italic
      geom_tiplab(align = TRUE, linetype = "dotted", size = 4, offset = 0.1)
    
  } else {
    
    message(ts(), "Adding custom tip labels as specified by 'plotting_name' columns of tree_metadata")
    
    # Check table is okay
    if ( ("plotting_name" %in% colnames(metadata_with_plotting_name)) == FALSE ) {
      stop("You set the '--genome_plotting_names' flag but did not provide a 'plotting_name' column in your tree metadata table. Cannot continue. Exiting...")
    }
    
    # Reduce the metadata down to just 'subject_name' and 'plotting_name'
    required_cols <- c("subject_name", "plotting_name")
    name_metadata <- dplyr::select(metadata_with_plotting_name, required_cols)
    
    # Add the labels
    # N.B. Adding the 'name_metadata' on the first line binds the metadata table onto the tree's data table for aesthetics later
    tree_plot <- tree_plot  %<+% name_metadata +
      # Align tips to far right
      # TODO - make the size a function of total plot size? Or is fixed okay? But need to get it to match the heatmap font on the x-axis
      # TODO - make the offset a function of total plot size
      # TODO - ideally, make italic
      geom_tiplab(aes(label = plotting_name), align = TRUE, linetype = "dotted", size = 4, offset = 0.1)
      
  }
  
  return(tree_plot)
  
}


# Function: decorates an existing ggtree plot with custom fill/colour based on user-supplied metadata
# Inputs: 'tree_plot' - ggtree plot; '
          # 'metadata' - data frame of user-supplied information with row names matching the tree's tip labels;
          # 'tree_decorator_colname' - character vector (length 1) of the name of the column in 'metadata' that you want to overlay on the plot
# Return: ggtree plot
decorate_tree_plot <- function(tree_plot, metadata, tree_decorator_colname) {
  # For decorating with metadata, see https://bioconductor.org/packages/devel/bioc/vignettes/ggtree/inst/doc/treeAnnotation.html#tree-annotation-with-user-specified-annotation (accessed Sept. 15, 2018)
  
  # TODO - add more options for the user (e.g., fill size, or maybe even colouring text instead?)
  
  # Confirm the tree_decorator_colname is real
  if ( (tree_decorator_colname %in% colnames(metadata)) == FALSE) {
    stop("Cannot find column named '", tree_decorator_colname, "' in the supplied metadata table for mapping onto the tree. Exiting...")
  }
  
  # Add the decorator onto the tree
  # N.B. Adding the 'metadata' on the first line binds the metadata table onto the tree's data table for aesthetics later
  tree_plot <- tree_plot  %<+% metadata +
    # TODO - check: does aes_string work?
    geom_tippoint(aes_string(fill = tree_decorator_colname), shape = 21, alpha = 0.8, size = 5) +
    # geom_tiplab(align = TRUE, linetype = NULL, size = 2) + # aes(colour = metabolism) , linetype = "dotted", size = 2
    theme(legend.position = "left")
  
  # Extract a vector of the decorator column to decide on appropriate colour scale
  decorators <- metadata[,match(tree_decorator_colname, colnames(metadata))]
  
  # Add colour scale
  if ( is.numeric(decorators) == TRUE ) {
    # Continuous data
    colour_palette <- rev(RColorBrewer::brewer.pal(name = "YlGnBu", n = 5))
    
    tree_plot <- tree_plot + scale_fill_gradientn(colours = colour_palette, na.value = "grey50")
    
  } else if ( is.character(decorators) == TRUE ) {
    # Discrete data
    # Determine the fill colours for the tree based on the number of unique entries
    num_factors <- length(unique(decorators))
    colour_palette <- choose_discrete_colour_scale(num_factors)
    
    tree_plot <- tree_plot + scale_fill_manual(values = colour_palette)
    
  } else {
    # Discrete data?
    warning("The decorator_column in the metadata table is not numeric or character; fishy! Watch for possible errors...")
    # Determine the fill colours for the tree based on the number of unique entries
    num_factors <- length(unique(decorators))
    colour_palette <- choose_discrete_colour_scale(num_factors)
    
    tree_plot <- tree_plot + scale_fill_manual(values = colour_palette)
    
  }
  
  return(tree_plot)
  
}


# Function: master function to load and plot the phylogenetic tree
# Inputs: described above
# Return: list of two: 'phylo_tree' - ggtree unplotted object ; 'phylo_tree_fig' - ggtree figure
load_and_plot_phylogenetic_tree <- function(input_phylogenetic_tree_filepath, root_name, bootstrap_cutoff, 
                                   tree_metadata_filename, tree_decorator_colname) {
  # Read tree
  message(ts(), "Reading input phylogenetic tree")
  phylo_tree <- read.tree(input_phylogenetic_tree_filepath)
  
  # Optionally re-root tree
  if (is.na(root_name) == FALSE) {
    message(ts(), "Re-rooting tree to '", root_name, "'")
    phylo_tree <- reroot_tree(phylo_tree, root_name)
  }
  
  # Set cutoff for bootstraps externally, to be overlaid onto the tree figure later
  # No cutoff is applied if bootstrap_cutoff is set to 'NA'
  message(ts(), "Generating bootstrap labels")
  bootstrap_label_data <- generate_bootstrap_labels(phylo_tree, bootstrap_cutoff)
  
  
  # Generate tree plot
  message(ts(), "Generating ggtree plot")
  phylo_tree_fig <- make_tree_plot(phylo_tree, bootstrap_label_data)
  
  # Add tip labels (either defaults or based on plotting_name)
  if ( is.na(tree_metadata_filename) == FALSE && is.na(map_genome_names) == FALSE ) {
    # Load metadata
    metadata <- load_tree_metadata(tree_metadata_filename, phylo_tree)
    
    phylo_tree_fig <- add_tip_labels_to_tree_plot(phylo_tree_fig, metadata)
  } else {
    phylo_tree_fig <- add_tip_labels_to_tree_plot(phylo_tree_fig)
  }
  
  # Optionally decorate with metadata
  if ( is.na(tree_metadata_filename) == FALSE && is.na(tree_decorator_colname) == FALSE ) {
    # Load metadata
    metadata <- load_tree_metadata(tree_metadata_filename, phylo_tree)
    
    # Decorate
    message(ts(), "Decorating ggtree plot with metadata column '", tree_decorator_colname, "'")
    phylo_tree_fig <- decorate_tree_plot(phylo_tree_fig, metadata, tree_decorator_colname)
  }
  
  # Make list to return to user
  tree_list <- list(phylo_tree, phylo_tree_fig)
  names(tree_list) <- c("phylo_tree", "phylo_tree_fig")
  
  return(tree_list)
  
}

# Function: loads BLAST table and checks for expected column names (HARD-CODED in function)
# Inputs: 'input_blast_table_filename'
# Return: blast_table (data frame)
read_blast_table <- function(input_blast_table_filename) {
  
  # Load the data
  blast_table <- read.table(input_blast_table_filename, sep = ",", header = TRUE, stringsAsFactors = FALSE)
  
  # HARD-CODED expected header names
  expected_header_names <- c("subject_name", "qseqid", "sseqid", "pident", "evalue", "qcovhsp", "bitscore")
  
  # Confirm the columns look okay
  if ( identical(colnames(blast_table), expected_header_names) == FALSE ) {
    
    stop("BLAST data table did not have expected column names (", glue::glue_collapse(expected_header_names, sep = ", "),
         "). Instead, had: ", glue::glue_collapse(colnames(blast_table, sep = ", ")), ". Exiting...")
    
  }
  
  return(blast_table)
  
}


# Function: changes the order of the subject_name in the BLAST table to match that of the ggtree tips. Checks for perfect match.
# Inputs: 'blast_table' data frame; 'tip_order' - character vector of the ggtree tips in order
# Return: blast_table (data frame) with subject_name as an ordered factor
order_blast_table_subjects <- function(blast_table, tip_order) {
  # Check that the tip_order labels match the blast table's subject_name labels
  if ( identical(sort(tip_order), sort(unique(blast_table$subject_name))) == FALSE) {
    stop(ts(), "Tree tip labels do not perfectly match the subject_name entries in the BLAST table. Exiting...")
  }
  
  # Make subjects the same order as in the tree
  blast_table$subject_name <- factor(blast_table$subject_name, levels = rev(tip_order), ordered = TRUE)
  
  return(blast_table)
}


# Function: reads the optional gene_naming_table and then overlays the names and order onto the BLAST table
# Inputs: 'blast_table' data frame; 'gene_naming_table_filename' - character vector (length 1) specifying the filename
# Return: blast_table (data frame) with qseqid as an ordered factor with human-readable names
overlay_gene_naming <- function(blast_table, gene_naming_table_filename) {
  
  # Load gene naming table
  gene_naming_table <- read.table(gene_naming_table_filename, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  
  # Check that the column names match expected (for the first two columns; doesn't matter after that)
  # HARD-CODED
  gene_table_expected_headers <- c("qseqid", "gene_name")
  
  if ( identical(colnames(gene_naming_table)[1:2], gene_table_expected_headers) == FALSE ) {
    stop(ts(), "The first two columns of the gene table should be: ", glue::glue_collapse(gene_table_expected_headers, sep = ", "),
         ". However, you provided something else: ",  glue::glue_collapse(colnames(gene_naming_table)[1:2], sep = ", "), ". Exiting...")
  }
  
  # Check that the query seq IDs in the gene table include those provided as queries for BLAST (can contain extra - that's fine)
  gene_table_queries <- sort(gene_naming_table$qseqid)
  blast_table_queries <- sort(unique(blast_table$qseqid))
  if ( length(unique(blast_table_queries %in% gene_table_queries)) > 1 ) {
    stop(ts(), "Query sequence IDs from the BLAST table are not contained in the gene naming table provided. Exiting...")
  } else if ( unique(blast_table_queries %in% gene_table_queries) == FALSE ) {
    stop(ts(), "None of the query sequence IDs from the BLAST table match the gene naming table provided. Exiting...")
  }
  
  # Change qseqid to gene_name and order according to the gene_naming_table
  blast_table$qseqid <- plyr::mapvalues(x = blast_table$qseqid, from = gene_naming_table$qseqid, to = gene_naming_table$gene_name)
  blast_table$qseqid <- factor(blast_table$qseqid, levels = gene_naming_table$gene_name, ordered = TRUE)
  
  return(blast_table)
  
}


# Function: plots the BLAST table as a heatmap in ggplot
# Inputs: 'blast_table' data frame
# Return: ggplot heatmap
plot_blast_heatmap <- function(blast_table) {
  
  blast_heatmap <- ggplot(blast_table, aes(y = subject_name, x = qseqid)) +
    geom_tile(aes(fill = pident)) +
    theme_bw() +
    theme(panel.grid = element_blank(), axis.title = element_text(size = 12),
          strip.text = element_text(size = 7), strip.background = element_rect(fill = "#e6e6e6"),
          panel.border = element_rect(colour = "black", size = 1),
          axis.text = element_text(size = 10, colour = "black"), 
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          axis.text.y = element_blank(),
          axis.ticks = element_line(size = 0.5), axis.line = element_line(colour = "black", size = 0.5),
          legend.text = element_text(size = 10, colour = "black"), legend.title = element_blank(),
          legend.key = element_blank(), legend.key.size = unit(4.5, "mm")) +
    xlab("Gene") +
    ylab(NULL)
  
  return(blast_heatmap)
}


# Function: master function to load and plot the BLAST table and gene_naming_table to produce a heatmap
# Inputs: see above
# Return: list of two: 'blast_table' data frame; 'blast_heatmap' ggplot object
load_and_plot_blast_table <- function(input_blast_table_filename, tip_order, gene_naming_table_filename) {
  
  # Load the blast table
  message(ts(), "Loading the BLAST table")
  blast_table <- read_blast_table(input_blast_table_filename)
  
  # Order the BLAST table subject names to match the ggtree
  message(ts(), "Aligning BLAST table's subject names to match order of the ggtree")
  blast_table <- order_blast_table_subjects(blast_table, tip_order)
  
  # Overlay gene names and gene naming order, if provided
  if ( is.na(gene_naming_table_filename) == FALSE ) {
    message(ts(), "Overlaying gene naming and ordering onto the BLAST table")
    blast_table <- overlay_gene_naming(blast_table, gene_naming_table_filename)
  }
  
  # Create the heatmap
  message(ts(), "Plotting BLAST heatmap")
  blast_heatmap <- plot_blast_heatmap(blast_table)
  
  # Make output list
  blast_table_list <- list(blast_table, blast_heatmap)
  names(blast_table_list) <- c("blast_table", "blast_heatmap")
  
  return(blast_table_list)
  
}


main <- function() {
  # Run command line version if requested
  if (RUN_COMMAND_LINE == TRUE) {
    parse_command_line_input()
  }
  
  # Startup messages
  message(ts(), "Running generate_BackBLAST_heatmap.R")
  message(ts(), "Input phylogenetic tree filepath: ", input_phylogenetic_tree_filepath)
  message(ts(), "Input BLAST table filename: ", input_blast_table_filename)
  message(ts(), "Output PDF filename: ", output_pdf_name)
  message(ts(), "Bootstrap display cutoff (%; ignored if 'NA'): ", bootstrap_cutoff)
  message(ts(), "Root name (ignored if 'NA'): ", root_name)
  message(ts(), "Tree metadata filename (ignored if 'NA'): ", tree_metadata_filename)
  message(ts(), "Decorator column name in tree metadata (ignored if 'NA'): ", tree_decorator_colname)
  message(ts(), "Input gene naming table filename (ignored if 'NA'): ", gene_naming_table_filename)
  
  # Load and plot the tree
  # TODO - after reading metadata table, CHECK that the first column names correspond to the tree tip labels
  phylo_tree_list <- load_and_plot_phylogenetic_tree(input_phylogenetic_tree_filepath, root_name, bootstrap_cutoff, 
                                     tree_metadata_filename, tree_decorator_colname)
  
  # Get tip order of the tree, to match with heatmap later
  # Based on https://groups.google.com/forum/#!topic/bioc-ggtree/LqRDK78m3U4 (accessed Sept. 15, 2018)
  message(ts(), "Exporting tip order of tree to correspond with heatmap")
  tip_order <- dplyr::filter(ggtree(phylo_tree_list[[1]])$data, isTip == TRUE)
  tip_order <- tip_order[order(tip_order$y, decreasing = TRUE),]$label # plotting_name

  
  # Load and plot the BLAST table as a heatmap
  blast_table_list <- load_and_plot_blast_table(input_blast_table_filename, tip_order, gene_naming_table_filename)
  
  # Combine the tree and heatmap
  # Got ggarrange ideas from https://cran.r-project.org/web/packages/egg/vignettes/Ecosystem.html (accessed Sept. 15, 2018)
  # TODO - make the dimensions a function of the gene number and the subject number
  message(ts(), "Combining the ggtree and the heatmap")
  combined_plot <- ggarrange(phylo_tree_list[[2]], blast_table_list[[2]], nrow = 1, widths = c(5, 8), heights = c(5), padding = unit(0, "mm"))
  # N.B., set debug = TRUE to see grid lines
  
  # Print a PDF of the combined plot
  # TODO - set dimensions as a function of the gene number and the subject number
  message(ts(), "Saving to PDF")
  pdf(file = output_pdf_name, width = 17, height = 7)
  print(combined_plot)
  dev.off()
  
  message(ts(), "generate_BackBLAST_heatmap.R: done.")
  
}

main()

