#!/usr/bin/env Rscript
# generate_heatmap.R
# Copyright Lee H. Bergstrand and Jackson M. Tsuji, 2019
# Plots a newick treefile and BLAST table together as a phylogenetic tree and heatmap
# Part of the BackBLAST pipeline

# Load libraries
# Note: warn.conflicts: Normally, when a library is loaded that has a function with identical 
  # name to another function (e.g., setdiff() in dplyr), a warning is given during package 
  # load. The warnings are disabled here to prevent excessive messages when the script is run.
  # Long-term, once possible to switch to R 3.6.0, options(conflicts.policy(list(warn = FALSE)))
  # can be used.
library(argparser)
library(futile.logger)
library(tools)
library(glue, warn.conflicts = FALSE)
library(plyr, warn.conflicts = FALSE)
library(dplyr, warn.conflicts = FALSE)
library(tibble, warn.conflicts = FALSE)
library(reshape2, warn.conflicts = FALSE)
library(RColorBrewer, warn.conflicts = FALSE)
library(ggplot2, warn.conflicts = FALSE)
library(ape, warn.conflicts = FALSE)
library(maps, warn.conflicts = FALSE)
library(phytools, warn.conflicts = FALSE)
library(tidytree, warn.conflicts = FALSE)
library(treeio, warn.conflicts = FALSE)
suppressPackageStartupMessages(library(ggtree))
library(gridExtra, warn.conflicts = FALSE)
library(egg, warn.conflicts = FALSE)

#' Reads a data table as a tibble with several default parameters. All parameters below are the same as read.table()
#' 
#' @param file Filepath to data table
#' @param sep Field separator character
#' @param header Logical; does the table have headers as the first row?
#' @param stringsAsFactors Logical; should character vectors be converted to factors?
#' @return Tibble of the data table
#' @export
read_tibble <- function(file, sep = "\t", header = TRUE, stringsAsFactors = FALSE) {
  data_table <- read.table(file, sep = sep, header = header, stringsAsFactors = stringsAsFactors) %>%
    tibble::as_tibble()
  return(data_table)
}

#' Convert character "NA" (from command line) into constant NA
#' 
#' @param entry single-length vector of any type
#' @return single-length vector; if it was "NA", it will now be NA; otherwise it will be the same as input
#' @export
convert_to_constant_NA <- function(entry) {
  conversion <- entry
  if (!is.na(entry)) {
    if (entry == "NA") {
      futile.logger::flog.debug("Converting 'NA' input to constant NA value")
      conversion <- NA
    }
  }
  return(conversion)
}

#' Re-root a ggtree
#' 
#' @param phylo_tree ggtree-format tree
#' @param root_name Character (length 1 vector); the exact name of the to-be root
#' @return ggtree-format tree, rooted
#' @export
reroot_ggtree <- function(phylo_tree, root_name) {

  # Midpoint root the tree if requested; otherwise root to a name
  if (root_name == "midpoint") {
    futile.logger::flog.info("Midpoint rooting the tree")
    tree_rooted <- phytools::midpoint.root(phylo_tree)
  } else {
    futile.logger::flog.info(glue::glue("Re-rooting tree to '", root_name, "'"))
    # Extract the data from the tree in tabular format
    tree_data <- ggtree::ggtree(phylo_tree)$data

    # Check that the root_name exists
    if ( !(root_name %in% tree_data$label) ) {
      futile.logger::flog.error(glue::glue("Could not find the root_name '", root_name,
                            "' in the provided tree. Cannot re-root. Exiting..."))
      quit(save = "no", status = 1)
    }

    # Re-root
    tree_rooted <- treeio::root(phylo_tree, outgroup = root_name)
  }

  return(tree_rooted)
}

#' Extracts tree label (bootstrap) data and applies bootstrap cutoff if desired
#' 
#' @param phylo_tree ggtree-format tree
#' @param bootstrap_cutoff Numeric (length 1 vector); the minimum bootstrap cutoff, in percent
#' @return extracted tree data as a data frame with labels adjusted to only contain bootstraps (numeric) (optionally) above bootstrap cutoff for branches.
#' To be mapped onto the main tree during plotting.
#' @export
generate_bootstrap_labels <- function(phylo_tree, bootstrap_cutoff) {
  # For adding bootstrap labels like this, see https://guangchuangyu.github.io/software/ggtree/faq/ (accessed Sept 14, 2018)
  # Note that 'labels' encompasses tip labels and node labels. You can to extract the tip labels alone.
  
  # Extract info from the tree into data frame format
  bootstrap_label_data <- ggtree::ggtree(phylo_tree)$data
  
  # Set the labels to numeric (instead of character; will make non character go 'NA', like perhaps tip labels)
  # This longer lapply function tries to avoid warnings from coersion while still allowing for transparent error reporting otherwise
  bootstrap_label_data$label <- unlist(lapply(bootstrap_label_data$label, function(label) {
    if( !is.na(suppressWarnings(as.numeric(label))) ) {
      return(as.numeric(label))
    } else {
      return(NA)}
  }))
  
  # Filter out tip labels and 'NA' labels
  bootstrap_label_data <- dplyr::filter(bootstrap_label_data, !is.na(label) & isTip == FALSE)
  
  # Set bootstrap cutoff if desired
  if ( !is.na(bootstrap_cutoff) ) {
    futile.logger::flog.info(glue::glue("Applying bootstrap display cutoff of ", bootstrap_cutoff))
    bootstrap_label_data <- dplyr::filter(bootstrap_label_data, label >= bootstrap_cutoff)
  }
  
  # Return whole tree data frame (helps for mapping the labels onto the existing tree later)
  return(bootstrap_label_data)
}

#' Makes an initial plot of the tree with overlaid bootstrap labels
#' 
#' @param phylo_tree ggtree-format tree
#' @param bootstrap_label_data data frame extracted from tree and possibly modified, generated in 'generate_bootstrap_labels'
#' @return ggtree plot
#' @export
plot_ggtree <- function(phylo_tree, bootstrap_label_data) {
  
  tree_plot <- ggtree::ggtree(phylo_tree, size = 1, colour = "black", ladderize = TRUE,
                      branch.length = 0.1) +
    
    # Add the bootstrap labels from the external file
    geom_text(data = bootstrap_label_data, aes(label = label), nudge_x = -0.035, nudge_y = 0.35, size = 3) +
    
    # Add the dotted lines from the tree tips
    geom_tiplab(align = TRUE, linetype = "dotted", size = 0.1, offset = 0.1) +
    
    # Add scale bar
    geom_treescale(x = 0, y = 1.5, linesize = 1, fontsize = 3, offset = 0.1) +
    
    # Manually tune the y-axis boundaries to match the heatmap: https://stackoverflow.com/a/18718196 (accessed Sept. 15, 2018)
    # TODO - this might need to be a function of the number of entries
    scale_y_discrete(expand = c(0,0.6)) +
  
    ## Manually set the righthand cutoff for the tree
    xlim_tree(1)
    # xlim etc.: https://guangchuangyu.github.io/2016/10/xlim_tree-set-x-axis-limits-for-only-tree-panel/ (accessed Sept 15, 2018)
  
  return(tree_plot)
}

#' Master function to load and plot the phylogenetic tree
#' 
#' @param input_phylogenetic_tree_filepath Character (length 1 vector); the filepath to the phylogenetic tree
#' @param root_name Character (length 1 vector); the exact name of the to-be root
#' @param bootstrap_cutoff Numeric (length 1 vector); the minimum bootstrap cutoff, in percent
#' @return list of two: 'phylo_tree' - ggtree unplotted object ; 'phylo_tree_fig' - ggtree figure
#' @export
load_and_plot_phylogenetic_tree <- function(input_phylogenetic_tree_filepath, root_name, bootstrap_cutoff) {
  # Read tree
  futile.logger::flog.info("Reading input phylogenetic tree")
  phylo_tree <- ape::read.tree(input_phylogenetic_tree_filepath)
  
  # Optionally re-root tree
  if ( !is.na(root_name) ) {
    phylo_tree <- reroot_ggtree(phylo_tree, root_name)
  }
  
  # Set cutoff for bootstraps externally, to be overlaid onto the tree figure later
  # No cutoff is applied if bootstrap_cutoff is set to 'NA'
  futile.logger::flog.debug("Generating bootstrap labels")
  bootstrap_label_data <- generate_bootstrap_labels(phylo_tree, bootstrap_cutoff)
  
  # Generate tree plot
  futile.logger::flog.info("Generating ggtree plot")
  phylo_tree_fig <- plot_ggtree(phylo_tree, bootstrap_label_data)
  
  # Make list to return to user
  tree_list <- list(phylo_tree, phylo_tree_fig)
  names(tree_list) <- c("phylo_tree", "phylo_tree_fig")
  
  return(tree_list)
}

#' Loads BLAST table and checks for expected column names (HARD-CODED in function)
#' 
#' @param input_blast_table_filepath Character (length 1 vector); the filepath to the BLAST table (comma-separated)
#' @return tibble of BLAST data
#' @export
read_blast_results <- function(input_blast_table_filepath) {
  
  # Load the data
  blast_results <- read_tibble(input_blast_table_filepath, sep = ",")
  
  # HARD-CODED expected header names
  expected_header_names <- c("subject_name", "qseqid", "sseqid", "pident", "evalue", "qcovhsp", "bitscore")
  
  # Confirm the columns look okay
  if ( !identical(colnames(blast_results), expected_header_names) ) {
    futile.logger::flog.error(glue::glue("BLAST data table did not have expected column names ('", 
                          glue::glue_collapse(expected_header_names, sep = "; "), 
                          "'). Instead, had: '", glue::glue_collapse(colnames(blast_results), sep = "; "), 
                          "'. Exiting..."))
    quit(save = "no", status = 1)
  }
  
  return(blast_results)
}

#' Changes the order of the subject_name in the BLAST table to match that of the ggtree tips
#' 
#' @param blast_results Tibble output of read_blast_results
#' @param tip_order Character vector of the exact order of the tip labels in the phylogenetic tree
#' @return tibble of BLAST data with subject_name as an ordered factor and reduced to match phylogenetic tree names (if needed)
#' @export
order_blast_subjects <- function(blast_results, tip_order) {
  
  # If entries are not identical, try filtering down to just what is in the phylogenetic tree
  if ( !identical(sort(unique(blast_results$subject_name)), sort(tip_order)) ) {
    
    # Make sure that the tree tips are contained in the blast table; if so, everything is okay
    # The length should be zero if all tree tips are contained in the blast table
    if(length(setdiff(tip_order, unique(blast_results$subject_name))) > 0) {
      
      # But if the length is > 0, it means some entries in the tree are MISSING in the blast table, a major issue
      missing_blast_table_entries <- setdiff(tip_order, unique(blast_results$subject_name))
      flog.error(glue::glue("The provided BLAST table is missing some entries in the phylogenetic tree: '",
                            glue::glue_collapse(missing_blast_table_entries, sep = ", "),
                            "'. Cannot continue -- exiting..."))
      quit(save = "no", status = 1)
    }
    
    # Report to the user which entries in the BLAST table are to be ignored
    extra_blast_table_entries <- setdiff(unique(blast_results$subject_name), tip_order)
    
    futile.logger::flog.info(glue::glue("Some entries in the BLAST table are missing in the phylogenetic tree ",
                                        "and will be removed in plotting: '", 
                                        glue::glue_collapse(extra_blast_table_entries, sep = ", "),  "'."))
    
    # Filter down the BLAST table to have the same subject_name's as the tree
    blast_results <- dplyr::filter(blast_results, subject_name %in% tip_order)
  }
  
  # Make subjects the same order as in the tree
  blast_results$subject_name <- factor(blast_results$subject_name, levels = rev(tip_order), ordered = TRUE)
  
  return(blast_results)
}

#' Overlays user-given genome names for the heatmap y-axis
#' 
#' @param blast_results Tibble output of read_blast_results
#' @param genome_metadata_filepath Character (length 1 vector); the filepath of the tab-separated genome metadata file.
#' Must at least have the columns 'subject_name' and 'plotting_name' as the first and second columns, respectively.
#' @return tibble of BLAST data with subject_name renamed with the user-desired values
#' @export
overlay_genome_naming <- function(blast_results, genome_metadata_filepath) {
  
  # Load the metadata table
  futile.logger::flog.info("Loading genome metadata")
  genome_metadata_table <- read_tibble(genome_metadata_filepath)
  
  # Check that the column names match expected (for the first two columns; doesn't matter after that)
  # HARD-CODED
  genome_metadata_table_expected_headers <- c("subject_name", "plotting_name")
  
  if ( !identical(colnames(genome_metadata_table)[1:2], genome_metadata_table_expected_headers) ) {
    futile.logger::flog.error(glue::glue("The first two columns of the genome metadata table should be: ", 
                                         glue::glue_collapse(genome_metadata_table_expected_headers, sep = ", "), 
                                         ". However, you provided something else: ", 
                                         glue::glue_collapse(colnames(genome_metadata_table)[1:2], sep = ", "), ". Exiting..."))
    quit(save = "no", status = 1)
  }
  
  # Reduce the metadata down to just the expected headers
  genome_metadata_table <- dplyr::select(genome_metadata_table, genome_metadata_table_expected_headers)
  
  # If entries are not identical, try filtering down to just what is in the BLAST table
  if ( !identical(sort(unique(blast_results$subject_name)), sort(genome_metadata_table$subject_name)) ) {
    
    # Make sure that the subject_name's are contained in the metadata; if so, everything is okay
    # The length should be zero if all tree tips are contained in the metadata table
    if(length(setdiff(unique(blast_results$subject_name), genome_metadata_table$subject_name)) > 0) {
      # But if the length is > 0, it means some entries are MISSING in the metadata, a major issue
      missing_metadata_entries <- setdiff(unique(blast_results$subject_name), genome_metadata_table$subject_name)
      flog.error(glue::glue("The provided gene metadata file is missing some entries in the BLAST table: '",
                            glue::glue_collapse(missing_metadata_entries, sep = ", "),
                            "'. Cannot continue -- exiting..."))
      quit(save = "no", status = 1)
    }
    
    # Report to the user which entries in the metadata are extra and will be ignored
    extra_metadata_entries <- setdiff(genome_metadata_table$subject_name, unique(blast_results$subject_name))
    
    futile.logger::flog.info(glue::glue("Some entries in the genome metadata are missing in the BLAST table ",
                                        "and will be removed in plotting: '", 
                                        glue::glue_collapse(extra_metadata_entries, sep = ", "),  "'."))
    
    # Filter down the metadata to match the BLAST table entries
    genome_metadata_table <- dplyr::filter(genome_metadata_table, subject_name %in% unique(blast_results$subject_name))
  }
  
  # Change the subject_name to be the plotting_name in the BLAST table
  blast_results$subject_name <- plyr::mapvalues(blast_results$subject_name, from = genome_metadata_table$subject_name, 
                                               to = genome_metadata_table$plotting_name, warn_missing = TRUE)
  
  return(blast_results)
}

#' Overlays user-given gene names for the heatmap x-axis
#' 
#' @param blast_results Tibble output of read_blast_results
#' @param gene_metadata_filepath Character (length 1 vector); the filepath of the tab-separated gene metadata file.
#' Must at least have the columns 'qseqid' and 'gene_name' as the first and second columns, respectively.
#' @return tibble of BLAST data with qseqid renamed with the user-desired values
#' @export
overlay_gene_naming <- function(blast_results, gene_metadata_filepath) {
  
  # Load gene naming table
  gene_metadata_table <- read_tibble(gene_metadata_filepath)
  
  # Check that the column names match expected (for the first two columns; doesn't matter after that)
  # HARD-CODED
  gene_metadata_table_expected_headers <- c("qseqid", "gene_name")
  
  if ( !identical(colnames(gene_metadata_table)[1:2], gene_metadata_table_expected_headers) ) {
    futile.logger::flog.error(glue::glue("The first two columns of the gene table should be: ", 
                          glue::glue_collapse(gene_metadata_table_expected_headers, sep = ", "), 
                          ". However, you provided something else: ", 
                          glue::glue_collapse(colnames(gene_metadata_table)[1:2], sep = ", "), ". Exiting..."))
    quit(save = "no", status = 1)
  }
  
  # If entries are not identical, try filtering down to just what is in the BLAST table
  if ( !identical(sort(unique(blast_results$qseqid)), sort(gene_metadata_table$qseqid)) ) {
    # Remove any entries from the BLAST table that are missing in the metadata
    if(length(setdiff(unique(blast_results$qseqid), gene_metadata_table$qseqid)) > 0) {
      missing_metadata_entries <- setdiff(unique(blast_results$qseqid), gene_metadata_table$qseqid)
      flog.warn(glue::glue("The provided gene metadata file is missing some entries in the BLAST table: '",
                            glue::glue_collapse(missing_metadata_entries, sep = ", "),
                            "'. These entries will be REMOVED from the BLAST table when plotting."))
      blast_results <- dplyr::filter(blast_results, !(qseqid %in% (missing_metadata_entries)))
    }
    
    # Remove any entries in the metadata that are missing in the BLAST table
    if (length(setdiff(gene_metadata_table$qseqid, unique(blast_results$qseqid))) > 0) {
      extra_metadata_entries <- dplyr::filter(gene_metadata_table,
                                              qseqid %in% setdiff(gene_metadata_table$qseqid, 
                                                                  unique(blast_results$qseqid)))
      extra_metadata_entries$user_reporting <- paste(extra_metadata_entries$qseqid, " (",
                                                     extra_metadata_entries$gene_name, ")",
                                                     sep = "")
      futile.logger::flog.info(glue::glue("Some entries in the gene metadata are missing in the BLAST table ",
                                          "and will be removed in plotting: '", 
                                          glue::glue_collapse(extra_metadata_entries$user_reporting, 
                                                              sep = ", "),  "'."))
      gene_metadata_table <- dplyr::filter(gene_metadata_table, qseqid %in% unique(blast_results$qseqid))
    }
  }
  
  # Change qseqid to gene_name and order according to the gene_metadata_table
  blast_results$qseqid <- plyr::mapvalues(x = blast_results$qseqid, from = gene_metadata_table$qseqid, 
                                        to = gene_metadata_table$gene_name)
  blast_results$qseqid <- factor(blast_results$qseqid, levels = gene_metadata_table$gene_name, ordered = TRUE)
  
  return(blast_results)
}

#' Plots the BLAST table as a heatmap in ggplot
#' 
#' @param blast_results Tibble output of read_blast_results
#' @return ggplot heatmap
#' @export
plot_blast_heatmap <- function(blast_results) {
  
  # Add NA values for missing grid values so that grid lines will appear in the final plot
  # TODO - A bit hacky
  blast_results <- reshape2::dcast(blast_results, subject_name ~ qseqid, value.var = "pident") %>%
    reshape2::melt(na.rm = FALSE, id.vars = c("subject_name"), variable.name = "qseqid",
                   value.name = "pident") %>%
    tibble::as_tibble()
  
  blast_heatmap <- ggplot2::ggplot(blast_results, aes(y = subject_name, x = qseqid)) +
    geom_tile(aes(fill = pident), colour = "black") +
    theme_bw() +
    theme(panel.grid = element_blank(), axis.title = element_text(size = 12),
          panel.border = element_rect(colour = "black", size = 1),
          axis.text = element_text(size = 10, colour = "black"),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          axis.ticks = element_line(size = 0.5), axis.line = element_line(colour = "black", size = 0.5),
          legend.text = element_text(size = 10, colour = "black"),
          legend.title = element_text(size = 10, face = "bold"),
          legend.key = element_rect(colour = "transparent")) +
    guides(fill = guide_legend(title = "Amino acid \nidentity (%)")) +
    scale_fill_gradientn(colours = RColorBrewer::brewer.pal(name = "Blues", n = 5), 
                         na.value = "transparent") +
    xlab(NULL) +
    ylab(NULL)
  
  return(blast_heatmap)
}

#' Master function to load and plot the BLAST table metadata to produce a heatmap
#' 
#' @param input_blast_table_filepath Character (length 1 vector); the filepath to the BLAST table (comma-separated)
#' @param tip_order Character vector of the exact order of the tip labels in the phylogenetic tree; provide 'NA' to skip alignment
#' @param genome_metadata_filepath Character (length 1 vector); the filepath of the tab-separated genome metadata file. 
#' Must at least have the columns 'subject_name' and 'plotting_name' as the first and second columns, respectively.
#' @param gene_metadata_filepath Character (length 1 vector); the filepath of the tab-separated gene metadata file. 
#' Must at least have the columns 'qseqid' and 'gene_name' as the first and second columns, respectively.
#' @return list of two: 'blast_results' data frame; 'blast_heatmap' ggplot object
#' @export
load_and_plot_blast_results <- function(input_blast_table_filepath, tip_order = NA, gene_metadata_filepath,
                                       genome_metadata_filepath) {
  # Load the blast table
  futile.logger::flog.info("Loading the BLAST table")
  blast_results <- read_blast_results(input_blast_table_filepath)
  
  # Order the BLAST table subject names to match the ggtree
  if (!is.na(tip_order[1])) {
    futile.logger::flog.debug("Aligning BLAST table's subject names to match order of the ggtree")
    blast_results <- order_blast_subjects(blast_results, tip_order)
  } else {
    futile.logger::flog.debug("Skipping alignment of BLAST table to phylogenetic tree")
  }
  
  # Overlay genome names and genome naming order, if provided
  if ( !is.na(genome_metadata_filepath) ) {
    futile.logger::flog.info("Overlaying genome naming and ordering onto the BLAST table")
    blast_results <- overlay_genome_naming(blast_results, genome_metadata_filepath)
  }
  
  # Overlay gene names and gene naming order, if provided
  if ( !is.na(gene_metadata_filepath) ) {
    futile.logger::flog.info("Overlaying gene naming and ordering onto the BLAST table")
    blast_results <- overlay_gene_naming(blast_results, gene_metadata_filepath)
  }
  
  # Create the heatmap
  futile.logger::flog.info("Plotting BLAST heatmap")
  blast_heatmap <- plot_blast_heatmap(blast_results)
  
  # Make output list
  blast_results_list <- list(blast_results, blast_heatmap)
  names(blast_results_list) <- c("blast_results", "blast_heatmap")
  
  return(blast_results_list)
}

main <- function(params) {
  # Startup messages
  futile.logger::flog.info("Running generate_heatmap.R")
  futile.logger::flog.info("######### Settings #########")
  futile.logger::flog.info(glue::glue("Input phylogenetic tree filepath (ignored if 'NA'): ", params$input_phylogenetic_tree_filepath))
  futile.logger::flog.info(glue::glue("Input BLAST table filepath: ", params$input_blast_table_filepath))
  futile.logger::flog.info(glue::glue("Output PDF filepath: ", params$output_pdf_filepath))
  futile.logger::flog.info(glue::glue("Bootstrap display cutoff (%; ignored if 'NA'): ", params$bootstrap_cutoff))
  futile.logger::flog.info(glue::glue("Root name (ignored if 'NA'): ", params$root_name))
  futile.logger::flog.info(glue::glue("Input genome metadata filepath (ignored if 'NA'): ", 
                                      params$genome_metadata_filepath))
  futile.logger::flog.info(glue::glue("Input gene metadata filepath (ignored if 'NA'): ", 
                                      params$gene_metadata_filepath))
  futile.logger::flog.info(glue::glue("Plot width (mm): ", params$plot_width))
  futile.logger::flog.info(glue::glue("Plot height (mm): ", params$plot_height))
  futile.logger::flog.info(glue::glue("Write data table: ", params$write_data))
  futile.logger::flog.info("############################")
  
  # Convert character "NA" (from command line) into true NA
  params <- lapply(params, convert_to_constant_NA)
  
  # Load and plot the tree
  if (!is.na(params$input_phylogenetic_tree_filepath)) {
    phylo_tree_list <- load_and_plot_phylogenetic_tree(params$input_phylogenetic_tree_filepath, 
                                                       params$root_name, params$bootstrap_cutoff)
    
    # Get tip order of the tree, to match with heatmap later
    # Based on https://groups.google.com/forum/#!topic/bioc-ggtree/LqRDK78m3U4 (accessed Sept. 15, 2018)
    futile.logger::flog.debug("Exporting tip order of tree to correspond with heatmap")
    tip_order <- dplyr::filter(ggtree::ggtree(phylo_tree_list[[1]])$data, isTip == TRUE)
    tip_order <- tip_order[order(tip_order$y, decreasing = TRUE),]$label # plotting_name
  } else {
    flog.info("Skipping plotting phylogenetic tree")
    tip_order <- NA
  }
  
  # Load and plot the BLAST table as a heatmap
  blast_results_list <- load_and_plot_blast_results(params$input_blast_table_filepath, 
                                                tip_order, params$gene_metadata_filepath,
                                                params$genome_metadata_filepath)

  # Save the heatmap data
  if (params$write_data == TRUE) {
    futile.logger::flog.info("Saving raw heatmap data to file")
    output_table_filepath = paste(tools::file_path_sans_ext(params$output_pdf_filepath), ".tsv", sep = "")
    write.table(blast_results_list[[1]], file = output_table_filepath, sep = "\t",
              col.names = TRUE, row.names = FALSE, quote = FALSE)
  }

  # Save the plot
  if (!is.na(params$input_phylogenetic_tree_filepath)) {
    # Combine the tree and heatmap
    # Got ggarrange ideas from https://cran.r-project.org/web/packages/egg/vignettes/Ecosystem.html (accessed Sept. 15, 2018)
    futile.logger::flog.info("Combining the ggtree and the heatmap")
    combined_plot <- egg::ggarrange(phylo_tree_list[[2]], blast_results_list[[2]], 
                                    nrow = 1, widths = c(1, 1.5), heights = c(1), padding = unit(0, "mm"))
    
    # Print a PDF of the combined plot
    # N.B., dimensions need to be input in inches (25.4 mm per inch)
    futile.logger::flog.info("Saving to PDF")
    pdf(file = params$output_pdf_filepath, width = params$plot_width / 25.4, 
        height = params$plot_height / 25.4)
    print(combined_plot)
    dev.off()
  } else {
    # Save heatmap alone
    futile.logger::flog.info("Saving heatmap to PDF")
    
    pdf(file = params$output_pdf_filepath, width = params$plot_width / 25.4, 
        height = params$plot_height / 25.4)
    print(blast_results_list[[2]])
    dev.off()
  }
  futile.logger::flog.info("generate_heatmap.R: done.")
}

if ( !interactive() ) {
  parser <- argparser::arg_parser(
      description = glue::glue("generate_heatmap.R: Binds a phylogenetic tree to a BLAST table heatmap.
                                Copyright Lee H. Bergstrand and Jackson M. Tsuji, 2019."))
  
  # Add required args
  parser <- argparser::add_argument(parser = parser, arg = "input_phylogenetic_tree_filepath", 
                                    help = "Input phylogenetic tree filepath (set to 'NA' to plot heatmap only)", 
                                    type = "character", default = NULL)
  parser <- argparser::add_argument(parser = parser, arg = "input_blast_table_filepath", 
                                    help = "Input BLAST table filepath", 
                                    type = "character", default = NULL)
  parser <- argparser::add_argument(parser = parser, arg = "output_pdf_filepath", 
                                    help = "Output PDF filepath", 
                                    type = "character", default = NULL)
  
  # Add optional args (set to 'NA' to ignore)
  parser <- argparser::add_argument(parser = parser, arg = "--genome_metadata_filepath", short = "-m",
                                    help = "Genome metadata filepath", 
                                    type = "character", default = NA)
  parser <- argparser::add_argument(parser = parser, arg = "--gene_metadata_filepath", short = "-g",
                                    help = "Gene metadata filepath",
                                    type = "character", default = NA)
  parser <- argparser::add_argument(parser = parser, arg = "--bootstrap_cutoff", short = "-b",
                                    help = "Bootstrap cutoff value",
                                    type = "numeric", default = NA)
  parser <- argparser::add_argument(parser = parser, arg = "--root_name", short = "-r",
                                    help = "Root name ('midpoint' to midpoint root or NA to keep the existing root; default NA)",
                                    type = "character", default = NA)
  parser <- argparser::add_argument(parser = parser, arg = "--plot_width", short = "-w",
                                    help = "Plot width (mm)",
                                    type = "numeric", default = 400)
  parser <- argparser::add_argument(parser = parser, arg = "--plot_height", short = "-z",
                                    help = "Plot height (mm)",
                                    type = "numeric", default = 200)
  parser <- argparser::add_argument(parser = parser, arg = "--write_data", short = "-d",
                                    help = "Write raw plotting data to disk (same basepath as the PDF, but as a .tsv file)",
                                    flag = TRUE)
  
  params <- argparser::parse_args(parser)
  
  main(params)
}
