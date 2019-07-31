#!/usr/bin/env Rscript
# generate_BackBLAST_heatmap.R
# Copyright Lee Bergstrand and Jackson M. Tsuji, 2018
# Plots a newick treefile and BLAST table together as a phylogenetic tree and heatmap
# See required R packages below.

# Load libraries
library(argparser, quietly = TRUE)
library(futile.logger)
library(glue, warn.conflicts = FALSE)
library(plyr, warn.conflicts = FALSE)
library(dplyr, warn.conflicts = FALSE)
library(tibble, warn.conflicts = FALSE)
library(reshape2, warn.conflicts = FALSE)
suppressMessages(library(ggplot2, warn.conflicts = FALSE))
library(RColorBrewer, warn.conflicts = FALSE)
suppressMessages(library(ggtree))
library(ape, warn.conflicts = FALSE)
library(maps, warn.conflicts = FALSE)
library(phytools, warn.conflicts = FALSE)
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

#' Writes a data table to a file with several default parameters. All params below are the same as in write.table()
#' 
#' @param x Data frame or tibble
#' @param file Output filepath
#' @param sep Field separator character
#' @param row.names Logical; display row names in the output table?
#' @param col.names Logical; display column names in the output table?
#' @param quote Logical; add quotation marks around fields?
#' @export
write_table <- function(x, file, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE) {
  write.table(x = x, file = file, sep = sep, row.names = row.names, col.names = col.names, quote = quote)
}

#' Chooses a nice discrete colour scale of the desired length
#' 
#' @param length Numeric; length of the desired colour scale vector
#' @return Character vector of HTML colour codes
#' @export
choose_discrete_colour_scale <- function(length) {
  
  # Choose the best colour scale based on the number of entries to be plotted
  if ( length == 2 ) {
    colour_palette <- RColorBrewer::brewer.pal(n = 3, name = "Dark2")[c(1,2)]
  } else if ( length <= 8 ) {
    colour_palette <- RColorBrewer::brewer.pal(n = length, name = "Dark2")
  } else if ( length <= 12 ) {
    colour_palette <- RColorBrewer::brewer.pal(n = length, name = "Set3")
  } else if ( length > 12 ) {
    colour_palette <- scales::hue_pal(h = c(20,290))(length)
  } else {
    futile.logger::flog.error(glue::glue("Something is wrong with the provided length ('", length, 
                          "'). Is it non-numeric? Exiting..."))
    quit(status = 1, save = "no")
  }
  
  return(colour_palette)
  
}

#' Re-root a ggtree
#' 
#' @param phylo_tree ggtree-format tree
#' @param root_name Character (length 1 vector); the exact name of the to-be root
#' @return ggtree-format tree, rooted
#' @export
reroot_ggtree <- function(phylo_tree, root_name) {
  
  # Extract the data from the tree in tabular format
  tree_data <- ggtree::ggtree(phylo_tree)$data
  
  # Check that the root_name exists
  if ( (root_name %in% tree_data$label) == FALSE ) {
    futile.logger::flog.error(glue::glue("Could not find the root_name '", root_name, 
                          "' in the provided tree. Cannot re-root. Exiting..."))
    quit(status = 1, save = "no")
  }
  
  # Get the ggtree ID # of the to-be root
  tip_label_index <- match(x = root_name, table = phylo_tree$tip.label)
  
  # Check that only one matching parent node exists
  if (length(tip_label_index) != 1) {
    futile.logger::flog.error(glue::glue("The provided root_name '", root_name, 
                          "' appears to have more than one associated node. Cannot re-root. Exiting..."))
    quit(status = 1, save = "no")
  }
  
  # Re-root
  # TODO - consider suppressMessages()
  tree_rooted <- ggtree::reroot(phylo_tree, node = tip_label_index)
  
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
  bootstrap_label_data$label <- as.numeric(bootstrap_label_data$label)
  # TODO - make this more elegant to avoid warning
  
  # Filter out tip labels and 'NA' labels
  bootstrap_label_data <- dplyr::filter(bootstrap_label_data, is.na(label) == FALSE & isTip == FALSE)
  
  # Set bootstrap cutoff if desired
  if ( is.na(bootstrap_cutoff) == FALSE ) {
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
    # TODO - fine-tune the 'nudge' value so that it works more generically
    geom_text(data = bootstrap_label_data, aes(label = label), nudge_x = -0.035, nudge_y = 0.35, size = 3) +
    # TODO - ALT CODE
    # geom_text2(aes(subset = (grepl(pattern = "^[0-9]+$", x = label) & !(isTip) & as.numeric(label) > bootstrap_cutoff), 
    #                label = as.numeric(label)),
    #            nudge_x = -0.03, nudge_y = 0.4, size = 2.5) +

    # Add the dotted lines from the tree tips
    # TODO - check if 'size' is okay
    geom_tiplab(align = TRUE, linetype = "dotted", size = 0.1, offset = 0.1) +
    
    # Add scale bar
    # TODO - set y based on tree topography?
    geom_treescale(x = 0, y = 8, linesize = 1, fontsize = 3, offset = 0.5) +
    
    # Manually tune the y-axis boundaries to match the heatmap: https://stackoverflow.com/a/18718196 (accessed Sept. 15, 2018)
    # TODO - change this to be a function of the number of entries
    scale_y_discrete(expand = c(0,0.6)) +
  
    ## Manually set the righthand cutoff for the tree
    # TODO - change this to be a function of the tree branch length?? Or is 1 always okay?
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
  if ( is.na(root_name) == FALSE && root_name != "NA" ) {
    futile.logger::flog.info(glue::glue("Re-rooting tree to '", root_name, "'"))
    phylo_tree <- reroot_ggtree(phylo_tree, root_name)
  }
  
  # Set cutoff for bootstraps externally, to be overlaid onto the tree figure later
  # No cutoff is applied if bootstrap_cutoff is set to 'NA'
  futile.logger::flog.info("Generating bootstrap labels")
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
read_blast_tibble <- function(input_blast_table_filepath) {
  
  # Load the data
  blast_tibble <- read_tibble(input_blast_table_filepath, sep = ",")
  
  # HARD-CODED expected header names
  expected_header_names <- c("subject_name", "qseqid", "sseqid", "pident", "evalue", "qcovhsp", "bitscore")
  
  # Confirm the columns look okay
  if ( identical(colnames(blast_tibble), expected_header_names) == FALSE ) {
    
    futile.logger::flog.error(glue::glue("BLAST data table did not have expected column names ('", 
                          glue::glue_collapse(expected_header_names, sep = "; "), 
                          "'). Instead, had: '", glue::glue_collapse(colnames(blast_tibble), sep = "; "), 
                          "'. Exiting..."))
    quit(save = "no", status = 1)
    
  }
  
  return(blast_tibble)
  
}


#' Changes the order of the subject_name in the BLAST table to match that of the ggtree tips
#' 
#' @param blast_tibble Tibble output of read_blast_tibble
#' @param tip_order Character vector of the exact order of the tip labels in the phylogenetic tree
#' @return tibble of BLAST data with subject_name as an ordered factor and reduced to match phylogenetic tree names (if needed)
#' @export
order_blast_tibble_subjects <- function(blast_tibble, tip_order) {
  
  # If entries are not identical, try filtering down to just what is in the phylogenetic tree
  if ( identical(sort(unique(blast_tibble$subject_name)), 
                 sort(tip_order)) == FALSE ) {
    
    # Make sure that the tree tips are contained in the blast tibble; if so, everything is okay
    # The length should be zero if all tree tips are contained in the blast tibble
    if(length(setdiff(tip_order, blast_tibble$subject_name)) > 0) {
      
      # But if the length is > 0, it means some entries in the tree are MISSING in the blast tibble, a major issue
      missing_blast_tibble_entries <- setdiff(blast_tibble$subject_name, tip_order)
      flog.error(glue::glue("The provided BLAST tibble is missing some entries in the phylogenetic tree: '",
                            glue::glue_collapse(missing_blast_tibble_entries, sep = ", "),
                            "'. Cannot continue -- exiting..."))
      quit(save = "no", status = 1)
      
    }
    
    # Report to the user which entries in the BLAST tibble are to be ignored
    extra_blast_tibble_entries <- !(unique(blast_tibble$subject_name) %in% tip_order)
    extra_blast_tibble_entries <- unique(blast_tibble$subject_name)[extra_blast_tibble_entries]
    
    futile.logger::flog.info(glue::glue("Some entries in the BLAST table are missing in the phylogenetic tree ",
                                        "and will be removed in plotting: '", 
                                        glue::glue_collapse(extra_blast_tibble_entries, sep = ", "),  "'."))
    
    # Filter down the BLAST tibble to have the same subject_name's as the tree
    blast_tibble <- dplyr::filter(blast_tibble, subject_name %in% tip_order)
    
  }
  
  # Make subjects the same order as in the tree
  blast_tibble$subject_name <- factor(blast_tibble$subject_name, levels = rev(tip_order), ordered = TRUE)
  
  return(blast_tibble)
}


#' Overlays user-given genome names for the heatmap y-axis
#' 
#' @param blast_tibble Tibble output of read_blast_tibble
#' @param genome_metadata_filepath Character (length 1 vector); the filepath of the tab-separated genome metadata file.
#' Must at least have the columns 'subject_name' and 'plotting_name' as the first and second columns, respectively.
#' @return tibble of BLAST data with subject_name renamed with the user-desired values
#' @export
overlay_genome_naming <- function(blast_tibble, genome_metadata_filepath) {
  
  # Load the metadata tibble
  futile.logger::flog.info("Loading genome metadata")
  genome_metadata <- read_tibble(tree_metadata_filename)
  
  # Check that the column names match expected (for the first two columns; doesn't matter after that)
  # HARD-CODED
  genome_metadata_expected_headers <- c("subject_name", "plotting_name")
  
  if ( identical(colnames(genome_metadata)[1:2], genome_metadata) == FALSE ) {
    futile.logger::flog.error(glue::glue("The first two columns of the genome metadata table should be: ", 
                                         glue::glue_collapse(genome_metadata_expected_headers, sep = ", "), 
                                         ". However, you provided something else: ", 
                                         glue::glue_collapse(colnames(genome_metadata)[1:2], sep = ", "), ". Exiting..."))
    quit(save = "no", status = 1)
  }
  
  # Reduce the metadata down to just the expected headers
  genome_metadata <- dplyr::select(genome_metadata, genome_metadata_expected_headers)
  
  # Reduce the metadata down to entries contained within the BLAST table
  # TODO - consider reporting to user what was omitted
  genome_metadata <- dplyr::filter(genome_metadata, subject_name %in% unique(blast_tibble$subject_name))
  
  # Change the subject_name to be the plotting_name in the BLAST table
  # TODO - change the warning here to a complete error out
  blast_tibble$subject_name <- plyr::mapvalues(blast_tibble$subject_name, from = genome_metadata$subject_name, 
                                               to = genome_metadata$plotting_name, warn_missing = TRUE)
  
  return(blast_tibble)
  
}


#' Overlays user-given gene names for the heatmap x-axis
#' 
#' @param blast_tibble Tibble output of read_blast_tibble
#' @param gene_metadata_filepath Character (length 1 vector); the filepath of the tab-separated gene metadata file.
#' Must at least have the columns 'qseqid' and 'gene_name' as the first and second columns, respectively.
#' @return tibble of BLAST data with qseqid renamed with the user-desired values
#' @export
overlay_gene_naming <- function(blast_tibble, gene_metadata_filepath) {
  
  # Load gene naming table
  gene_naming_tibble <- read_tibble(gene_metadata_filepath)
  
  # Check that the column names match expected (for the first two columns; doesn't matter after that)
  # HARD-CODED
  gene_tibble_expected_headers <- c("qseqid", "gene_name")
  
  if ( identical(colnames(gene_naming_tibble)[1:2], gene_tibble_expected_headers) == FALSE ) {
    futile.logger::flog.error(glue::glue("The first two columns of the gene table should be: ", 
                          glue::glue_collapse(gene_tibble_expected_headers, sep = ", "), 
                          ". However, you provided something else: ", 
                          glue::glue_collapse(colnames(gene_naming_tibble)[1:2], sep = ", "), ". Exiting..."))
    quit(save = "no", status = 1)
  }
  
  # Check that the query seq IDs in the gene table include those provided as queries for BLAST (can contain extra - that's fine)
  gene_table_queries <- sort(gene_naming_table$qseqid)
  blast_tibble_queries <- sort(unique(blast_tibble$qseqid))
  if ( length(unique(blast_tibble_queries %in% gene_table_queries)) > 1 ) {
    
    # Determine missing entries
    missing_entries <- blast_tibble_queries[!(blast_tibble_queries %in% gene_table_queries)]
    good_entries <- blast_tibble_queries[(blast_tibble_queries %in% gene_table_queries)]
    
    # Warn user
    futile.logger::flog.warn(glue::glue("Some query sequence IDs from the BLAST table are not contained ",
                                        "in the gene naming table provided. These will be REMOVED from the ",
                                        "heatmap plot! The following will be removed: '", 
                                        glue::glue_collapse(missing_entries, sep = "; "), "'."))
    
    # Remove entries from BLAST table and proceed
    blast_tibble <- dplyr::filter(blast_tibble, qseqid %in% good_entries)
    
  } else if ( unique(blast_tibble_queries %in% gene_table_queries) == FALSE ) {
    flog.error("None of the query sequence IDs from the BLAST table match the gene naming table provided. Exiting...")
    quit(save = "no", status = 1)
  }
  
  # Change qseqid to gene_name and order according to the gene_naming_table
  blast_tibble$qseqid <- plyr::mapvalues(x = blast_tibble$qseqid, from = gene_naming_table$qseqid, 
                                        to = gene_naming_table$gene_name)
  blast_tibble$qseqid <- factor(blast_tibble$qseqid, levels = gene_naming_table$gene_name, ordered = TRUE)
  
  return(blast_tibble)
  
}


#' Plots the BLAST table as a heatmap in ggplot
#' 
#' @param blast_tibble Tibble output of read_blast_tibble
#' @return ggplot heatmap
#' @export
plot_blast_heatmap <- function(blast_tibble) {
  
  # Add NA values for missing grid values so that grid lines will appear in the final plot
  # A bit hacky
  blast_tibble <- reshape2::dcast(blast_tibble, subject_name ~ qseqid, value.var = "pident") %>%
    reshape2::melt(na.rm = FALSE, id.vars = c("subject_name"), variable.name = "qseqid",
                   value.name = "pident") %>%
    tibble::as_tibble()
  
  blast_heatmap <- ggplot2::ggplot(blast_tibble, aes(y = subject_name, x = qseqid)) +
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


#' Master function to load and plot the BLAST table and gene_naming_table to produce a heatmap
#' 
#' @param input_blast_table_filepath Character (length 1 vector); the filepath to the BLAST table (comma-separated)
#' @param tip_order Character vector of the exact order of the tip labels in the phylogenetic tree
#' @param genome_metadata_filepath Character (length 1 vector); the filepath of the tab-separated genome metadata file. 
#' Must at least have the columns 'subject_name' and 'plotting_name' as the first and second columns, respectively.
#' @param gene_metadata_filepath Character (length 1 vector); the filepath of the tab-separated gene metadata file. 
#' Must at least have the columns 'qseqid' and 'gene_name' as the first and second columns, respectively.
#' @return list of two: 'blast_tibble' data frame; 'blast_heatmap' ggplot object
#' @export
load_and_plot_blast_tibble <- function(input_blast_table_filepath, tip_order, gene_metadata_filepath,
                                       genome_metadata_filepath) {
  
  # Load the blast table
  futile.logger::flog.info("Loading the BLAST table")
  blast_tibble <- read_blast_tibble(input_blast_table_filepath)
  
  # Order the BLAST table subject names to match the ggtree
  futile.logger::flog.info("Aligning BLAST table's subject names to match order of the ggtree")
  blast_tibble <- order_blast_tibble_subjects(blast_tibble, tip_order)
  
  # Overlay gene names and gene naming order, if provided
  if ( is.na(genome_metadata_filepath) == FALSE ) {
    futile.logger::flog.info("Overlaying genome naming and ordering onto the BLAST table")
    blast_tibble <- overlay_genome_naming(blast_tibble, genome_metadata_filepath)
  }
  
  # Overlay gene names and gene naming order, if provided
  if ( is.na(gene_metadata_filepath) == FALSE ) {
    futile.logger::flog.info("Overlaying gene naming and ordering onto the BLAST table")
    blast_tibble <- overlay_gene_naming(blast_tibble, gene_metadata_filepath)
  }
  
  # Create the heatmap
  futile.logger::flog.info("Plotting BLAST heatmap")
  blast_heatmap <- plot_blast_heatmap(blast_tibble)
  
  # Make output list
  blast_tibble_list <- list(blast_tibble, blast_heatmap)
  names(blast_tibble_list) <- c("blast_tibble", "blast_heatmap")
  
  return(blast_tibble_list)
  
}


main <- function(params) {
  # Startup messages
  futile.logger::flog.info("Running generate_BackBLAST_heatmap.R")
  futile.logger::flog.info(glue::glue("Input phylogenetic tree filepath: ", params$input_phylogenetic_tree_filepath))
  futile.logger::flog.info(glue::glue("Input BLAST table filepath: ", params$input_blast_table_filepath))
  futile.logger::flog.info(glue::glue("Output PDF filepath: ", params$output_pdf_filepath))
  futile.logger::flog.info(glue::glue("Bootstrap display cutoff (%; ignored if 'NA'): ", params$bootstrap_cutoff))
  futile.logger::flog.info(glue::glue("Root name (ignored if 'NA'): ", params$root_name))
  futile.logger::flog.info(glue::glue("Genome metadata filepath (ignored if 'NA'): ", params$genome_metadata_filepath))
  futile.logger::flog.info(glue::glue("Input gene naming table filepath (ignored if 'NA'): ", 
                                      params$gene_metadata_filepath))
  
  # Load and plot the tree
  # TODO - after reading metadata table, CHECK that the first column names correspond to the tree tip labels
  phylo_tree_list <- load_and_plot_phylogenetic_tree(params$input_phylogenetic_tree_filepath, 
                                                     params$root_name, params$bootstrap_cutoff)
  
  # Get tip order of the tree, to match with heatmap later
  # Based on https://groups.google.com/forum/#!topic/bioc-ggtree/LqRDK78m3U4 (accessed Sept. 15, 2018)
  futile.logger::flog.info("Exporting tip order of tree to correspond with heatmap")
  tip_order <- dplyr::filter(ggtree::ggtree(phylo_tree_list[[1]])$data, isTip == TRUE)
  tip_order <- tip_order[order(tip_order$y, decreasing = TRUE),]$label # plotting_name

  
  # Load and plot the BLAST table as a heatmap
  blast_tibble_list <- load_and_plot_blast_tibble(params$input_blast_table_filepath, 
                                                tip_order, params$gene_metadata_filepath,
                                                params$genome_metadata_filepath)
  
  # Combine the tree and heatmap
  # Got ggarrange ideas from https://cran.r-project.org/web/packages/egg/vignettes/Ecosystem.html (accessed Sept. 15, 2018)
  # TODO - make the dimensions a function of the gene number and the subject number
  futile.logger::flog.info("Combining the ggtree and the heatmap")
  combined_plot <- egg::ggarrange(phylo_tree_list[[2]], blast_tibble_list[[2]], 
                                  nrow = 1, widths = c(5, 8), heights = c(5), padding = unit(0, "mm"))
  # N.B., set debug = TRUE to see grid lines
  
  # Print a PDF of the combined plot
  # TODO - set dimensions as a function of the gene number and the subject number
  futile.logger::flog.info("Saving to PDF")
  pdf(file = params$output_pdf_filepath, width = 17, height = 7)
  print(combined_plot)
  dev.off()
  
  futile.logger::flog.info("generate_BackBLAST_heatmap.R: done.")
  
}

if (interactive() == FALSE) {
  
  parser <- argparser::arg_parser(description = glue::glue("generate_BackBLAST_heatmap.R: Binds a phylogenetic tree to a BLAST table heatmap.
                                             Copyright Lee Bergstrand and Jackson M. Tsuji, 2019."))
  
  # Add required args
  parser <- argparser::add_argument(parser = parser, arg = "input_phylogenetic_tree_filepath", 
                                    help = "Input phylogenetic tree filepath", 
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
                                    help = "Root name",
                                    type = "character", default = NA)
  
  params <- argparser::parse_args(parser)
  
  main(params)
  
} else {
  
  # Manually set when running in RStudio for development
  params <- list()
  setwd("/home/jmtsuji/Research_General/Bioinformatics/02_git/BackBLAST_Reciprocal_BLAST/test")
  
  # Required inputs
  params$input_phylogenetic_tree_filepath <- "GSB_riboproteins_tree_vs3.treefile"
  params$input_blast_table_filepath <- "GSB_pathway_vs7_e40_best_hits_MOD3.csv"
  params$output_pdf_filepath <- "test4.pdf"
  
  # Optional inputs (set to 'NA' to ignore)
  params$genome_metadata_filepath <- NA
  params$gene_metadata_filepath <- NA
  params$bootstrap_cutoff <- 90
  params$root_name <- "Ignavibacterium_album_JCM_16511_NC_017464.1" # Optional; set to NA if you want to use the tree as-is.
  
  main(params)
  
}

