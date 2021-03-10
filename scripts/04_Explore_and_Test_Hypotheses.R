# -----------------------------------------------------------------------------#
# Microbiome analysis workshop
# Exploring cleaned data using phyloseq and corncob packages
# Author: Geoffrey Zahn
# Software versions:  R v 3.6.3
#                     tidyverse v 1.3.0
#                     phyloseq v 1.30.0
#                     purrr v 0.3.4
#                     Biostrings v 2.54.0
#                     corncob v 0.1.0
#                     vegan v 2.6.0
# -----------------------------------------------------------------------------#

start.time <- Sys.time()

# PACKAGES, SCRIPTS, AND SETUP ####
library(tidyverse); packageVersion("tidyverse")
library(phyloseq); packageVersion("phyloseq")
library(purrr); packageVersion("purrr")
library(Biostrings); packageVersion("Biostrings")
library(corncob); packageVersion("corncob")
library(vegan); packageVersion("vegan")


#################################################################################
#                               Main workflow                                   #
#  Explore alpha and beta diversity, visualize data set, test hypotheses,       #
#  search for differentially abundant taxa                                      #
#                                                                               #
#################################################################################


# Load cleaned phyloseq object ####

ps <- readRDS("./output/clean_phyloseq_object.RDS")


# Alpha diversity metrics ####

# richness

# simpson / shannon

# transform raw counts to relative abundance ####

# Alpha diversity by grouped data

# merge samples (using raw counts)

# plot_bar

# Beta-diversity ####

# Ordination

# permanova



# Network analysis





