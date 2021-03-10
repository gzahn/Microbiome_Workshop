# -----------------------------------------------------------------------------#
# Microbiome analysis workshop
# Cleaning up the phyloseq object
# Author: Geoffrey Zahn
# Software versions:  R v 3.6.3
#                     tidyverse v 1.3.0
#                     phyloseq v 1.30.0
# -----------------------------------------------------------------------------#

start.time <- Sys.time()

# PACKAGES, SCRIPTS, AND SETUP ####
library(tidyverse); packageVersion("tidyverse")
library(phyloseq); packageVersion("phyloseq")

#################################################################################
#                               Main workflow                                   #
#        Remove non-bacteria, chloroplast and mitochondrial sequences           #
#                                                                               #
#################################################################################



# load phyloseq object with phylogenetic tree ####

ps <- readRDS("./output/ps_not-cleaned_w_tree.RDS")


# Find and remove non-bacteria ####
ps_nonbact <- subset_taxa(ps, Kingdom != "Bacteria")

# quick plot to look at kingdom-level taxonomy
ps %>% transform_sample_counts(function(x){x/sum(x)}) %>%
  plot_bar(fill="Kingdom")
ggsave("./output/figs/16S_Kingdom-Level_Taxonomic_Proportions.png",dpi=300) # save figure for later use

# same plot, but non-bacteria, for sanity check
ps_nonbact %>% 
  transform_sample_counts(function(x){x/sum(x)}) %>%
  plot_bar(fill="Kingdom")


# REMOVE NON-BACTERIA, CHLOROPLASTS, MITOCHONDRIA, and empty samples/taxa ####
ps <- subset_taxa(ps, Kingdom == "Bacteria")
ps <- subset_taxa(ps,Class != "Chloroplast")
ps <- subset_taxa(ps, taxa_sums(ps) > 0)
ps <- subset_samples(ps, sample_sums(ps) > 0)

# Save DNA sequences apart from rownames (from subsetted ps object)
seqs <- taxa_names(ps)
seqs <- DNAStringSet(seqs)
saveRDS(seqs,"./output/16S_ASV_reference_sequences.RDS")

# Save RDS object for cleaned up Phyloseq object
saveRDS(ps, file = "./output/clean_phyloseq_object.RDS")

end.time <- Sys.time()

end.time - start.time