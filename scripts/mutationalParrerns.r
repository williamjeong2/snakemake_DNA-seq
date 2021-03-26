list.of.packages <- c("optparse", "data.table", "gridExtra", "R.utils")
list.of.bio.packages <- c("MutationalPatterns", "BSgenome", "ggplot2")

not.installed.bio.packages <- list.of.bio.packages[!(list.of.bio.packages %in% installed.packages()[, "Package"])]
if(length(not.installed.bio.packages)){
  if(!requireNamespace("BiocManager", qqietly = TRUE))
    install.packages("BiocManager")
    BiocManager::install(not.installed.bio.packages, suppressUpdates = TRUE)
}
lapply(list.of.bio.packages, require, character.only = TRUE)
lapply(list.of.packages, require, character.only = TRUE)

# arguments to provide
option_list = list(
  # default params
  make_option("--input", type="character", metavar="character"),
  make_option("--ref", type="character", default="BSgenome.Hsapiens.UCSC.hg38",metavar="character"), # set reference 
  make_option("--outdir", type="character", metavar="character"), #dir to save
)

# parse the command-line arguments and pass them to a list called 'opt'
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

detectGroups <- function(x){
  tem <- gsub("[0-9]*$", "", x)
  tem <- gsub("_$", "", tem)
  return(tem)
}

savePlot <- function(path, plot){
  for(i in c("png", "svg")) {
    ggsave(filename = paste0(opt$outdir, path, ".", i),
    plot = plot,
    device = i, scale = 2,
    width = 7, height = 7, units = "in",
    dpi = 320)
  }
}
# Data

library(BSgenome)
if (grepl("homo|human|sapiens", opt$ref, ignore.case = TRUE)){
  ref_genome <- "BSgenome.Hsapiens.UCSC.hg38"
}
if (grepl("mouse|mus|musculur", opt$ref, ignore.case = TRUE)){
  ref_genome <- "BSgenome.Mmusculus.UCSC.mm10"
}

BiocManager::install(pkgs = ref_genome)
library(ref_genome, character.only = TRUE)

vcf_files <- fread(opt$input)
colName <- colnames(vcf_files)[10:ncol(vcf_files)]
sample_names <- colnames(colName)

# comparison <- read_excel(opt$metadata, skip = 2, sheet = "comparisons")
# comp = as.data.frame(paste(paste(comparison$ctrl_cellType,comparison$ctrl_group, sep = "_"), paste(comparison$treat_cellType,comparison$treat_group, sep = "_"), sep = "-"))
# colnames(comp) <- "comparison_list"

# for(i in 1:nrow(comp)){
#   ifelse(!dir.exists(paste0(opt$outdir, comp[i, ])), dir.create(paste0(opt$outdir, comp[i, ]), showWarnings = FALSE, recursive = T), print("pass"))
# }

grl <- read_vcfs_as_granges(vcf_files, sample_names, ref_genome)
indel_grl <- read_vcfs_as_granges(vcf_files, sample_name, ref_genome, type = "indel")
dbs_grl <- read_vcfs_as_granges(vcf_files, sample_name, ref_genome, type = "dbs")
mbs_grl <- read_vcfs_as_granges(vcf_files, sample_name, ref_genome, type = "bms")
tissue <- detectGroups(sample_names)

# 4. Mutation characteristics
# 4.1 SNVS
# 4.1.1 Base substitution types
muts <- mutations_from_vcf(grl[[1]])
types <- mut_type(grl[[1]])
context <- mut_context(grl[[1]], ref_genome)

type_context <- type_context(grl[[1]], ref_genome)
type_occurrences <- mut_type_occurrences(grl, ref_genome)

# 4.1.2 Muatation spectrum
library("gridExtra")

p1 <- plot_spectrum(type_occurrences)
p2 <- plot_spectrum(type_occurrences, CT = TRUE)
p3 <- plot_spectrum(type_occurrences, CT = TRUE, 
                    indv_points = TRUE, legend = FALSE)
p <- grid.arrange(p1, p2, p3, ncol = 3, widths = c(3, 3, 1.75))
savePlot("1_Mutation_spectrum_1", p)

p4 <- plot_spectrum(type_occurrences, by = tissue, CT = TRUE, legend = TRUE)
p5 <- plot_spectrum(type_occurrences, CT = TRUE, legend = TRUE, error_bars = "stdev")
p <- grid.arrange(p4, p5, ncol = 2, widths = c(4, 2.3))
savePlot("1_Mutation_spectrum_2", p)

# 4.1.3 96 mutational profile
mut_mat <- mut_matrix(vcf_list = grl, ref_genome = ref_genome)
p <- plot_96_profile(mut_mat)
savePlot("2_96_mutational_profile", p)

# 4.2 Indels
indel_grl <- get_indel_context(indel_grl, ref_genome)
indel_counts <- count_indel_contexts(indel_grl)

savePlot("3_indel_spectra_1", plot_indel_contexts(indel_counts, condensed = TRUE))
savePlot("3_indel_spectra_2", plot_main_indel_contexts(indel_counts))

# 4.3 DBSs
dbs_grl <- get_dbs_context(dbs_grl)
dbs_counts <- count_dbs_contexts(dbs_grl)

savePlot("4_DBSs_contexts", plot_dbs_contexts(dbs_counts, same_y = TRUE))
savePlot("4_DBSs_based_on_reference_bases", plot_main_dbs_contexts(dbs_counts, same_y = TRUE))

# 4.4 MBSs
mbs_counts <- count_mbs_contexts(mbs_grl)

savePlot("5_MBSs_contexts", plot_mbs_contexts(mbs_counts, sampe_y = TRUE))


# 5.2 Signature refitting
# 5.2.1 Find mathematically optimal contribution of COSMIC signatures
signatures = get_known_signatures()
fir_res <- fir_to_signatures(mut_mat, signatures)

p <- plot_contribution(fir_res$contribution,
  coord_flip = FALSE,
  mode = "absolute")
savePlot("6_contribution_of_COSMIC_signatures", p)

p <- plot_original_vs_reconstructed(mut_mat, fit_res$reconstructed, 
  y_intercept = 0.95)
savePlot("6_cosine_similarity_with_reconstructed_profiles", p)

# 7. Genomic distribution
# 7.1 Rainfall plot
# Define autosomal chromosomes
chromosomes <- seqnames(get(ref_genome))[1:22]

# Make a rainfall plot
p <- plot_rainfall(grl[[1]],
  title = names(grl[1]),
  chromosomes = chromosomes, cex = 1.5, ylim = 1e+09
)
savePlot("7_rainfall_plot", p)