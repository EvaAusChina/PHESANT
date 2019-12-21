library(data.table)
source("../summarise_phenotypes.r")

# Biomarkers for first 100K run
# n_chunks <- 1
# filename_root <- "biomarkers/pharma_exomes_biomarkers_parsed_output_100k_chunk."
# variable_info <- "variable-info/outcome_info_final_round3.tsv"
# coding_info_file <- "variable-info/data-coding-ordinal-info.txt"

# Note: coding info file for the first 100K run was the following file, rather than the recently updated file referred to below.
# coding_info_file <- "variable-info/data-coding-ordinal-info.txt"

n_exomes <- '200k'
n_chunks <- 10

outcome_info_final <- fread("G:/COURSERA/GITHUB/PHESANT/variable-info/outcome_info_final_pharma_nov2019.tsv")
variable_info <- outcome_info_final
# variable_info <- "variable-info/outcome_info_final_pharma_nov2019.tsv"
filename_root <- paste0("../../pharma_exomes_parsed_output_", n_exomes, "_chunk.")
coding_info_file <- "../variable-info/data-coding-ordinal-info-nov2019-update.txt"

# Read in all the data table codings
to_read <- paste("../WAS/codings/", dir("../WAS/codings"), sep="")
codings_tables <- list()

for(i in to_read) {
	name <- gsub("../WAS/codings/coding(.*).tsv", "\\1", i)
	codings_tables[[name]] <- read.table(i, header=TRUE, sep='\t', quote="", stringsAsFactors=FALSE)
}

for (i in 1:n_chunks)
{
	# Both sexes
	hist_filename <- paste0(filename_root, i, "_hist")
	pheno_summary <- paste0(filename_root, i, "_phenosummary.tsv")

	filename <- paste0(filename_root, i)
	tsv_filename <- paste0(filename, ".tsv")
	log_filename <- paste0(filename_root, i)
	log_file <- paste0(log_filename, ".log")
	tsv_data <- read.table(tsv_filename, header=TRUE, sep='\t')
	names(tsv_data)[1] <- "userId"

	outcome_info <- read.table(variable_info, sep='\t', quote="", comment.char="", header=TRUE)
	print(paste0("both sexes ", i))
	summary_file <-  get_hists_and_notes(hist_filename, tsv_data, log_file, outcome_info, codings_tables, samples_for_inclusion=TRUE, check=FALSE, start_column=4)
	write.table(summary_file, file=pheno_summary, col.names=TRUE, row.names=TRUE, sep='\t', quote=FALSE)
}
