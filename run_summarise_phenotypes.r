source("summarise_phenotypes.r")

print("I am here")
## VERNERI'S APPLICATION.

# Running Verneri's application through.
hist_filename <- "~/results/ukb1859_hist"
pheno_summary <- "~/results/ukb1859_phenosummary.tsv"

filename <- "~/results/ukb7127_output"
tsv_filename <- paste(filename, ".tsv", sep="")
log_file <- paste(filename, ".log", sep="")
tsv_data <- read.table(tsv_filename, header=TRUE, sep='\t')

qc_data <- read.table("~/results/ukb1859_qc.tsv", header=TRUE)

outcome_info <- read.table("~/Repositories/PHESANT/variable-info/outcome_info_final.tsv",
					   sep='\t', quote="", comment.char="", header=TRUE)

samples_for_removal <- as.character(read.table("~/results/w1859_20170726_participantwithdrawallist.csv")$V1)

notes_for_manny_1859 <-  get_hists_and_notes(hist_filename, tsv_data, log_file, outcome_info, codings_tables, qc_data, samples_for_removal)

write.table(notes_for_manny_1859, file=pheno_summary, col.names=TRUE, row.names=TRUE, sep='\t', quote=FALSE)

pheno_summary_PHESANT <- "~/results/ukb1859_phenosummary_final.tsv"
notes_for_manny_1859_PHESANT_codings_included <- include_PHESANT_reassignment_names(pheno_summary, outcome_info)
write.table(notes_for_manny_1859_PHESANT_codings_included, file=pheno_summary_PHESANT, col.names=TRUE, row.names=TRUE, sep='\t', quote=FALSE)


# Check against the samples you get if using Liam's collection.
hist_filename <- "~/results/ukb1859_hist_check"
pheno_summary <- "~/results/ukb1859_phenosummary_check.tsv"

filename <- "~/results/ukb7127_output"
tsv_filename <- paste(filename, ".tsv", sep="")
log_file <- paste(filename, ".log", sep="")
tsv_data <- read.table(tsv_filename, header=TRUE, sep='\t')

qc_data <- ""

outcome_info <- read.table("~/Repositories/PHESANT/variable-info/outcome_info_final.tsv",
					   sep='\t', quote="", comment.char="", header=TRUE)

samples_for_inclusion <- as.character(read.table("~/results/ukb1859_qc.sample_list", header=TRUE)$sample)
check <- TRUE

notes_for_manny_1859 <-  get_hists_and_notes(hist_filename, tsv_data, log_file, outcome_info, codings_tables, qc_data, samples_for_removal, samples_for_inclusion, check)

write.table(notes_for_manny_1859, file=pheno_summary, col.names=TRUE, row.names=TRUE, sep='\t', quote=FALSE)

pheno_summary_PHESANT <- "~/results/ukb1859_phenosummary_check_final.tsv"
notes_for_manny_1859_PHESANT_codings_included <- include_PHESANT_reassignment_names(pheno_summary, outcome_info)
write.table(notes_for_manny_1859_PHESANT_codings_included, file=pheno_summary_PHESANT, col.names=TRUE, row.names=TRUE, sep='\t', quote=FALSE)

# Summary of ALL the samples in Verneri's application.
hist_filename <- "~/results/ukb1859_hist_all"
pheno_summary <- "~/results/ukb1859_phenosummary_all.tsv"

filename <- "~/results/ukb7127_output"
tsv_filename <- paste(filename, ".tsv", sep="")
log_file <- paste(filename, ".log", sep="")
tsv_data <- read.table(tsv_filename, header=TRUE, sep='\t')

qc_data <- ""

outcome_info <- read.table("~/Repositories/PHESANT/variable-info/outcome_info_final.tsv",
					   sep='\t', quote="", comment.char="", header=TRUE)

samples_for_removal <- ""
samples_for_inclusion <- tsv_data$userId
check <- FALSE

notes_for_manny_1859 <-  get_hists_and_notes(hist_filename, tsv_data, log_file, outcome_info, codings_tables, qc_data, samples_for_removal, samples_for_inclusion, check)

write.table(notes_for_manny_1859, file=pheno_summary, col.names=TRUE, row.names=TRUE, sep='\t', quote=FALSE)

pheno_summary_PHESANT <- "~/results/ukb1859_phenosummary_all_final.tsv"
notes_for_manny_1859_PHESANT_codings_included <- include_PHESANT_reassignment_names(pheno_summary, outcome_info)
write.table(notes_for_manny_1859_PHESANT_codings_included, file=pheno_summary_PHESANT, col.names=TRUE, row.names=TRUE, sep='\t', quote=FALSE)


## JOEL'S APPLICATION.

# Running Joel's application through.
hist_filename <- "~/results/ukb1189_hist"
pheno_summary <- "~/results/ukb1189_phenosummary.tsv"

filename <- "~/results/ukb1189_output"
tsv_filename <- paste(filename, ".tsv", sep="")
log_file <- paste(filename, ".log", sep="")
tsv_data <- read.table(tsv_filename, header=TRUE, sep='\t')

qc_data <- read.table("~/results/ukb1189_qc.tsv", header=TRUE)

outcome_info <- read.table("~/Repositories/PHESANT/variable-info/outcome_info_final.tsv",
					   sep='\t', quote="", comment.char="", header=TRUE)

samples_for_removal <- as.character(read.table("~/results/w1189_20170726_participantwithdrawallist.csv")$V1)

notes_for_manny_1189 <-  get_hists_and_notes(hist_filename, tsv_data, log_file, outcome_info, codings_tables, qc_data, samples_for_removal)

write.table(notes_for_manny_1189, file=pheno_summary, col.names=TRUE, row.names=TRUE, sep='\t', quote=FALSE)

pheno_summary_PHESANT <- "~/results/ukb1189_phenosummary_final.tsv"
notes_for_manny_1189_PHESANT_codings_included <- include_PHESANT_reassignment_names(pheno_summary, outcome_info)
write.table(notes_for_manny_1189_PHESANT_codings_included, file=pheno_summary_PHESANT, col.names=TRUE, row.names=TRUE, sep='\t', quote=FALSE)

# Check against the samples you get if using Liam's collection.
hist_filename <- "~/results/ukb1189_hist_check"
pheno_summary <- "~/results/ukb1189_phenosummary_check.tsv"

filename <- "~/results/ukb1189_output"
tsv_filename <- paste(filename, ".tsv", sep="")
log_file <- paste(filename, ".log", sep="")
tsv_data <- read.table(tsv_filename, header=TRUE, sep='\t')

qc_data <- ""

outcome_info <- read.table("~/Repositories/PHESANT/variable-info/outcome_info_final.tsv",
					   sep='\t', quote="", comment.char="", header=TRUE)

samples_for_removal <- ""
samples_for_inclusion <- as.character(read.table("~/results/ukb1189_qc.sample_list", header=TRUE)$sample)

notes_for_manny_1189 <-  get_hists_and_notes(hist_filename, tsv_data, log_file, outcome_info, codings_tables, qc_data, samples_for_removal, samples_for_inclusion)

write.table(notes_for_manny_1189, file=pheno_summary, col.names=TRUE, row.names=TRUE, sep='\t', quote=FALSE)

pheno_summary_PHESANT <- "~/results/ukb1189_phenosummary_check_final.tsv"
notes_for_manny_1189_PHESANT_codings_included <- include_PHESANT_reassignment_names(pheno_summary, outcome_info)
write.table(notes_for_manny_1189_PHESANT_codings_included, file=pheno_summary_PHESANT, col.names=TRUE, row.names=TRUE, sep='\t', quote=FALSE)

# Summary of ALL the samples in Joel's application.
hist_filename <- "~/results/ukb1189_hist_all"
pheno_summary <- "~/results/ukb1189_phenosummary_all.tsv"

filename <- "~/results/ukb1189_output"
tsv_filename <- paste(filename, ".tsv", sep="")
log_file <- paste(filename, ".log", sep="")
tsv_data <- read.table(tsv_filename, header=TRUE, sep='\t')

qc_data <- ""

outcome_info <- read.table("~/Repositories/PHESANT/variable-info/outcome_info_final.tsv",
					   sep='\t', quote="", comment.char="", header=TRUE)

samples_for_removal <- ""
samples_for_inclusion <- as.character(read.table("~/results/ukb1189_qc.sample_list", header=TRUE)$sample)

notes_for_manny_1189 <-  get_hists_and_notes(hist_filename, tsv_data, log_file, outcome_info, codings_tables, qc_data, samples_for_removal, samples_for_inclusion)

write.table(notes_for_manny_1189, file=pheno_summary, col.names=TRUE, row.names=TRUE, sep='\t', quote=FALSE)

pheno_summary_PHESANT <- "~/results/ukb1189_phenosummary_all_final.tsv"
notes_for_manny_1189_PHESANT_codings_included <- include_PHESANT_reassignment_names(pheno_summary, outcome_info)
write.table(notes_for_manny_1189_PHESANT_codings_included, file=pheno_summary_PHESANT, col.names=TRUE, row.names=TRUE, sep='\t', quote=FALSE)

# Run the ICD10 codes in Joel's application through.
hist_filename <- "~/results/ukb1189_icd10_hist"
pheno_summary <- "~/results/ukb1189_icd10_phenosummary.tsv"

filename <- "~/results/ukb1189_icd10_flags.tsv"
log_file <- FALSE
tsv_data <- read.table(filename, header=TRUE, sep='\t')
names(tsv_data)[1] <- "userId"
names(tsv_data)[-1] <- paste("X41202", names(tsv_data)[-1], sep="_")

qc_data <- read.table("~/results/ukb1189_qc.tsv", header=TRUE)

outcome_info <- read.table("~/Repositories/PHESANT/variable-info/outcome_info_final.tsv",
					   sep='\t', quote="", comment.char="", header=TRUE)

samples_for_removal <- as.character(read.table("~/results/w1189_20170726_participantwithdrawallist.csv")$V1)

notes_for_manny_1189_icd10 <- get_hists_and_notes(hist_filename, tsv_data, log_file, outcome_info, codings_tables, qc_data, samples_for_removal, start_column=2)
write.table(notes_for_manny_1189_icd10, file=pheno_summary, col.names=TRUE, row.names=TRUE, sep='\t', quote=FALSE)

## ANDREA'S APPLICATION

# Running Andrea's application through.
# Can only do 'all', as we don't have the QC information and set of redacted/individuals that have removed consent.
hist_filename <- "~/results/ukb9438_hist"
pheno_summary <- "~/results/ukb9438_phenosummary.tsv"

filename <- "~/results/ukb9438_output"
tsv_filename <- paste(filename, ".tsv", sep="")
log_file <- paste(filename, ".log", sep="")
tsv_data <- read.table(tsv_filename, header=TRUE, sep='\t')

qc_data <- ""

outcome_info <- read.table("~/Repositories/PHESANT/variable-info/outcome_info_final.tsv",
					   sep='\t', quote="", comment.char="", header=TRUE)

samples_for_removal <- ""
samples_for_inclusion <- tsv_data$userId
check <- FALSE

notes_for_manny_9438 <-  get_hists_and_notes(hist_filename, tsv_data, log_file, outcome_info, codings_tables, qc_data, samples_for_removal, samples_for_inclusion, check)
write.table(notes_for_manny_9438, file=pheno_summary, col.names=TRUE, row.names=TRUE, sep='\t', quote=FALSE)

pheno_summary_PHESANT <- "~/results/ukb9438_phenosummary_final.tsv"
notes_for_manny_9438_PHESANT_codings_included <- include_PHESANT_reassignment_names(pheno_summary, outcome_info)
write.table(notes_for_manny_9438_PHESANT_codings_included, file=pheno_summary_PHESANT, col.names=TRUE, row.names=TRUE, sep='\t', quote=FALSE)

