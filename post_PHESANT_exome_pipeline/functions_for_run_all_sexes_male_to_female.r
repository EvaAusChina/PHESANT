library(data.table)
library(dplyr)

minimum_bin <- function(df_column, minimum=100) {
	
	if(any(df_column=="", na.rm=TRUE)){
		df_column <- df_column[-which(df_column=="")]
		if(length(df_column) == 0) return(FALSE)
	}
	
	if(sum(is.na(df_column)) == length(df_column)) return(FALSE)
	result <- table(df_column)

	if(length(result) == 1) {
		# This is the case when everyone is in a single bin.
		return(FALSE)
	} else {
		# Otherwise, ask the size of the smallest bin.
		return(min(result) > minimum)
	}
}

keep_cat_ordered <- function(df_column, minimum=5000) {
	
	if(any(df_column=="", na.rm=TRUE)){
		df_column <- df_column[-which(df_column=="")]
		if(length(df_column) == 0) return(FALSE)
	}

	if(sum(!is.na(df_column)) >= 5000) {
		return(TRUE)
	} else {
		return(FALSE)
	}
}
get_raw_cts(trybd$f.eid, trybd)
trybd[1:5,1:5]

load("G:/phenotype estimation/F__phenotype estimation_trybd.Rdata")
get_raw_cts <- function(cts_variable, file_header, file=paste0(QCed_input_file_and_output_name , ".tsv"), single=TRUE)
{	
	# Get matches in the header file.
	match <- grep(paste0('x', cts_variable, '_'), header)

	if(single) {
		if(length(match)==1) {
			return(match)
		} else {
			return(c())
		}
	} else {
		if(length(match)==1) {
			return(c())
		} else {
			return(match)
		}
	}
}
library(stringr)
colnames(trybd) <- str_extract(colnames(trybd),"[0-9]+")
head(colnames(trybd))
coding_info <- reassign_all_values(colnames(trybd)[-1], variable_info)
reassign_all_values <- function(cts_variables, phenofile) {
	where <- phenofile$FieldID %in% cts_variables
	changes <- cbind(phenofile$FieldID[where], phenofile$DATA_CODING[where])
	changes <- changes[!is.na(changes[,2]),]
	# create a list of variables for each coding
	codings <- list()
	for(i in names(table(changes[,2]))) {
		codings[[i]] <- changes[which(changes[,2]==strtoi(i)),1]
	}
	return(codings)
}

get_reassignment <- function(reassignment) {
	matrix(as.integer(unlist(strsplit(strsplit(reassignment, split='\\|')[[1]], split='='))), ncol=2, byrow=TRUE)
}

make_the_changes <- function(reassignment_matrix, cts_variables, data_frame)
{
	for(i in 1:nrow(reassignment_matrix)) {
		matches <- data_frame[cts_variables] == reassignment_matrix[i,1]
		if(length(matches) > 0)
			data_frame[cts_variables][matches] <- reassignment_matrix[i,2]
	}
	return(data_frame)
}

change_values <- function(codings, data_frame, coding_info_file) {
	# Need to read in the encoding and determine the changes to make
	reassignments <- fread(coding_info_file, sep=',', header=TRUE, data.table=FALSE)
	reassignments <- reassignments[reassignments$dataCode %in% names(codings),]
	reassignment_matrices <- sapply(reassignments$reassignments, get_reassignment)
	names(reassignment_matrices) <- reassignments$dataCode
	codings <- codings[order(names(codings))]
	reassignment_matrices <- reassignment_matrices[order(names(reassignment_matrices))]

	for(i in 1:length(reassignment_matrices)) {
		data_frame <- make_the_changes(reassignment_matrices[[i]], as.character(codings[[i]]), data_frame)
	}

	return(data_frame)
}


### take average for each variables and integrate them together
average_over_cts_multi <- function(cts_variable, data_frame)
{
	where <- grep(cts_variable, names(data_frame))
	column <- rowMeans(data_frame[,where], na.rm=TRUE)
	column[is.nan(column)] <- NA
	return(column)
}

irnt <- function(cts_variable) {
    set.seed(1234) # This is the same as was used by PHESANT - for checking.
    n_cts <- length(which(!is.na(cts_variable)))
    quantile_cts <- (rank(cts_variable, na.last = "keep", ties.method = "random") - 0.5) / n_cts
    # use the above to check, but also use frank for the real thing
    cts_IRNT <- qnorm(quantile_cts)	
    return(cts_IRNT)
}

look_for_logical <- function(column) {
	return("TRUE" %in% column | "FALSE" %in% column)
}