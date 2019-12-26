# Performs binary logistic regression on the phenotype stored in this data 
# and stores result in 'results-logistic-binary' results file.
binaryLogisticRegression <- function(varName, varType, thisdata, varlogfile)
{   
    to_save <- thisdata[, phenoStartIdx]
    phenoFactor <- factor(thisdata[, phenoStartIdx])
    facLevels <- levels(phenoFactor)

    # Assert variable has exactly two distinct values
    if (length(facLevels) != 2) {
        cat("BINARY-NOT2LEVELS- (", length(facLevels), ") || ",
            sep="", file=varlogfile, append=TRUE)
        incrementCounter("binary.nottwolevels")
    }

    idxTrue <- which(phenoFactor == facLevels[2])
    idxFalse <- which(phenoFactor == facLevels[1])
    numNotNA <- length(which(!is.na(phenoFactor)))

    # Added to ensure that TRUE/FALSE is returned for these binary variables.
    phenoFactor_tmp <- rep.int(NA, length(phenoFactor))
    phenoFactor_tmp[idxTrue] <- TRUE
    phenoFactor_tmp[idxFalse] <- FALSE
    phenoFactor <- factor(phenoFactor_tmp)

    idxTrue <- length(idxTrue)
    idxFalse <- length(idxFalse)

    if (idxTrue < opt$bintruecutoff || idxFalse < opt$bintruecutoff) {
        cat("BINARY-LOGISTIC-SKIP-", opt$bintruecutoff, " (", idxTrue, "/", idxFalse, ") || ",
            sep="", file=varlogfile, append=TRUE)
        incrementCounter(paste("binary.", opt$bintruecutoff, sep=""))
        return(NULL)
    } else if (numNotNA < 500) {    
        cat("BINARY-LOGISTIC-SKIP-", opt$binnacutoff, " (", numNotNA, ") || ",
            sep="", file=varlogfile, append=TRUE)
        incrementCounter(paste("binary.", opt$binnacutoff, sep=""))
        return(NULL)
    } else {
        cat("sample ", idxTrue, "/", idxFalse, "(", numNotNA, ") || ",
            sep="", file=varlogfile, append=TRUE)

        tryCatch(
        { 
            incrementCounter("success.binary")
            return(list(phenoFactor, varName))
            # END TRYCATCH
        }, error = function(e) {
            print(paste("ERROR:", varName, gsub("[\r\n]", "", e)))
            incrementCounter("binary.error")
            return(NULL)
        })
    }
}
