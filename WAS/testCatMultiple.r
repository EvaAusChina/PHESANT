# Performs preprocessing of categorical (multiple) fields, namely:
# 1) Reassigning values as specified in data coding file
# 2) Generating binary variable for each category in field, restricting to correct set of participants as specified
# in CAT_MULT_INDICATOR_FIELDS field of variable info file (either NO_NAN, ALL or a field ID)
# 3) Checking derived variable has at least catmultcutoff cases in each group
# 4) Calling binaryLogisticRegression function for this derived binary variable
variablelistfile <- fread("variable-info/outcome_info_final_pharma_nov2019.tsv")
datacodingfile <- read.table("variable-info/data-coding-ordinal-info-nov2019-update.txt")
load("G:/phenotype estimation/F__phenotype estimation_trybd.Rdata")
colnames(trybd) <- gsub("[.]", "_", colnames(trybd))
colnames(trybd) <- gsub("^f_", "x", colnames(trybd))
colnames(trybd)[1] <- "userId"


testCategoricalMultiple(colnames(trybd)[-1], datacodingfile$V1,trybd , variablelistfile)
head(variablelistfile)
testCategoricalMultiple <- function(varName, varType, thisdata, varlogfile)
{
    data_to_add_mat <- matrix(nrow=nrow(thisdata), ncol=0)
    data_to_add_names <- c()

    cat("CAT-MULTIPLE || ", file=varlogfile, append=TRUE)

    pheno <- thisdata[,phenoStartIdx:ncol(thisdata), drop=FALSE]
    pheno <- reassignValue(pheno, varName, varlogfile)

    # Get unique values from all columns of this variable
    # uniqueValues <- unique(na.omit(c(pheno)))
    numRows <- nrow(pheno)
    numCols <- ncol(pheno)

    uniqueValues = unique(na.omit(pheno[,1]));
    numCols = ncol(pheno);
    numRows = nrow(pheno);
    if (numCols>1) {
        for (num in 2:numCols) {
            u = unique(na.omit(pheno[,num]))
            uniqueValues = union(uniqueValues,u);
        }
    }

    print(uniqueValues)

    # For each value create a binary variable and test this
    if (length(uniqueValues) > 1) {
        for (variableVal in uniqueValues) {
            # Numeric negative values we assume are missing.
            if (is.numeric(variableVal) & variableVal < 0) {
                cat("SKIP_val:", variableVal," < 0 || ", sep="",
                    file=varlogfile, append=TRUE)
                next
            }

             # Make variable for this value - make sure that this is just the sample (the rows).
            idxForVar <- which(pheno == variableVal, arr.ind=TRUE)[,1]

            cat("CAT-MUL-BINARY-VAR ", variableVal, " || ", sep="", file=varlogfile, append=TRUE)
            incrementCounter("catMul.binary")

            # Make zero vector and set 1s for those with this variable value
            varBinary <- rep.int(FALSE, numRows)
            varBinary[idxForVar] <- TRUE
            varBinaryFactor <- factor(varBinary)

            # Data for this new binary variable
            newthisdata <- cbind.data.frame(thisdata[,1:numPreceedingCols], varBinaryFactor)
            newthisdata_to_save <- cbind.data.frame(thisdata[,1:numPreceedingCols], varBinary)
            
            # One of 3 ways to decide which examples are negative
            idxsToRemove <- restrictSample(varName, pheno, variableVal, varlogfile)
            
            # Create an ids to save vector, and plug these in NAs at the right positions. 
        	if (!is.null(idxsToRemove)) {
                newthisdata_to_save <- newthisdata
                newthisdata_to_save[idxsToRemove,] <- NA
                newthisdata <- newthisdata[-idxsToRemove,]
            }

            facLevels <- levels(newthisdata_to_save[,phenoStartIdx])
            idxTrue <- length(which(newthisdata_to_save[,phenoStartIdx] == TRUE))
            idxFalse <- length(which(newthisdata_to_save[,phenoStartIdx] == FALSE))
                    
            if (idxTrue < opt$catmultcutoff || idxFalse < opt$catmultcutoff) {
                cat("CAT-MULT-SKIP-", opt$catmultcutoff, " (", idxTrue, " vs ", idxFalse, ") || ",
                    sep="", file=varlogfile, append=TRUE)
                incrementCounter(paste("catMul.", opt$catmultcutoff, sep=""))
            } else {
                incrementCounter(paste("catMul.over", opt$catmultcutoff, sep=""))
            	# Binary so logistic regression
                data_to_add <- binaryLogisticRegression(paste(varName, variableVal, sep="_"),
                    varType, newthisdata_to_save, varlogfile)
                data_to_add_mat <- cbind(data_to_add_mat, as.logical(data_to_add[[1]]))
                data_to_add_names <- c(data_to_add_names, data_to_add[[2]])
        	}
        }
    }
    return(list(data_to_add_mat, data_to_add_names))
}

# Restricts sample based on value in CAT_MULT_INDICATOR_FIELDS column of variable info file,
# either NO_NAN, ALL or a field ID.
# Returns idx's that should be removed from the sample.
restrictSample <- function(varName, pheno, variableVal, varlogfile)
{
    # Get definition for sample for this variable either NO_NAN, ALL or a variable ID
    varIndicator <- vl$phenoInfo$CAT_MULT_INDICATOR_FIELDS[which(vl$phenoInfo$FieldID == varName)]
    return(restrictSample2(varName, pheno, varIndicator, variableVal, varlogfile))
}

restrictSample2 <- function(varName,pheno, varIndicator, variableVal, varlogfile)
{
    if (varIndicator=="NO_NAN") { 
        # Remove NAs (remove all people with no value for this variable)

        # Row indexes with NA in all columns of this cat mult field		
        ind <- apply(pheno, 1, function(x) all(is.na(x)))
        naIdxs <- which(ind == TRUE)
        cat("NO_NAN Remove NA participants", length(naIdxs), "|| ",
            file=varlogfile, append=TRUE)
    } else if (varIndicator == "ALL") {
        # Use all people (no missing assumed) so return empty vector
        # e.g. hospital data and death registry
        naIdxs <- c()
        cat("ALL || ", file=varlogfile, append=TRUE)
    } else if (varIndicator != "") {
        # Remove people who have no value for indicator variable
        print(varIndicator)
        indName <- grep(paste0("x", varIndicator,"_0_"), names(data), value=TRUE)
        print(indName)
        cat("Indicator name ", indName, " || ", sep="", file=varlogfile, append=TRUE)
        indicatorVar <- data[,indName, drop=FALSE]

        # Remove participants with all NAs in this related field
        where_change <- which(sapply(indicatorVar, is.nan) | indicatorVar == "NaN" | indicatorVar == "", arr.ind=TRUE)
        
        if(nrow(where_change) > 0) {
            indicatorVar[which(sapply(indicatorVar, is.nan) | indicatorVar == "NaN" | indicatorVar == "", arr.ind=TRUE)] <- NA
        }

        naIdxs <- which(apply(indicatorVar, 1, allNAs))

        cat("Remove indicator var NAs:", length(naIdxs), "|| ", file=varlogfile, append=TRUE)

        if (is.numeric(as.matrix(indicatorVar))) {
            # Remove participants with value <0 in this related field - assumed missing indicators
            lessZero <- which(indicatorVar < 0, arr.ind=TRUE)[,1]
            naIdxs <- union(naIdxs, lessZero)
            cat("Remove indicator var <0:", length(lessZero), "|| ", file=varlogfile, append=TRUE)
        }
    } else {
        stop("Categorical multiples variables need a value for CAT_MULT_INDICATOR_FIELDS", call.=FALSE)
    }

    # Remove people with pheno < 0 if they aren't a positive example for this 
    # variable indicator because we can't know if they are a negative example or not.
    if (is.numeric(as.matrix(pheno))) {

        idxForVar <- which(pheno == variableVal, arr.ind = TRUE)[,1]
        idxMissing <- which(pheno < 0, arr.ind = TRUE)[,1]

        # All people with < 0 value and not variableVal
        naMissing <- setdiff(idxMissing, idxForVar)
        # Add these people with unknowns to set to remove from sample
        naIdxs <- union(naIdxs, naMissing)

        cat("Removed", length(naMissing) , "examples !=",
            variableVal, "but with missing value (<0) || ",
            file=varlogfile, append=TRUE)
    } else {
        cat("Not numeric || ", file=varlogfile, append=TRUE)
    }
    return(naIdxs)
}
