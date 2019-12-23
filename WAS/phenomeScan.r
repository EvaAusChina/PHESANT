# The MIT License (MIT)
# Copyright (c) 2017 Louise AC Millard, MRC Integrative Epidemiology Unit, University of Bristol
#
# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated 
# documentation files (the "Software"), to deal in the Software without restriction, including without 
# limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of 
# the Software, and to permit persons to whom the Software is furnished to do so, subject to the following 
# conditions:
#
# The above copyright notice and this permission notice shall be included in all copies or substantial portions 
# of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED 
# TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL 
# THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF 
# CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
# DEALINGS IN THE SOFTWARE.

library("optparse")
library(readr)
library("data.table")

option_list <- list(
  make_option(c("-f", "--phenofile"), type="character", default=NULL,
    help="phenotype dataset file name", metavar="character"),
  make_option(c("-v", "--variablelistfile"), type="character", default=NULL,
    help="variablelistfile file name (should be tab separated)", metavar="character"),
  make_option(c("-d", "--datacodingfile"), type="character", default=NULL,
    help="datacodingfile file name (should be comma separated)", metavar="character"),
  make_option(c("-r", "--resDir"), type="character", default=NULL,
    help="resDir option should specify directory where results files should be stored", metavar="character"),
  make_option(c("-u", "--userId"), type="character", default="userId",
    help=paste("userId option should specify user ID column in trait of interest and phenotype files [default = %default]"), 
    metavar="character"),
  make_option(c("-t", "--test"), action="store_true", default=FALSE,
    help="run test phenome scan on test data (see test subfolder) [default = %default]"),
  make_option(c("-a", "--partIdx"), type="integer", default=NULL,
    help="part index of phenotype (used to parallelise)"),
  make_option(c("-b", "--numParts"), type="integer", default=NULL,
    help="number of phenotype parts (used to parellelise)"),
  # make_option(c("-l", "--log"), type="character", default='log',
  #   help="name of the logfile [default = %default]"),
  make_option(c("-o", "--out"), type="character", default='output',
    help="name of the output .tsv file containing the parsed columns in the provided phenofile [default = %default]"),
  make_option(c("-c", "--catmultcutoff"), type="integer", default=50,
    help="The cutoff for exclusion when creating dichotomous variables for CAT-MULTIPLE."),
  make_option(c("-z", "--catordnacutoff"), type="integer", default=5000,
    help="The cutoff for exclusion for number of non-NAs in ordered categorical variables."),
  make_option(c("-w", "--catunordnacutoff"), type="integer", default=5000,
    help="The cutoff for exclusion for number of non-NAs in unordered categorical variables."),
  make_option(c("-x", "--contnacutoff"), type="integer", default=5000,
    help="The cutoff for exclusion for number of non-NAs in continuous variables."),
  make_option(c("-n", "--binnacutoff"), type="integer", default=5000,
    help="The cutoff for exclusion for number of non-NAs in binary-variables."),
  make_option(c("-i", "--bintruecutoff"), type="integer", default=100,
    help="The cutoff for exclusion for numbers of members of a category in binary-variables."),
  make_option(c("-m", "--mincategorysize"), type="integer", default=10,
    help="The minimum number of samples in a category for categorical single, integer, and continous variables."),
  make_option(c("-e", "--maxunorderedcategories"), type="integer", default=1000,
    help="The maximum number of categories in an unordered categorical variable."),
  make_option(c("-g", "--propforcontinuous"), type="double", default=0.2,
    help="The cutoff for proportion of samples with the same value for the variable to not be considered continuous."),
  make_option(c("-j", "--inttocontcutoff"), type="integer", default=10,
    help="The cutoff for the number of distinct integer values to send an integer variable to a continuous variable.")
)


load("G:/phenotype estimation/F__phenotype estimation_trybd.Rdata")
write.csv(trybd, file = "WAS/phenotypes.csv")

variablelistfile <- fread("variable-info/outcome_info_final_pharma_nov2019.tsv")
datacodingfile <- read.table("variable-info/data-coding-ordinal-info-nov2019-update.txt")
testDir = "../testWAS/"
phenofile="${testDir}data/phenotypes.csv"


datadir <- '../WAS/'
opt$phenofile <<-  paste(datadir,'phenotypes.csv', sep="")
opt$variablelistfile <- './testWAS/variable-lists/outcome-info.tsv'
opt$datacodingfile <<- '../testWAS/variable-lists/data-coding-ordinal-info.txt'
opt$traitofinterest <<- 'exposure'
opt$resDir <<- '../testWAS/results/'
opt$userId <<- 'userId'

opt$variablelistfile <- variablelistfile
opt$phenofile <- phenofile
opt$userId <- "userId"
opt$datacodingfile <- datacodingfile
opt$resDir =  'testWAS/results'

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

varlogfile <- paste(opt$resDir, '/', opt$out, '.', opt$partIdx, '.log', sep='')
outputfile <- paste(opt$resDir, '/', opt$out, '.', opt$partIdx, '.tsv', sep='')

source("WAS/processArgs.r")
opt <<- opt
print(opt)
processArgs(phenofile)

source("WAS/initFunctions.r")
loadSource()

# Load the files we write to and use
counters <- initCounters()
vl <- initVariableLists()
cat("LOADING...\n")

# Load data
d <- loadData()
data <- d$datax
confounders <- d$confounders
numPreceedingCols <- ncol(confounders)+1
phenoStartIdx <- numPreceedingCols+1
print(paste0('pheno start Idx: ', phenoStartIdx))
cat("LOADING DONE.\n")

phenoVars <- colnames(data)

# First and second columns are the id and snpScore, respectively, 
# as determined in loadData.r
phenoVars <- phenoVars[-c(1,2)] 

# This decides on the start and end idxs of phenotypes that we test,
# so that we can parallelise into multiple jobs
if (opt$varTypeArg != "all") {
    
    partSize <- ceiling(length(phenoVars)/opt$numParts)
    partStart <- (opt$partIdx-1)*partSize + 1

    if (opt$partIdx == opt$numParts) {
        partEnd <- length(phenoVars)
    } else {
        partEnd <- partStart + partSize - 1
    }

} else {
    partStart <- 1
    partEnd <- length(phenoVars)
}

cat(partStart, '-', partEnd, '\n')

currentVar <- ""
currentVarShort <- ""
first <- TRUE

# Zero because then the idx is the position of the previous variable, i.e. the var in currentVar.
phenoIdx <- 0 
i <- 1

data_to_store <- matrix(nrow = nrow(data), ncol = 0)
data_to_store <- as.data.frame(data_to_store)
data_to_store_var <- c()

for (var in phenoVars) { 

    varx <- gsub("^x", "", var)
    varx <- gsub("_[0-9]+$", "", varx)
    varxShort <- gsub("^x", "", var)
    varxShort <- gsub("_[0-9]+_[0-9]+$", "", varxShort)

    # Test this variable
    cat("The current variable is", currentVar, '\n')
    cat(var, '\n')

    if (currentVar == varx) {
        thisCol <- data[,eval(var)]
        thisCol <- replaceNaN(thisCol)
        currentVarValues <- cbind.data.frame(currentVarValues, thisCol)
    } else if (currentVarShort == varxShort) {
    	# Different time point of this variable so skip
    } else {
        # New variable so run test for previous (we have collected all the columns now)
        if (first == FALSE) {
            thisdata <- cbind.data.frame(data$geno, confounders, currentVarValues)
            colnames(thisdata)[1] <- "geno"
            # Only start new variable processing if last column of it is within 
            # the idx range for this part
            if (phenoIdx >= partStart && phenoIdx <= partEnd) {
                next_test <- testAssociations(currentVar, currentVarShort, thisdata, varlogfile)

                if (is.null(next_test) == FALSE) {
                    data_to_store <- cbind(data_to_store, next_test[[1]])
                    data_to_store_var <- c(data_to_store_var, next_test[[2]])
                }
            }
        }

        first <- FALSE
        # New variable so set values
        currentVar <- varx
        currentVarShort <- varxShort
        currentVarValues <- data[,eval(var)]
        currentVarValues <- replaceNaN(currentVarValues)
    }
    phenoIdx <- phenoIdx + 1
    i <- i+1
}

# Last variable so test association
thisdata <- cbind.data.frame(data$geno, confounders, currentVarValues)
colnames(thisdata)[1] <- "geno"

if (phenoIdx >= partStart && phenoIdx <= partEnd) {
    next_test <- testAssociations(currentVar, currentVarShort, thisdata, varlogfile)
    if (is.null(next_test) == FALSE) {
        data_to_store <- cbind(data_to_store, next_test[[1]])
        data_to_store_var <- c(data_to_store_var, next_test[[2]])
    }
}

colnames(data_to_store) <- data_to_store_var
data_to_store <- cbind.data.frame(data[opt$userId], confounders, data_to_store)
colnames(data_to_store)[1] <- "userId"

fwrite(data_to_store, sep='\t', quote=TRUE, row.names=FALSE, file=outputfile)

# Save counters of each path in variable flow
saveCounts()
warnings()
