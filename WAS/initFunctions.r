# Load the required R files
loadSource <- function() {
    source("WAS/loadData.r")
    source("WAS/reassignValue.r")
    source("WAS/validateInput.r")
    source("WAS/testNumExamples.r")
    source("WAS/binaryLogisticRegression.r")
    source("WAS/equalSizedBins.r")
    source("WAS/fixOddFieldsToCatMul.r")
    source("WAS/replaceMissingandNaN.r")
    source("WAS/testAssociations.r")
    source("WAS/testCatMultiple.r")
    source("WAS/testCatSingle.r")
    source("WAS/testContinuous.r")
    source("WAS/testInteger.r")
    source("WAS/testCategoricalOrdered.r")
    source("WAS/testCategoricalUnordered.r")
    source("WAS/counts.r")
}

# Init the counters used to determine how many variables took each path 
# in the variable processing flow.
initCounters <- function() {
	counters <- data.frame(name=character(),
						   countValue=integer(),
                           stringsAsFactors=FALSE)
	return(counters)
}

# Load the variable information and data code information files
initVariableLists <- function()
{
	phenoInfo <- read.table(opt$variablelistfile, sep="\t", header=1, comment.char="", quote="")
	dataCodeInfo <- read.table(opt$datacodingfile,sep=",", header=1)
	vars <- list(phenoInfo=phenoInfo, dataCodeInfo=dataCodeInfo)
	return(vars)
}
