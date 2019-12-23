# Parse the arguments input by the user
# if argument 'test' is used then run test phenome scan
processArgs()

opt$phenofile <<-  paste(datadir,'phenotypes.csv', sep="")
opt$variablelistfile <- './testWAS/variable-lists/outcome-info.tsv'
opt$datacodingfile <<- '../testWAS/variable-lists/data-coding-ordinal-info.txt'
opt$traitofinterest <<- 'exposure'
opt$resDir <<- '../testWAS/results/'
opt$userId <<- 'userId'
processArgs <- function()
{
    if (opt$test==TRUE) {
        # Set up the test phenome scan settings
        datadir <- '../testWAS/data/'
        opt$phenofile <<-  paste(datadir,'phenotypes.csv', sep="")
        opt$variablelistfile <<- '../testWAS/variable-lists/outcome-info.tsv'
        opt$datacodingfile <<- '../testWAS/variable-lists/data-coding-ordinal-info.txt'
        opt$traitofinterest <<- 'exposure'
        opt$resDir <<- '../testWAS/results/'
        opt$userId <<- 'userId'

        processParts(opt$partIdx, opt$numParts)	
    } else {
        # Check arguments are supplied correctly
        if (is.null(opt$phenofile)){
            print_help(opt_parser)
            stop("phenofile argument must be supplied", call.=FALSE)
        } else if (!file.exists(opt$phenofile)) {
            stop(paste("phenotype data file phenofile=", opt$phenofile,
                " does not exist", sep=""), call.=FALSE)
        }

        if (is.null(opt$variablelistfile)){
            print_help(opt_parser)
            stop("variablelistfile argument must be supplied", call.=FALSE)
        } else if (!file.exists(opt$variablelistfile)) {
            stop(paste("variable listing file variablelistfile=",
                opt$variablelistfile, " does not exist", sep=""), call.=FALSE)
        }

        if (is.null(opt$datacodingfile)){
            print_help(opt_parser)
            stop("datacodingfile argument must be supplied", call.=FALSE)
        } else if (!file.exists(opt$datacodingfile)) {
            stop(paste("data coding file datacodingfile=",
                opt$datacodingfile, " does not exist", sep=""), call.=FALSE)
        }

        if (is.null(opt$resDir)){
            print_help(opt_parser)
            stop("resDir argument must be supplied", call.=FALSE)
        } else if (!file.exists(opt$resDir)) {
            stop(paste("results directory resDir=",
                opt$resDir, " does not exist", sep=""), call.=FALSE)
        }
        processParts(opt$partIdx, opt$numParts)
    }
    
    # Just some information to the user
    print("Adjusting for age and sex")
}


# Parse the 'part' arguments and check they are valid
processParts <- function(pIdx, nParts)
{
    if (is.null(pIdx) && is.null(nParts)) {
        opt$varTypeArg <<- "all"
        print(paste("Running with all traits in phenotype file:", opt$phenofile))
    } else if (is.null(pIdx)) {
        print_help(opt_parser)
        stop("pIdx argument must be supplied when nParts argument is supplied", call.=FALSE)
    } else if (is.null(nParts)) {
        print_help(opt_parser)
        stop("nParts argument must be supplied when pIdx argument is supplied", call.=FALSE)
    } else if (pIdx < 1 || pIdx > nParts) {
        print_help(opt_parser)
        stop("pIdx arguments must be between 1 and nParts inclusive", call.=FALSE)
    } else {
        opt$varTypeArg <<- paste(pIdx, "-", nParts, sep="")
        print(paste("Running with part",pIdx,"of",nParts," in phenotype file:", opt$phenofile))
    }
}
