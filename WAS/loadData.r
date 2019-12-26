# Loads phenotype and trait of interest data files
# Creates phenotype / trait of interest data frame
# Creates confounder data frame
# Returns an object holding these two data frames
loadData <- function()
{
    # Load phenotype
    cat("Loading phenotypes...\n")

    if (length(grep('.tsv', opt$phenofile)) == 1 | length(grep('.tab', opt$phenofile)) == 1) {
        phenotype <- fread(opt$phenofile, header=TRUE, sep='\t', data.table=FALSE)
    } else if (length(grep('.csv', opt$phenofile)) == 1) {
            phenotype <- fread(opt$phenofile, header=TRUE, sep=',', data.table=FALSE)
    } else if (length(grep('.Rdata', opt$phenofile)) == 1) {
        load(opt$phenofile)
        phenotype <- phenotypes
        rm("phenotypes")
        if (!exists("phenotype")) {
            stop("Error: phenotype not found in .Rdata file")
        }

    } else {
        stop("Cannot detect the format of the phenotype file")
    }
    if (ncol(phenotype) == 1)
        stop("Number of columns of the phenotype file is 1: is it tab separated?")

    validatePhenotypeInput(phenotype)

    # Load SNPs.
    # cat("Loading trait of interest file...\n")
    # snpScores <- read.table(opt$traitofinterestfile, sep=",", header=1)
    # validateTraitInput(snpScores)

    # # Keep only the userID and exposure variable
    # idx2 <- which(names(snpScores) == opt$traitofinterest)

    # # Remove all rows with no trait of interest
    # idxNotEmpty <- which(!is.na(snpScores[,idx2]))
    # cat("Trait of interest has", nrow(snpScores), "rows with",
    #     length(idxNotEmpty), "not NA.")
    # snpScores <- snpScores[idxNotEmpty,]

    # DEV: Hack - set the userId in the trait of interest file to be 
    # the same as that in the phenotype file.
    snpScores <- cbind.data.frame(phenotype[opt$userId], 0)#snpScores[,idx2])
    colnames(snpScores)[1] <- opt$userId
    colnames(snpScores)[2] <- "geno"

    print("Merging trait of interest and phenotype data")
    # Merge to one matrix
    datax <- merge(snpScores, phenotype, by=opt$userId, all=FALSE)
    print("checking size of datax:")
    print(dim(datax))

    # datax <- phenotype[,-1]
    # 
    
    if (nrow(datax)==0)
    	stop("No examples with row in both trait of interest and phenotype files", call.=FALSE)

    datax <- fixOddFieldsToCatMul(datax)
    print("checking size of datax:")
    print(dim(datax))
    
    # Variables we adjust for.
    age <- datax[,"x21022_0_0"]
    sex <- datax[,"x31_0_0"]
    confounders <- cbind.data.frame(age,sex)

    d <- list(datax=datax, confounders=confounders)

    return(d)

}
