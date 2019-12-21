# 1. Run PHESANT on the cloud using neale_lab_parsed_and_restricted_to_QCed_samples.tsv as the phenotype file, see parsing_raw_pheno_initiation_script.

# 2. Run 01_run_summarise_phenotypes_initial_PHESANT_run.r to get the summary file for the PHESANT run on the cloud (after downloading the resultant .tsv and .log files to the current directory)

# 3. Run 02_run_all_sexes_to_male_female.r. This restricts the PHESANT run to males and females
# and remove phenotypes with too few phenotypes. Note that this creates IRNT - here a new renormalised version for each sex. To get the pheno on the same scale as both_sexes, you need to restrict those phenotypes to men/women.

# 4. Run 03_run_summarise_phenotypes_cloud.r. This then generates summary files for the sex-specific phenotype files that have been parsed by PHESANT, and restricts the summary files to the sex-specific.

# 5. Run 04_run_combine_PHESANT_output.r to combine the categorical (1-4), the IRNT cts, and the raw cts variables.

print("look! I am here")