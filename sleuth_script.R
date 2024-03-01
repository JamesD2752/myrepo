library(data.table)
library(dplyr)
library(sleuth)

# Reads the data
stab <- read.table("sample_to_covariates.txt", header = TRUE, stringsAsFactors = FALSE, sep = '\t')

# Checks if 'sample' and 'path' columns exist in your data, if not, it creates them
if (!("sample" %in% colnames(stab))) {
  stab$sample <- 1:nrow(stab)
}

if (!("path" %in% colnames(stab))) {
  stab$path <- 1:nrow(stab)
}

# Assigns column names to stab dataframe
colnames(stab) <- c("sample", "path", colnames(stab)[3:length(colnames(stab))])

# Prints head of stab dataframe to check if data loading and column renaming if it is successful
print(head(stab))

# Now, you can proceed with sleuth_prep
so <- sleuth_prep(stab, ~sample)

# Check if sleuth_prep is successful
if (is.null(so)) {
  stop("Error: Sleuth preparation failed.")
}

# a differential expression analysis
# fits a model comparing conditions
so <- sleuth_fit(so, ~condition, 'full')

# Checks if sleuth_fit is successful
if (is.null(so)) {
  stop("Error: Sleuth fitting failed.")
}

# a reduced model
so <- sleuth_fit(so, ~1, 'reduced')

# Checks if sleuth_fit for reduced model is successful
if (is.null(so)) {
  stop("Error: Sleuth fitting for reduced model failed.")
}

# likelihood ratio test for differential expression
so <- sleuth_lrt(so, 'reduced', 'full')

# Checks if likelihood ratio test is successful
if (is.null(so)) {
  stop("Error: Likelihood ratio test failed.")
}

# extracts results
sleuth_table <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)

# Checks if sleuth_results is successful
if (is.null(sleuth_table)) {
  stop("Error: Failed to extract results.")
}

# Prints head of sleuth_table to verify if it's created successfully
print(head(sleuth_table))

# filters significant results
sleuth_significant <- filter(sleuth_table, qval <= 0.05) %>% arrange(pval)
sig_sleuth <- sleuth_significant %>% select(target_id, test_stat, pval, qval)

# Prints head of sig_sleuth to verify if significant results are filtered correctly
print(head(sig_sleuth))

# writes target id, test stat, pval and qval for significant transcript
# includes header, tab-delimit
write.table(sig_sleuth, file = "topten.txt", quote = FALSE, row.names = FALSE)
