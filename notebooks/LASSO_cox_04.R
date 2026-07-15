library(glmnet) # https://www.serdarbalci.com/jsurvival/articles/09-lassocox-comprehensive.html#example-1-lung-cancer-prognostic-model
library(survival) # https://www.statology.org/lasso-regression-in-r/, https://www.youtube.com/watch?v=5GZ5BHOugBQ
library(readr)

df <- read.csv("/home/oliver/oligodendroglioma_project/data/TCGA/rnaseq_survival_ready.csv")

X <- as.matrix(df[, !colnames(df) %in% c("sampleID", "PFI", "PFI.time")]) # Requires numeric values only for this part, Input matrix X
Y <- Surv(time = df$PFI.time, event = df$PFI) # this is where the 'clinical' data comes in (PFI, etc), response matrix Y https://glmnet.stanford.edu/articles/glmnet.html

set.seed(42) # Answer to the Universe...
# https://glmnet.stanford.edu/articles/glmnet.html
cv_fit <- cv.glmnet(X, Y, family = "cox", alpha = 1) # alpha=1 to run LASSO, family = "cox" for survival
# Cross Validation, automatically doing 10 CV (9 train 1 test loop)
plot(cv_fit)
coefs <- coef(cv_fit, s = "lambda.min") # We use lambda.min to extract the most amount of genes so we dont lose any good data, we can always change it later
# This is the penalizing ^
# I can skip all the lambda work by just selecting lambda.min and seeing how many genes 'make the cut', 
# I see that relatively few do, so this is perfect and I can skip all the lambda formatting.
# Further, I do not want to accidentally lose any genes so min is good.


active_coefficients <- as.matrix(coefs) # https://www.rdocumentation.org/packages/data.table/versions/1.18.4/topics/as.matrix
surviving_genes <- data.frame(
  Gene = rownames(active_coefficients),
  Coefficient = active_coefficients[, 1] # takes out only the first row and all the columns, forms a clean dataframe from a sparse matrix
)
final_signature <- surviving_genes[surviving_genes$Coefficient != 0, ] # filters out genes that were shrunk to 0 to not clog
print(final_signature)



# Onto basic cox in R
# https://metricgate.com/blogs/how-to-run-cox-regression-in-r/, https://www.youtube.com/watch?v=CUmyZms9tvQ for R

selected_gene_names <- final_signature$Gene

final_cox_data <- df[, c("PFI", "PFI.time", selected_gene_names)]

cox_formula <- as.formula(paste("Surv(PFI.time, PFI) ~", paste(selected_gene_names, collapse = " + ")))

final_cox_model <- coxph(cox_formula, data = final_cox_data)

print(summary(final_cox_model))
# Some gemini help for converting python to R.