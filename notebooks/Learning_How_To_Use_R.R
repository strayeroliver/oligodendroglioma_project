# Test R- learning how to use it

library(glmnet) # https://www.serdarbalci.com/jsurvival/articles/09-lassocox-comprehensive.html#example-1-lung-cancer-prognostic-model
library(survival) # https://www.statology.org/lasso-regression-in-r/, https://www.youtube.com/watch?v=5GZ5BHOugBQ
library(readr)

# This is similar to pd.read_csv essentially just creates a basic dataframe
df <- read.csv("/home/oliver/oligodendroglioma_project/data/TCGA/rnaseq_survival_ready.csv")

# Here, from the stanford resource:
# X stands for the input matrix, where its purely numbers
# Y stands for the response vector Y with the results that come from the specific genes.
X <- as.matrix(df[, !colnames(df) %in% c("sampleID", "PFI", "PFI.time")]) # Requires numeric values only for this part, Input matrix X
Y <- Surv(time = df$PFI.time, event = df$PFI) # this is where the 'clinical' data comes in (PFI, etc), response matrix Y https://glmnet.stanford.edu/articles/glmnet.html

fit <- glmnet(X, Y, family = "cox")
# This fits a basic model in glmnet, which is essentially visualizing what LASSO does.
# It tests a wide spread of lambdas and shows you which variables get minimized and by how much.
# This helps me visualize what glmnet and lambdas are actually doing to my data:

plot(fit)
