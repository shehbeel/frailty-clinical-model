# Building ML Model for Frailty prediction
# Author: Shehbeel Arif
# NASA GeneLab Multi-Omics AWG

## LOAD LIBRARIES
# Library data manipulation
library(dplyr)
# Set seed
set.seed(1234)

################################################################################
## SET DIRECTORIES
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
analysis_dir <- file.path(root_dir, "analyses", "02-frailty-modeling")


# Set output directories
results_dir <- file.path(analysis_dir, "results")
plots_dir <- file.path(analysis_dir, "plots")

# Make output directories if they don't exist
if (!dir.exists(results_dir)) {
  dir.create(results_dir)
}
if (!dir.exists(plots_dir)) {
  dir.create(plots_dir)
}

# Set input directory
data_dir <- file.path(root_dir, "data")

# Declare input file paths
metadata_file <- file.path(root_dir, "data", "GSE144304_meta.txt")
data_file <- file.path(root_dir, "analyses", "01-convert-ensembl-ids", "results", "GSE144304_raw_counts_with_genes.csv")

#######
# Import metadata and data 
metadata <- readr::read_delim(metadata_file)
expression_df <- readr::read_csv(data_file)

#######
## PREPROCESS THE DATA
expression_df <- expression_df %>%
  distinct(gene_symbol, .keep_all = TRUE) %>% # Remove duplicate genes
  tibble::column_to_rownames("gene_symbol")

# Make the data in the order of the metadata
expression_df <- expression_df %>%
  dplyr::select(metadata$sample_id)

# Check if the expression and metadata are in the same order
all.equal(colnames(expression_df), metadata$sample_id)


# Make cluster a factor and set the levels appropriately
metadata <- metadata %>%
  dplyr::mutate(
    # Here we define the values our factor variable can have and their order.
    condition = factor(condition, levels = c("frail", "fit", "young"))
  )

# Define a minimum counts cutoff and filter the data to include
# only rows (genes) that have total counts above the cutoff
filtered_expression_df <- expression_df %>%
  dplyr::filter(rowSums(.) >= 10)

gene_matrix <- filtered_expression_df %>%
  t() %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "sample_id")

# Merge meta and gene matrix data
df <- merge(metadata, gene_matrix, by = "sample_id")
df <- df %>%
  select(-c("sample_id", "condition"))
#df$gender <- unclass(df$gender)
df$gender<-c(male=1,female=0)[df$gender]
df$frailty<-c(1=frail,female=normal)[df$frailty]

# Move condition column to as first column
df <- df %>%
  select(frailty, everything())

# Convert to matrix
#df <- df %>% as.matrix()
#df <- as.matrix(sapply(df, as.numeric))  

# look at the class of each column
sapply(df, class)

#df <- df %>% mutate_at(2:20299, as.numeric)
#df[,2:20]

################################################################################
# Build MSGL model
library(msgl)


fit_cv <- msgl::cv(as.matrix(sapply(df[,2:20299], as.numeric)), df$frailty, 
                   fold = 10,alpha = 0.5, lambda = 0.1, use_parallel = TRUE, 
                   standardize = FALSE)
fit_cv

fit <- msgl::fit(as.matrix(sapply(df[,2:20299], as.numeric)), df$frailty, 
                 alpha = 0.5, lambda = 1, standardize = FALSE)
fit

#find the best model index from cv
features(fit)[[best_model(fit_cv)]]
parameters(fit)[[best_model(fit_cv)]]
coef(fit, best_model(fit_cv))
coef = coef(fit, best_model(fit_cv))
coef = as.matrix(coef)
coef[which(coef == 0)] = NA
coef = data.frame(coef)
coef = coef[,-1]
coef$condition = rownames(coef)

library(reshape2)
coef = melt(coef)
coef$variable = as.character(coef$variable)

options(scipen=10000)
library(corrplot)
library(viridis)
library(ggplot2)
library(RColorBrewer)
p1 = ggplot(coef, aes(variable, condition, fill = value > 0)) + 
  geom_tile() + 
  geom_text(aes(label = formatC(value, format = "e", digits = 2))) + 
  scale_fill_discrete(labels = c("Coefficient < 0", "Coefficient > 0", "Not selected")) + 
  labs(x = "", y = "", fill = "") + 
  theme(text = element_text(family = "Times New Roman"),
        axis.text.x = element_text(angle = 30, vjust = 0.5, hjust=0.5, color = "grey20", size = 12, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 12, face = "plain"),  
        axis.title.x = element_text(color = "grey20", size = 12, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 12, face = "plain"),
        strip.text.x = element_text(size = 12))

p2 = ggplot(df, aes(x = condition, y = `ACTA1`, group = condition, color = condition)) + 
  geom_boxplot(outlier.shape = NA) + geom_jitter(shape=16, size = 6,position=position_jitter(0.2)) + 
  labs(x = "", y = "ACTA1", color = "") + theme_bw() +
  theme(legend.position = "right",
        text = element_text(family = "Times New Roman"),
        axis.text.x = element_text(color = "grey20", size = 12, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 12, face = "plain"),  
        axis.title.x = element_text(color = "grey20", size = 12, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 12, face = "plain")) + 
  stat_summary(fun=mean, geom="point", shape=20, size=6, color="red", fill="red") + 
  scale_color_manual(values = viridis(6))

p3 = ggplot(df, aes(x = condition, y = `MYH2`, group = condition, color = condition)) + 
  geom_boxplot(outlier.shape = NA) + geom_jitter(shape=16, size = 6,position=position_jitter(0.2)) + 
  labs(x = "", y = "MYH2", color = "") + theme_bw() +
  theme(legend.position = "right",
        text = element_text(family = "Times New Roman"),
        axis.text.x = element_text(color = "grey20", size = 12, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 12, face = "plain"),  
        axis.title.x = element_text(color = "grey20", size = 12, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 12, face = "plain")) + 
  stat_summary(fun=mean, geom="point", shape=20, size=6, color="red", fill="red") + 
  scale_color_manual(values = viridis(6))

library(ggpubr)
library(grid)
ggarrange(p1,
          labels = c("A"))
ggarrange(p2, p3,
          labels = c("B","C"))
####
require(nnet)
multinom_fit <- multinom(condition ~ `gender` + `ACTA1` + `MYH2`, data = df)
summary(multinom_fit)
round((1 - pnorm(abs(summary(multinom_fit)$coefficients/summary(multinom_fit)$standard.errors), 0, 1)) * 2, digits = 3)
exp(coef(multinom_fit))

head(pp <- fitted(multinom_fit))
library(sjPlot)
tab_model(multinom_fit, digits = 6)

prediction = predict(multinom_fit, newdata = df, "class")
sum(prediction == df$condition)/nrow(df)


result1 = predict(multinom_fit, newdata = df, type='probs')
library(HandTill2001)
auc(multcap(
  response = factor(df$condition),
  predicted = as.matrix(result1)
))



##################
# Other Genes
#################
require(nnet)
multinom_fit <- multinom(frailty ~ `gender` + `RNF115` + `STRADA` + `NCOA4` + 
                           `FXYD1` + `ILF3-DT` + `ZNF285` + `MAGIX` + `TMEM185A` +
                           `EGLN2` + `CAHM` + `TAF15` + `KMT2B` + `POM121C` +
                           `ATP6V1FNB` + `PAXIP1-DT` + `ZNHIT3` + `CCL14` + 
                           `GGNBP2` + `EBLN3P` + `LOC101929130`, data = df, model = TRUE)
summary(multinom_fit)
round((1 - pnorm(abs(summary(multinom_fit)$coefficients/summary(multinom_fit)$standard.errors), 0, 1)) * 2, digits = 3)
exp(coef(multinom_fit))

head(pp <- fitted(multinom_fit))
library(sjPlot)
tab_model(multinom_fit, digits = 6)

prediction = predict(multinom_fit, newdata = df, "class")
sum(prediction == df$frailty)/nrow(df)


result1 = predict(multinom_fit, newdata = df)
library(HandTill2001)
auc(multcap(
  response = factor(df$frailty),
  predicted = as.matrix(result1)
))

# Load ROC-AUC library
library(pROC)
# create roc curve
roc_object <- roc(df$frailty, as.numeric(prediction), plot=TRUE, 
                  legacy.axes=TRUE, percent=TRUE, xlab="False Positive (%)",
                  ylab="True Positive (%)", col="#377eb8", lwd=4, print.auc=TRUE,
                  auc.polygon=TRUE)

# calculate area under curve
auc(roc_object)

plot.roc(df$frailty, as.numeric(prediction), percent=TRUE, col="#4daf4a", lwd=4,
         print.auc=TRUE, add=TRUE, print.auc.y=40)

# Goodness-of-fit Test
# https://peopleanalytics-regression-book.org/multinomial-logistic-regression-for-nominal-category-outcomes.html
DescTools::PseudoR2(multinom_fit, 
                    which = c("McFadden", "CoxSnell", "Nagelkerke"))
# The most approachable method to assess model confidence is the Hosmer-Lemeshow test


