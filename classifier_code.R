#Code for classifiers for "Distinct Host and Microbial Biology Distinguishes Lower Respiratory Tract Infection from Incidental Pathogen Carriage"
#Includes code for generating the relevant panels of Figures 5

#loading necessary packages
library(dplyr)
library(limma)
library(ggplot2)
library(ggpubr)
library(glmnet)
library(pROC)
library(DESeq2)
library(tidyr)
library(ggbeeswarm)

#importing metadata and gene counts files
metadata <- read.csv("/Users/emilylydon/Library/CloudStorage/Box-Box/VAP_IvC_project/Code:Source Data File/metadata.csv", check.names = FALSE)
counts <- read.csv("/Users/emilylydon/Library/CloudStorage/Box-Box/VAP_IvC_project/Code:Source Data File/gene_counts.csv", check.names = FALSE, row.names=1)

#extracting key of ENSEMBL IDs/gene symbols
gene_symbols <- counts[, "gene_symbol", drop = FALSE] #extracting gene names for later
counts <- counts[, -1] #getting rid of gene symbol to have a clean gene counts table

metadata$group <- factor(
  metadata$group,
  levels = c("LRTI", "IPC", "CTRL")
)

#VST transformation of gene counts
dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts, 
                                      colData = metadata, 
                                      design = ~1)
vsd <- DESeq2::varianceStabilizingTransformation(dds)
counts_VST <- SummarizedExperiment::assay(vsd) %>%
  round(.,digits=2)
counts_VST <- as.data.frame(counts_VST)

#extracting FABP4 (the most DE gene in the prior DE analyses) and adding to the metadata file
FABP4.id <- rownames(gene.symbol %>% subset(gene_symbol=="FABP4")) 
metadata$FABP4 <- as.numeric(t(counts_VST[FABP4.id, , drop = FALSE]))

#generating folds for cross-validation that are consistent across all classifiers
set.seed(123)
min.IPC <- 13 # minimum number of IPC sample per fold (70 IPC across 5 folds)
min.CTRL <- 9 # minimum number of ctrls sample per fold (49 CTRL across 5 folds)
while (TRUE) {
  cv.folds <- metadata %>%
    dplyr::select(Patient, group) %>%
    mutate(fold=sample(rep(1:5, length.out=nrow(.)))) #randomly generates the folds
  
  cv.folds.table <- cv.folds %>%
    group_by(fold) %>%
    dplyr::count(group) #counts the number per group
  
  print("Generated CV folds with following counts:")
  print(cv.folds.table)
  
  if ((min(cv.folds.table[cv.folds.table$group=="IPC","n"]) < min.IPC) | min(cv.folds.table[cv.folds.table$group=="CTRL", "n"]) < min.ctrl) {
    print("At least one fold has too few IPC or CTRL samples. Regenerating CV folds...")
  } else {
    break
  }
}
metadata$fold <- cv.folds$fold #adding the fold to the metadata df for 5-fold CV
metadata <- metadata %>% mutate(Infection_status = ifelse(group == "LRTI", 1, 0)) #generating new column, splitting cases into Infection or no infection for classifier

###########################################
#Classifier for FABP4 - LRTI vs [IPC+CTRL]
###########################################
fold.roc <- list() # For storing ROC results
predictions <- list() #for storing predictions

for (k in c(1:5)) {
  test.fold <- metadata[metadata$fold==k,]
  train.folds <- metadata[metadata$fold!=k,]
  
  glm_model <- glm(
    Infection_status ~ FABP4,
    data = train.folds,
    family="binomial"
  )
  
  #predicting on test set and storing predictions for each fold
  pred_probs <- predict(glm_model, newdata = test.fold, type = "response")
  predictions[[k]] <- data.frame(
    Patient = test.fold$Patient,
    true_label = test.fold$Infection_status,
    predicted_prob = pred_probs,
    fold = k
  )
  
  #generating ROCs for each fold
  fold.roc[[k]] <- roc(
    response = test.fold$Infection_status,
    predictor = pred_probs,
    quiet = TRUE
  )
  
  cat(sprintf("Fold %d AUC: %.3f\n", k, fold.roc[[k]]$auc))
}

predictions.df <- do.call(rbind, predictions) #creating a dataframe of all of the out-of-fold predictions

overall_roc <- roc(response = predictions.df$true_label,
                   predictor = predictions.df$predicted_prob,
                   quiet = TRUE) #generating a ROC curve with the out-of-fold predictions

fold_auc <- sapply(fold.roc, function(x) auc(x)) #calculating mean AUC value by averaging the AUCs of each of folds
mean_auc_FABP4 <- mean(fold_auc)
cat(sprintf("Mean AUC across folds: %.3f\n", mean_auc_FABP4))

n_boot <- 1000 #generating bootstrapped values for violin plot/AUC 95% CI
auc_bootstrap <- numeric(n_boot)

for (i in 1:n_boot) {
  boot_sample <- predictions.df %>% 
    slice_sample(n = nrow(predictions.df), replace = TRUE)   # Resample the predictions dataframe with replacement
  
  roc_boot <- roc(response = boot_sample$true_label,
                  predictor = boot_sample$predicted_prob,
                  quiet = TRUE)   # Compute AUC for this bootstrap sample
  
  auc_bootstrap[i] <- auc(roc_boot) # store the AUC (x1000)
}

ci_95 <- quantile(auc_bootstrap, c(0.025, 0.975)) ##calculation of the 95% CI
cat("95% CI for AUC:", round(ci_95[1], 2), "-", round(ci_95[2], 2), "\n")

fabp4_auc_df <- data.frame(
  Classifier = "FABP4",
  AUC = auc_bootstrap
) #storing the AUCs for the violin plot (fig 5a)

###########################################
#Classifier for SDI - LRTI vs [IPC+CTRL]
###########################################
fold.roc <- list() # For storing ROC results
predictions <- list() #for storing predictions

for (k in c(1:5)) {
  test.fold <- metadata[metadata$fold==k,]
  train.folds <- metadata[metadata$fold!=k,]
  
  glm_model <- glm(
    Infection_status ~ SDI,
    data = train.folds,
    family="binomial"
  )
  
  #predicting on test set and storing predictions for each fold
  pred_probs <- predict(glm_model, newdata = test.fold, type = "response")
  predictions[[k]] <- data.frame(
    Patient = test.fold$Patient,
    true_label = test.fold$Infection_status,
    predicted_prob = pred_probs,
    fold = k
  )
  
  #generating ROCs for each fold
  fold.roc[[k]] <- roc(
    response = test.fold$Infection_status,
    predictor = pred_probs,
    quiet = TRUE
  )
  
  cat(sprintf("Fold %d AUC: %.3f\n", k, fold.roc[[k]]$auc))
}
predictions.df <- do.call(rbind, predictions) #collapsing into a dataframe with out-of-fold predictions

fold_auc <- sapply(fold.roc, function(x) auc(x)) #calculating mean AUC value by average of each fold
mean_auc_SDI <- mean(fold_auc)
cat(sprintf("Mean AUC across folds: %.3f\n", mean_auc_SDI))

#generating bootstrapped values for violin plot/AUC 95% CI
n_boot <- 1000
auc_bootstrap <- numeric(n_boot)

for (i in 1:n_boot) {
  boot_sample <- predictions.df %>% 
    slice_sample(n = nrow(predictions.df), replace = TRUE)   # Resample the predictions dataframe with replacement
  
  roc_boot <- roc(response = boot_sample$true_label,
                  predictor = boot_sample$predicted_prob,
                  quiet = TRUE)   # Compute AUC for this bootstrap sample
  
  auc_bootstrap[i] <- auc(roc_boot)
}

# 95% CI
ci_95 <- quantile(auc_bootstrap, c(0.025, 0.975))
cat("95% CI for AUC:", round(ci_95[1], 2), "-", round(ci_95[2], 2), "\n")

SDI_auc_df <- data.frame(
  Classifier = "Alpha diversity",
  AUC = auc_bootstrap
) #storing for later violin plot

###########################################
#LASSO classifier - LRTI vs [IPC+CTRL]
###########################################
fold.roc <- list()
predictions <- list()
coefficients <- list() #for storing model coefficients from each fold
set.seed(123)

for (k in c(1:5)) {
  test.fold <- metadata[metadata$fold==k,]
  rownames(test.fold) <- test.fold$Patient
  train.folds <- metadata[metadata$fold!=k,]
  
  counts.train <-counts_VST[,train.folds$Patient] #train counts
  counts.test <- counts_VST[,test.fold$Patient] #test counts
  
  # Generating model
  model_LASSO <- glmnet::cv.glmnet(t(counts.train), 
                                   train.folds$Infection_status,
                                   family = "binomial")
  
  #extracting model coefficients and extracting gene symbol
  coefficients_fold <- as.data.frame(coef(model_LASSO, s='lambda.1se', gamma=c("gamma.1se"))[,1] %>% 
                                       .[. != 0])
  colnames(coefficients_fold) <- "coefficient"
  coefficients_fold$ENSEMBL_ID <- rownames(coefficients_fold)
  coefficients_fold <- merge(coefficients_fold, gene.symbol, by.x = "ENSEMBL_ID", by.y = "row.names", all.x = TRUE)
  coefficients_fold <- coefficients_fold[, c("ENSEMBL_ID", "gene_symbol", "coefficient")]
  coefficients_fold$fold <- k
  coefficients[[k]] <- coefficients_fold
  
  #predicting on test set and storing predictions for each fold
  pred_probs <- as.data.frame(predict(model_LASSO, t(counts.test), type = "response"))
  predictions[[k]] <- data.frame(
    Patient = test.fold$Patient,
    true_label = test.fold$Infection_status,
    predicted_prob = pred_probs$lambda.1se,
    fold = k
  )
  
  #generating ROCs for each fold
  fold.roc[[k]] <- roc(
    response = test.fold$Infection_status,
    predictor = pred_probs$lambda.1se,
    quiet = TRUE)
  print(sprintf(
    "Fold %d AUC: %.3f", k, fold.roc[[k]]$auc))
}
coefficients_all <- do.call(rbind, coefficients) #collapsing coefficients
predictions.df <- do.call(rbind, predictions) #collapsing out-of-fold predictions

fold_auc <- sapply(fold.roc, function(x) auc(x)) #calculating mean AUC value for annotation
mean_auc_LASSO <- mean(fold_auc)
cat(sprintf("Mean AUC across folds: %.3f\n", mean_auc_LASSO))

#generating bootstrapped values for violin plot/AUC 95% CI
n_boot <- 1000
auc_bootstrap <- numeric(n_boot)

# Bootstrap loop
for (i in 1:n_boot) {
  boot_sample <- predictions.df %>% 
    slice_sample(n = nrow(predictions.df), replace = TRUE) # Resample the predictions dataframe with replacement

  roc_boot <- roc(response = boot_sample$true_label,
                  predictor = boot_sample$predicted_prob,
                  quiet = TRUE) # Compute AUC for this bootstrap sample
  
  auc_bootstrap[i] <- auc(roc_boot)
}

# 95% CI
ci_95 <- quantile(auc_bootstrap, c(0.025, 0.975))
cat("95% CI for AUC:", round(ci_95[1], 2), "-", round(ci_95[2], 2), "\n")

LASSO_auc_df <- data.frame(
  Classifier = "LASSO",
  AUC = auc_bootstrap
)

###########################################
#Classifier for FABP4+SDI - LRTI vs [IPC+CTRL]
###########################################
set.seed(123)
fold.roc <- list()
predictions <- list()

for (k in c(1:5)) {
  test.fold <- metadata[metadata$fold==k,]
  rownames(test.fold) <- test.fold$Patient
  train.folds <- metadata[metadata$fold!=k,]
  
  glm_model <- glm(
    Infection_status ~ SDI + FABP4,
    data = train.folds,
    family="binomial"
  )
  
  pred <- predict(
    glm_model,
    newdata=test.fold,
    type="response"
  )
  
  #predicting on test set and storing predictions for each fold
  pred_probs <- predict(glm_model, newdata = test.fold, type = "response")
  predictions[[k]] <- data.frame(
    Patient = test.fold$Patient,
    true_label = test.fold$Infection_status,
    predicted_prob = pred_probs,
    fold = k
  )
  
  #storing ROC curves for each fold
  fold.roc[[k]] <- roc(
    response = test.fold$Infection_status,
    predictor = pred_probs,
    quiet = TRUE)
  print(sprintf(
    "Fold %d AUC: %.3f", k, fold.roc[[k]]$auc))
}
predictions.df <- do.call(rbind, predictions)

fold_auc <- sapply(fold.roc, function(x) auc(x))
mean_auc_FABP4_SDI <- mean(fold_auc) #calculating AUC value for annotation
cat(sprintf("Mean AUC across folds: %.3f\n", mean_auc_FABP4_SDI))

#generating bootstrapped values for violin plot/AUC 95% CI
n_boot <- 1000
auc_bootstrap <- numeric(n_boot)

for (i in 1:n_boot) {
  boot_sample <- predictions.df %>% 
    slice_sample(n = nrow(predictions.df), replace = TRUE)   # Resample the predictions dataframe with replacement
  
  # Compute AUC for this bootstrap sample
  roc_boot <- roc(response = boot_sample$true_label,
                  predictor = boot_sample$predicted_prob,
                  quiet = TRUE)
  
  auc_bootstrap[i] <- auc(roc_boot)
}

# 95% CI
ci_95 <- quantile(auc_bootstrap, c(0.025, 0.975))
cat("95% CI for AUC:", round(ci_95[1], 2), "-", round(ci_95[2], 2), "\n")

FABP4_diversity_auc_df <- data.frame(
  Classifier = "FABP4 + alpha diversity",
  AUC = auc_bootstrap
)

###########################################
#Classifier for SDI + LASSO - LRTI vs [IPC+CTRL]
###########################################
set.seed(123)
fold.roc <- list()
predictions <- list()
coefficients <- list()

for (k in c(1:5)) {
  test.fold <- metadata[metadata$fold==k,]
  rownames(test.fold) <- test.fold$Patient
  train.folds <- metadata[metadata$fold!=k,]
  
  # Train data, adding in diversity
  counts.train <- counts_VST[,train.folds$Patient]
  diversity.train <- train.folds$SDI
  counts.div.train <- cbind(t(counts.train), diversity=diversity.train)
  
  # Test data, adding in diversity
  counts.test <- counts_VST[,test.fold$Patient]
  diversity.test <- test.fold$SDI
  counts.div.test <- cbind(t(counts.test), diversity=diversity.test)
  
  # Generating model
  model <- glmnet::cv.glmnet(counts.div.train, 
                             train.folds$Infection_status,
                             family = "binomial",
  )
  
  #extracting model coefficients
  coefficients_fold <- as.data.frame(coef(model, s='lambda.1se', gamma=c("gamma.1se"))[,1] %>% 
                                       .[. != 0])
  colnames(coefficients_fold) <- "coefficient"
  coefficients_fold$ENSEMBL_ID <- rownames(coefficients_fold)
  coefficients_fold <- merge(coefficients_fold, gene.symbol, by.x = "ENSEMBL_ID", by.y = "row.names", all.x = TRUE)
  coefficients_fold <- coefficients_fold[, c("ENSEMBL_ID", "gene_symbol", "coefficient")]
  coefficients_fold$fold <- k
  coefficients[[k]] <- coefficients_fold
  
  #predicting on test set and storing predictions for each fold
  pred_probs <- as.data.frame(predict(model, counts.div.test, type = "response"))
  predictions[[k]] <- data.frame(
    Patient = test.fold$Patient,
    true_label = test.fold$Infection_status,
    predicted_prob = pred_probs$lambda.1se,
    fold = k
  )
  
  #storing ROC curves for each fold
  fold.roc[[k]] <- roc(
    response = test.fold$Infection_status,
    predictor = pred_probs$lambda.1se,
    quiet = TRUE)
  print(sprintf(
    "Fold %d AUC: %.3f", k, fold.roc[[k]]$auc))
}

coefficients.df <- do.call(rbind, coefficients) #collapsing coefficients
predictions.df <- do.call(rbind, predictions) #collapsing the out-of-fld predictions

fold_auc <- sapply(fold.roc, function(x) auc(x))
mean_auc_LASSO_SDI <- mean(fold_auc) # AUC value for annotation
cat(sprintf("Mean AUC across folds: %.3f\n", mean_auc_LASSO_SDI))

#generating bootstrapped values for violin plot/AUC 95% CI
n_boot <- 1000
auc_bootstrap <- numeric(n_boot)

for (i in 1:n_boot) {
  boot_sample <- predictions.df %>% 
    slice_sample(n = nrow(predictions.df), replace = TRUE) # Resample the predictions dataframe with replacement
  
  roc_boot <- roc(response = boot_sample$true_label,
                  predictor = boot_sample$predicted_prob,
                  quiet = TRUE)   # Compute AUC for this bootstrap sample
  
  auc_bootstrap[i] <- auc(roc_boot)
}

# 95% CI
ci_95 <- quantile(auc_bootstrap, c(0.025, 0.975))
cat("95% CI for AUC:", round(ci_95[1], 2), "-", round(ci_95[2], 2), "\n")

LASSO_diversity_auc_df <- data.frame(
  Classifier = "LASSO + alpha diversity",
  AUC = auc_bootstrap
)

#############
#Plotting
#############

#For best performing classifier (LASSO + SDI), plotting ROC curves for each of the 5 folds with an average of the 5 folkds
roc.approx <- data.frame(
  fpr.out=seq(0, 1, length.out=100),
  tpr1=0,
  tpr2=0,
  tpr3=0,
  tpr4=0,
  tpr5=0)

for (k in 1:5) { #computing the mean of the curves
  fpr <- rev(1-fold.roc[[k]]$specificities)
  tpr <- rev(fold.roc[[k]]$sensitivities)
  tpr.out <- approx(fpr, tpr,
                    xout=roc.approx$fpr.out,
                    method="linear", ties="ordered")
  roc.approx[,k+1] <- tpr.out$y
  plot(fpr, tpr, main=k)
  points(tpr.out$x, tpr.out$y, col="red")
}
roc.approx$tpr.mean <- rowMeans(roc.approx[,2:6])
roc.approx[1,"tpr.mean"] <- 0 # Force the mean ROC to start at (0,0)
plot(roc.approx$fpr.out, roc.approx$tpr.mean)

#ROC curve w/ all folds shown separately
p <- ggplot()
for (k in 1:5) {
  p <- p + geom_line(
    data=data.frame(
      x=rev(1-fold.roc[[k]]$specificities),
      y=rev(fold.roc[[k]]$sensitivities)),
    aes(x=x,y=y), color="#70A494FF", alpha=0.4, linewidth=0.4)
}
p <- p + geom_line(
  data=roc.approx,
  aes(x=fpr.out, y=tpr.mean),
  col="#70A494FF", linewidth=1)

# Add the diagonal line
p <- p +
  geom_segment(aes(x=0,xend=1,y=0,yend=1),
               linetype="dashed", col="black", linewidth=0.4) +
  labs(x="1 - Specificity",
       y="Sensitivity") +
  xlim(-0.01,1.01) +
  ylim(-0.01,1.01) +
  expand_limits(x = 0, y = 0) +
  annotate("text", x = 0.75, y = 0.2,
           label = paste0("AUC = ", round(mean_auc_LASSO_SDI, 2), "\n", "(", round(ci_95[1], 2), " - ", round(ci_95[2], 2), ")"),
           size = 4, family = "Arial") +
  theme_minimal(base_family = "Arial", base_size = 12) + 
  theme(
    plot.title = element_text(hjust = 0.5, size = 12, family = "Arial"),
    axis.text.x = element_text(size = 12, family = "Arial"),
    axis.text.y = element_text(size = 12, family = "Arial"),
    axis.title.y = element_text(size = 12, family = "Arial"),
    text = element_text(size=12, family="Arial"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line.x = element_line(color = "black", linewidth = 0.4),
    axis.line.y = element_line(color = "black", linewidth = 0.4),
    axis.ticks = element_line(color = "black", linewidth = 0.4)
  ) +
  ggtitle("Multi-gene + SDI performance")
p

#Violin plot for each of the diagnostic classifiers

#first accumulating all of the bootstrapping values for each classifier
diversity_auc_df$Classifier <- "SDI"
fabp4_auc_df$Classifier <- "FABP4"
FABP4_diversity_auc_df$Classifier <- "FABP4 + SDI"
LASSO_auc_df$Classifier <- "Multi-gene"
LASSO_diversity_auc_df$Classifier <- "Multi-gene + SDI"

#all the mean AUCs based on average of 5 folds for each classifier
all_auc_df <- rbind(
  diversity_auc_df,
  fabp4_auc_df,
  FABP4_diversity_auc_df,
  LASSO_auc_df,
  LASSO_diversity_auc_df
)

all_auc_df$Classifier <- factor(all_auc_df$Classifier,
                                levels = c("SDI", "FABP4", "FABP4 + SDI", "Multi-gene", "Multi-gene + SDI")) #setting thte order

classifier_colors <- c(
  "SDI" = "#DE8A5AFF",
  "FABP4" = "#EDBB8AFF",
  "FABP4 + SDI" = "#F6EDBDFF",
  "Multi-gene" = "#B4C8A8FF",
  "Multi-gene + SDI" = "#70A494FF"
)

auc_summary <- all_auc_df %>%
  group_by(Classifier) %>%
  summarise(
    auc_lower = round(quantile(AUC, 0.025), 2),
    auc_upper = round(quantile(AUC, 0.975), 2)
  ) %>%
  mutate(
    auc_mean = case_when(
      Classifier == "SDI" ~ formatC(mean_auc_SDI, format = "f", digits = 2),
      Classifier == "FABP4" ~ formatC(mean_auc_FABP4, format = "f", digits = 2),
      Classifier == "FABP4 + SDI" ~ formatC(mean_auc_FABP4_SDI, format = "f", digits = 2),
      Classifier == "Multi-gene" ~ formatC(mean_auc_LASSO, format = "f", digits = 2),
      Classifier == "Multi-gene + SDI" ~ formatC(mean_auc_LASSO_SDI, format = "f", digits = 2) #formatting of the labels on the plot above each violin
    ),
    label = paste0("AUC = ", auc_mean, "\n(", auc_lower, "–", auc_upper, ")")
  )

violin_classifiers <- ggplot(all_auc_df, aes(x = Classifier, y = AUC, fill = Classifier)) +
  geom_violin(alpha = 0.8, color = "black") +
  geom_boxplot(width = 0.1, outlier.shape = NA, color = "black") +
  geom_text(data = auc_summary, aes(x = Classifier, y = 0.97, label = label),
            inherit.aes = FALSE, family = "Arial", size = 4) +
  scale_fill_manual(values = classifier_colors)+
  scale_y_continuous(limits = c(0.67, 1.0), breaks = seq(0.6, 1.0, by = 0.1)) +
  labs(
    title = "Diagnostic performance of LRTI classifiers",
    x = NULL,
    y = "AUC"
  ) +
  theme_minimal(base_family = "Arial", base_size = 12) +
  scale_x_discrete(labels = c(
    "SDI" = "SDI",
    "FABP4" = expression(italic("FABP4")),
    "FABP4 + SDI" = expression(italic("FABP4")~"+"~SDI),
    "Multi-gene" = "Multi-gene",
    "Multi-gene + SDI" = expression("Multi-gene + SDI")
  )) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 12, family = "Arial"),
    axis.text.x = element_text(size = 12, family = "Arial"),
    axis.text.y = element_text(size = 12, family = "Arial"),
    axis.title.y = element_text(size = 12, family = "Arial"),
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line.x = element_line(color = "black", linewidth = 0.4),
    axis.line.y = element_line(color = "black", linewidth = 0.4),
    axis.ticks = element_line(color = "black", linewidth = 0.4),
  )
violin_classifiers

#############
#Proteomics: FABP4 from SomaScan data
#############

proteomics_data <- read.csv("/Users/emilylydon/Library/CloudStorage/Box-Box/VAP_IvC_project/Code:Source Data File/proteomics_data.csv")

#generating folds for cross-validation (3 folds this time, since total n is lower)
set.seed(123)
min.IPC <- 3 # minimum number of IC sample per fold (10 IPC in 3 folds)
min.CTRL <- 7 # minimum number of ctrls sample per fold (22 CTRL in 3 folds)
while (TRUE) {
  # Generate 5 fold
  cv.folds <- proteomics_data %>%
    dplyr::select(Patient, group) %>%
    mutate(fold=sample(rep(1:3, length.out=nrow(.))))
  
  cv.folds.table <- cv.folds %>%
    group_by(fold) %>%
    dplyr::count(group) #counts breakdown in each fold
  
  print("Generated CV folds with following counts:")
  print(cv.folds.table)
  
  if ((min(cv.folds.table[cv.folds.table$group=="IPC","n"]) < min.IPC) | min(cv.folds.table[cv.folds.table$group=="CTRL", "n"]) < min.CTRL) {
    print("At least one fold has too few colonizer or control samples. Regenerating CV folds...")
  } else {
    break
  }
}
proteomics_data$fold <- cv.folds$fold ##adding the fold to the metadata df for 5-fold CV

proteomics_data$infection_status <- ifelse(proteomics_data$group == "LRTI", 1, 0)  #binary for classifier

#now doing 3-fold CV and saving results
fold.roc <- list() # For storing ROC results
predictions <- list() #for storing predictions
fold_auc <- numeric(3)  

for (k in c(1:3)) {
  test.fold <- proteomics_data[proteomics_data$fold==k,]
  train.folds <- proteomics_data[proteomics_data$fold!=k,]
  
  glm_model <- glm(
    infection_status ~ FABP4_RFU_somascan,
    data = train.folds,
    family="binomial"
  )
  
  #predicting on test set and storing predictions for each fold
  pred_probs <- predict(glm_model, newdata = test.fold, type = "response")
  predictions[[k]] <- data.frame(
    Patient = test.fold$Patient,
    true_label = test.fold$infection_status,
    predicted_prob = pred_probs,
    fold = k
  )
  
  #generating ROCs for each fold
  fold.roc[[k]] <- roc(
    response = test.fold$infection_status,
    predictor = pred_probs,
    quiet = TRUE
  )
  
  fold_auc[k] <- auc(fold.roc[[k]])
  cat(sprintf("Fold %d AUC: %.3f\n", k, fold.roc[[k]]$auc))
}

fold_auc <- sapply(fold.roc, function(x) auc(x)) #calculating mean AUC of the 3 folds
mean_auc_FABP4_protein <- mean(fold_auc)
cat(sprintf("Mean AUC across folds: %.3f\n", mean_auc_FABP4_protein))

predictions.df <- do.call(rbind, predictions) #collapsing out-of-fold predictions

overall_roc <- roc(response = predictions.df$true_label,
                   predictor = predictions.df$predicted_prob,
                   quiet = TRUE)

roc_df <- data.frame(
  fpr = rev(1 - overall_roc$specificities),  # False Positive Rate = 1 - specificity
  tpr = rev(overall_roc$sensitivities))       # True Positive Rate = sensitivity

#generating bootstrapped values for AUC 95% CI
n_boot <- 1000
auc_bootstrap <- numeric(n_boot)

for (i in 1:n_boot) {
  # Resample the predictions dataframe with replacement
  boot_sample <- predictions.df %>% 
    slice_sample(n = nrow(predictions.df), replace = TRUE)
  
  # Compute AUC for this bootstrap sample
  roc_boot <- roc(response = boot_sample$true_label,
                  predictor = boot_sample$predicted_prob,
                  quiet = TRUE)
  
  auc_bootstrap[i] <- auc(roc_boot)
}

# 95% CI
ci_95 <- quantile(auc_bootstrap, c(0.025, 0.975))
cat("95% CI for AUC:", round(ci_95[1], 2), "-", round(ci_95[2], 2), "\n")

# creating the overlapping ROC curve, like we did for panel 5b
roc.approx <- data.frame(
  fpr.out=seq(0, 1, length.out=100),
  tpr1=0,
  tpr2=0,
  tpr3=0)

for (k in 1:3) {
  fpr <- rev(1-fold.roc[[k]]$specificities)
  tpr <- rev(fold.roc[[k]]$sensitivities)
  tpr.out <- approx(fpr, tpr,
                    xout=roc.approx$fpr.out,
                    method="linear", ties="ordered")
  roc.approx[,k+1] <- tpr.out$y
  plot(fpr, tpr, main=k)
  points(tpr.out$x, tpr.out$y, col="red")
}
roc.approx$tpr.mean <- rowMeans(roc.approx[,2:4])
roc.approx[1,"tpr.mean"] <- 0 # Force the mean ROC to start at (0,0)
plot(roc.approx$fpr.out, roc.approx$tpr.mean)

# ROC curve w/ all folds shown separately
FABP4_protein_classifier <- ggplot()
for (k in 1:3) {
  FABP4_protein_classifier <- FABP4_protein_classifier + geom_line(
    data=data.frame(
      x=rev(1-fold.roc[[k]]$specificities),
      y=rev(fold.roc[[k]]$sensitivities)),
    aes(x=x,y=y), color="#70A494FF", alpha=0.4, linewidth=0.4)
}
FABP4_protein_classifier <- FABP4_protein_classifier + geom_line(
  data=roc.approx,
  aes(x=fpr.out, y=tpr.mean),
  col="#70A494FF", linewidth=1)

# Add the diagonal line
FABP4_protein_classifier <- FABP4_protein_classifier +
  geom_segment(aes(x=0,xend=1,y=0,yend=1),
               linetype="dashed", col="black", linewidth=0.4) +
  labs(x="1 - Specificity",
       y="Sensitivity") +
  xlim(-0.01,1.01) +
  ylim(-0.01,1.01) +
  expand_limits(x = 0, y = 0) +
  annotate("text", x = 0.75, y = 0.2,
           label = paste0("AUC = ", sprintf("%.2f", mean_auc_FABP4_protein), "\n", "(", round(ci_95[1], 2), " - ", round(ci_95[2], 2), ")"),
           size = 4, family = "Arial") +
  theme_minimal(base_family = "Arial", base_size = 12) + 
  theme(
    plot.title = element_text(hjust = 0.5, size = 12, family = "Arial"),
    axis.text.x = element_text(size = 12, family = "Arial"),
    axis.text.y = element_text(size = 12, family = "Arial"),
    axis.title.y = element_text(size = 12, family = "Arial"),
    text = element_text(size=12, family="Arial"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line.x = element_line(color = "black", linewidth = 0.4),
    axis.line.y = element_line(color = "black", linewidth = 0.4),
    axis.ticks = element_line(color = "black", linewidth = 0.4)
  ) +
  ggtitle("FABP4 protein diagnostic performance")
FABP4_protein_classifier
