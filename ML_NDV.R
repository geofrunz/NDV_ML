# Install required packages if not already installed:
# install.packages(c("ape", "phangorn", "ips", "ggplot2", "lattice", "caret", "randomForest", "neuralnet", "kernlab", "mboost", "ipred"))

# Load libraries
library(ape)
library(seqinr)
library(phangorn)
library(ips)
library(ggplot2)
library(lattice)
library(caret)
library(randomForest)
library(neuralnet)
library(kernlab)
library(mboost)
library(ipred)
library(reshape2)
library(grid)
library(gridExtra)

#########################
# FUNCTIONS
#########################

# Function to compute the absolute difference between amino acid scale values 
# for the vaccine and challenge sequences.
aa_dif <- function(aa, data) {
  comparison <- matrix(NA, nrow = nrow(data), ncol = length(aa[[data[1, 1]]]))
  for (i in 1:nrow(data)) {
    Vaccine   <- aa[[data[i, 1]]]
    Challenge <- aa[[data[i, 2]]]
    for (pos in 1:length(Vaccine)) {
      comparison[i, pos] <- abs(Vaccine[pos] - Challenge[pos])
    }
  }
  return(comparison)
}

# Function to calculate DNA pairwise distances using the "raw" model.
dna_pair_distance <- function(data, dna) {
  x <- NA
  for (i in 1:nrow(data)) {
    Vaccine   <- data[i, 1]
    Challenge <- data[i, 2]
    pair <- as.DNAbin(dna[c(Vaccine, Challenge), ])
    x[i] <- dist.dna(pair, model = "raw")
  }
  return(x)
}

# Function to remove columns that have identical (zero) values.
remove_identical_columns <- function(data) {
  pos <- 1
  diverse <- c()
  for (i in colnames(data)) {
    ifelse(all(data[, i] == 0),
           no = diverse[pos] <- i,
           yes = pos <- pos)
    pos <- pos + 1
  }
  return(data[, na.exclude(diverse)])
}

# Function to convert amino acid sequences into a numeric scale.
# 'seq' is the amino acid sequence, 'index' is the column in the scale to use,
# and 'scale' is the provided scale (e.g. Atchley_scale).
convertScale <- function(seq, index, scale) {
  for (i in unique(seq)) {
    ifelse(test = i == "X",
           no  = seq[seq == i] <- scale[scale$Amino.acid == i, index],
           yes = seq[seq == i] <- 0)
  }
  return(seq)
}

# Function to apply the conversion to all sequences in a list.
convertListAA <- function(seq, index, scale) {
  aa_scaled <<- list()
  for (aa in labels(seq)) {
    aa_scaled[[aa]] <- as.numeric(convertScale(seq[[aa]], index, scale))
  }
  return(aa_scaled)
}

#########################
# SECTION 1: F Protein Analysis
#########################

# Set working directory (change path as needed)
setwd("XXX")

# Read Atchley scale values
Atchley_scale <- read.csv(file = "~/Atchley_scale.csv", header = TRUE, sep = ";")

# Read DNA sequences for F protein (in FASTA format)
dna <- read.dna(file = "F_sequences.fas", format = "fasta", as.character = TRUE)

# (Optionally, one may select columns based on guidance scores;
# here we simply use the entire alignment.)
dna_selected <- dna

# Translate DNA sequences into amino acids
aa <- list()
for (name in rownames(dna_selected)) {
  aa[[name]] <- translate(dna_selected[name, ])
}

# Read cross-neutralization data
data <- read.csv("CrossNeutralization_data.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Compute scaled amino acid values for various properties
aa_scaled_Polarity          <- convertListAA(aa, "Polarity.index", Atchley_scale)
aa_scaled_Structure         <- convertListAA(aa, "Secondary.structure.factor", Atchley_scale)
aa_scaled_Volume            <- convertListAA(aa, "Volume", Atchley_scale)
aa_scaled_Refractivity_Heat <- convertListAA(aa, "Refractivity_Heat.Capacity", Atchley_scale)
aa_scaled_Charge           <- convertListAA(aa, "Charge_Iso.electric.point", Atchley_scale)

# Check ordering between the DNA sequences and the vaccine accession numbers
cbind(sort(rownames(dna)), sort(levels(factor(data$Vaccine.Acc))))

# Build the working data frame for F protein including DNA distance and AA differences
data_work <- data.frame(
  data,
  distDNA = dna_pair_distance(data, dna_selected),
  aa_dif(aa_scaled_Polarity, data),
  aa_dif(aa_scaled_Structure, data),
  aa_dif(aa_scaled_Volume, data),
  aa_dif(aa_scaled_Refractivity_Heat, data),
  aa_dif(aa_scaled_Charge, data)
)
# (Optional: remove columns with no variation)
# data_work <- remove_identical_columns(data_work)

# Remove first two columns (if these are non-predictors) and check dimensions
data_work <- data_work[, 3:ncol(data_work)]
dim(data_work)

# --- Linear Model for F Protein ---
plot(x = data_work$distDNA, y = data_work$Protection, main = "F Protein: DNA Distance vs Protection")
dist_lm <- lm(Protection ~ distDNA, data = data_work)
summary(dist_lm)
plot(dist_lm)

# Remove near-zero variance predictors
nzv <- nearZeroVar(data_work, saveMetrics = TRUE)
nzv[nzv$nzv, ][1:10, ]
nzv <- nearZeroVar(data_work)
data_work <- data_work[, -nzv]
dim(data_work)

# Remove highly correlated predictors (cutoff > 0.7)
descrCor <- cor(data_work)
summary(descrCor[upper.tri(descrCor)])
highlyCorDescr <- findCorrelation(descrCor, cutoff = 0.7)
data_work <- data_work[, -highlyCorDescr]
descrCor2 <- cor(data_work)
summary(descrCor2[upper.tri(descrCor2)])
dim(data_work)
# Reassign the DNA distance (if removed during correlation filtering)
data_work$distDNA <- dna_pair_distance(data, dna)

# --- Splitting Data for Modeling (F Protein) ---
set.seed(3456)
trainIndex <- createDataPartition(data_work$Protection, p = 0.6, list = FALSE, times = 1)
TestData  <- data_work[-trainIndex, ]
data_work <- data_work[ trainIndex, ]

# --- Linear Model Evaluation with caret ---
ctrl <- trainControl(method = "repeatedcv", number = 5, repeats = 10, selectionFunction = "oneSE")
lm_evaluationF <- train(Protection ~ ., data = data_work[, 1:2], method = "lm", metric = "RMSE", trControl = ctrl)
lm_evaluationF$finalModel
predict_values <- predict(lm_evaluationF, newdata = TestData)
plot(predict_values, TestData$Protection, main = "LM Prediction (F Protein)")
abline(a = 0, b = 1)

# --- Random Forest Modeling ---
ctrl <- trainControl(method = "repeatedcv", number = 5, repeats = 10, selectionFunction = "oneSE")
grid <- expand.grid(.mtry = seq(from = 10, to = ncol(data_work) / 3, by = 2))
random_forest <- train(Protection ~ ., data = data_work, method = "rf", metric = "RMSE", 
                       trControl = ctrl, tuneGrid = grid, ntree = 500, na.action = "na.omit")
random_forest$finalModel
plot(random_forest)

# Optimize random forest tuning
grid <- expand.grid(.mtry = seq(20, 80, 5))
random_forest_optF <- train(Protection ~ ., data = data_work, method = "rf", metric = "RMSE", 
                            trControl = ctrl, tuneGrid = grid, ntree = 500, na.action = "na.omit")
random_forest_optF$finalModel
plot(random_forest_optF)
predict_values <- predict(random_forest_optF, newdata = TestData)
plot(predict_values, TestData$Protection, main = "RF Prediction (F Protein)")
abline(a = 0, b = 1)

# --- Treebag Modeling ---
treebagF <- train(Protection ~ ., data = data_work, method = "treebag", metric = "RMSE", trControl = ctrl)
predict_values <- predict(treebagF, newdata = TestData)
plot(predict_values, TestData$Protection, main = "Treebag Prediction (F Protein)")
abline(a = 0, b = 1)

# --- Adaboost (gamboost) Modeling ---
grid <- expand.grid(.mstop = 500, .prune = TRUE)
gamboost <- train(Protection ~ ., data = data_work, method = "gamboost", metric = "RMSE", 
                  trControl = ctrl, tuneGrid = grid)
gamboost

# --- Neural Network Modeling ---
grid <- expand.grid(.layer1 = 1:3, .layer2 = 0:1, .layer3 = 0)
annF <- train(Protection ~ ., data = data_work, method = "neuralnet", metric = "RMSE", 
              trControl = ctrl, tuneGrid = grid, stepmax = 1e+06, preProcess = "range")
annF$bestTune
plot(annF)  # Plot neural network model

# --- SVM Modeling (Linear and Radial) ---

# SVM Linear
grid <- expand.grid(.C = c(seq(0.01, 0.1, 0.01), 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 2, 4, 6, 10, 20, 30, 40, 50, 60, 70))
svm_linear <- train(Protection ~ ., data = data_work, method = "svmLinear", metric = "RMSE", 
                    trControl = ctrl, tuneGrid = grid, na.action = "na.omit")
svm_linear$finalModel
plot(svm_linear)
predict_values <- predict(svm_linear, newdata = data_work)
plot(predict_values, data_work$Protection, main = "SVM Linear Prediction (F Protein)")
abline(a = 0, b = 1)

# SVM Linear with Centering/Scaling
svm_linear_scaled <- train(Protection ~ ., data = data_work, method = "svmLinear", metric = "RMSE", 
                           trControl = ctrl, tuneGrid = grid, na.action = "na.omit", 
                           preProcess = c("center", "scale"))
svm_linear_scaled$finalModel
plot(svm_linear_scaled)

# SVM Linear Optimization
grid <- expand.grid(.C = seq(0.001, 0.05, 0.001))
svm_linear_optF <- train(Protection ~ ., data = data_work, method = "svmLinear", metric = "RMSE", 
                         trControl = ctrl, tuneGrid = grid, na.action = "na.omit")
svm_linear_optF$finalModel
plot(svm_linear_optF)
gbmImp <- varImp(svm_linear_optF, scale = FALSE)
gbmImp
plot(gbmImp, top = 20)

# SVM Linear Bagging (using caret's bag method)
bagctrl <- bagControl(fit = svmBag$fit, predict = svmBag$pred, aggregate = svmBag$aggregate)
grid <- expand.grid(.C = seq(0.001, 0.05, 0.005))
svm_linear_opt_bag <- train(Protection ~ ., data = data_work, method = "bag", metric = "RMSE", 
                            trControl = data.frame(C = seq(0.001, 0.05, 0.005)),
                            tuneGrid = grid, na.action = "na.omit", bagControl = bagctrl, B = 10)
svm_linear_opt_bag$finalModel
plot(svm_linear_opt_bag)

# SVM Radial
grid <- expand.grid(.C = c(seq(0.01, 0.1, 0.01), 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 2, 4, 6, 10, 20, 30, 40, 50, 60, 70),
                    .sigma = c(seq(0.01, 0.1, 0.01), 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 2, 4, 6, 10))
svm_radial <- train(Protection ~ ., data = data_work, method = "svmRadial", metric = "RMSE", 
                    trControl = ctrl, tuneGrid = grid, na.action = "na.omit")
svm_radial$finalModel
plot(svm_radial)

# SVM Radial Optimization
grid <- expand.grid(.C = seq(0.1, 5, 0.1), .sigma = seq(0.0001, 0.01, 0.005))
svm_radial_optF <- train(Protection ~ ., data = data_work, method = "svmRadial", metric = "RMSE", 
                         trControl = ctrl, tuneGrid = grid, na.action = "na.omit")
svm_radial_optF$finalModel
plot(svm_radial_optF)

# --- Model Evaluation for F Protein ---

# Evaluate several models using resamples
resamps <- resamples(list(RF = random_forest_optF, "SVM linear" = svm_linear_optF,
                          "SVM radial" = svm_radial_optF, LM = lm_evaluationF,
                          TB = treebagF, ANN = annF), metric = c("Rsquared", "RMSE"))
print(resamps)
summary(resamps)
trellis.par.set(theme1)
FMetric <- bwplot(resamps, layout = c(1, 3))

# Compute differences between models
difValues <- diff(resamps)
print(difValues)
summary(difValues)
FRMSE <- dotplot(difValues, metric = "RMSE")
FMAE  <- dotplot(difValues, metric = "MAE")
FR2   <- dotplot(difValues, metric = "Rsquared")
diffF <- grid.arrange(FRMSE, FMAE, FR2, ncol = 1)
grid.arrange(FMetric, diffF, ncol = 2)

# Edit boxplots for performance metrics using ggplot2
data_long <- melt(data = resamps$values, id.vars = "Resample", 
                  measure.vars = c("RF~RMSE", "RF~Rsquared", "SVM~RMSE", "SVM~Rsquared",  
                                   "LM~RMSE", "LM~Rsquared", "TB~RMSE", "TB~Rsquared",  
                                   "ANN~RMSE", "ANN~Rsquared"), na.rm = TRUE)
colnames(data_long) <- c("Resample", "Method", "Value")
Rsquare <- grep("Rsquared", x = data_long$Method)
RMSE    <- grep("RMSE", x = data_long$Method)
data_long$Measure <- "Rsquare"
data_long[RMSE, "Measure"] <- "RMSE"

plot1 <- ggplot(data = data_long[data_long$Measure == "Rsquare", ], aes(x = Method, y = Value)) +
  geom_boxplot() +
  labs(title = "Rsquare")
plot2 <- ggplot(data = data_long[data_long$Measure == "RMSE", ], aes(x = Method, y = Value)) +
  geom_boxplot() +
  labs(title = "RMSE")
grid.newpage()
vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
pushViewport(viewport(layout = grid.layout(2, 1)))
print(plot1, vp = vplayout(1, 1))
print(plot2, vp = vplayout(2, 1))

# Test predictions using RF and SVM models
RF_predict <- predict(random_forest_optF, data_work)
plot(RF_predict, data_work$Protection, main = "RF Predictions (F Protein)")
abline(a = 0, b = 1)
RF_predict <- predict(svm_linear_optF, data_work)
plot(RF_predict, data_work$Protection, main = "SVM Predictions (F Protein)")
abline(a = 0, b = 1)

#########################
# SECTION 2: HN Protein Analysis
#########################

# Read DNA sequences for HN protein
dnaHN <- read.dna(file = "HN_sequences.fas", format = "fasta", as.character = TRUE)
dna_selected_dnaHN <- dnaHN

# Translate HN DNA sequences into amino acids
aaHN <- list()
for (name in rownames(dna_selected_dnaHN)) {
  aaHN[[name]] <- translate(dna_selected_dnaHN[name, ])
}
aaHN

# Read cross-neutralization data (re-read to ensure consistency)
data <- read.csv("CrossNeutralization.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Compute scaled amino acid values for HN protein
aa_scaled_PolarityHN          <- convertListAA(aaHN, "Polarity.index", Atchley_scale)
aa_scaled_StructureHN         <- convertListAA(aaHN, "Secondary.structure.factor", Atchley_scale)
aa_scaled_VolumeHN            <- convertListAA(aaHN, "Volume", Atchley_scale)
aa_scaled_Refractivity_HeatHN <- convertListAA(aaHN, "Refractivity_Heat.Capacity", Atchley_scale)
aa_scaled_ChargeHN           <- convertListAA(aaHN, "Charge_Iso.electric.point", Atchley_scale)

# (Optional) Check ordering/alignment of sequence names and vaccine/challenge accessions
d <- data.frame(cbind(sort(rownames(dna_selected_dnaHN)), sort(levels(factor(data$Vaccine.Acc)))))
print(d)
print(d$X1 == d$X2)

d <- data.frame(cbind(sort(rownames(dna_selected_dnaHN)), sort(rownames(dna_selected)), 
                      sort(levels(factor(data$Vaccine.Acc))), sort(levels(factor(data$Challenge.Acc)))))
print(d)
print(d$X1 == d$X2)

# Build working data frame for HN protein
data_workHN <- data.frame(
  data,
  distDNA = dna_pair_distance(data, dna_selected_dnaHN),
  aa_dif(aa_scaled_PolarityHN, data),
  aa_dif(aa_scaled_StructureHN, data),
  aa_dif(aa_scaled_VolumeHN, data),
  aa_dif(aa_scaled_Refractivity_HeatHN, data),
  aa_dif(aa_scaled_ChargeHN, data)
)
data_workHN <- data_workHN[, 3:ncol(data_workHN)]
data_workHN[1:10, 1:20]

# --- Linear Model for HN Protein ---
plot(x = data_workHN$distDNA, y = data_workHN$Protection, main = "HN Protein: DNA Distance vs Protection")
dist_lm <- lm(Protection ~ distDNA, data = data_workHN)
summary(dist_lm)
plot(dist_lm)

# Remove near-zero variance predictors for HN data
nzv <- nearZeroVar(data_workHN, saveMetrics = TRUE)
nzv[nzv$nzv, ][1:10, ]
nzv <- nearZeroVar(data_workHN)
data_workHN <- data_workHN[, -nzv]
dim(data_workHN)

# Remove highly correlated predictors (cutoff > 0.7)
descrCor <- cor(data_workHN)
summary(descrCor[upper.tri(descrCor)])
highlyCorDescr <- findCorrelation(descrCor, cutoff = 0.7)
data_workHN <- data_workHN[, -highlyCorDescr]
descrCor2 <- cor(data_workHN)
summary(descrCor2[upper.tri(descrCor2)])
dim(data_workHN)
# Reassign DNA distance if necessary
data_workHN$distDNA <- dna_pair_distance(data, dna)

# --- Splitting Data for HN Protein Modeling ---
set.seed(3456)
trainIndex <- createDataPartition(data_workHN$Protection, p = 0.6, list = FALSE, times = 1)
TestDataHN <- data_workHN[-trainIndex, ]
data_workHN <- data_workHN[trainIndex, ]

# --- Linear Model Evaluation (HN Protein) ---
ctrl <- trainControl(method = "repeatedcv", number = 5, repeats = 10, selectionFunction = "oneSE")
lm_evaluationHN <- train(Protection ~ ., data = data_workHN[, 1:2], method = "lm", metric = "RMSE", trControl = ctrl)
lm_evaluationHN$finalModel
predict_values <- predict(lm_evaluationHN, newdata = TestDataHN)
plot(predict_values, TestDataHN$Protection, main = "LM Prediction (HN Protein)")
abline(a = 0, b = 1)
postResample(pred = predict_values, obs = TestDataHN$Protection)

# --- Random Forest Modeling for HN Protein ---
grid <- expand.grid(.mtry = seq(from = 10, to = ncol(data_workHN)/3, by = 2))
random_forest <- train(Protection ~ ., data = data_workHN, method = "rf", metric = "RMSE", 
                       trControl = ctrl, tuneGrid = grid, ntree = 500, na.action = "na.omit")
random_forest$finalModel
plot(random_forest)
grid <- expand.grid(.mtry = seq(2, 80, 5))
random_forest_optHN <- train(Protection ~ ., data = data_workHN, method = "rf", metric = "RMSE", 
                             trControl = ctrl, tuneGrid = grid, ntree = 500, na.action = "na.omit")
random_forest_optHN$finalModel
plot(random_forest_optHN)
predict_values <- predict(random_forest_optHN, newdata = TestDataHN)
plot(predict_values, TestDataHN$Protection, main = "RF Prediction (HN Protein)")
abline(a = 0, b = 1)
postResample(pred = predict_values, obs = TestDataHN$Protection)

# Optimize with centering and scaling
grid <- expand.grid(.mtry = seq(2, 60, 3))
random_forest_optHN_CS <- train(Protection ~ ., data = data_workHN, method = "rf", metric = "RMSE", 
                                trControl = ctrl, tuneGrid = grid, ntree = 500, na.action = "na.omit", 
                                preProcess = c("center", "scale"))
random_forest_optHN_CS$finalModel
plot(random_forest_optHN_CS)
predict_values <- predict(random_forest_optHN_CS, newdata = TestDataHN)
plot(predict_values, TestDataHN$Protection, main = "RF (Center/Scale) Prediction (HN Protein)")
abline(a = 0, b = 1)

# --- Treebag Modeling (HN Protein) ---
treebagHN <- train(Protection ~ ., data = data_workHN, method = "treebag", metric = "RMSE", trControl = ctrl)
predict_values <- predict(treebagHN, newdata = TestDataHN)
plot(predict_values, TestDataHN$Protection, main = "Treebag Prediction (HN Protein)")
abline(a = 0, b = 1)

# --- Adaboost (gamboost) Modeling (HN Protein) ---
grid <- expand.grid(.mstop = 500, .prune = TRUE)
gamboost <- train(Protection ~ ., data = data_workHN, method = "gamboost", metric = "RMSE", 
                  trControl = ctrl, tuneGrid = grid)
gamboost

# --- Neural Network Modeling (HN Protein) ---
grid <- expand.grid(.layer1 = 1:5, .layer2 = 0:2, .layer3 = 0)
annHN <- train(Protection ~ ., data = data_workHN, method = "neuralnet", metric = "RMSE", 
               trControl = ctrl, tuneGrid = grid, stepmax = 1e+06, preProcess = "range")
annHN$bestTune
plot(annHN)

# --- SVM Modeling (HN Protein) ---

# SVM Linear
grid <- expand.grid(.C = c(seq(0.01, 0.1, 0.01), 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 2, 4, 6, 10, 20, 30, 40, 50, 60, 70))
svm_linear <- train(Protection ~ ., data = data_workHN, method = "svmLinear", metric = "RMSE", 
                    trControl = ctrl, tuneGrid = grid, na.action = "na.omit")
svm_linear$finalModel
plot(svm_linear)
predict_values <- predict(svm_linear, newdata = data_workHN)
plot(predict_values, data_workHN$Protection, main = "SVM Linear Prediction (HN Protein)")
abline(a = 0, b = 1)

# SVM Linear (Scaled)
svm_linear_scaled <- train(Protection ~ ., data = data_workHN, method = "svmLinear", metric = "RMSE", 
                           trControl = ctrl, tuneGrid = grid, na.action = "na.omit", 
                           preProcess = c("center", "scale"))
svm_linear_scaled$finalModel
plot(svm_linear_scaled)

# SVM Linear Optimization
grid <- expand.grid(.C = seq(0.001, 0.05, 0.001))
svm_linear_optHN <- train(Protection ~ ., data = data_workHN, method = "svmLinear", metric = "RMSE", 
                          trControl = ctrl, tuneGrid = grid, na.action = "na.omit")
svm_linear_optHN$finalModel
plot(svm_linear_optHN)
gbmImp <- varImp(svm_linear_optHN, scale = FALSE)
gbmImp
plot(gbmImp, top = 20)

# SVM Linear Optimization with Scaling
grid <- expand.grid(.C = seq(0.001, 0.05, 0.001))
svm_linear_optHN_SC <- train(Protection ~ ., data = data_workHN, method = "svmLinear", metric = "RMSE", 
                             trControl = ctrl, tuneGrid = grid, na.action = "na.omit", 
                             preProcess = c("center", "scale"))
svm_linear_optHN_SC$finalModel
plot(svm_linear_optHN_SC)
gbmImp <- varImp(svm_linear_optHN, scale = FALSE)
gbmImp
plot(gbmImp, top = 20)

# SVM Radial
grid <- expand.grid(.C = c(seq(0.01, 0.1, 0.01), 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 2, 4, 6, 10, 20, 30, 40, 50, 60, 70),
                    .sigma = c(seq(0.01, 0.1, 0.01), 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 2, 4, 6, 10))
svm_radial <- train(Protection ~ ., data = data_workHN, method = "svmRadial", metric = "RMSE", 
                    trControl = ctrl, tuneGrid = grid, na.action = "na.omit")
svm_radial$finalModel
plot(svm_radial)

# SVM Radial Optimization
grid <- expand.grid(.C = seq(0.1, 7, 0.1), .sigma = seq(0.0001, 0.02, 0.001))
svm_radial_optHN <- train(Protection ~ ., data = data_workHN, method = "svmRadial", metric = "RMSE", 
                          trControl = ctrl, tuneGrid = grid, na.action = "na.omit")
svm_radial_optHN$finalModel
plot(svm_radial_optHN)

# --- Model Evaluation for HN Protein ---
resamps <- resamples(list(RF = random_forest_optHN, "SVM linear" = svm_linear_optHN,
                          "SVM radial" = svm_radial_optHN, LM = lm_evaluationHN,
                          TB = treebagHN, ANN = annHN), metric = c("Rsquared", "RMSE"))
print(resamps)
summary(resamps)
trellis.par.set(theme1)
HNMetric <- bwplot(resamps, layout = c(1, 3))
difValues <- diff(resamps)
print(difValues)
summary(difValues)
HNRMSE <- dotplot(difValues, metric = "RMSE")
HNMAE  <- dotplot(difValues, metric = "MAE")
HNR2   <- dotplot(difValues, metric = "Rsquared")
diffHN <- grid.arrange(HNRMSE, HNMAE, HNR2, ncol = 1)
grid.arrange(HNMetric, diffHN, ncol = 2)

# Edit boxplots for HN performance metrics using ggplot2
data_long <- melt(data = resamps$values, id.vars = "Resample", 
                  measure.vars = c("RF~RMSE", "RF~Rsquared", "SVM~RMSE", "SVM~Rsquared",  
                                   "LM~RMSE", "LM~Rsquared", "TB~RMSE", "TB~Rsquared",  
                                   "ANN~RMSE", "ANN~Rsquared"), na.rm = TRUE)
colnames(data_long) <- c("Resample", "Method", "Value")
Rsquare <- grep("Rsquared", x = data_long$Method)
RMSE    <- grep("RMSE", x = data_long$Method)
data_long$Measure <- "Rsquare"
data_long[RMSE, "Measure"] <- "RMSE"

plot1 <- ggplot(data = data_long[data_long$Measure == "Rsquare", ], aes(x = Method, y = Value)) +
  geom_boxplot() +
  labs(title = "Rsquare")
plot2 <- ggplot(data = data_long[data_long$Measure == "RMSE", ], aes(x = Method, y = Value)) +
  geom_boxplot() +
  labs(title = "RMSE")
grid.newpage()
vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
pushViewport(viewport(layout = grid.layout(2, 1)))
print(plot1, vp = vplayout(1, 1))
print(plot2, vp = vplayout(2, 1))

# Test predictions using RF and SVM models (HN Protein)
RF_predict <- predict(random_forest_optHN, data_workHN)
plot(RF_predict, data_workHN$Protection, main = "RF Predictions (HN Protein)")
abline(a = 0, b = 1)
RF_predict <- predict(svm_linear_optHN, data_workHN)
plot(RF_predict, data_workHN$Protection, main = "SVM Predictions (HN Protein)")
abline(a = 0, b = 1)

#########################
# SECTION 3: Merging F and HN Datasets
#########################

# Build F protein dataset
data_work <- data.frame(
  data,
  distDNA = dna_pair_distance(data, dna_selected),
  aa_dif(aa_scaled_Polarity, data),
  aa_dif(aa_scaled_Structure, data),
  aa_dif(aa_scaled_Volume, data),
  aa_dif(aa_scaled_Refractivity_Heat, data),
  aa_dif(aa_scaled_Charge, data)
)
data_work <- data_work[, 3:ncol(data_work)]
data_work[1:10, 1:20]

# Build HN protein dataset
data_workHN <- data.frame(
  data,
  distDNA = dna_pair_distance(data, dna_selected_dnaHN),
  aa_dif(aa_scaled_PolarityHN, data),
  aa_dif(aa_scaled_StructureHN, data),
  aa_dif(aa_scaled_VolumeHN, data),
  aa_dif(aa_scaled_Refractivity_HeatHN, data),
  aa_dif(aa_scaled_ChargeHN, data)
)
data_workHN <- data_workHN[, 3:ncol(data_workHN)]
data_workHN[1:10, 1:20]

# Merge datasets by combining common predictors from F and HN proteins
dataComplete <- cbind(data_work, data_workHN[, 4:ncol(data_workHN)])

# Remove near-zero variance predictors
nzv <- nearZeroVar(dataComplete, saveMetrics = TRUE)
nzv[nzv$nzv, ][1:10, ]
nzv <- nearZeroVar(dataComplete)
dataComplete <- dataComplete[, -nzv]
dim(dataComplete)

# Remove highly correlated predictors (cutoff > 0.7)
descrCor <- cor(dataComplete)
summary(descrCor[upper.tri(descrCor)])
highlyCorDescr <- findCorrelation(descrCor, cutoff = 0.7)
dataComplete <- dataComplete[, -highlyCorDescr]
descrCor2 <- cor(dataComplete)
summary(descrCor2[upper.tri(descrCor2)])
dim(dataComplete)
# Add DNA distance values for F and HN proteins separately
dataComplete$distDNAF  <- dna_pair_distance(data, dna_selected)
dataComplete$distDNAHN <- dna_pair_distance(data, dna_selected_dnaHN)

# --- Splitting Data for Merged Dataset ---
set.seed(3456)
trainIndex <- createDataPartition(dataComplete$Protection, p = 0.6, list = FALSE, times = 1)
TestDataComplete <- dataComplete[-trainIndex, ]
dataComplete <- dataComplete[trainIndex, ]

# --- Linear Model Evaluation (Merged Dataset) ---
ctrl <- trainControl(method = "repeatedcv", number = 5, repeats = 10, selectionFunction = "oneSE")
lm_evaluationComplete <- train(Protection ~ ., data = dataComplete[, 1:2], method = "lm", metric = "RMSE", trControl = ctrl)
lm_evaluationComplete$finalModel
predict_values <- predict(lm_evaluationComplete, newdata = TestDataComplete)
plot(predict_values, TestDataComplete$Protection, main = "LM Prediction (Merged)")
abline(a = 0, b = 1)
postResample(pred = predict_values, obs = TestDataComplete$Protection)

# --- Random Forest Modeling (Merged Dataset) ---
grid <- expand.grid(.mtry = seq(from = 10, to = ncol(dataComplete)/3, by = 2))
random_forest <- train(Protection ~ ., data = dataComplete, method = "rf", metric = "RMSE", 
                       trControl = ctrl, tuneGrid = grid, ntree = 500, na.action = "na.omit")
random_forest$finalModel
plot(random_forest)
grid <- expand.grid(.mtry = seq(40, 120, 5))
random_forest_optComplete <- train(Protection ~ ., data = dataComplete, method = "rf", metric = "RMSE", 
                                   trControl = ctrl, tuneGrid = grid, ntree = 500, na.action = "na.omit")
random_forest_optComplete$finalModel
plot(random_forest_optComplete)
predict_values <- predict(random_forest_optComplete, newdata = TestDataComplete)
plot(predict_values, TestDataComplete$Protection, main = "RF Prediction (Merged)")
abline(a = 0, b = 1)
postResample(pred = predict_values, obs = TestDataComplete$Protection)

# --- Treebag Modeling (Merged Dataset) ---
treebagComplete <- train(Protection ~ ., data = dataComplete, method = "treebag", metric = "RMSE", trControl = ctrl)
predict_values <- predict(treebagComplete, newdata = TestDataComplete)
plot(predict_values, TestDataComplete$Protection, main = "Treebag Prediction (Merged)")
abline(a = 0, b = 1)

# --- Adaboost (gamboost) Modeling (Merged Dataset) ---
grid <- expand.grid(.mstop = 500, .prune = TRUE)
gamboost <- train(Protection ~ ., data = dataComplete, method = "gamboost", metric = "RMSE", 
                  trControl = ctrl, tuneGrid = grid)
gamboost

# --- Neural Network Modeling (Merged Dataset) ---
grid <- expand.grid(.layer1 = 1:5, .layer2 = 0:2, .layer3 = 0)
annComplete <- train(Protection ~ ., data = dataComplete, method = "neuralnet", metric = "RMSE", 
                     trControl = ctrl, tuneGrid = grid, stepmax = 1e+06, preProcess = "range")
annComplete$bestTune
plot(annComplete)

# --- SVM Modeling (Merged Dataset) ---

# SVM Linear
grid <- expand.grid(.C = c(seq(0.01, 0.1, 0.01), 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 2, 4, 6, 10, 20, 30, 40, 50, 60, 70))
svm_linear <- train(Protection ~ ., data = dataComplete, method = "svmLinear", metric = "RMSE", 
                    trControl = ctrl, tuneGrid = grid, na.action = "na.omit")
svm_linear$finalModel
plot(svm_linear)
predict_values <- predict(svm_linear, newdata = dataComplete)
plot(predict_values, dataComplete$Protection, main = "SVM Linear Prediction (Merged)")
abline(a = 0, b = 1)

# SVM Linear with Centering/Scaling (Merged)
svm_linear_scaled <- train(Protection ~ ., data = dataComplete, method = "svmLinear", metric = "RMSE", 
                           trControl = ctrl, tuneGrid = grid, na.action = "na.omit", 
                           preProcess = c("center", "scale"))
svm_linear_scaled$finalModel
plot(svm_linear_scaled)

# SVM Linear Optimization (Merged)
grid <- expand.grid(.C = seq(0.001, 0.05, 0.001))
svm_linear_optComplete <- train(Protection ~ ., data = dataComplete, method = "svmLinear", metric = "RMSE", 
                                trControl = ctrl, tuneGrid = grid, na.action = "na.omit")
svm_linear_optComplete$finalModel
plot(svm_linear_optComplete)
gbmImp <- varImp(svm_linear_optComplete, scale = FALSE)
gbmImp
plot(gbmImp, top = 20)

# SVM Linear Bagging (Merged)
bagctrl <- bagControl(fit = svmBag$fit, predict = svmBag$pred, aggregate = svmBag$aggregate)
grid <- expand.grid(.C = seq(0.001, 0.05, 0.005))
svm_linear_opt_bag <- train(Protection ~ ., data = dataComplete, method = "bag", metric = "RMSE", 
                            trControl = data.frame(C = seq(0.001, 0.05, 0.005)),
                            tuneGrid = grid, na.action = "na.omit", bagControl = bagctrl, B = 10)
svm_linear_opt_bag$finalModel
plot(svm_linear_opt_bag)

# SVM Radial
grid <- expand.grid(.C = c(seq(0.01, 0.1, 0.01), 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 2, 4, 6, 10, 20, 30, 40, 50, 60, 70),
                    .sigma = c(seq(0.01, 0.1, 0.01), 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 2, 4, 6, 10))
svm_radial <- train(Protection ~ ., data = dataComplete, method = "svmRadial", metric = "RMSE", 
                    trControl = ctrl, tuneGrid = grid, na.action = "na.omit")
svm_radial$finalModel
plot(svm_radial)

# SVM Radial Optimization (Merged)
grid <- expand.grid(.C = seq(0.1, 7, 1), .sigma = seq(0.0001, 0.02, 0.001))
svm_radial_optComplete <- train(Protection ~ ., data = dataComplete, method = "svmRadial", metric = "RMSE", 
                                trControl = ctrl, tuneGrid = grid, na.action = "na.omit")
svm_radial_optComplete$finalModel
plot(svm_radial_optComplete)

# --- Model Evaluation for Merged Dataset ---

resamps <- resamples(list("RF Merged" = random_forest_optComplete,
                          "SVM linear Merged" = svm_linear_optComplete,
                          "SVM radial Merged" = svm_radial_optComplete,
                          "LM Complete" = lm_evaluationComplete,
                          "TB Complete" = treebagComplete,
                          "ANN Complete" = annComplete,
                          "RF F" = random_forest_optF,
                          "SVM linear F" = svm_linear_optF,
                          "SVM radial F" = svm_radial_optF,
                          "LM F" = lm_evaluationF,
                          "TB F" = treebagF,
                          "ANN F" = annF,
                          "RF HN" = random_forest_optHN,
                          "SVM linear HN" = svm_linear_optHN,
                          "SVM radial HN" = svm_radial_optHN,
                          "LM HN" = lm_evaluationHN,
                          "TB HN" = treebagHN,
                          "ANN HN" = annHN),
                     metric = c("Rsquared", "RMSE"))
print(resamps)
summary(resamps)
trellis.par.set(theme1)
AllMetric <- bwplot(resamps, layout = c(3, 1))

# Compute differences between all models and plot
difValues <- diff(resamps)
print(difValues)
summary(difValues)
bwplot(difValues, layout = c(3, 1), scales = list(tck = c(1, 0), x = list(cex = 1), y = list(cex = 0.5)))
dotplot(difValues, metric = c("MAE", "RMSE", "Rsquared"))
AllRMSE <- dotplot(difValues, metric = "RMSE", scales = list(tck = c(1, 0), x = list(cex = 1), y = list(cex = 0.3)), cex = 0.5)
AllMAE  <- dotplot(difValues, metric = "MAE", scales = list(tck = c(1, 0), x = list(cex = 1), y = list(cex = 0.3)), cex = 0.5)
AllR2   <- dotplot(difValues, metric = "Rsquared", scales = list(tck = c(1, 0), x = list(cex = 1), y = list(cex = 0.3)), cex = 0.5)
grid.arrange(AllRMSE, AllMAE, AllR2, ncol = 3)

# (Optional) Plot variable importance for the best RF model in the merged dataset
gbmImp <- varImp(random_forest_optComplete, scale = FALSE)
gbmImp
plot(gbmImp, top = 20)

# --- Save / Load Workspace ---
# Uncomment the following line to save the workspace image:
# save.image(file = "ML_NDV.RData")
#load(file = "ML_NDV.RData")
