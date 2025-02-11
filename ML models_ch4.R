# ===============================================================
# Methane Emissions Modeling Framework v2.0
# Author: Phoo Pye Zone
# Last Update: 2025-02-11
# Dataset: ch4_simulated (duplication NOT removed, missing values 
#          from country averages)
# ===============================================================


processed_path <- "/Users/phoopyezone/Documents/2_Research Projects/IRRI_MSI/data/processed"
ch4_simulated <- read.csv(file.path(processed_path, "CH4 emission_Synthetic.csv"))


# --- Dependencies --- #
# Core modeling packages
needs.these <- c("caret", "e1071", "randomForest", "gbm", "nnet", "kernlab")
invisible(lapply(needs.these, require, character.only = TRUE))

# @note: For missing packages, run this:
# install.packages(needs.these)

#====================== DATA PREP SECTION ======================#
# > Target variable: CH4_seasonal
# > Predictors: Environmental + Management factors
#==========================================================#

# !important - These vars showed significance in pilot study
key_predictors <- c(
  # Climate factors
  "temperature",    # °C, daily mean
  "precipitation",  # mm/day
  
  # Soil characteristics 
  "SOC",           # Soil organic carbon (%)
  "pH",            # Soil pH
  "bulk_density",  # g/cm³
  "total_N",       # Total nitrogen (%)
  "Clay",          # Clay content (%)
  "Silt",          # Silt content (%)
  
  # Management 
  "N_amount",      # N fertilizer (kg/ha)
  "organic_fert_amount",
  "water_regime"   # *categorical - water management type
)

# --- Data cleaning & subset creation --- #
model_ready_data <- ch4_simulated[, c("CH4_seasonal", key_predictors)] %>%
  na.omit()  # Removes incomplete cases

# Pro-tip: Always check data structure after cleaning
# str(model_ready_data)    # Uncomment to inspect

#====================== MODEL SETUP ======================#
# Using 10-fold CV for robust performance estimation
#=====================================================#

# --- Cross-validation configuration --- #
set.seed(42)  # My lucky number - ensures reproducibility
cv_specs <- trainControl(
  method = "cv",
  number = 10,           # 10 folds chosen based on sample size
  savePredictions = TRUE,
  allowParallel = TRUE   # Speeds up training on multi-core systems
)

# Common formula for all models
target_formula <- as.formula("CH4_seasonal ~ .")

#====================== MODEL ENSEMBLE ======================#


# --- 1. Classic Linear Model --- #
# @rationale: Baseline model, checks linear relationships
mlr_model <- train(
  form = target_formula,
  data = model_ready_data,
  method = "lm",
  trControl = cv_specs,
  metric = "RMSE"
)

# --- 2. Random Forest --- #
# @rationale: Handles non-linear relationships & interactions
rf_model <- train(
  form = target_formula,
  data = model_ready_data,
  method = "rf",
  trControl = cv_specs,
  metric = "RMSE",
  tuneLength = 5,    # !Adjust if more tuning needed
  importance = TRUE  # Enables variable importance calculation
)

# --- 3. Gradient Boosting Machine --- #
# @rationale: Sequential learning, often good with noisy data
gbm_model <- train(
  form = target_formula,
  data = model_ready_data,
  method = "gbm",
  trControl = cv_specs,
  metric = "RMSE",
  verbose = FALSE,
  tuneLength = 5
)

# --- 4. k-NN Regressor --- #
# @rationale: Local patterns in feature space
# @warning: Sensitive to feature scaling
knn_model <- train(
  form = target_formula,
  data = model_ready_data,
  method = "knn",
  trControl = cv_specs,
  metric = "RMSE",
  tuneLength = 5  # Tries different k values
)

# --- 5. SVR (Radial Kernel) --- #
# @rationale: Good for non-linear patterns with moderate n
# @note: Can be slow with large datasets
svr_model <- train(
  form = target_formula,
  data = model_ready_data,
  method = "svmRadial",
  trControl = cv_specs,
  metric = "RMSE",
  tuneLength = 5
)

# --- 6. Neural Net (Single Hidden Layer) --- #
# @limitation: Basic architecture, consider deep learning for complex patterns
# @todo: Explore keras for deeper architectures if needed
nnet_model <- train(
  form = target_formula,
  data = model_ready_data,
  method = "nnet",
  trControl = cv_specs,
  metric = "RMSE",
  linout = TRUE,
  trace = FALSE,
  MaxNWts = 10000,
  tuneLength = 5
)

#====================== RESULTS SECTION ======================#
# Performance comparison across all models
#=======================================================#

cat("\n=== CH4 Emission Model Performance Summary ===\n")
cat("RMSE Comparison (lower is better):\n")
cat("-------------------------------------------\n")

# Collect and display results
model_performance <- data.frame(
  Model = c("MLR", "RF", "GBM", "kNN", "SVR", "NNet"),
  RMSE = c(
    min(mlr_model$results$RMSE),
    min(rf_model$results$RMSE),
    min(gbm_model$results$RMSE),
    min(knn_model$results$RMSE),
    min(svr_model$results$RMSE),
    min(nnet_model$results$RMSE)
  )
)

# Display formatted results
print(model_performance[order(model_performance$RMSE),])

#====================== DETAILED METRICS ======================#
# Comprehensive model evaluation across multiple metrics
#=======================================================#

# --- Performance Metrics Function --- #
get_model_metrics <- function(model, data) {
  predictions <- predict(model, data)
  actual <- data$CH4_seasonal
  
  # Calculate metrics
  rmse <- sqrt(mean((predictions - actual)^2))
  mae <- mean(abs(predictions - actual))
  r2 <- cor(predictions, actual)^2
  mape <- mean(abs((actual - predictions)/actual)) * 100
  
  # Calculate residuals
  residuals <- actual - predictions
  
  # Additional metrics
  bias <- mean(residuals)
  residual_se <- sd(residuals)
  
  return(list(
    RMSE = rmse,
    MAE = mae,
    R2 = r2,
    MAPE = mape,
    Bias = bias,
    ResidualSE = residual_se
  ))
}

# --- Compute Metrics for All Models --- #
models_list <- list(
  MLR = mlr_model,
  RF = rf_model,
  GBM = gbm_model,
  KNN = knn_model,
  SVR = svr_model,
  NNET = nnet_model
)

# Calculate metrics for each model
detailed_metrics <- lapply(models_list, get_model_metrics, data = model_ready_data)

# Convert to data frame for nice display
metrics_df <- do.call(rbind, lapply(names(detailed_metrics), function(model_name) {
  metrics <- detailed_metrics[[model_name]]
  data.frame(
    Model = model_name,
    RMSE = round(metrics$RMSE, 2),
    MAE = round(metrics$MAE, 2),
    R2 = round(metrics$R2, 3),
    MAPE = round(metrics$MAPE, 2),
    Bias = round(metrics$Bias, 2),
    ResidualSE = round(metrics$ResidualSE, 2)
  )
}))

# --- Display Comprehensive Results --- #
cat("\n=== Comprehensive Model Performance Analysis ===\n")
cat("All metrics (sorted by R² performance):\n")
cat("-------------------------------------------\n")
print(metrics_df[order(-metrics_df$R2),])

# --- Model Interpretation --- #
cat("\n=== Variable Importance (Random Forest) ===\n")
cat("-------------------------------------------\n")
rf_importance <- varImp(rf_model)
print(rf_importance)

# --- Cross-Model Performance Visualization --- #

library(ggplot2)
metrics_long <- reshape2::melt(metrics_df, id.vars = "Model", 
                             variable.name = "Metric", value.name = "Value")
ggplot(metrics_long, aes(x = Model, y = Value)) +
     geom_bar(stat = "identity") +
     facet_wrap(~Metric, scales = "free_y") +
     theme_minimal() +
     theme(axis.text.x = element_text(angle = 45, hjust = 1))

# --- Prediction Template --- #
# new_data <- read.csv("your_new_data.csv")
# predictions <- predict(rf_model, new_data)  # Using RF as example

