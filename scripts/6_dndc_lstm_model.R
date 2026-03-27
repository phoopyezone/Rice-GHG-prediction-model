##### 6_lstm_ch4.R #####
# LSTM model to predict daily CH4-flux from DNDC daily panel
# Backend: torch + luz (no Python/TensorFlow dependency)
# Target: CH4-flux (kg C/ha/day), log1p-transformed
###############################################################################

library(data.table)
library(torch)
library(luz)

# ── 0. Paths ─────────────────────────────────────────────────────────────────
output_dir <- "C:/DNDC/simulate_dndc/dndc_outputs"
data       <- readRDS(file.path(output_dir, "daily_panel_allcells.rds"))

# ── 1. Features ───────────────────────────────────────────────────────────────
# Selected based on DNDC fermentation submodel drivers:
# - Hydrology:   Ponding, IniSoilWater, Precipitation, Irrigation_mm
# - Temperature: SoilT_1cm, SoilT_10cm
# - C substrate: DOC, SOC, Humads, Microbe
# - Plant:       RootC, LAI, GrainC (aerenchyma transport)
# - Management:  Fertilizer_kgN_ha
features <- c("Ponding","SoilT_1cm","SoilT_10cm",
              "DOC","SOC","Humads","Microbe",
              "RootC","LAI","GrainC",
              "Precipitation","IniSoilWater",
              "Fertilizer_kgN_ha")
target   <- "CH4-flux"
lookback <- 7     # 7-day sequence window
n_feat   <- length(features)

# ── 2. Prepare data ───────────────────────────────────────────────────────────
setorder(data, cell_id, year, Day)
dt <- data[, c("cell_id","year","Day", features, target), with=FALSE]
dt <- dt[complete.cases(dt)]
dt[, ch4_log := log1p(get(target))]

# Normalise features (save scaler for inference)
scaler <- data.table(
  feature = features,
  mean    = sapply(features, function(f) mean(dt[[f]], na.rm=TRUE)),
  sd      = sapply(features, function(f) sd(dt[[f]],   na.rm=TRUE))
)
for (f in features) {
  m <- scaler[feature == f, mean]
  s <- scaler[feature == f, sd]
  if (s > 0) dt[, (f) := (get(f) - m) / s]
}

# ── 3. Build sequences per cell-year ─────────────────────────────────────────
make_sequences <- function(mat_x, vec_y, lb) {
  n <- nrow(mat_x)
  if (n <= lb) return(NULL)
  idx <- seq(lb + 1, n)
  X   <- array(NA_real_, c(length(idx), lb, ncol(mat_x)))
  for (i in seq_along(idx)) X[i,,] <- as.matrix(mat_x[(idx[i]-lb):(idx[i]-1), ])
  list(X=X, y=vec_y[idx])
}

groups <- unique(dt[, .(cell_id, year)])
X_list <- vector("list", nrow(groups))
y_list <- vector("list", nrow(groups))

cat("Building sequences...\n")
for (k in seq_len(nrow(groups))) {
  sub <- dt[cell_id == groups$cell_id[k] & year == groups$year[k]]
  res <- make_sequences(sub[, features, with=FALSE], sub$ch4_log, lookback)
  if (!is.null(res)) { X_list[[k]] <- res$X; y_list[[k]] <- res$y }
}
X_list <- Filter(Negate(is.null), X_list)
y_list <- Filter(Negate(is.null), y_list)

X_all <- do.call(abind::abind, c(X_list, list(along=1)))  # [N, lookback, features]
y_all <- unlist(y_list)
cat(sprintf("Total sequences: %s\n", format(length(y_all), big.mark=",")))

# ── 4. Train/test split (80/20 chronological) ────────────────────────────────
n      <- length(y_all)
tr_idx <- seq_len(floor(0.8 * n))
X_tr   <- X_all[tr_idx,,];  y_tr <- y_all[tr_idx]
X_te   <- X_all[-tr_idx,,]; y_te <- y_all[-tr_idx]

# ── 5. torch Dataset ─────────────────────────────────────────────────────────
ch4_dataset <- dataset(
  name       = "ch4_dataset",
  initialize = function(X, y) {
    self$X <- torch_tensor(X, dtype=torch_float())
    self$y <- torch_tensor(matrix(y, ncol=1), dtype=torch_float())
  },
  .getitem = function(i) list(x=self$X[i,,], y=self$y[i,]),
  .length  = function()  nrow(self$X)
)

train_ds <- ch4_dataset(X_tr, y_tr)
test_ds  <- ch4_dataset(X_te, y_te)

# ── 6. LSTM model ─────────────────────────────────────────────────────────────
ch4_lstm <- nn_module(
  "ch4_lstm",
  initialize = function(n_features, hidden1=64, hidden2=32, dropout=0.2) {
    self$lstm1 <- nn_lstm(n_features, hidden1, batch_first=TRUE)
    self$drop1 <- nn_dropout(dropout)
    self$lstm2 <- nn_lstm(hidden1, hidden2, batch_first=TRUE)
    self$drop2 <- nn_dropout(dropout)
    self$fc    <- nn_linear(hidden2, 1)
  },
  forward = function(x) {
    # x: [batch, seq, features]
    c(out, .) %<-% self$lstm1(x)
    out <- self$drop1(out)
    c(out, .) %<-% self$lstm2(out)
    out <- self$drop2(out)
    out <- out[, dim(out)[2], ]   # last timestep: [batch, hidden2]
    self$fc(out)
  }
)

# ── 7. Train with luz ────────────────────────────────────────────────────────
cat("Training LSTM...\n")
fitted <- ch4_lstm |>
  setup(
    loss      = nn_mse_loss(),
    optimizer = optim_adam,
    metrics   = list(luz_metric_mae())
  ) |>
  set_hparams(n_features=n_feat, hidden1=64, hidden2=32, dropout=0.2) |>
  set_opt_hparams(lr=1e-3) |>
  fit(
    train_ds,
    epochs     = 30,
    valid_data = test_ds,
    dataloader_options = list(batch_size=512, shuffle=TRUE),
    callbacks  = list(luz_callback_early_stopping(monitor="valid_loss", patience=5))
  )

# ── 8. Evaluate ──────────────────────────────────────────────────────────────
cat("Evaluating...\n")
preds <- predict(fitted, test_ds)
y_pred    <- expm1(as.numeric(preds))
y_actual  <- expm1(y_te)

rmse <- sqrt(mean((y_pred - y_actual)^2))
mae  <- mean(abs(y_pred - y_actual))
r2   <- 1 - sum((y_actual - y_pred)^2) / sum((y_actual - mean(y_actual))^2)
cat(sprintf("Test RMSE : %.5f kg C/ha/day\n", rmse))
cat(sprintf("Test MAE  : %.5f kg C/ha/day\n", mae))
cat(sprintf("Test R²   : %.4f\n", r2))

# ── 9. Save model + scaler ───────────────────────────────────────────────────
luz_save(fitted, file.path(output_dir, "lstm_ch4.pt"))
saveRDS(scaler,  file.path(output_dir, "lstm_ch4_scaler.rds"))
cat("Saved to:", output_dir, "\n")