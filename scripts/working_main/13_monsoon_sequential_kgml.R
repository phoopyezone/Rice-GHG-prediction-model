###############################################################################
# Script 13 — Monsoon-Only Sequential KGML Pipeline
#             LSTM + Transformer Encoder vs RF / XGBoost / KNN
#             Myanmar Rice Paddy CH4 Emission | Carbon Credit MRV
#
# Input   : C:/DNDC/simulate_dndc/dndc_outputs/enriched_daily_panel.rds
# Models  : LSTM, Transformer, Random Forest, XGBoost, KNN
# Season  : MONSOON ONLY (Day 166–309). Summer left entirely untouched.
#
# Changes in v1.1
# ─────────────────────────────────────────────────────────────────────────
#  - Added coro to package list (required by torch dataloader loops).
#  - Fixed flat feature matrix cbind: explicit as.matrix() on both sides;
#    replaced !is.finite() on data.frame with column-wise loop on matrix.
#  - Expanded STATIC_MGMT with full management factors + inline comments.
#
# Feature Policy
# ─────────────────────────────────────────────────────────────────────────
#  Sequence features (daily, observable):
#    Precipitation, Irrigation, Ponding
#
#  Static features (farm-observable or public data):
#    Soil   : clay, bulk_density, pH, SOC_static, porosity, sand, silt, elevation
#    Mgmt   : N fertiliser (total/urea/ammonium/n_apps/method),
#             flooding (total/FL1/FL2 days, AWD flag, irrigation code),
#             tillage (method code, n_passes),
#             residue (monsoon fraction, summer fraction)
#    Spatial: lat, lon
#
#  FORBIDDEN (runtime hard-stop — DNDC internal states/process variables):
#    CH4_flux/prod/oxid/pool, SOC/DOC/Microbe/Humads/Humus, NEE, NPP,
#    GrainC/LAI/TotalCropN/LeafC/StemC/RootC, Water_stress, N_stress,
#    Evaporation, Transpiration, IniSoilWater, EndSoilWater,
#    Leaching, Runoff, SoilHetResp, dSOC, LitterC, ManureC, DOC_leach
#
# Split Strategy
# ─────────────────────────────────────────────────────────────────────────
#  Spatial holdout : 20% cell_id → test (never seen during training)
#  Temporal split  : train 2005–2020 | val 2021–2024 (within train cells)
###############################################################################

cat("\n", strrep("=", 72), "\n")
cat("  Script 13 v1.1 — Monsoon Sequential KGML\n")
cat("  LSTM + Transformer vs RF / XGBoost / KNN\n")
cat(strrep("=", 72), "\n\n")

# =============================================================================
# 0. PACKAGES & CONFIG
# =============================================================================

pkgs_cran <- c("data.table","ggplot2","patchwork","scales",
               "ranger","xgboost","FNN","openxlsx","coro")
for (p in pkgs_cran) {
  if (!requireNamespace(p, quietly = TRUE)) {
    message(sprintf("  Installing %s ...", p))
    install.packages(p, repos = "https://cloud.r-project.org")
  }
}

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork)
  library(scales); library(ranger); library(xgboost)
  library(FNN); library(openxlsx); library(coro)
})

if (!requireNamespace("torch", quietly = TRUE))
  install.packages("torch", repos = "https://cloud.r-project.org")
if (!torch::torch_is_installed()) {
  message("  Installing torch backend (~500 MB, one-time only) ...")
  torch::install_torch()
}
library(torch)

set.seed(42)
torch::torch_manual_seed(42)

# ── Paths ─────────────────────────────────────────────────────────────────────
PANEL_RDS <- "C:/DNDC/simulate_dndc/dndc_outputs/enriched_daily_panel.rds"
MODEL_DIR <- "C:/DNDC/simulate_dndc/models"
FIG_DIR   <- "C:/DNDC/simulate_dndc/figures/models"
OUT_DIR   <- "C:/DNDC/simulate_dndc/dndc_outputs"
for (d in c(MODEL_DIR, FIG_DIR, OUT_DIR)) dir.create(d, FALSE, TRUE)

# ── Season window ─────────────────────────────────────────────────────────────
MONSOON_DOY <- c(166L, 309L)
SEQ_LEN     <- MONSOON_DOY[2] - MONSOON_DOY[1] + 1L   # 144 days

# ── Split config ──────────────────────────────────────────────────────────────
TEST_CELL_FRAC <- 0.20
TRAIN_YEARS    <- 2005:2020
VAL_YEARS      <- 2021:2024

# ── LSTM hyperparameters ──────────────────────────────────────────────────────
LSTM_HIDDEN   <- 64L
LSTM_LAYERS   <- 2L
LSTM_DROPOUT  <- 0.25
LSTM_EPOCHS   <- 100L
LSTM_LR       <- 5e-4
LSTM_BATCH    <- 512L
LSTM_PATIENCE <- 15L
LSTM_CLIP     <- 1.0

# ── Transformer hyperparameters ───────────────────────────────────────────────
TF_D_MODEL  <- 64L      # must be divisible by TF_NHEAD
TF_NHEAD    <- 4L
TF_LAYERS   <- 2L
TF_FF_DIM   <- 256L
TF_DROPOUT  <- 0.20
TF_EPOCHS   <- 100L
TF_LR       <- 5e-4
TF_BATCH    <- 512L
TF_PATIENCE <- 15L
TF_CLIP     <- 1.0

device <- if (cuda_is_available()) "cuda" else "cpu"
message(sprintf("  Torch device: %s", device))


# =============================================================================
# 1. LOAD DAILY PANEL — MONSOON DAYS ONLY
# =============================================================================
cat("\n[1] Loading daily panel ...\n")

panel <- readRDS(PANEL_RDS)
setDT(panel)

.rename_safe <- function(dt, old, new) {
  present <- old %in% names(dt)
  if (any(present)) setnames(dt, old[present], new[present], skip_absent = TRUE)
}
.rename_safe(panel,
             old = c("CH4-flux","CH4-prod.","CH4-oxid.","CH4-pool",
                     "Litter-C","Manure-C","DOC-leach","Soil-heterotrophic-respiration"),
             new = c("CH4_flux","CH4_prod","CH4_oxid","CH4_pool",
                     "LitterC","ManureC","DOC_leach","SoilHetResp"))

message(sprintf("  Panel: %s rows × %d cols",
                format(nrow(panel), big.mark=","), ncol(panel)))
message("  Cols: ", paste(head(names(panel), 25), collapse=", "), " ...")

# Filter monsoon only — summer never enters any model object
monsoon_dt <- panel[Day >= MONSOON_DOY[1] & Day <= MONSOON_DOY[2]]
rm(panel); invisible(gc())

message(sprintf("  Monsoon subset: %s rows (Day %d–%d = %d days/year)",
                format(nrow(monsoon_dt), big.mark=","),
                MONSOON_DOY[1], MONSOON_DOY[2], SEQ_LEN))


# =============================================================================
# 2. FORBIDDEN VARIABLE GUARD
# =============================================================================
cat("\n[2] Enforcing feature policy ...\n")

FORBIDDEN_VARS <- c(
  # Target + CH4 process intermediates (direct leakage)
  "CH4_flux","CH4_prod","CH4_oxid","CH4_pool",
  # Carbon pool dynamics — DNDC internal states
  "SOC","dSOC","DOC","Microbe","Humads","Humus",
  "LitterC","ManureC","DOC_leach","SoilHetResp",
  # Ecosystem fluxes — DNDC computed
  "NEE","NPP",
  # Crop physiological state — simulated by DNDC crop growth module
  "GrainC","LAI","TotalCropN","LeafC","StemC","RootC",
  "Water_stress","N_stress",
  # Internal water balance — DNDC computed fluxes
  "Evaporation","Transpiration","IniSoilWater","EndSoilWater",
  "Leaching","Runoff"
)

.enforce_feature_policy <- function(feature_vec, context="") {
  hits <- intersect(feature_vec, FORBIDDEN_VARS)
  if (length(hits) > 0)
    stop(sprintf(
      "[FORBIDDEN FEATURE — %s] DNDC internal variables not permitted:\n  %s\n",
      context, paste(hits, collapse=", ")))
  invisible(TRUE)
}

message(sprintf("  %d forbidden variables registered.", length(FORBIDDEN_VARS)))


# =============================================================================
# 3. COLUMN RENAMING & FEATURE DECLARATION
# =============================================================================
cat("\n[3] Declaring features ...\n")

# Soil columns arrive hyphenated from cells.rds
SOIL_RENAME_MAP <- c(
  "clay_0-5cm"="clay_fraction", "bdod_0-5cm"="bulk_density",
  "pH_0-5cm"="soil_pH",         "soc_0-5cm"="SOC_static",
  "poros_0-5cm"="porosity",     "sand_0-5cm"="sand_fraction",
  "silt_0-5cm"="silt_fraction", "elevation"="elevation_m"
)
present_soil <- names(SOIL_RENAME_MAP)[names(SOIL_RENAME_MAP) %in% names(monsoon_dt)]
if (length(present_soil) > 0)
  setnames(monsoon_dt, present_soil, SOIL_RENAME_MAP[present_soil])

if (all(c("x","y") %in% names(monsoon_dt)) && !"lon" %in% names(monsoon_dt))
  setnames(monsoon_dt, c("x","y"), c("lon","lat"))

# ── Time-varying sequence features (daily, observable) ────────────────────────
SEQ_FEATURES <- c(
  "Precipitation",   # mm/day — NASA POWER rainfall (drives DNDC water balance)
  "Irrigation",      # mm/day — management schedule input from scene_registry
  "Ponding"          # mm     — water depth; set by flooding management + rainfall
)
.enforce_feature_policy(SEQ_FEATURES, "SEQ_FEATURES")

# ── Static soil features (SoilGrids 0–5 cm, public data) ─────────────────────
STATIC_SOIL <- c(
  "clay_fraction",   # g/g   — texture; controls CH4 diffusion and drainage
  "bulk_density",    # g/cm³ — affects porosity and aeration status
  "soil_pH",         # pH    — controls methanogen community activity
  "SOC_static",      # kg/kg — initial organic C substrate for methanogenesis
  "porosity",        # v/v   — saturated zone capacity; flooding persistence
  "sand_fraction",   # g/g   — drainage / water retention proxy
  "silt_fraction",   # g/g   — texture complement
  "elevation_m"      # m     — altitude proxy for temperature + drainage
)

# ── Static management features (from scene_registry; farm-observable) ─────────
STATIC_MGMT <- c(
  # ── Fertiliser ─────────────────────────────────────────────────────────────
  "fert_total_N_kgNha",    # kg N/ha  total N applied across all events
  "fert_urea_kgNha",       # kg N/ha  urea fraction of total N
  "fert_ammonium_kgNha",   # kg N/ha  ammonium sulphate fraction
  "fert_n_applications",   # integer  number of split application events
  "fert_method_code",      # code     0=broadcast, 1=incorporated, 2=injected
  # ── Water management ───────────────────────────────────────────────────────
  "total_flood_days",      # days     FL1 + FL2 combined flooding duration
  "FL1_flood_days",        # days     first flood event duration
  "FL2_flood_days",        # days     second flood event duration (if present)
  "awd_flag",              # binary   0=continuous flooding, 1=AWD
  "irrigation_control_code", # code   0=rainfed, 1=controlled, 2=continuous
  # ── Tillage ────────────────────────────────────────────────────────────────
  "till_method_code",      # code     0=no-till, 1=conventional, 2=deep plow
  "till_applications",     # integer  number of tillage passes per season
  # ── Residue management ─────────────────────────────────────────────────────
  "residue_monsoon_frac",  # fraction monsoon crop residue returned to soil
  "residue_summer_frac"    # fraction summer crop residue returned to soil
)

# ── Spatial context ───────────────────────────────────────────────────────────
STATIC_SPATIAL <- c(
  "lat",   # degrees N — temperature / solar radiation gradient
  "lon"    # degrees E — regional precipitation pattern proxy
)

STATIC_FEATURES <- c(STATIC_SOIL, STATIC_MGMT, STATIC_SPATIAL)
.enforce_feature_policy(STATIC_FEATURES, "STATIC_FEATURES")

TARGET <- "ann_CH4_kgCha"   # seasonal monsoon CH4 sum (kg C/ha)

# Resolve which declared columns are actually present
avail_seq    <- intersect(SEQ_FEATURES,    names(monsoon_dt))
avail_static <- intersect(STATIC_FEATURES, names(monsoon_dt))

if (length(avail_seq) == 0)
  stop("No sequence features found. Check Precipitation/Irrigation/Ponding in panel.")
if (!"CH4_flux" %in% names(monsoon_dt))
  stop("CH4_flux missing — run Script 8 to regenerate enriched_daily_panel.rds.")

message(sprintf("  Seq features   (%d/%d): %s",
                length(avail_seq), length(SEQ_FEATURES),
                paste(avail_seq, collapse=", ")))
message(sprintf("  Static features(%d/%d) — soil:%d  mgmt:%d  spatial:%d",
                length(avail_static), length(STATIC_FEATURES),
                sum(avail_static %in% STATIC_SOIL),
                sum(avail_static %in% STATIC_MGMT),
                sum(avail_static %in% STATIC_SPATIAL)))

missing_static <- setdiff(STATIC_FEATURES, names(monsoon_dt))
if (length(missing_static) > 0)
  message("  Not found (skipped): ", paste(missing_static, collapse=", "))


# =============================================================================
# 4. AGGREGATE TARGET + STATIC FEATURES
# =============================================================================
cat("\n[4] Aggregating target and static features ...\n")

monsoon_dt[, monsoon_day := Day - MONSOON_DOY[1] + 1L]

sample_agg <- monsoon_dt[, {
  ch4_sum     <- sum(CH4_flux, na.rm=TRUE)
  static_list <- lapply(avail_static, function(col) {
    v <- .SD[[col]]; if (length(v) > 0) v[[1L]] else NA_real_
  })
  names(static_list) <- avail_static
  c(list(ann_CH4_kgCha=ch4_sum), static_list)
}, by=.(cell_id, scene_id, cal_year),
.SDcols=c("CH4_flux", avail_static)]

setorder(sample_agg, cell_id, scene_id, cal_year)
sample_agg[, sample_idx := .I]
N_samples <- nrow(sample_agg)

message(sprintf("  Samples: %s (%d cells × %d scenes × %d years)",
                format(N_samples, big.mark=","),
                uniqueN(sample_agg$cell_id),
                uniqueN(sample_agg$scene_id),
                uniqueN(sample_agg$cal_year)))

monsoon_dt <- monsoon_dt[
  sample_agg[, .(cell_id, scene_id, cal_year, sample_idx)],
  on=.(cell_id, scene_id, cal_year)
]


# =============================================================================
# 5. BUILD 3D SEQUENCE ARRAY  [N, SEQ_LEN, n_seq_feat]
# =============================================================================
cat(sprintf("\n[5] Building sequence array [%d × %d × %d] ...\n",
            N_samples, SEQ_LEN, length(avail_seq)))

setorder(monsoon_dt, sample_idx, monsoon_day)

rows_check <- monsoon_dt[, .N, by=sample_idx]
n_short    <- sum(rows_check$N != SEQ_LEN)
if (n_short > 0)
  message(sprintf("  WARNING: %d samples have != %d days (padding with 0).",
                  n_short, SEQ_LEN))

n_seq_feat <- length(avail_seq)
X_seq_raw  <- array(0.0, dim=c(N_samples, SEQ_LEN, n_seq_feat),
                    dimnames=list(NULL, NULL, avail_seq))

n_total <- N_samples * SEQ_LEN
for (fi in seq_len(n_seq_feat)) {
  col      <- avail_seq[fi]
  raw_vals <- as.numeric(monsoon_dt[[col]])
  raw_vals[!is.finite(raw_vals)] <- 0.0
  n_have   <- length(raw_vals)
  if (n_have < n_total) raw_vals <- c(raw_vals, rep(0.0, n_total - n_have))
  X_seq_raw[,,fi] <- matrix(raw_vals[seq_len(n_total)],
                            nrow=N_samples, ncol=SEQ_LEN, byrow=TRUE)
}
message(sprintf("  Array: %.1f MB", object.size(X_seq_raw) / 1024^2))

# Static matrix [N, n_static] — explicitly numeric matrix
n_static_feat <- length(avail_static)
X_static_raw  <- matrix(0.0, nrow=N_samples, ncol=n_static_feat,
                        dimnames=list(NULL, avail_static))
for (col in avail_static) {
  vals <- as.numeric(sample_agg[[col]])
  vals[!is.finite(vals)] <- 0.0
  X_static_raw[, col] <- vals
}
message(sprintf("  Static matrix: %d × %d", N_samples, n_static_feat))

y_raw <- as.numeric(sample_agg[[TARGET]])
y_raw[!is.finite(y_raw)] <- 0.0


# =============================================================================
# 6. TRAIN / VAL / TEST SPLIT
# =============================================================================
cat("\n[6] Splitting data ...\n")

all_cells   <- sort(unique(sample_agg$cell_id))
n_test      <- max(1L, round(length(all_cells) * TEST_CELL_FRAC))
set.seed(42)
test_cells  <- sort(sample(all_cells, n_test))
train_cells <- setdiff(all_cells, test_cells)

idx_trainval <- which(sample_agg$cell_id %in% train_cells)
idx_test     <- which(sample_agg$cell_id %in% test_cells)
idx_train    <- idx_trainval[sample_agg$cal_year[idx_trainval] %in% TRAIN_YEARS]
idx_val      <- idx_trainval[sample_agg$cal_year[idx_trainval] %in% VAL_YEARS]

stopifnot(!any(sample_agg$cell_id[idx_train] %in% test_cells))
stopifnot(!any(sample_agg$cell_id[idx_val]   %in% test_cells))
stopifnot(length(intersect(idx_train, idx_val)) == 0L)

message(sprintf("  Cells  : total=%d | train=%d | test=%d (%.0f%% holdout)",
                length(all_cells), length(train_cells), length(test_cells),
                TEST_CELL_FRAC*100))
message(sprintf("  Samples: train=%d | val=%d | test=%d",
                length(idx_train), length(idx_val), length(idx_test)))


# =============================================================================
# 7. FEATURE SCALING  (fit ONLY on training data)
# =============================================================================
cat("\n[7] Fitting scalers on training set ...\n")

# Target: log1p + z-score
y_log      <- log1p(pmax(0.0, y_raw))
y_log_mean <- mean(y_log[idx_train])
y_log_sd   <- max(sd(y_log[idx_train]), 1e-8)
y_scaled   <- (y_log - y_log_mean) / y_log_sd
back_transform <- function(y_sc) expm1(y_sc * y_log_sd + y_log_mean)

# Sequence scaler (per-feature standardisation over training time steps)
seq_means <- vapply(seq_len(n_seq_feat),
                    function(fi) mean(X_seq_raw[idx_train,,fi]), numeric(1))
seq_sds   <- vapply(seq_len(n_seq_feat),
                    function(fi) max(sd(as.vector(X_seq_raw[idx_train,,fi])), 1e-8), numeric(1))
names(seq_means) <- names(seq_sds) <- avail_seq

.scale_seq <- function(arr, m, s) {
  for (fi in seq_along(m)) arr[,,fi] <- (arr[,,fi] - m[fi]) / s[fi]
  arr
}
X_seq_sc <- .scale_seq(X_seq_raw, seq_means, seq_sds)

# Static scaler
static_means <- colMeans(X_static_raw[idx_train,, drop=FALSE], na.rm=TRUE)
static_sds   <- pmax(apply(X_static_raw[idx_train,, drop=FALSE], 2, sd, na.rm=TRUE), 1e-8)
.scale_static <- function(mat, m, s) sweep(sweep(mat, 2, m, "-"), 2, s, "/")
X_static_sc <- .scale_static(X_static_raw, static_means, static_sds)
X_static_sc[!is.finite(X_static_sc)] <- 0.0

message("  Seq means: ", paste(names(seq_means),"=",round(seq_means,2), collapse=" | "))
message("  Seq sds  : ", paste(round(seq_sds,2), collapse=" | "))


# =============================================================================
# 8. FLAT FEATURE MATRIX  (for RF / XGBoost / KNN)
# =============================================================================
cat("\n[8] Building flat feature matrix for ML baselines ...\n")

# Summarise each daily sequence into 6 statistics
flat_blocks <- lapply(seq_len(n_seq_feat), function(fi) {
  mat  <- X_seq_raw[,,fi]          # [N, T] raw values
  feat <- avail_seq[fi]
  out  <- cbind(
    rowSums(mat,  na.rm=TRUE),
    rowMeans(mat, na.rm=TRUE),
    apply(mat, 1, max,      na.rm=TRUE),
    apply(mat, 1, sd,       na.rm=TRUE),
    apply(mat, 1, quantile, probs=0.75, na.rm=TRUE),
    mat[, SEQ_LEN]
  )
  colnames(out) <- paste0(feat, c("_sum","_mean","_max","_sd","_q75","_last"))
  out
})
# Both sides are explicitly matrices — no list coercion
flat_seq_mat <- do.call(cbind, flat_blocks)               # [N, 6 × n_seq_feat]
flat_all_raw <- cbind(flat_seq_mat, X_static_raw)         # [N, seq_stats + static]

# Replace non-finite values column by column (avoids type error on data.frames)
for (j in seq_len(ncol(flat_all_raw)))
  flat_all_raw[!is.finite(flat_all_raw[, j]), j] <- 0.0

.enforce_feature_policy(colnames(flat_all_raw), "ML flat feature matrix")

flat_means <- colMeans(flat_all_raw[idx_train,, drop=FALSE])
flat_sds   <- pmax(apply(flat_all_raw[idx_train,, drop=FALSE], 2, sd), 1e-8)
flat_sc    <- sweep(sweep(flat_all_raw, 2, flat_means, "-"), 2, flat_sds, "/")
for (j in seq_len(ncol(flat_sc)))
  flat_sc[!is.finite(flat_sc[, j]), j] <- 0.0

X_flat_train <- flat_sc[idx_train,, drop=FALSE]
X_flat_val   <- flat_sc[idx_val,,   drop=FALSE]
X_flat_test  <- flat_sc[idx_test,,  drop=FALSE]
y_train_raw  <- y_raw[idx_train]
y_val_raw    <- y_raw[idx_val]
y_test_raw   <- y_raw[idx_test]

n_feat_flat <- ncol(flat_sc)
message(sprintf("  Flat features: %d (%d seq-stats + %d static) | train:%d val:%d test:%d",
                n_feat_flat, ncol(flat_seq_mat), n_static_feat,
                length(idx_train), length(idx_val), length(idx_test)))


# =============================================================================
# 9. EVALUATION HELPER
# =============================================================================

evaluate_model <- function(y_actual, y_pred, model_name) {
  ya <- as.numeric(y_actual); yp <- as.numeric(y_pred)
  ok <- is.finite(ya) & is.finite(yp)
  ya <- ya[ok]; yp <- yp[ok]; n <- length(ya)
  ss_res <- sum((ya-yp)^2); ss_tot <- sum((ya-mean(ya))^2)
  r2    <- if (ss_tot < 1e-12) NA_real_ else 1 - ss_res/ss_tot
  rmse  <- sqrt(mean((ya-yp)^2))
  mae   <- mean(abs(ya-yp))
  mbe   <- mean(yp-ya)
  nrmse <- if (mean(ya) > 1e-6) rmse/mean(ya)*100 else NA_real_
  rpd   <- sd(ya) / max(rmse, 1e-12)
  data.table(Model=model_name,
             R2=round(r2,4), RMSE=round(rmse,3), MAE=round(mae,3),
             MBE=round(mbe,3), NRMSE_pct=round(nrmse,2),
             RPD=round(rpd,3), N=n)
}

results_list <- list()
pred_list    <- list()


# =============================================================================
# 10. TORCH DATASET & DATALOADERS
# =============================================================================
cat("\n[10] Building torch datasets ...\n")

.enforce_feature_policy(avail_seq,    "DL seq features")
.enforce_feature_policy(avail_static, "DL static features")

ch4_seq_dataset <- dataset(
  name = "CH4SeqDataset",
  initialize = function(seq_arr, static_mat, y_vec) {
    self$x_seq    <- torch_tensor(seq_arr,    dtype=torch_float())$to(device=device)
    self$x_static <- torch_tensor(static_mat, dtype=torch_float())$to(device=device)
    self$y        <- torch_tensor(y_vec,       dtype=torch_float())$to(device=device)
  },
  .getitem = function(i) list(
    x_seq    = self$x_seq[i,,],
    x_static = self$x_static[i,],
    y        = self$y[i]$unsqueeze(1L)
  ),
  .length = function() self$x_seq$size(1L)
)

.mk_ds <- function(idx)
  ch4_seq_dataset(
    seq_arr    = X_seq_sc[idx,,,  drop=FALSE],
    static_mat = X_static_sc[idx,, drop=FALSE],
    y_vec      = y_scaled[idx]
  )

train_ds <- .mk_ds(idx_train)
val_ds   <- .mk_ds(idx_val)
test_ds  <- .mk_ds(idx_test)

.mk_dl <- function(ds, batch, shuffle)
  dataloader(ds, batch_size=batch, shuffle=shuffle, drop_last=FALSE)


# =============================================================================
# 11. DL UTILITIES
# =============================================================================

# Sinusoidal PE → R matrix [T, d_model]
.make_sinusoidal_pe <- function(T, d_model) {
  pe <- matrix(0.0, nrow=T, ncol=d_model)
  for (pos in seq_len(T))
    for (i in seq(1, d_model-1, 2)) {
      freq       <- (pos-1) / (10000^((i-1L)/d_model))
      pe[pos, i] <- sin(freq)
      if (i+1L <= d_model) pe[pos, i+1L] <- cos(freq)
    }
  pe
}

.train_epoch <- function(model, dl, optimizer, clip) {
  model$train(); total <- 0.0; n <- 0L
  coro::loop(for (batch in dl) {
    optimizer$zero_grad()
    pred <- model(batch$x_seq, batch$x_static)
    loss <- nnf_mse_loss(pred, batch$y)
    loss$backward()
    nn_utils_clip_grad_norm_(model$parameters, max_norm=clip)
    optimizer$step()
    total <- total + loss$item(); n <- n + 1L
  })
  total / max(n, 1L)
}

.eval_epoch <- function(model, dl) {
  model$eval(); total <- 0.0; n <- 0L
  with_no_grad({
    coro::loop(for (batch in dl) {
      loss  <- nnf_mse_loss(model(batch$x_seq, batch$x_static), batch$y)
      total <- total + loss$item(); n <- n + 1L
    })
  })
  total / max(n, 1L)
}

.get_dl_preds <- function(model, dl) {
  model$eval(); preds <- numeric(0)
  with_no_grad({
    coro::loop(for (batch in dl) {
      out   <- model(batch$x_seq, batch$x_static)
      preds <- c(preds, as.numeric(out$squeeze(2L)$cpu()))
    })
  })
  preds
}

.fit_dl_model <- function(model, train_dl, val_dl, optimizer, scheduler,
                          n_epochs, patience, clip, model_name) {
  best_val <- Inf; pat_ctr <- 0L; best_state <- NULL
  t_hist   <- numeric(n_epochs); v_hist <- numeric(n_epochs); ep_run <- 0L
  
  for (ep in seq_len(n_epochs)) {
    t_loss <- .train_epoch(model, train_dl, optimizer, clip)
    v_loss <- .eval_epoch(model, val_dl)
    scheduler$step()
    t_hist[ep] <- t_loss; v_hist[ep] <- v_loss; ep_run <- ep
    
    if (ep %% 10L == 0L || ep == 1L)
      message(sprintf("    [%s] Ep%3d | train=%.5f  val=%.5f",
                      model_name, ep, t_loss, v_loss))
    
    if (v_loss < best_val - 1e-6) {
      best_val <- v_loss; best_state <- model$state_dict(); pat_ctr <- 0L
    } else {
      pat_ctr <- pat_ctr + 1L
      if (pat_ctr >= patience) {
        message(sprintf("    [%s] Early stop ep%d", model_name, ep)); break
      }
    }
  }
  if (!is.null(best_state)) model$load_state_dict(best_state)
  message(sprintf("    [%s] Best val=%.5f", model_name, best_val))
  list(model=model, train_hist=t_hist[seq_len(ep_run)],
       val_hist=v_hist[seq_len(ep_run)], n_epochs=ep_run)
}


# =============================================================================
# 12. LSTM MODEL
# =============================================================================
cat("\n[12] === LSTM ===\n")

lstm_net <- nn_module(
  "CH4_LSTM",
  initialize = function(n_seq, n_static, hidden, n_layers, dropout) {
    self$lstm <- nn_lstm(
      input_size=n_seq, hidden_size=hidden, num_layers=n_layers,
      batch_first=TRUE, dropout=if(n_layers>1L) dropout else 0.0
    )
    self$seq_norm  <- nn_layer_norm(hidden)
    self$seq_drop  <- nn_dropout(dropout)
    self$static_fc <- nn_sequential(nn_linear(n_static,32L), nn_relu(), nn_dropout(dropout))
    self$head      <- nn_sequential(nn_linear(hidden+32L,64L), nn_relu(),
                                    nn_dropout(dropout), nn_linear(64L,1L))
  },
  forward = function(x_seq, x_static) {
    out    <- self$lstm(x_seq)
    h_last <- out[[1]][,-1L,]         # last time step [B, hidden]
    h_last <- self$seq_drop(self$seq_norm(h_last))
    s      <- self$static_fc(x_static)
    self$head(torch_cat(list(h_last, s), dim=2L))
  }
)

lstm_model <- lstm_net(n_seq_feat, n_static_feat,
                       LSTM_HIDDEN, LSTM_LAYERS, LSTM_DROPOUT)$to(device=device)
lstm_opt   <- optim_adam(lstm_model$parameters, lr=LSTM_LR, weight_decay=1e-5)
lstm_sch   <- lr_step(lstm_opt, step_size=20L, gamma=0.5)

train_dl_lstm <- .mk_dl(train_ds, LSTM_BATCH, TRUE)
val_dl_lstm   <- .mk_dl(val_ds,   LSTM_BATCH, FALSE)
test_dl_lstm  <- .mk_dl(test_ds,  LSTM_BATCH, FALSE)

t0_lstm  <- proc.time()
lstm_fit <- .fit_dl_model(lstm_model, train_dl_lstm, val_dl_lstm,
                          lstm_opt, lstm_sch,
                          LSTM_EPOCHS, LSTM_PATIENCE, LSTM_CLIP, "LSTM")
elapsed_lstm <- (proc.time()-t0_lstm)[3]
message(sprintf("  LSTM: %.1f s | %d epochs", elapsed_lstm, lstm_fit$n_epochs))

lstm_pred_test  <- back_transform(.get_dl_preds(lstm_fit$model, test_dl_lstm))
lstm_pred_train <- back_transform(.get_dl_preds(lstm_fit$model,
                                                .mk_dl(train_ds, LSTM_BATCH, FALSE)))
results_list[["LSTM_train"]] <- evaluate_model(y_train_raw, lstm_pred_train, "LSTM (train)")
results_list[["LSTM_test"]]  <- evaluate_model(y_test_raw,  lstm_pred_test,  "LSTM (test)")
pred_list[["LSTM"]] <- data.table(actual=y_test_raw, predicted=lstm_pred_test,
                                  model="LSTM", residual=lstm_pred_test-y_test_raw)
torch_save(lstm_fit$model$state_dict(), file.path(MODEL_DIR,"lstm_monsoon_state.pt"))
message("  Saved: lstm_monsoon_state.pt")


# =============================================================================
# 13. TRANSFORMER ENCODER MODEL
# =============================================================================
cat("\n[13] === Transformer Encoder ===\n")

tf_has_batch_first <- tryCatch({
  tl <- nn_transformer_encoder_layer(d_model=8L, nhead=2L, batch_first=TRUE)
  rm(tl); TRUE
}, error=function(e) FALSE)
message(sprintf("  batch_first=%s (torch >= 0.10 required)", tf_has_batch_first))

pe_mat <- .make_sinusoidal_pe(SEQ_LEN, TF_D_MODEL)   # [T, d_model]

transformer_net <- nn_module(
  "CH4_Transformer",
  initialize = function(n_seq, n_static, d_model, nhead, n_layers,
                        ff_dim, dropout, pe_matrix) {
    self$input_proj <- nn_linear(n_seq, d_model)
    # Store PE as plain [T, d_model] tensor; unsqueeze to [1, T, d_model] in
    # forward() so it broadcasts correctly against [B, T, d_model].
    # Pre-permuting at init time produces wrong strides that cause a
    # "size of tensor a (B) must match b (T)" error at runtime.
    self$register_buffer(
      "pe",
      torch_tensor(pe_matrix, dtype=torch_float())   # [T, d_model]
    )
    enc_layer <- nn_transformer_encoder_layer(
      d_model=d_model, nhead=nhead, dim_feedforward=ff_dim,
      dropout=dropout, batch_first=TRUE
    )
    self$encoder   <- nn_transformer_encoder(enc_layer, num_layers=n_layers)
    self$static_fc <- nn_sequential(nn_linear(n_static,32L), nn_relu(), nn_dropout(dropout))
    self$pool_drop <- nn_dropout(dropout)
    self$head      <- nn_sequential(nn_linear(d_model+32L,64L), nn_relu(),
                                    nn_dropout(dropout), nn_linear(64L,1L))
  },
  forward = function(x_seq, x_static) {
    # x_seq   : [B, T, F_seq]
    # self$pe : [T, d_model] -> unsqueeze(1L) -> [1, T, d_model]
    #           broadcasts against [B, T, d_model] automatically
    h_proj <- self$input_proj(x_seq)                         # [B, T, d_model]
    h_enc  <- self$encoder(h_proj + self$pe$unsqueeze(1L))   # [B, T, d_model]
    h_pool <- self$pool_drop(h_enc$mean(dim=2L))             # [B, d_model]
    self$head(torch_cat(list(h_pool, self$static_fc(x_static)), dim=2L))
  }
)

tf_model <- transformer_net(n_seq_feat, n_static_feat,
                            TF_D_MODEL, TF_NHEAD, TF_LAYERS, TF_FF_DIM,
                            TF_DROPOUT, pe_matrix=pe_mat)$to(device=device)
tf_opt   <- optim_adam(tf_model$parameters, lr=TF_LR, weight_decay=1e-5)
tf_sch   <- lr_step(tf_opt, step_size=20L, gamma=0.5)

train_dl_tf <- .mk_dl(train_ds, TF_BATCH, TRUE)
val_dl_tf   <- .mk_dl(val_ds,   TF_BATCH, FALSE)
test_dl_tf  <- .mk_dl(test_ds,  TF_BATCH, FALSE)

t0_tf  <- proc.time()
tf_fit <- .fit_dl_model(tf_model, train_dl_tf, val_dl_tf,
                        tf_opt, tf_sch,
                        TF_EPOCHS, TF_PATIENCE, TF_CLIP, "Transformer")
elapsed_tf <- (proc.time()-t0_tf)[3]
message(sprintf("  Transformer: %.1f s | %d epochs", elapsed_tf, tf_fit$n_epochs))

tf_pred_test  <- back_transform(.get_dl_preds(tf_fit$model, test_dl_tf))
tf_pred_train <- back_transform(.get_dl_preds(tf_fit$model,
                                              .mk_dl(train_ds, TF_BATCH, FALSE)))
results_list[["TF_train"]] <- evaluate_model(y_train_raw, tf_pred_train, "Transformer (train)")
results_list[["TF_test"]]  <- evaluate_model(y_test_raw,  tf_pred_test,  "Transformer (test)")
pred_list[["Transformer"]] <- data.table(actual=y_test_raw, predicted=tf_pred_test,
                                         model="Transformer",
                                         residual=tf_pred_test-y_test_raw)
torch_save(tf_fit$model$state_dict(), file.path(MODEL_DIR,"transformer_monsoon_state.pt"))
message("  Saved: transformer_monsoon_state.pt")


# =============================================================================
# 14. RANDOM FOREST (ranger)
# =============================================================================
cat("\n[14] === Random Forest ===\n")
t0 <- proc.time()

rf_model <- ranger::ranger(
  x=X_flat_train, y=y_train_raw,
  num.trees=500L, mtry=max(1L, floor(sqrt(n_feat_flat))),
  min.node.size=5L, importance="impurity", seed=42L, num.threads=4L
)
rf_pred_train <- rf_model$predictions
rf_pred_test  <- predict(rf_model, data=X_flat_test)$predictions
message(sprintf("  %.1f s | OOB R2=%.4f", (proc.time()-t0)[3], rf_model$r.squared))

results_list[["RF_train"]] <- evaluate_model(y_train_raw, rf_pred_train, "RF (train OOB)")
results_list[["RF_test"]]  <- evaluate_model(y_test_raw,  rf_pred_test,  "RF (test)")
pred_list[["RF"]] <- data.table(actual=y_test_raw, predicted=rf_pred_test,
                                model="Random Forest", residual=rf_pred_test-y_test_raw)
rf_imp <- sort(rf_model$variable.importance, decreasing=TRUE)
saveRDS(rf_model, file.path(MODEL_DIR,"rf_monsoon.rds"))
message("  Saved: rf_monsoon.rds")


# =============================================================================
# 15. GRADIENT BOOSTING (XGBoost)
# =============================================================================
cat("\n[15] === Gradient Boosting (XGBoost) ===\n")
t0 <- proc.time()

xgb_dm_train <- xgb.DMatrix(data=X_flat_train, label=y_train_raw)
xgb_dm_val   <- xgb.DMatrix(data=X_flat_val,   label=y_val_raw)
xgb_dm_test  <- xgb.DMatrix(data=X_flat_test,  label=y_test_raw)

xgb_params <- list(
  booster="gbtree", objective="reg:squarederror",
  eta=0.05, max_depth=6L, subsample=0.80, colsample_bytree=0.80,
  min_child_weight=5L, lambda=1.0, alpha=0.1, eval_metric="rmse"
)

xgb_cv <- xgb.cv(params=xgb_params, data=xgb_dm_train,
                 nrounds=600L, nfold=5L, early_stopping_rounds=30L,
                 verbose=1L, print_every_n=100L, seed=42L)

best_rounds <- xgb_cv$best_iteration
if (is.null(best_rounds) || is.na(best_rounds) || best_rounds < 1L)
  best_rounds <- which.min(xgb_cv$evaluation_log$test_rmse_mean)
best_rounds <- as.integer(best_rounds)
message(sprintf("  Best rounds: %d", best_rounds))

xgb_model <- xgb.train(
  params=xgb_params, data=xgb_dm_train, nrounds=best_rounds,
  watchlist=list(train=xgb_dm_train, val=xgb_dm_val),
  print_every_n=100L, verbose=1L
)
xgb_pred_train <- predict(xgb_model, xgb_dm_train)
xgb_pred_test  <- predict(xgb_model, xgb_dm_test)
message(sprintf("  %.1f s", (proc.time()-t0)[3]))

results_list[["GBM_train"]] <- evaluate_model(y_train_raw, xgb_pred_train, "GBM (train CV)")
results_list[["GBM_test"]]  <- evaluate_model(y_test_raw,  xgb_pred_test,  "GBM (test)")
pred_list[["GBM"]] <- data.table(actual=y_test_raw, predicted=xgb_pred_test,
                                 model="Gradient Boosting",
                                 residual=xgb_pred_test-y_test_raw)
xgb_imp_raw <- xgb.importance(feature_names=colnames(X_flat_train), model=xgb_model)
xgb.save(xgb_model, file.path(MODEL_DIR,"xgb_monsoon.bin"))
message("  Saved: xgb_monsoon.bin")


# =============================================================================
# 16. K-NEAREST NEIGHBOURS (FNN)
# =============================================================================
cat("\n[16] === KNN ===\n")
t0 <- proc.time()

k_cands  <- c(3L,5L,7L,10L,15L,20L)
fold_ids <- sample(rep(1:5, length.out=nrow(X_flat_train)))
knn_cv_rmse <- vapply(k_cands, function(k) {
  mean(vapply(1:5, function(f) {
    Xtr <- X_flat_train[fold_ids!=f,, drop=FALSE]
    Xva <- X_flat_train[fold_ids==f,, drop=FALSE]
    ytr <- y_train_raw[fold_ids!=f]; yva <- y_train_raw[fold_ids==f]
    sqrt(mean((yva - FNN::knn.reg(Xtr,Xva,ytr,k)$pred)^2))
  }, numeric(1)))
}, numeric(1))

best_k <- k_cands[which.min(knn_cv_rmse)]
message(sprintf("  Best k=%d (CV RMSE=%.4f)", best_k, min(knn_cv_rmse)))

knn_pred_train <- FNN::knn.reg(X_flat_train, X_flat_train, y_train_raw, best_k)$pred
knn_pred_test  <- FNN::knn.reg(X_flat_train, X_flat_test,  y_train_raw, best_k)$pred
message(sprintf("  %.1f s", (proc.time()-t0)[3]))

knn_label <- sprintf("KNN k=%d", best_k)
results_list[["KNN_train"]] <- evaluate_model(y_train_raw, knn_pred_train,
                                              paste0(knn_label," (train)"))
results_list[["KNN_test"]]  <- evaluate_model(y_test_raw,  knn_pred_test,
                                              paste0(knn_label," (test)"))
pred_list[["KNN"]] <- data.table(actual=y_test_raw, predicted=knn_pred_test,
                                 model=knn_label, residual=knn_pred_test-y_test_raw)
saveRDS(list(X_train=X_flat_train, y_train=y_train_raw, k=best_k,
             flat_means=flat_means, flat_sds=flat_sds),
        file.path(MODEL_DIR,"knn_monsoon.rds"))
message("  Saved: knn_monsoon.rds")


# =============================================================================
# 17. COMPILE & RANK RESULTS
# =============================================================================
cat("\n[17] Compiling results ...\n")

results_all  <- rbindlist(results_list)
results_all[, Split     := fifelse(grepl("train|OOB|CV", Model, ignore.case=TRUE),
                                   "Train/OOB", "Test")]
results_all[, ModelName := sub(" \\(.*","", Model)]
test_results <- results_all[Split=="Test"][order(-R2)]

message("\nTEST SET PERFORMANCE:")
print(test_results[, .(Model, R2, RMSE, MAE, MBE, NRMSE_pct, RPD, N)])

fwrite(results_all,  file.path(OUT_DIR,"S13_metrics_all.csv"))
fwrite(test_results, file.path(OUT_DIR,"S13_metrics_test.csv"))

all_preds_dt <- rbindlist(pred_list, use.names=TRUE)
all_preds_dt[, scene_id := rep(sample_agg$scene_id[idx_test], length(pred_list))]
all_preds_dt[, cal_year := rep(sample_agg$cal_year[idx_test],  length(pred_list))]
all_preds_dt[, cell_id  := rep(sample_agg$cell_id[idx_test],   length(pred_list))]
fwrite(all_preds_dt, file.path(OUT_DIR,"S13_predictions_test.csv"))
message("  CSVs written.")


# =============================================================================
# 18. FIGURES (PNG ONLY)
# =============================================================================
cat("\n[18] Generating figures ...\n")

.theme_ml <- function(base=11)
  theme_bw(base_size=base) +
  theme(plot.title=element_text(face="bold"), panel.grid.minor=element_blank(),
        strip.background=element_rect(fill="grey92"), strip.text=element_text(face="bold"))

.save_png <- function(p, fname, w=12, h=8, dpi=300) {
  ggsave(file.path(FIG_DIR,fname), plot=p, width=w, height=h, dpi=dpi, device="png")
  message(sprintf("  Saved: %s", fname))
}

MODEL_COLS <- c("LSTM"="#C62828","Transformer"="#6A1B9A",
                "Random Forest"="#1565C0","Gradient Boosting"="#2E7D32")
MODEL_COLS[[knn_label]] <- "#E65100"

# ── F1: Performance bar chart ──────────────────────────────────────────────────
perf_long <- melt(test_results[, .(ModelName,R2,RMSE,MAE,NRMSE_pct)],
                  id.vars="ModelName", variable.name="Metric", value.name="Value")
perf_long[, Metric := factor(Metric, levels=c("R2","RMSE","MAE","NRMSE_pct"),
                             labels=c("R²","RMSE (kg C/ha)","MAE (kg C/ha)","NRMSE (%)"))]

p_perf <- ggplot(perf_long, aes(x=ModelName, y=Value, fill=ModelName)) +
  geom_col(alpha=0.85, width=0.65, colour="white", linewidth=0.3) +
  geom_text(aes(label=sprintf("%.3f",Value)), vjust=-0.4, size=3.0, fontface="bold") +
  scale_fill_brewer(palette="Set2", guide="none") +
  facet_wrap(~Metric, scales="free_y", ncol=4) +
  labs(title=sprintf("Model Performance — Monsoon CH4 (Test Set, %d spatial holdout cells)",
                     length(test_cells)),
       x=NULL, y="Value",
       caption="R²: higher better | RMSE/MAE/NRMSE: lower better | Day 166–309") +
  .theme_ml() + theme(axis.text.x=element_text(angle=30,hjust=1))
.save_png(p_perf, "S13_F1_performance_comparison.png", w=15, h=6)

# ── F2: Predicted vs Actual ────────────────────────────────────────────────────
all_preds_plot <- copy(all_preds_dt)
all_preds_plot[, model_f := factor(model,
                                   levels=c("LSTM","Transformer","Random Forest","Gradient Boosting",knn_label))]
lim_max <- max(c(all_preds_plot$actual, all_preds_plot$predicted), na.rm=TRUE)*1.05

r2_ann <- data.table(
  model_f=factor(c("LSTM","Transformer","Random Forest","Gradient Boosting",knn_label),
                 levels=c("LSTM","Transformer","Random Forest","Gradient Boosting",knn_label)),
  R2_v   =c(results_list[["LSTM_test"]]$R2, results_list[["TF_test"]]$R2,
            results_list[["RF_test"]]$R2,   results_list[["GBM_test"]]$R2,
            results_list[["KNN_test"]]$R2),
  RMSE_v =c(results_list[["LSTM_test"]]$RMSE, results_list[["TF_test"]]$RMSE,
            results_list[["RF_test"]]$RMSE,   results_list[["GBM_test"]]$RMSE,
            results_list[["KNN_test"]]$RMSE)
)
r2_ann[, label := sprintf("R²=%.3f\nRMSE=%.1f", R2_v, RMSE_v)]

p_scatter <- ggplot(all_preds_plot, aes(x=actual, y=predicted, colour=model_f)) +
  geom_abline(slope=1, intercept=0, colour="grey40", linewidth=0.7, linetype="dashed") +
  geom_point(alpha=0.20, size=0.7) +
  geom_smooth(method="lm", se=FALSE, linewidth=1.0, formula=y~x) +
  geom_text(data=r2_ann, aes(x=lim_max*0.04,y=lim_max*0.96,label=label),
            inherit.aes=FALSE, hjust=0, vjust=1, size=3.0, fontface="bold", colour="black") +
  scale_colour_manual(values=MODEL_COLS, guide="none") +
  coord_fixed(xlim=c(0,lim_max), ylim=c(0,lim_max)) +
  facet_wrap(~model_f, ncol=3L) +
  labs(title="Predicted vs Actual — Monsoon CH4",
       subtitle="Spatial holdout | Dashed=1:1 | Solid=OLS fit",
       x="Actual CH4 (kg C/ha/season)", y="Predicted CH4 (kg C/ha/season)",
       caption="DNDC 9.5 | Myanmar 2005–2024") +
  .theme_ml()
.save_png(p_scatter, "S13_F2_predicted_vs_actual.png", w=14, h=9)

# ── F3: Residual distributions ────────────────────────────────────────────────
p_resid <- ggplot(all_preds_plot, aes(x=residual, fill=model_f, colour=model_f)) +
  geom_histogram(aes(y=after_stat(density)), bins=50, alpha=0.40, linewidth=0.15) +
  geom_density(alpha=0, linewidth=0.9) +
  geom_vline(xintercept=0, colour="grey25", linewidth=0.7, linetype="dashed") +
  scale_fill_manual(values=MODEL_COLS, guide="none") +
  scale_colour_manual(values=MODEL_COLS, guide="none") +
  facet_wrap(~model_f, ncol=3L, scales="free_y") +
  labs(title="Residual Distributions — All Models",
       subtitle="Residual = Predicted – Actual | Centred at 0 = unbiased",
       x="Residual (kg C/ha/season)", y="Density", caption="Test set | Monsoon only") +
  .theme_ml()
.save_png(p_resid, "S13_F3_residual_distributions.png", w=14, h=8)

# ── F4a: LSTM training curves ─────────────────────────────────────────────────
p_lstm_loss <- ggplot(
  melt(data.table(epoch=seq_along(lstm_fit$train_hist),
                  Train=lstm_fit$train_hist, Val=lstm_fit$val_hist),
       id.vars="epoch", variable.name="Split", value.name="Loss"),
  aes(x=epoch, y=Loss, colour=Split)) +
  geom_line(linewidth=0.9) +
  scale_colour_manual(values=c(Train="#C62828",Val="#1565C0"), name="") +
  labs(title="LSTM Training Curves",
       subtitle=sprintf("hidden=%d | layers=%d | lr=%.0e | batch=%d | pat=%d",
                        LSTM_HIDDEN,LSTM_LAYERS,LSTM_LR,LSTM_BATCH,LSTM_PATIENCE),
       x="Epoch", y="MSE Loss (scaled log-target)") +
  .theme_ml()
.save_png(p_lstm_loss, "S13_F4a_LSTM_training_curves.png", w=10, h=5)

# ── F4b: Transformer training curves ──────────────────────────────────────────
p_tf_loss <- ggplot(
  melt(data.table(epoch=seq_along(tf_fit$train_hist),
                  Train=tf_fit$train_hist, Val=tf_fit$val_hist),
       id.vars="epoch", variable.name="Split", value.name="Loss"),
  aes(x=epoch, y=Loss, colour=Split)) +
  geom_line(linewidth=0.9) +
  scale_colour_manual(values=c(Train="#6A1B9A",Val="#AB47BC"), name="") +
  labs(title="Transformer Training Curves",
       subtitle=sprintf("d_model=%d | heads=%d | layers=%d | lr=%.0e | batch=%d",
                        TF_D_MODEL,TF_NHEAD,TF_LAYERS,TF_LR,TF_BATCH),
       x="Epoch", y="MSE Loss (scaled log-target)") +
  .theme_ml()
.save_png(p_tf_loss, "S13_F4b_Transformer_training_curves.png", w=10, h=5)

# ── F5: Feature importance ────────────────────────────────────────────────────
n_imp     <- min(15L, length(rf_imp))
rf_imp_dt <- data.table(Feature=names(rf_imp)[1:n_imp],
                        Importance=as.numeric(rf_imp)[1:n_imp],
                        Model="Random Forest")
xgb_imp_dt <- xgb_imp_raw[seq_len(min(15L,nrow(xgb_imp_raw))),
                          .(Feature, Importance=Gain, Model="Gradient Boosting")]
imp_all <- rbind(rf_imp_dt, xgb_imp_dt)
imp_all[, Feature := reorder(Feature, Importance, FUN=mean)]

p_imp <- ggplot(imp_all, aes(x=Feature, y=Importance, fill=Model)) +
  geom_col(position="dodge", alpha=0.85, width=0.7) +
  scale_fill_manual(values=c("Random Forest"="#1565C0","Gradient Boosting"="#2E7D32")) +
  coord_flip() +
  labs(title="Feature Importance — RF & XGBoost (Top 15)",
       subtitle="RF: impurity | XGBoost: gain | Flat features: {sum/mean/max/sd/q75/last} × seq + static",
       x=NULL, y="Importance") +
  .theme_ml(base=10)
.save_png(p_imp, "S13_F5_feature_importance.png", w=13, h=8)

# ── F6: Residuals by scenario ─────────────────────────────────────────────────
p_scene <- ggplot(all_preds_plot,
                  aes(x=scene_id, y=residual, fill=model_f, colour=model_f)) +
  geom_hline(yintercept=0, colour="grey40", linewidth=0.5, linetype="dashed") +
  geom_boxplot(alpha=0.40, outlier.size=0.5, linewidth=0.45) +
  scale_fill_manual(values=MODEL_COLS, name="Model") +
  scale_colour_manual(values=MODEL_COLS, name="Model") +
  labs(title="Residuals by Scenario (Test Set)",
       subtitle="Zero = perfect. Spread shows scenario-specific bias.",
       x="Scenario ID", y="Residual (kg C/ha/season)", caption="Monsoon spatial holdout") +
  .theme_ml() +
  theme(axis.text.x=element_text(angle=30,hjust=1), legend.position="bottom")
.save_png(p_scene, "S13_F6_residuals_by_scenario.png", w=16, h=7)

# ── F7: Calibration ───────────────────────────────────────────────────────────
calib_dt <- all_preds_plot[, {
  brks <- quantile(predicted, probs=seq(0,1,0.1), na.rm=TRUE)
  bin  <- as.integer(cut(predicted, breaks=brks, include.lowest=TRUE))
  tmp  <- data.table(predicted=predicted, actual=actual, bin=bin)
  tmp[!is.na(bin), .(mean_pred=mean(predicted), mean_actual=mean(actual), n=.N), by=bin]
}, by=model_f]

p_calib <- ggplot(calib_dt, aes(x=mean_actual, y=mean_pred, colour=model_f)) +
  geom_abline(slope=1, intercept=0, colour="grey40", linewidth=0.7, linetype="dashed") +
  geom_line(alpha=0.60, linewidth=0.7) +
  geom_point(aes(size=n), alpha=0.75) +
  scale_colour_manual(values=MODEL_COLS, name="Model") +
  scale_size_continuous(name="N", range=c(2,6)) +
  labs(title="Calibration — Mean Predicted vs Actual by Decile",
       subtitle="Dashed = perfect calibration",
       x="Mean Actual CH4 (kg C/ha)", y="Mean Predicted CH4 (kg C/ha)") +
  .theme_ml()
.save_png(p_calib, "S13_F7_calibration_decile.png", w=10, h=7)


# =============================================================================
# 19. EXCEL WORKBOOK (5 sheets)
# =============================================================================
cat("\n[19] Writing Excel workbook ...\n")

hparam_dt <- data.table(
  Model = c(rep("LSTM",7), rep("Transformer",7), rep("Random Forest",3),
            rep("Gradient Boosting",7), "KNN", rep("All",6)),
  Parameter = c(
    "hidden","n_layers","dropout","epochs_max","epochs_actual","lr","batch",
    "d_model","n_heads","n_layers","ff_dim","dropout","epochs_max","epochs_actual",
    "n_trees","mtry","min_node_size",
    "eta","max_depth","subsample","colsample_bytree","n_rounds_cv","n_rounds_final","min_child_weight",
    "k_neighbours",
    "seq_len","n_seq_feat","n_static_feat","n_flat_feat","test_cell_frac","train_years"
  ),
  Value = c(
    LSTM_HIDDEN,LSTM_LAYERS,LSTM_DROPOUT,LSTM_EPOCHS,lstm_fit$n_epochs,LSTM_LR,LSTM_BATCH,
    TF_D_MODEL,TF_NHEAD,TF_LAYERS,TF_FF_DIM,TF_DROPOUT,TF_EPOCHS,tf_fit$n_epochs,
    500,max(1L,floor(sqrt(n_feat_flat))),5,
    0.05,6,0.80,0.80,600,best_rounds,5,
    best_k,
    SEQ_LEN,n_seq_feat,n_static_feat,n_feat_flat,
    TEST_CELL_FRAC,paste0(min(TRAIN_YEARS),"-",max(TRAIN_YEARS))
  )
)

preds_export <- copy(all_preds_dt)
preds_export[, `:=`(actual=round(actual,4), predicted=round(predicted,4),
                    residual=round(residual,4))]

resid_summary <- all_preds_dt[, .(
  mean_res=round(mean(residual,na.rm=TRUE),3),
  sd_res  =round(sd(residual,  na.rm=TRUE),3),
  p5      =round(quantile(residual,0.05,na.rm=TRUE),3),
  p25     =round(quantile(residual,0.25,na.rm=TRUE),3),
  median  =round(median(residual,  na.rm=TRUE),3),
  p75     =round(quantile(residual,0.75,na.rm=TRUE),3),
  p95     =round(quantile(residual,0.95,na.rm=TRUE),3),
  pct_within_10pct=round(mean(abs(residual/(actual+1e-8))<0.10)*100,1)
), by=model]

wb <- createWorkbook()
addWorksheet(wb,"1_Overall_Metrics");  writeData(wb,1,results_all)
addWorksheet(wb,"2_Test_Metrics");     writeData(wb,2,test_results)
addWorksheet(wb,"3_Predictions");      writeData(wb,3,preds_export)
addWorksheet(wb,"4_Residual_Summary"); writeData(wb,4,resid_summary)
addWorksheet(wb,"5_Hyperparameters");  writeData(wb,5,hparam_dt)
conditionalFormatting(wb,2, cols=which(names(test_results)=="R2"),
                      rows=2:(nrow(test_results)+1),
                      style=c("#FCE4D6","#F4B8A4","#E26B5A"),
                      rule=c(0.5,0.8,0.95), type="colourScale")

xlsx_path <- file.path(OUT_DIR,"S13_model_results.xlsx")
saveWorkbook(wb, xlsx_path, overwrite=TRUE)
message(sprintf("  Saved: %s", basename(xlsx_path)))


# =============================================================================
# 20. SAVE SCALING OBJECTS & ARCH SPECS
# =============================================================================
cat("\n[20] Saving metadata ...\n")

saveRDS(list(
  y_log_mean=y_log_mean, y_log_sd=y_log_sd, back_transform=back_transform,
  seq_features=avail_seq, seq_means=seq_means, seq_sds=seq_sds,
  static_features=avail_static, static_means=static_means, static_sds=static_sds,
  flat_means=flat_means, flat_sds=flat_sds, flat_feature_names=colnames(X_flat_train),
  test_cells=test_cells, train_cells=train_cells,
  idx_train=idx_train, idx_val=idx_val, idx_test=idx_test,
  TRAIN_YEARS=TRAIN_YEARS, VAL_YEARS=VAL_YEARS,
  MONSOON_DOY=MONSOON_DOY, SEQ_LEN=SEQ_LEN
), file.path(MODEL_DIR,"S13_scaling_meta.rds"))

saveRDS(list(n_seq=n_seq_feat, n_static=n_static_feat,
             hidden=LSTM_HIDDEN, n_layers=LSTM_LAYERS, dropout=LSTM_DROPOUT),
        file.path(MODEL_DIR,"S13_lstm_arch.rds"))

saveRDS(list(n_seq=n_seq_feat, n_static=n_static_feat,
             d_model=TF_D_MODEL, nhead=TF_NHEAD, n_layers=TF_LAYERS,
             ff_dim=TF_FF_DIM, dropout=TF_DROPOUT, pe_matrix=pe_mat, seq_len=SEQ_LEN),
        file.path(MODEL_DIR,"S13_transformer_arch.rds"))

message("  Saved: S13_scaling_meta.rds | S13_lstm_arch.rds | S13_transformer_arch.rds")


# =============================================================================
# 21. CONSOLE SUMMARY
# =============================================================================
cat("\n", strrep("=", 72), "\n")
cat("  SCRIPT 13 v1.1 — COMPLETE\n")
cat(strrep("=", 72), "\n\n")

cat(sprintf("  Season   : monsoon only (Day %d–%d = %d days)\n",
            MONSOON_DOY[1], MONSOON_DOY[2], SEQ_LEN))
cat(sprintf("  Samples  : %s  (%d cells × %d scenes × %d years)\n",
            format(N_samples, big.mark=","),
            uniqueN(sample_agg$cell_id), uniqueN(sample_agg$scene_id),
            uniqueN(sample_agg$cal_year)))
cat(sprintf("  Split    : train=%d | val=%d | test=%d (%d cells, %.0f%% holdout)\n",
            length(idx_train), length(idx_val), length(idx_test),
            length(test_cells), TEST_CELL_FRAC*100))
cat(sprintf("  Features : %d seq × %d time-steps | %d static | %d flat\n",
            n_seq_feat, SEQ_LEN, n_static_feat, n_feat_flat))

cat("\n  TEST SET RESULTS (sorted by R²):\n")
for (i in seq_len(nrow(test_results))) {
  r <- test_results[i]
  cat(sprintf("  %-24s  R²=%6.4f  RMSE=%7.2f  MAE=%6.2f  NRMSE=%5.1f%%\n",
              r$Model, r$R2, r$RMSE, r$MAE, r$NRMSE_pct))
}
cat(sprintf("\n  LSTM       : %.1f s | %d epochs\n", elapsed_lstm, lstm_fit$n_epochs))
cat(sprintf("  Transformer: %.1f s | %d epochs\n",  elapsed_tf,   tf_fit$n_epochs))

cat("\n  MODELS:\n")
for (f in c("lstm_monsoon_state.pt","transformer_monsoon_state.pt",
            "rf_monsoon.rds","xgb_monsoon.bin","knn_monsoon.rds",
            "S13_scaling_meta.rds","S13_lstm_arch.rds","S13_transformer_arch.rds"))
  cat(sprintf("    %s/%s\n", MODEL_DIR, f))

cat("\n  FIGURES (PNG):\n")
for (f in c("S13_F1_performance_comparison.png","S13_F2_predicted_vs_actual.png",
            "S13_F3_residual_distributions.png","S13_F4a_LSTM_training_curves.png",
            "S13_F4b_Transformer_training_curves.png","S13_F5_feature_importance.png",
            "S13_F6_residuals_by_scenario.png","S13_F7_calibration_decile.png"))
  cat(sprintf("    %s/%s\n", FIG_DIR, f))

cat("\n  TABLES:\n")
for (f in c("S13_metrics_all.csv","S13_metrics_test.csv",
            "S13_predictions_test.csv","S13_model_results.xlsx"))
  cat(sprintf("    %s/%s\n", OUT_DIR, f))

cat("\n  Run: Rscript 13_monsoon_sequential_kgml.R\n")
cat(strrep("=", 72), "\n\n")