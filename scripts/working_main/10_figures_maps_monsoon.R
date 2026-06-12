###############################################################################
# Script 9+10 — Myanmar CH4 Figures & Spatial Maps | MONSOON ONLY
# DOY 166–309 (Jun 15 – Nov 5)  |  13 scenarios  |  263 grid cells
#
# INPUTS : enriched_daily_panel.rds
# OUTPUTS: PNG figures -> FIGURE_DIR (daily/flux) and MAP_DIR (spatial)
#
# FLUX FIGURES (F1–F7)
#   F1  All 13 scenarios — monsoon daily CH4, single panel
#   F2  Faceted by scenario category vs baseline
#   F3  CH4 components (prod / oxid / net) — Baseline, AWD, Hi-Residue
#   F4  Interannual variability 2005-2024 (3 key scenarios)
#   F5  Heatmap: scenario × DOY (7-day rolling mean)
#   F6  Flooding window overlay (mechanism view)
#   F7  Monsoon cumulative CH4 — scenarios ranked by season-end total
#
# SPATIAL MAPS (M1–M6)
#   M1  Baseline mean monsoon CH4
#   M2  All scenarios — small-multiple tile maps
#   M3  Delta vs baseline (diverging)
#   M4  AWD mitigation — absolute + % reduction
#   M5  Coefficient of variation (interannual variability)
#   M6  Optimal (lowest-emitting) scenario per cell
###############################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
})
###############################################################################
# Script 9+10 — Myanmar CH4 Figures & Spatial Maps | MONSOON ONLY
# DOY 166–309 (Jun 15 – Nov 5)  |  13 scenarios  |  263 grid cells
#
# INPUTS : enriched_daily_panel.rds
# OUTPUTS: PNG figures -> FIGURE_DIR (daily/flux) and MAP_DIR (spatial)
#
# FLUX FIGURES (F1–F7)
#   F1  All 13 scenarios — monsoon daily CH4, single panel
#   F2  Faceted by scenario category vs baseline
#   F3  CH4 components (prod / oxid / net) — Baseline, AWD, Hi-Residue
#   F4  Interannual variability 2005-2024 (3 key scenarios)
#   F5  Heatmap: scenario × DOY (7-day rolling mean)
#   F6  Flooding window overlay (mechanism view)
#   F7  Monsoon cumulative CH4 — scenarios ranked by season-end total
#
# SPATIAL MAPS (M1–M6)
#   M1  Baseline mean monsoon CH4
#   M2  All scenarios — small-multiple tile maps
#   M3  Delta vs baseline (diverging)
#   M4  AWD mitigation — absolute + % reduction
#   M5  Coefficient of variation (interannual variability)
#   M6  Optimal (lowest-emitting) scenario per cell
###############################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

# =============================================================================
# 0. PATHS & CONSTANTS
# =============================================================================

PANEL_RDS  <- "C:/DNDC/simulate_dndc/dndc_outputs/enriched_daily_panel.rds"
FIGURE_DIR <- "C:/DNDC/simulate_dndc/figures/monsoon"
MAP_DIR    <- "C:/DNDC/simulate_dndc/figures/monsoon/maps"
dir.create(FIGURE_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(MAP_DIR,    showWarnings = FALSE, recursive = TRUE)

MONSOON    <- c(166L, 309L)   # DOY window
TILE_SIZE  <- 0.48            # slightly < 0.5° to leave gaps between tiles

stopifnot("enriched_daily_panel.rds not found" = file.exists(PANEL_RDS))

# =============================================================================
# 1. METADATA
# =============================================================================

SCENE_ORDER <- c("S000","S001","S002","S003","S004","S005","S006",
                 "S008","S009","S010","S011","S012","S013")

SCENE_SHORT <- c(
  S000="Baseline",   S001="Inj.Fert",   S002="Med-Hi Fert",
  S003="Ext.Fert",   S004="Split App",  S005="Short Flood",
  S006="AWD",        S008="No-Till",    S009="Deep Plow",
  S010="Cont.Irrig", S011="Rainfed",    S012="Low Residue",
  S013="Hi Residue"
)

SCENE_COLS <- c(
  S000="#2C3E50", S001="#1565C0", S002="#42A5F5", S003="#00BCD4",
  S004="#80DEEA", S005="#2E7D32", S006="#66BB6A", S008="#E65100",
  S009="#FF8F00", S010="#AB47BC", S011="#CE93D8", S012="#C62828",
  S013="#EF9A9A"
)

CAT_COLS <- c(
  baseline="#2C3E50", fertilizer="#1565C0", water="#2E7D32",
  tillage="#E65100",  irrigation="#AB47BC",  residue="#C62828"
)

# Month x-axis within monsoon window
MON_BREAKS <- c(166, 182, 213, 244, 274, 305, 309)
MON_LABELS <- c("Jun 15","Jul","Aug","Sep","Oct","Nov","Nov 5")

# =============================================================================
# 2. SHARED HELPERS
# =============================================================================

theme_ch4 <- function(b = 11) {
  theme_bw(base_size = b) +
    theme(plot.title = element_text(face = "bold", size = b + 2),
          plot.subtitle = element_text(size = b - 1, colour = "grey40"),
          plot.caption  = element_text(size = b - 2, colour = "grey50", hjust = 0),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_line(colour = "grey92"),
          strip.background = element_rect(fill = "grey93", colour = "grey70"),
          strip.text = element_text(face = "bold", size = b - 1))
}

theme_map <- function(base = 10) {
  b <- base
  theme_bw(base_size = b) +
    theme(plot.title = element_text(face = "bold", size = b + 2),
          plot.subtitle = element_text(size = b - 1, colour = "grey40"),
          plot.caption  = element_text(size = b - 2, colour = "grey50", hjust = 0),
          panel.grid.major = element_line(colour = "grey88", linewidth = 0.2),
          panel.grid.minor = element_blank(),
          legend.key.width = unit(0.5, "cm"),
          axis.text  = element_text(size = b - 2),
          axis.title = element_blank(),
          strip.background = element_rect(fill = "grey90", colour = "grey70"),
          strip.text = element_text(face = "bold", size = b - 2))
}

save_png <- function(p, name, dir = FIGURE_DIR, w = 14, h = 8, dpi = 300) {
  path <- file.path(dir, name)
  ggsave(path, p, width = w, height = h, dpi = dpi, device = "png", bg = "white")
  message(sprintf("  [PNG] %s  (%d KB)", name,
                  round(file.info(path)$size / 1024)))
  invisible(path)
}

ch4_fill_scale <- function(name = "Monsoon CH4\n(kg C/ha)", ...) {
  scale_fill_gradientn(
    colours = c("#FFFDE7","#FFF176","#FFB300","#E65100","#B71C1C"),
    name = name, labels = label_number(accuracy = 1), ...)
}

mon_x <- function() {
  scale_x_continuous(breaks = MON_BREAKS, labels = MON_LABELS, expand = c(0.01, 0))
}

# =============================================================================
# 3. LOAD & AGGREGATE — MONSOON ONLY
# =============================================================================

message("Loading panel ...")
panel <- readRDS(PANEL_RDS)
setnames(panel,
         old = c("CH4-flux","CH4-prod.","CH4-oxid.","CH4-pool"),
         new = c("CH4_flux","CH4_prod","CH4_oxid","CH4_pool"),
         skip_absent = TRUE)

stopifnot("CH4_flux missing" = "CH4_flux" %in% names(panel),
          "'x' missing"      = "x"        %in% names(panel),
          "'y' missing"      = "y"        %in% names(panel))

panel <- panel[Day >= MONSOON[1] & Day <= MONSOON[2]]
message(sprintf("  Monsoon panel: %s rows", format(nrow(panel), big.mark = ",")))

# Optional CH4_prod / CH4_oxid guard
if (!"CH4_prod" %in% names(panel)) panel[, CH4_prod := NA_real_]
if (!"CH4_oxid" %in% names(panel)) panel[, CH4_oxid := NA_real_]

# ── Agg A: scene × DOY ────────────────────────────────────────────────────────
message("Aggregating: scene × Day ...")
agg <- panel[, .(
  ch4_med  = median(CH4_flux, na.rm = TRUE),
  ch4_mean = mean(  CH4_flux, na.rm = TRUE),
  ch4_q25  = quantile(CH4_flux, 0.25, na.rm = TRUE),
  ch4_q75  = quantile(CH4_flux, 0.75, na.rm = TRUE),
  ch4_q10  = quantile(CH4_flux, 0.10, na.rm = TRUE),
  ch4_q90  = quantile(CH4_flux, 0.90, na.rm = TRUE),
  prod_med  = median(CH4_prod, na.rm = TRUE),
  oxid_med  = median(CH4_oxid, na.rm = TRUE)
), by = .(scene_id, scene_category, Day)]

agg[, scene_id    := factor(scene_id, levels = SCENE_ORDER)]
agg[, scene_short := factor(SCENE_SHORT[as.character(scene_id)],
                            levels = SCENE_SHORT[SCENE_ORDER])]
agg[, scene_category := factor(scene_category,
                               levels = c("baseline","fertilizer","water",
                                          "tillage","irrigation","residue"))]

# ── Agg B: scene × year × DOY (3 focal scenarios only) ────────────────────────
message("Aggregating: scene × year × Day ...")
agg_yr <- panel[scene_id %in% c("S000","S006","S013"), .(
  ch4_med = median(CH4_flux, na.rm = TRUE),
  ch4_q25 = quantile(CH4_flux, 0.25, na.rm = TRUE),
  ch4_q75 = quantile(CH4_flux, 0.75, na.rm = TRUE)
), by = .(scene_id, cal_year, Day)]

# ── Agg C: cell-level seasonal totals (for maps) ──────────────────────────────
message("Aggregating: cell × scene × year totals ...")
cell_ann <- panel[, .(
  ann_CH4 = sum(CH4_flux, na.rm = TRUE)
), by = .(cell_id, x, y, scene_id, cal_year)]

map_data <- cell_ann[, .(
  mean_CH4 = mean(ann_CH4, na.rm = TRUE),
  sd_CH4   = sd(ann_CH4,   na.rm = TRUE),
  cv_CH4   = sd(ann_CH4, na.rm = TRUE) / mean(ann_CH4, na.rm = TRUE) * 100,
  n_years  = .N
), by = .(cell_id, x, y, scene_id)]

map_data[, scene_id := as.character(scene_id)]
map_data[, scene_short := factor(SCENE_SHORT[scene_id], levels = SCENE_SHORT[SCENE_ORDER])]

# Delta vs S000
base_map <- map_data[scene_id == "S000", .(cell_id, x, y, base_CH4 = mean_CH4)]
map_data  <- merge(map_data, base_map, by = c("cell_id","x","y"), all.x = TRUE)
map_data[, delta_CH4 := mean_CH4 - base_CH4]
map_data[, pct_change := delta_CH4 / base_CH4 * 100]

rm(panel, cell_ann); gc()
message("  Panel freed. Aggregation complete.\n")

# =============================================================================
# 4. MYANMAR BORDER
# =============================================================================

mmr_border <- tryCatch({
  if (!requireNamespace("maps", quietly = TRUE))
    install.packages("maps", repos = "https://cloud.r-project.org")
  library(maps)
  bd <- map_data("world", region = "Myanmar")
  if (!nrow(bd)) map_data("world", region = "Burma") else bd
}, error = function(e) {
  message("  WARN: Myanmar border unavailable — maps will lack outline")
  data.frame(long = numeric(0), lat = numeric(0), group = integer(0))
})

x_range <- range(map_data$x, na.rm = TRUE)
y_range <- range(map_data$y, na.rm = TRUE)
xlim    <- c(max(92.0, x_range[1] - 1.5), min(102.0, x_range[2] + 1.5))
ylim    <- c(max( 9.5, y_range[1] - 1.0), min( 29.0, y_range[2] + 1.0))

base_map_layers <- function() list(
  geom_polygon(data = mmr_border, aes(x = long, y = lat, group = group),
               fill = "grey95", colour = "grey50", linewidth = 0.35),
  coord_fixed(ratio = 1, xlim = xlim, ylim = ylim, expand = FALSE),
  scale_x_continuous(breaks = seq(92, 102, 2),
                     labels = paste0(seq(92, 102, 2), "°E")),
  scale_y_continuous(breaks = seq(10,  30,  4),
                     labels = paste0(seq(10,  30, 4), "°N"))
)

# =============================================================================
# F1 — All scenarios, monsoon daily CH4, single panel
# =============================================================================
message("F1: All-scenario monsoon flux ...")

p1 <- ggplot(agg, aes(x = Day, colour = scene_id, fill = scene_id)) +
  geom_ribbon(aes(ymin = ch4_q25, ymax = ch4_q75), alpha = 0.10, colour = NA) +
  geom_line(aes(y = ch4_med), linewidth = 0.75) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey60") +
  scale_colour_manual(values = SCENE_COLS, labels = SCENE_SHORT, name = "Scenario") +
  scale_fill_manual(  values = SCENE_COLS, labels = SCENE_SHORT, name = "Scenario") +
  mon_x() +
  scale_y_continuous(labels = label_number(accuracy = 0.01)) +
  labs(title    = "Monsoon Daily CH4 Flux — All 13 Scenarios",
       subtitle = "Median across 263 cells × 20 years. Ribbons = IQR (25th–75th pct). DOY 166–309.",
       x = "Day of year", y = "CH4 flux (kg C/ha/day)",
       caption = "Source: DNDC 9.5 batch simulation | Myanmar irrigated rice") +
  theme_ch4() + theme(legend.position = "right")

save_png(p1, "F1_monsoon_CH4_all_scenarios.png")

# =============================================================================
# F2 — Faceted by category, baseline reference in each facet
# =============================================================================
message("F2: Faceted by category ...")

baseline_ref <- agg[scene_id == "S000", .(Day, ch4_baseline = ch4_med)]

p2 <- ggplot(agg, aes(x = Day, colour = scene_id, fill = scene_id)) +
  geom_line(data = baseline_ref, aes(x = Day, y = ch4_baseline),
            inherit.aes = FALSE, colour = "grey50", linewidth = 0.5, linetype = "dashed") +
  geom_ribbon(aes(ymin = ch4_q25, ymax = ch4_q75), alpha = 0.18, colour = NA) +
  geom_line(aes(y = ch4_med), linewidth = 0.85) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey70") +
  scale_colour_manual(values = SCENE_COLS, labels = SCENE_SHORT, name = "Scenario") +
  scale_fill_manual(  values = SCENE_COLS, labels = SCENE_SHORT, name = "Scenario") +
  mon_x() +
  facet_wrap(~ scene_category, ncol = 3, scales = "free_y") +
  labs(title    = "Monsoon CH4 Flux by Management Category",
       subtitle = "Dashed grey = Baseline (S000). Ribbons = IQR. Scales free within category.",
       x = "Day of year", y = "CH4 flux (kg C/ha/day)",
       caption = "Source: DNDC 9.5 batch simulation | Myanmar irrigated rice") +
  theme_ch4() + theme(legend.position = "bottom")

save_png(p2, "F2_monsoon_CH4_by_category.png", w = 15, h = 10)

# =============================================================================
# F3 — CH4 components: production / oxidation / net (Baseline, AWD, Hi-Residue)
# =============================================================================
message("F3: CH4 components ...")

focal <- c(S000 = "Baseline (S000)", S006 = "AWD (S006)", S013 = "High Residue (S013)")
comp_long <- melt(
  agg[scene_id %in% names(focal),
      .(scene_id, Day,
        "Net flux"   = ch4_med,
        "Production" = prod_med,
        "Oxidation"  = oxid_med)],
  id.vars = c("scene_id","Day"), variable.name = "component", value.name = "ch4"
)
comp_long[, scene_label := focal[as.character(scene_id)]]
comp_long[, component   := factor(component, levels = c("Production","Oxidation","Net flux"))]

p3 <- ggplot(comp_long, aes(x = Day, y = ch4, colour = component,
                            linetype = component, linewidth = component)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey60") +
  geom_line() +
  scale_colour_manual(  values = c(Production="#C62828", Oxidation="#1565C0", "Net flux"="#2C3E50"), name = NULL) +
  scale_linetype_manual(values = c(Production="solid",   Oxidation="dashed",  "Net flux"="solid"),  name = NULL) +
  scale_linewidth_manual(values = c(Production=0.9,      Oxidation=0.7,       "Net flux"=1.2),      name = NULL) +
  mon_x() +
  facet_wrap(~ scene_label, ncol = 3) +
  labs(title    = "CH4 Production, Oxidation, and Net Flux — Monsoon",
       subtitle = "Red solid = gross production | Blue dashed = oxidation | Dark = net flux",
       x = "Day of year", y = "CH4 (kg C/ha/day)",
       caption = "Source: DNDC 9.5 | 263 cells × 20 years median") +
  theme_ch4() + theme(legend.position = "bottom")

save_png(p3, "F3_monsoon_CH4_components.png", w = 15, h = 6)

# =============================================================================
# F4 — Interannual variability 2005-2024 (3 scenarios)
# =============================================================================
message("F4: Interannual variability ...")

YR_PAL <- colorRampPalette(c("#1565C0","#66BB6A","#C62828"))(20)
names(YR_PAL) <- as.character(2005:2024)
agg_yr[, scene_label := focal[as.character(scene_id)]]

p4 <- ggplot(agg_yr, aes(x = Day, colour = factor(cal_year))) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey70") +
  geom_line(aes(y = ch4_med), linewidth = 0.5, alpha = 0.85) +
  scale_colour_manual(values = YR_PAL, name = "Year") +
  mon_x() +
  facet_wrap(~ scene_label, ncol = 3, scales = "free_y") +
  labs(title    = "Interannual Variability of Monsoon CH4 (2005-2024)",
       subtitle = "Each line = one year's median across 263 cells. Blue→green→red: early→mid→late years.",
       x = "Day of year", y = "CH4 flux (kg C/ha/day)",
       caption = "Source: DNDC 9.5 | Baseline, AWD, High Residue") +
  theme_ch4() + theme(legend.key.height = unit(0.35, "cm"),
                      legend.text = element_text(size = 8))

save_png(p4, "F4_monsoon_CH4_interannual.png", w = 16, h = 6)

# =============================================================================
# F5 — Heatmap: scenario × DOY (7-day rolling mean)
# =============================================================================
message("F5: Heatmap scene × DOY ...")

agg_sm <- copy(agg)
agg_sm[, ch4_7d := frollmean(ch4_mean, 7, align = "center", na.rm = TRUE), by = scene_id]

# Order scenarios by peak monsoon emission
peak_ord <- agg[, .(peak = max(ch4_med, na.rm = TRUE)), by = scene_id][order(-peak)]$scene_id
agg_sm[, scene_ord := factor(
  scene_id,
  levels = rev(as.character(peak_ord)),
  labels = rev(paste0(as.character(peak_ord), "  ", SCENE_SHORT[as.character(peak_ord)]))
)]

p5 <- ggplot(agg_sm, aes(x = Day, y = scene_ord, fill = ch4_7d)) +
  geom_tile() +
  scale_fill_gradient2(low = "#2980B9", mid = "#F9E79F", high = "#C0392B",
                       midpoint = 0.05, na.value = "grey85",
                       name = "Mean CH4\n(kg C/ha/d)",
                       labels = label_number(accuracy = 0.01)) +
  mon_x() +
  labs(title    = "CH4 Flux Heatmap — Scenario × Monsoon Day",
       subtitle = "7-day rolling mean | Scenarios ranked by peak emission (high→low top→bottom)",
       x = "Day of year", y = NULL,
       caption = "Source: DNDC 9.5 | 263 cells × 20 years") +
  theme_ch4() +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(colour = "grey60", fill = NA),
        axis.text.y = element_text(size = 9))

save_png(p5, "F5_monsoon_CH4_heatmap.png", w = 13, h = 7)

# =============================================================================
# F6 — All scenarios with flooding window overlay
# =============================================================================
message("F6: CH4 flux with flooding windows ...")

flood_ann <- data.frame(
  xmin = c(185, 247), xmax = c(226, 290),
  label = c("Flood period 1", "Flood period 2")
)

p6 <- ggplot(agg, aes(x = Day, colour = scene_id, fill = scene_id)) +
  annotate("rect", xmin = flood_ann$xmin, xmax = flood_ann$xmax,
           ymin = -Inf, ymax = Inf, alpha = 0.08, fill = "#1565C0") +
  annotate("text", x = (flood_ann$xmin + flood_ann$xmax) / 2, y = Inf,
           label = flood_ann$label, vjust = 1.4, size = 2.8,
           colour = "#1565C0", fontface = "italic") +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey60") +
  geom_ribbon(aes(ymin = ch4_q10, ymax = ch4_q90), alpha = 0.05, colour = NA) +
  geom_ribbon(aes(ymin = ch4_q25, ymax = ch4_q75), alpha = 0.12, colour = NA) +
  geom_line(aes(y = ch4_med), linewidth = 0.75) +
  scale_colour_manual(values = SCENE_COLS, labels = SCENE_SHORT, name = "Scenario") +
  scale_fill_manual(  values = SCENE_COLS, labels = SCENE_SHORT, name = "Scenario") +
  mon_x() +
  labs(title    = "Monsoon CH4 Flux with Flooding Windows",
       subtitle = "Blue shading = monsoon flooding periods (S000 baseline). IQR + 10-90th pct ribbons.",
       x = "Day of year", y = "CH4 flux (kg C/ha/day)",
       caption = "Flooding activates methanogenesis; AWD reduces emission within same window | DNDC 9.5") +
  theme_ch4() + theme(legend.position = "right")

save_png(p6, "F6_monsoon_CH4_flooding_overlay.png", w = 16, h = 8)

# =============================================================================
# F7 — Cumulative monsoon CH4 by scenario (ranked)
# =============================================================================
message("F7: Cumulative monsoon CH4 ...")

agg_cum <- copy(agg)[order(scene_id, Day)]
agg_cum[, cum_ch4 := cumsum(ch4_med), by = scene_id]

# Final totals for ranking
finals <- agg_cum[, .(total = max(cum_ch4, na.rm = TRUE)), by = .(scene_id, scene_short)]
setorder(finals, -total)
agg_cum[, scene_short := factor(scene_short, levels = finals$scene_short)]

p7 <- ggplot(agg_cum, aes(x = Day, y = cum_ch4, colour = scene_id)) +
  geom_line(linewidth = 0.9) +
  geom_text(data = finals,
            aes(x = 310, y = total, label = sprintf("%s  %.0f", scene_short, total),
                colour = scene_id),
            hjust = 0, size = 2.8, inherit.aes = FALSE, show.legend = FALSE) +
  scale_colour_manual(values = SCENE_COLS, labels = SCENE_SHORT, name = "Scenario") +
  scale_x_continuous(breaks = MON_BREAKS, labels = MON_LABELS,
                     expand = c(0.01, 0), limits = c(166, 360)) +
  scale_y_continuous(labels = label_number(accuracy = 1)) +
  labs(title    = "Cumulative Monsoon CH4 Emission by Scenario",
       subtitle = "Season-end totals (kg C/ha) labelled at right. Scenarios diverge with management intensity.",
       x = "Day of year", y = "Cumulative CH4 (kg C/ha)",
       caption = "Median across 263 cells × 20 years | Source: DNDC 9.5") +
  theme_ch4() + theme(legend.position = "none")

save_png(p7, "F7_monsoon_CH4_cumulative_ranked.png")

# =============================================================================
# M1 — Baseline monsoon CH4 tile map
# =============================================================================
message("\nM1: Baseline monsoon CH4 map ...")

m1 <- ggplot() +
  base_map_layers() +
  geom_tile(data = map_data[scene_id == "S000"],
            aes(x = x, y = y, fill = mean_CH4),
            width = TILE_SIZE, height = TILE_SIZE) +
  ch4_fill_scale(name = "Mean monsoon\nCH4 (kg C/ha)") +
  geom_polygon(data = mmr_border, aes(x = long, y = lat, group = group),
               fill = NA, colour = "grey30", linewidth = 0.5) +
  labs(title    = "Mean Monsoon CH4 — Myanmar Rice Paddies",
       subtitle = "Baseline (S000) | 20-year mean (2005-2024) | DOY 166–309 | 0.5° resolution",
       caption  = "Source: DNDC 9.5 | SoilGrids | NASA POWER") +
  theme_map()

save_png(m1, "M1_monsoon_CH4_baseline.png", dir = MAP_DIR, w = 9, h = 12)

# =============================================================================
# M2 — All 13 scenarios small-multiple tile maps
# =============================================================================
message("M2: All-scenario monsoon maps ...")

ch4_lim <- range(map_data$mean_CH4, na.rm = TRUE)
m2_data <- map_data[scene_id %in% SCENE_ORDER]
m2_data[, scene_short := factor(SCENE_SHORT[scene_id], levels = SCENE_SHORT[SCENE_ORDER])]

m2 <- ggplot() +
  geom_polygon(data = mmr_border, aes(x = long, y = lat, group = group),
               fill = "grey93", colour = "grey55", linewidth = 0.25) +
  geom_tile(data = m2_data, aes(x = x, y = y, fill = mean_CH4),
            width = TILE_SIZE, height = TILE_SIZE) +
  ch4_fill_scale(limits = ch4_lim) +
  geom_polygon(data = mmr_border, aes(x = long, y = lat, group = group),
               fill = NA, colour = "grey30", linewidth = 0.3) +
  coord_fixed(ratio = 1, xlim = xlim, ylim = ylim, expand = FALSE) +
  facet_wrap(~ scene_short, ncol = 4) +
  scale_x_continuous(breaks = seq(94, 100, 3), labels = paste0(seq(94, 100, 3), "°E")) +
  scale_y_continuous(breaks = seq(12, 28,  4), labels = paste0(seq(12, 28,  4), "°N")) +
  labs(title    = "Mean Monsoon CH4 by Scenario — Myanmar Rice Paddies",
       subtitle = "20-year mean | DOY 166–309 | Same colour scale across all panels",
       caption  = "Source: DNDC 9.5 batch simulation") +
  theme_map(base = 8) +
  theme(legend.position = "right", strip.text = element_text(size = 7, face = "bold"),
        axis.text = element_text(size = 6))

save_png(m2, "M2_monsoon_CH4_all_scenarios.png", dir = MAP_DIR, w = 18, h = 20)

# =============================================================================
# M3 — Delta vs baseline (diverging)
# =============================================================================
message("M3: Monsoon delta vs baseline ...")

m3_data <- map_data[scene_id != "S000" & scene_id %in% SCENE_ORDER]
m3_data[, scene_short := factor(SCENE_SHORT[scene_id],
                                levels = SCENE_SHORT[SCENE_ORDER[SCENE_ORDER != "S000"]])]
dlim <- max(abs(m3_data$delta_CH4), na.rm = TRUE)

m3 <- ggplot() +
  geom_polygon(data = mmr_border, aes(x = long, y = lat, group = group),
               fill = "grey93", colour = "grey55", linewidth = 0.25) +
  geom_tile(data = m3_data, aes(x = x, y = y, fill = delta_CH4),
            width = TILE_SIZE, height = TILE_SIZE) +
  scale_fill_gradient2(low = "#1565C0", mid = "white", high = "#C62828",
                       midpoint = 0, limits = c(-dlim, dlim),
                       name = "ΔMonsoon CH4\n(kg C/ha)",
                       labels = label_number(accuracy = 1)) +
  geom_polygon(data = mmr_border, aes(x = long, y = lat, group = group),
               fill = NA, colour = "grey30", linewidth = 0.3) +
  coord_fixed(ratio = 1, xlim = xlim, ylim = ylim, expand = FALSE) +
  facet_wrap(~ scene_short, ncol = 4) +
  scale_x_continuous(breaks = seq(94, 100, 3), labels = paste0(seq(94, 100, 3), "°E")) +
  scale_y_continuous(breaks = seq(12, 28,  4), labels = paste0(seq(12, 28,  4), "°N")) +
  labs(title    = "Monsoon CH4 Change vs Baseline by Scenario",
       subtitle = "Blue = reduction | Red = increase | 20-year mean | DOY 166–309",
       caption  = "Source: DNDC 9.5 | Myanmar irrigated rice") +
  theme_map(base = 8) +
  theme(legend.position = "right", strip.text = element_text(size = 7, face = "bold"),
        axis.text = element_text(size = 6))

save_png(m3, "M3_monsoon_CH4_delta_baseline.png", dir = MAP_DIR, w = 18, h = 20)

# =============================================================================
# M4 — AWD mitigation potential (S000 - S006)
# =============================================================================
message("M4: AWD mitigation ...")

awd <- merge(
  map_data[scene_id == "S000", .(cell_id, x, y, ch4_base = mean_CH4)],
  map_data[scene_id == "S006", .(cell_id, ch4_awd = mean_CH4)],
  by = "cell_id"
)
awd[, abs_red := ch4_base - ch4_awd]
awd[, pct_red := abs_red / ch4_base * 100]

m4_abs <- ggplot() + base_map_layers() +
  geom_tile(data = awd, aes(x = x, y = y, fill = abs_red),
            width = TILE_SIZE, height = TILE_SIZE) +
  scale_fill_gradientn(colours = c("#E3F2FD","#64B5F6","#1565C0","#0D47A1"),
                       name = "CH4 reduction\n(kg C/ha)",
                       labels = label_number(accuracy = 1)) +
  geom_polygon(data = mmr_border, aes(x = long, y = lat, group = group),
               fill = NA, colour = "grey30", linewidth = 0.5) +
  labs(title = "AWD — Absolute CH4 Reduction", subtitle = "S000 − S006 | Monsoon") +
  theme_map()

m4_pct <- ggplot() + base_map_layers() +
  geom_tile(data = awd, aes(x = x, y = y, fill = pct_red),
            width = TILE_SIZE, height = TILE_SIZE) +
  scale_fill_gradientn(colours = c("#E8F5E9","#81C784","#2E7D32","#1B5E20"),
                       name = "CH4 reduction\n(%)",
                       labels = label_percent(scale = 1, accuracy = 1)) +
  geom_polygon(data = mmr_border, aes(x = long, y = lat, group = group),
               fill = NA, colour = "grey30", linewidth = 0.5) +
  labs(title = "AWD — % CH4 Reduction", subtitle = "Darker = greater reduction") +
  theme_map()

m4 <- m4_abs + m4_pct +
  plot_annotation(
    title   = "AWD Mitigation Potential — Monsoon Season",
    caption = "Left: absolute (kg C/ha) | Right: % vs baseline | Source: DNDC 9.5",
    theme   = theme(plot.title = element_text(face = "bold", size = 14))
  )

save_png(m4, "M4_AWD_mitigation_monsoon.png", dir = MAP_DIR, w = 18, h = 12)

# =============================================================================
# M5 — Coefficient of variation (interannual variability) — baseline
# =============================================================================
message("M5: Interannual CV map ...")

m5 <- ggplot() + base_map_layers() +
  geom_tile(data = map_data[scene_id == "S000"],
            aes(x = x, y = y, fill = cv_CH4),
            width = TILE_SIZE, height = TILE_SIZE) +
  scale_fill_gradientn(
    colours = c("#F3E5F5","#CE93D8","#7B1FA2","#4A148C"),
    name = "CV (%)\ninterannual",
    labels = label_number(accuracy = 1)
  ) +
  geom_polygon(data = mmr_border, aes(x = long, y = lat, group = group),
               fill = NA, colour = "grey30", linewidth = 0.5) +
  labs(title    = "Interannual Variability of Monsoon CH4 — Baseline",
       subtitle = "Coefficient of variation (SD/mean × 100%) | 2005-2024 | Higher = more variable",
       caption  = "Source: DNDC 9.5 | Variability driven by climate variability") +
  theme_map()

save_png(m5, "M5_monsoon_CH4_CV_baseline.png", dir = MAP_DIR, w = 9, h = 12)

# =============================================================================
# M6 — Optimal (lowest-emitting) scenario per cell
# =============================================================================
message("M6: Optimal scenario map ...")

opt <- map_data[, .(
  best_scene = scene_id[which.min(mean_CH4)],
  best_CH4   = min(mean_CH4, na.rm = TRUE),
  base_CH4   = mean_CH4[scene_id == "S000"],
  max_red_pct= (mean_CH4[scene_id == "S000"] - min(mean_CH4, na.rm = TRUE)) /
    mean_CH4[scene_id == "S000"] * 100
), by = .(cell_id, x, y)]

opt[, scene_cat := fcase(
  best_scene == "S000", "baseline",
  best_scene %in% c("S001","S002","S003","S004"), "fertilizer",
  best_scene %in% c("S005","S006"), "water",
  best_scene %in% c("S008","S009"), "tillage",
  best_scene %in% c("S010","S011"), "irrigation",
  best_scene %in% c("S012","S013"), "residue",
  default = "baseline"
)]
opt[, scene_cat := factor(scene_cat, levels = names(CAT_COLS))]

m6_cat <- ggplot() + base_map_layers() +
  geom_tile(data = opt, aes(x = x, y = y, fill = scene_cat),
            width = TILE_SIZE, height = TILE_SIZE) +
  scale_fill_manual(values = CAT_COLS, name = "Best category", drop = FALSE) +
  geom_polygon(data = mmr_border, aes(x = long, y = lat, group = group),
               fill = NA, colour = "grey30", linewidth = 0.5) +
  labs(title = "Best Management Category", subtitle = "Category of lowest-emitting scenario per cell") +
  theme_map()

m6_pct <- ggplot() + base_map_layers() +
  geom_tile(data = opt, aes(x = x, y = y, fill = max_red_pct),
            width = TILE_SIZE, height = TILE_SIZE) +
  scale_fill_gradientn(
    colours = c("#FFF8E1","#FFCC02","#E65100","#B71C1C"),
    name = "Max reduction\n(%)", labels = label_number(accuracy = 1, suffix = "%")
  ) +
  geom_polygon(data = mmr_border, aes(x = long, y = lat, group = group),
               fill = NA, colour = "grey30", linewidth = 0.5) +
  labs(title = "Maximum Achievable CH4 Reduction", subtitle = "% vs baseline | Best scenario per cell") +
  theme_map()

m6 <- m6_cat + m6_pct +
  plot_annotation(
    title   = "Optimal Management per Rice Cell — Monsoon Season",
    caption = "Left: category | Right: % reduction vs S000 | Source: DNDC 9.5",
    theme   = theme(plot.title = element_text(face = "bold", size = 14),
                    plot.caption = element_text(size = 9, colour = "grey50"))
  )

save_png(m6, "M6_monsoon_optimal_scenario.png", dir = MAP_DIR, w = 18, h = 12)

# =============================================================================
# SUMMARY
# =============================================================================

all_png <- c(
  list.files(FIGURE_DIR, pattern = "\\.png$", full.names = TRUE),
  list.files(MAP_DIR,    pattern = "\\.png$", full.names = TRUE)
)

message("\n", strrep("=", 65))
message("COMPLETE — Monsoon figures saved")
message(strrep("=", 65))
message(sprintf("Run command    : Rscript 9_10_monsoon_figures.R"))
message(sprintf("Figures total  : %d PNG files\n", length(all_png)))
message("Flux figures -> ", FIGURE_DIR)
for (f in list.files(FIGURE_DIR, pattern = "\\.png$")) message(sprintf("  %s", f))
message("Spatial maps -> ", MAP_DIR)
for (f in list.files(MAP_DIR,    pattern = "\\.png$")) message(sprintf("  %s", f))
message(strrep("=", 65))