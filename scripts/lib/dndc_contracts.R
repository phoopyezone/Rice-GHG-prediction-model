# Shared contracts and helpers for Myanmar/global DNDC pipeline.

dndc_get_sim_root <- function() {
  Sys.getenv(
    "DNDC_SIM_DIR",
    unset = "G:/My Drive/Research/simulation/main_dndc/simulate_dndc"
  )
}

dndc_is_leap_year <- function(year) {
  (year %% 4 == 0 & year %% 100 != 0) | (year %% 400 == 0)
}

dndc_make_cell_id <- function(lon, lat, tile_id = 0L, prefix = "CELL") {
  lon_tag <- sprintf("%06d", as.integer(round((lon + 180) * 100)))
  lat_tag <- sprintf("%05d", as.integer(round((lat + 90) * 100)))
  sprintf("%s_t%03d_%s_%s", prefix, as.integer(tile_id), lon_tag, lat_tag)
}

dndc_locked_climate_schema <- function() {
  data.frame(
    field_order = 1:7,
    field_name = c("julian_day", "tmax_c", "tmin_c", "precip_cm", "wind_mps", "radiation_mj_m2_day", "humidity_pct"),
    dndc_format5_column = c("Julian day", "Tmax", "Tmin", "Precip", "Wind", "Radiation", "Humidity"),
    required = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE),
    units = c("day", "degC", "degC", "cm/day", "m/s", "MJ/m2/day", "percent"),
    stringsAsFactors = FALSE
  )
}

dndc_locked_conversion_map <- function() {
  data.frame(
    variable_id = c(
      "Daily_max_air_temperature",
      "Daily_min_air_temperature",
      "Daily_precipitation",
      "Daily_avg_wind_speed",
      "Daily_solar_radiation",
      "Daily_relative_humidity",
      "Clay_fraction",
      "Top_layer_SOC",
      "Bulk_density",
      "pH"
    ),
    source_units = c("degC", "degC", "mm/day", "m/s", "kWh/m2/day", "percent", "percent", "g/kg", "kg/dm3", "unitless"),
    target_units = c("degC", "degC", "cm/day", "m/s", "MJ/m2/day", "percent", "fraction", "fraction", "g/cm3", "unitless"),
    formula = c("x", "x", "x/10", "x", "x*3.6", "x", "x/100", "x/1000", "x", "x"),
    stringsAsFactors = FALSE
  )
}

dndc_manifest_schema <- function() {
  data.frame(
    manifest_name = c(
      rep("manifest_climate", 8),
      rep("manifest_soil", 7),
      rep("manifest_cells", 7)
    ),
    column_name = c(
      "scope", "tile_id", "year", "variable", "path", "rows", "created_utc", "source",
      "scope", "tile_id", "variable", "path", "rows", "created_utc", "source",
      "scope", "tile_id", "cell_id", "lon", "lat", "path", "created_utc"
    ),
    column_type = c(
      "character", "integer", "integer", "character", "character", "integer", "character", "character",
      "character", "integer", "character", "character", "integer", "character", "character",
      "character", "integer", "character", "numeric", "numeric", "character", "character"
    ),
    required = c(
      TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
      TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
      TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE
    ),
    stringsAsFactors = FALSE
  )
}

dndc_apply_conversion <- function(variable_id, x) {
  if (variable_id == "Daily_precipitation") return(x / 10)
  if (variable_id == "Daily_solar_radiation") return(x * 3.6)
  if (variable_id == "Clay_fraction") return(x / 100)
  if (variable_id == "Top_layer_SOC") return(x / 1000)
  x
}

dndc_write_table <- function(df, path_without_ext) {
  if (requireNamespace("arrow", quietly = TRUE)) {
    out <- paste0(path_without_ext, ".parquet")
    arrow::write_parquet(df, out)
    return(out)
  }
  out <- paste0(path_without_ext, ".csv")
  utils::write.csv(df, out, row.names = FALSE, na = "")
  out
}

dndc_read_table <- function(path_without_ext) {
  pq <- paste0(path_without_ext, ".parquet")
  cs <- paste0(path_without_ext, ".csv")
  if (file.exists(pq) && requireNamespace("arrow", quietly = TRUE)) {
    return(as.data.frame(arrow::read_parquet(pq)))
  }
  if (file.exists(cs)) {
    return(utils::read.csv(cs, stringsAsFactors = FALSE))
  }
  stop("Table not found for base path: ", path_without_ext)
}

dndc_write_locked_contracts <- function(sim_root = dndc_get_sim_root(),
                                        registry_csv = "G:/My Drive/Research/simulation/main_dndc/scripts/ai agent/dndc_variable_gpt.csv",
                                        guide_txt = "G:/My Drive/Research/simulation/main_dndc/scripts/ai agent/GuideDNDC95.txt") {
  contracts_dir <- file.path(sim_root, "contracts")
  dir.create(contracts_dir, recursive = TRUE, showWarnings = FALSE)

  climate_schema <- dndc_locked_climate_schema()
  conversion_map <- dndc_locked_conversion_map()
  manifest_schema <- dndc_manifest_schema()

  utils::write.csv(climate_schema, file.path(contracts_dir, "dndc_climate_format5_schema.csv"), row.names = FALSE)
  utils::write.csv(conversion_map, file.path(contracts_dir, "dndc_conversion_map.csv"), row.names = FALSE)
  utils::write.csv(manifest_schema, file.path(contracts_dir, "dndc_manifest_schema.csv"), row.names = FALSE)

  registry_rows <- data.frame()
  if (file.exists(registry_csv)) {
    reg <- utils::read.csv(registry_csv, stringsAsFactors = FALSE, fileEncoding = "UTF-8")
    if ("Selection" %in% names(reg)) {
      keep <- tolower(trimws(reg$Selection)) == "yes"
      registry_rows <- reg[keep, , drop = FALSE]
      utils::write.csv(registry_rows, file.path(contracts_dir, "registry_selected_yes.csv"), row.names = FALSE)
    }
  }

  guide_digest <- if (file.exists(guide_txt)) tools::md5sum(guide_txt) else NA_character_
  registry_digest <- if (file.exists(registry_csv)) tools::md5sum(registry_csv) else NA_character_
  lock_meta <- data.frame(
    lock_timestamp_utc = format(as.POSIXct(Sys.time(), tz = "UTC"), "%Y-%m-%dT%H:%M:%SZ"),
    guide_file = normalizePath(guide_txt, winslash = "/", mustWork = FALSE),
    guide_md5 = as.character(guide_digest),
    registry_file = normalizePath(registry_csv, winslash = "/", mustWork = FALSE),
    registry_md5 = as.character(registry_digest),
    selected_registry_rows = nrow(registry_rows),
    stringsAsFactors = FALSE
  )
  utils::write.csv(lock_meta, file.path(contracts_dir, "lock_metadata.csv"), row.names = FALSE)

  invisible(list(
    climate_schema = climate_schema,
    conversion_map = conversion_map,
    manifest_schema = manifest_schema,
    lock_metadata = lock_meta
  ))
}

dndc_now_utc <- function() {
  format(as.POSIXct(Sys.time(), tz = "UTC"), "%Y-%m-%dT%H:%M:%SZ")
}
