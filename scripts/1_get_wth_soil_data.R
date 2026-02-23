##### 1_dndc_get_wth_soil_data.R - DNDC VERSION #####
# Purpose : Download NASA POWER climate data + SoilGrids soil data for Myanmar,
#           then identify rice-growing grid cells for DNDC simulation input.
# Outputs : Per-variable NetCDF climate files, elevation.tif, soil_agg.tif, cells.rds
# Fix     : Rice-cell selection uses touches=FALSE + a minimum rice area threshold
#           (rice_min_ha) instead of the too-permissive original >0 approach.
# Resume  : All file writes are skipped if the output already exists — safe to rerun.

suppressPackageStartupMessages({library(nasapower);library(terra);library(dplyr);library(geodata);library(httr)})

# --- Paths ---
path <- "G:/My Drive/Research/simulation/main_dndc/simulate_dndc"
dir.create(path,FALSE,TRUE); setwd(path)
raw_dir <- file.path("data","raw"); power_dir <- file.path(raw_dir,"weather","power")
dir.create(power_dir,recursive=TRUE,showWarnings=FALSE)

# --- Config ---
years <- 2005:2024
ext   <- terra::ext(91.5,101.5,8,29)          # Myanmar bounding box (xmin,xmax,ymin,ymax)
rice_min_ha <- 1000   # minimum SPAM rice area (ha) required per 0.5° cell — key fix
touches_aoi <- FALSE  # exclude cells that only touch (not overlap) the Myanmar boundary — key fix

# Target grid: NASA POWER delivers data on a native 0.5° grid.
# lon_seq/lat_seq_desc define the expected row/column order used when filling raster values.
lon_seq      <- seq(xmin(ext),xmax(ext),by=0.5)
lat_seq_desc <- sort(seq(ymin(ext),ymax(ext),by=0.5),decreasing=TRUE) # N→S (raster row order)

# Helper: convert one day's POWER data.frame (LAT, LON, value) into a raster layer.
# Rows must be sorted N→S then W→E to match terra's fill order.
.make_layer <- function(df,vc){
  r <- terra::rast(ncols=length(lon_seq),nrows=length(lat_seq_desc),
                   xmin=xmin(ext),xmax=xmax(ext),ymin=ymin(ext),ymax=ymax(ext),crs="EPSG:4326")
  terra::values(r) <- (df|>dplyr::arrange(dplyr::desc(LAT),LON))[[vc]]; r}

# Helper: fetch one full calendar year for one POWER parameter using 4°×4° tiles.
# NASA POWER limits requests to ≤100 grid points AND ≥2° spatial extent per call.
# At 0.5° resolution, a 4°×4° tile = 8×8 = 64 points — safely within the 100-point cap.
# The last tile is dropped if it would be narrower than 2° (the API minimum).
# Retry logic: up to 5 attempts per tile, with escalating waits for timeouts / rate limits.
.fetch_year <- function(parm,yr){
  d1 <- paste0(yr,"-01-01"); d2 <- paste0(yr,"-12-31")
  httr::set_config(httr::timeout(300)); on.exit(httr::reset_config()) # 5-min HTTP timeout
  ts <- 4.0; ms <- 2.0                                               # tile size; minimum size
  lon0s <- seq(xmin(ext),xmax(ext),by=ts); lat0s <- seq(ymin(ext),ymax(ext),by=ts)
  # Drop last tile start if the remaining extent would fall below the 2° API minimum
  if(length(lon0s)>1&&(xmax(ext)-lon0s[length(lon0s)])<ms) lon0s <- lon0s[-length(lon0s)]
  if(length(lat0s)>1&&(ymax(ext)-lat0s[length(lat0s)])<ms) lat0s <- lat0s[-length(lat0s)]
  out <- list()
  for(lon0 in lon0s) for(lat0 in lat0s){
    lon1 <- min(lon0+ts,xmax(ext)); lat1 <- min(lat0+ts,ymax(ext))
    if((lon1-lon0)<ms||(lat1-lat0)<ms) next  # skip tiles narrower than API minimum
    for(k in 1:5){
      ok <- tryCatch({
        out[[length(out)+1]] <- nasapower::get_power(community="ag",lonlat=c(lon0,lat0,lon1,lat1),
                                                     pars=parm,dates=c(d1,d2),temporal_api="daily")
        Sys.sleep(3); TRUE  # polite delay between tiles to avoid rate limiting
      },error=function(e){msg <- as.character(e)
        if(grepl("Timeout|timed out",msg,ignore.case=TRUE)){Sys.sleep(30*k);FALSE}   # escalating wait
        else if(grepl("429|rate limit",msg,ignore.case=TRUE)){Sys.sleep(120*k);FALSE} # 2/4/6… min
        else{Sys.sleep(20);FALSE}})
      if(ok) break
    }
  }
  if(!length(out)) stop("No tiles downloaded for ",parm," ",yr)
  # Combine tiles, deduplicate overlapping edge points, then stack into a daily raster time-series
  dat <- dplyr::bind_rows(out)|>dplyr::rename(DATE=YYYYMMDD)
  dat$DATE <- as.Date(dat$DATE); dat <- dat[!duplicated(dat[,c("LAT","LON","DATE")]),]
  days <- sort(unique(dat$DATE))
  r <- terra::rast(lapply(seq_along(days),\(i).make_layer(dat[dat$DATE==days[i],c("LAT","LON",parm)],parm)))
  terra::time(r) <- days; names(r) <- paste0(parm,"_",format(days,"%Y%j")); r}

# Helper: load a variable's multi-year NetCDF if it already exists, otherwise download and save it.
# The consistent filename encodes parameter, year range, and spatial extent for traceability.
.get_or_load_var <- function(parm){
  fn <- file.path(power_dir,sprintf("%s-2005_2024-91.5x101.5x8x29.nc",parm))
  if(file.exists(fn)){message("Found existing ",basename(fn),"; loading..."); return(terra::rast(fn))}
  message("Downloading ",parm," ... (resume-safe; rerun if interrupted)")
  R <- terra::rast(lapply(years,\(yy).fetch_year(parm,yy))) # fetch year-by-year, then stack
  terra::writeCDF(R,fn,varname=parm,longname=parm,overwrite=TRUE); R}

# Helper: write a NetCDF only if it does not yet exist (idempotent), then report status.
write_if_missing <- function(x,fn,vn,ln,u){
  if(!file.exists(fn)) terra::writeCDF(x,fn,varname=vn,longname=ln,unit=u,overwrite=TRUE)
  message("✓ ",basename(fn),if(file.exists(fn))" (exists)"else" (written)")}

# =============================================================================
# 1) Download / load raw POWER climate variables
#    Downloads 20 years × 7 variables from NASA POWER (ag community).
#    Each variable is saved as a single multi-year NetCDF (~7 files total).
# =============================================================================
message("=== Downloading (or loading) DNDC Climate Variables ===")
R_TMAX<-.get_or_load_var("T2M_MAX")          # daily max air temperature (°C)
R_TMIN<-.get_or_load_var("T2M_MIN")          # daily min air temperature (°C)
R_TAVG<-.get_or_load_var("T2M")              # daily mean air temperature (°C)
R_PREC<-.get_or_load_var("PRECTOTCORR")      # precipitation (mm day⁻¹)
R_WIND<-.get_or_load_var("WS2M")             # wind speed at 2 m (m s⁻¹)
R_ALLSKY<-.get_or_load_var("ALLSKY_SFC_SW_DWN") # downward shortwave radiation (kWh m⁻² day⁻¹)
R_RH<-.get_or_load_var("RH2M")              # relative humidity at 2 m (%)

# =============================================================================
# 2) Write DNDC-unit NetCDFs (skipped if file already exists)
#    DNDC expects specific units that differ from raw POWER outputs.
#    Unit conversions applied: precip mm→cm (/10); radiation kWh→MJ (*3.6).
#    Temperature and wind are already in DNDC-compatible units.
# =============================================================================
message("=== Converting to DNDC Units ===")
write_if_missing(R_TMAX,      file.path(power_dir,"tmax-2005_2024-91.5x101.5x8x29.nc"),"tmax","Daily Max Temperature","C")
write_if_missing(R_TMIN,      file.path(power_dir,"tmin-2005_2024-91.5x101.5x8x29.nc"),"tmin","Daily Min Temperature","C")
write_if_missing(R_TAVG,      file.path(power_dir,"tavg-2005_2024-91.5x101.5x8x29.nc"),"tavg","Daily Avg Temperature","C")
write_if_missing(R_PREC/10,   file.path(power_dir,"prec-2005_2024-91.5x101.5x8x29.nc"),"prec","Daily Precipitation","cm")   # mm → cm
write_if_missing(R_WIND,      file.path(power_dir,"wind-2005_2024-91.5x101.5x8x29.nc"),"wind","Wind Speed","m/s")
write_if_missing(R_ALLSKY*3.6,file.path(power_dir,"srad-2005_2024-91.5x101.5x8x29.nc"),"srad","Solar Radiation","MJ/m2/day") # kWh → MJ
write_if_missing(R_RH,        file.path(power_dir,"rhum-2005_2024-91.5x101.5x8x29.nc"),"rhum","Relative Humidity","%")

# =============================================================================
# 3) Elevation
#    SRTM 30-arc-second DEM from geodata, resampled (averaged) to the 0.5° POWER
#    grid so every climate cell has a matching elevation value for DNDC input.
# =============================================================================
message("=== Processing Elevation ===")
elv_tif <- file.path(raw_dir,"elevation.tif")
if(!file.exists(elv_tif)){
  ref <- terra::rast(file.path(power_dir,"T2M_MAX-2005_2024-91.5x101.5x8x29.nc"))[[1]] # use as grid template
  elv <- geodata::elevation_30s("Myanmar",path=raw_dir)
  terra::resample(elv,ref,"average",filename=elv_tif,overwrite=TRUE); rm(ref,elv)
  message("✓ elevation.tif written")
} else message("✓ elevation.tif exists")

# =============================================================================
# 4) Soil (SoilGrids via geodata)
#    Variables: bulk density, clay/sand/silt fractions, SOC, pH at 3 depths.
#    Unit conversions (SoilGrids raw → physical units):
#      bdod : cg cm⁻³  → g cm⁻³  (/100)
#      clay/sand/silt: g kg⁻¹ → fraction (/1000)
#      soc  : dg kg⁻¹  → fraction (/1000)
#      pH   : pH×10    → pH     (/10)
#    Derived: porosity = 1 - (bdod / 2.65)  [particle density = 2.65 g cm⁻³]
#    USDA texture class (integer code) is computed from 0-5 cm clay/sand/silt %.
#    All layers are resampled (averaged) to the 0.5° elevation grid.
# =============================================================================
message("=== Processing Soil Data ===")
soil_tif <- file.path(raw_dir,"soil.tif"); soil_agg_tif <- file.path(raw_dir,"soil_agg.tif")
if(!file.exists(soil_agg_tif)){
  elv_ref <- terra::rast(elv_tif) # target grid for resampling
  if(!file.exists(soil_tif)){
    # Download SoilGrids mean values for 6 properties at depths 0-5, 5-15, 15-30 cm
    soil <- geodata::soil_world(c("bdod","clay","sand","silt","soc","phh2o"),c(5,15,30),stat="mean",vsi=TRUE)
    soil <- terra::crop(soil,ext,filename=soil_tif)
  } else soil <- terra::rast(soil_tif)
  # Robust layer lookup: handles both plain names and names with "_mean" suffix
  getL <- function(r,p){n<-names(r); if(p%in%n) return(r[[p]]); if(paste0(p,"_mean")%in%n) return(r[[paste0(p,"_mean")]]); m<-grep(paste0("^",p),n,value=TRUE); if(length(m)) return(r[[m[1]]]); stop("Missing layer: ",p)}
  # Convert units in-place (overwrites raw layers with physical-unit values)
  soil[["bdod_0-5cm"]]  <- getL(soil,"bdod_0-5cm")/100;   soil[["bdod_5-15cm"]]  <- getL(soil,"bdod_5-15cm")/100;  soil[["bdod_15-30cm"]] <- getL(soil,"bdod_15-30cm")/100
  soil[["clay_0-5cm"]]  <- getL(soil,"clay_0-5cm")/1000;  soil[["sand_0-5cm"]]   <- getL(soil,"sand_0-5cm")/1000;  soil[["silt_0-5cm"]]   <- getL(soil,"silt_0-5cm")/1000
  soil[["soc_0-5cm"]]   <- getL(soil,"soc_0-5cm")/1000;   soil[["pH_0-5cm"]]     <- getL(soil,"phh2o_0-5cm")/10;  soil[["poros_0-5cm"]]  <- 1-(soil[["bdod_0-5cm"]]/2.65)
  # USDA texture classification from clay/sand/silt percentages (0-5 cm layer)
  # Codes: 1=Sand, 2=Loamy sand, 3=Sandy loam, 4=Loam, 5=Silt loam, 6=Silt,
  #        7=Sandy clay loam, 8=Clay loam, 9=Silty clay loam,
  #        10=Sandy clay, 11=Silty clay, 12=Clay
  cl <- terra::values(soil[["clay_0-5cm"]])*100; sa <- terra::values(soil[["sand_0-5cm"]])*100; si <- terra::values(soil[["silt_0-5cm"]])*100
  tex <- rep(4L,length(cl))  # default: Loam
  tex[cl>40&sa>45]<-10L; tex[cl>40&si>40]<-11L; tex[cl>40&sa<=45&si<=40]<-12L
  tex[cl>=27&cl<=40&sa>20&sa<45]<-7L; tex[cl>=27&cl<=40&si>40]<-9L; tex[cl>=27&cl<=40&sa>=20&sa<=45&si<=40]<-8L
  tex[cl<27&cl>=7&sa<50&si>=28]<-5L; tex[cl<27&cl>=7&si>=80]<-6L; tex[cl<27&cl>=7&sa>=23&sa<52&si<50]<-4L
  tex[cl<27&cl>=7&sa>=50]<-3L; tex[cl<7&sa>=85]<-1L; tex[cl<7&sa>=70&sa<85]<-2L
  texture <- soil[["clay_0-5cm"]]; terra::values(texture) <- tex; names(texture) <- "texture_0-5cm"
  soil <- c(soil,texture); rm(cl,sa,si,tex,texture)
  # Resample all soil layers to the 0.5° climate grid using area-weighted average
  soil2 <- terra::resample(soil,elv_ref,"average",filename=soil_agg_tif,overwrite=TRUE); rm(soil,elv_ref)
  message("✓ soil_agg.tif written")
} else {soil2 <- terra::rast(soil_agg_tif); message("✓ soil_agg.tif exists")}

# =============================================================================
# 5) Rice cells table
#    Identifies which 0.5° grid cells within Myanmar have sufficient rice area.
#    Reference grid: the POWER climate grid (ensures cell IDs are consistent
#    with the climate files written above).
#    Rice area: SPAM 5-arc-min harvested area (ha), aggregated to 0.5° by sum
#    to preserve total area, then filtered by rice_min_ha threshold.
#    Output (cells.rds): data.frame with one row per rice cell containing
#    lon (x), lat (y), cell index, elevation, and all soil attributes.
#    This table drives Steps 2+ which loop over cells to write DNDC input files.
# =============================================================================
message("=== Creating Rice Cells Table (fixed) ===")
cells_rds <- file.path("data","cells.rds")
if(!file.exists(cells_rds)){
  aoi      <- geodata::gadm("Myanmar",level=1,path=raw_dir)
  ref_grid <- terra::rast(file.path(power_dir,"T2M_MAX-2005_2024-91.5x101.5x8x29.nc"))[[1]]
  # Assign sequential cell IDs on the POWER grid, masked to Myanmar (no boundary cells)
  cell_r   <- terra::mask(terra::init(ref_grid,"cell"),aoi,touches=touches_aoi)
  # SPAM rice: crop to Myanmar, sum to 0.5°, keep only cells meeting area threshold
  rice5m   <- geodata::crop_spam(crop="rice",var="area",raw_dir)|>terra::crop(aoi,mask=TRUE)
  rice_05  <- terra::resample(rice5m[[1]],ref_grid,method="sum") # area-preserving aggregation
  cell_r   <- terra::mask(cell_r,rice_05>=rice_min_ha,maskvalue=FALSE)
  # Build the cells table: coordinates + cell index + elevation + soil attributes
  cells    <- data.frame(cell_r,na.rm=TRUE)[,1]; xy <- data.frame(terra::xyFromCell(cell_r,cells)); xy$cell <- cells
  elv      <- terra::rast(elv_tif); xy$elevation <- round(elv[cells])
  xy       <- cbind(xy,terra::extract(soil2,cells)[,-1,drop=FALSE])
  saveRDS(xy,cells_rds)
  message(sprintf("✓ Created cells table with %d rice cells (min %.0f ha; touches=%s)",nrow(xy),rice_min_ha,touches_aoi))
} else {xy <- readRDS(cells_rds); message(sprintf("✓ Loaded existing cells table with %d rice cells",nrow(xy)))}

message("\n=== SUMMARY ===")
message(sprintf("Climate: 2005-2024 (%d years)",length(years)))
message("Vars: tmax,tmin,tavg,prec,wind,srad,rhum"); message("Grid: 0.5° x 0.5°")
message(sprintf("Rice-cell rule: rice_area >= %.0f ha; touches=%s",rice_min_ha,touches_aoi))
message(sprintf("Rice cells: %d",nrow(xy))); message(sprintf("Saved: %s",cells_rds))
message("\n✓ Step 1 complete. Ready for Step 2 (make_climate_files.R)")
