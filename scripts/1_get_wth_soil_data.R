##### 1_dndc_get_wth_soil_data.R - DNDC VERSION (Fix: touches=FALSE + rice_min_ha threshold) #####
suppressPackageStartupMessages({library(nasapower);library(terra);library(dplyr);library(geodata);library(httr)})
path <- "G:/My Drive/Research/simulation/main_dndc/simulate_dndc"
dir.create(path,FALSE,TRUE); setwd(path)
raw_dir <- file.path("data","raw"); power_dir <- file.path(raw_dir,"weather","power")
dir.create(power_dir,recursive=TRUE,showWarnings=FALSE)
years <- 2005:2024; ext <- terra::ext(91.5,101.5,8,29)
rice_min_ha <- 1000; touches_aoi <- FALSE  # key fix
lon_seq <- seq(xmin(ext),xmax(ext),by=0.5); lat_seq_desc <- sort(seq(ymin(ext),ymax(ext),by=0.5),decreasing=TRUE)
.make_layer <- function(df,vc){
  r <- terra::rast(ncols=length(lon_seq),nrows=length(lat_seq_desc),
                   xmin=xmin(ext),xmax=xmax(ext),ymin=ymin(ext),ymax=ymax(ext),crs="EPSG:4326")
  terra::values(r) <- (df|>dplyr::arrange(dplyr::desc(LAT),LON))[[vc]]; r}
.fetch_year <- function(parm,yr){
  d1 <- paste0(yr,"-01-01"); d2 <- paste0(yr,"-12-31")
  httr::set_config(httr::timeout(300)); on.exit(httr::reset_config())
  ts <- 4.0; ms <- 2.0
  lon0s <- seq(xmin(ext),xmax(ext),by=ts); lat0s <- seq(ymin(ext),ymax(ext),by=ts)
  if(length(lon0s)>1&&(xmax(ext)-lon0s[length(lon0s)])<ms) lon0s <- lon0s[-length(lon0s)]
  if(length(lat0s)>1&&(ymax(ext)-lat0s[length(lat0s)])<ms) lat0s <- lat0s[-length(lat0s)]
  out <- list()
  for(lon0 in lon0s) for(lat0 in lat0s){
    lon1 <- min(lon0+ts,xmax(ext)); lat1 <- min(lat0+ts,ymax(ext))
    if((lon1-lon0)<ms||(lat1-lat0)<ms) next
    for(k in 1:5){
      ok <- tryCatch({
        out[[length(out)+1]] <- nasapower::get_power(community="ag",lonlat=c(lon0,lat0,lon1,lat1),
                                                     pars=parm,dates=c(d1,d2),temporal_api="daily")
        Sys.sleep(3); TRUE
      },error=function(e){msg <- as.character(e)
        if(grepl("Timeout|timed out",msg,ignore.case=TRUE)){Sys.sleep(30*k);FALSE}
        else if(grepl("429|rate limit",msg,ignore.case=TRUE)){Sys.sleep(120*k);FALSE}
        else{Sys.sleep(20);FALSE}})
      if(ok) break
    }
  }
  if(!length(out)) stop("No tiles downloaded for ",parm," ",yr)
  dat <- dplyr::bind_rows(out)|>dplyr::rename(DATE=YYYYMMDD)
  dat$DATE <- as.Date(dat$DATE); dat <- dat[!duplicated(dat[,c("LAT","LON","DATE")]),]
  days <- sort(unique(dat$DATE))
  r <- terra::rast(lapply(seq_along(days),\(i).make_layer(dat[dat$DATE==days[i],c("LAT","LON",parm)],parm)))
  terra::time(r) <- days; names(r) <- paste0(parm,"_",format(days,"%Y%j")); r}
.get_or_load_var <- function(parm){
  fn <- file.path(power_dir,sprintf("%s-2005_2024-91.5x101.5x8x29.nc",parm))
  if(file.exists(fn)){message("Found existing ",basename(fn),"; loading..."); return(terra::rast(fn))}
  message("Downloading ",parm," ... (resume-safe; rerun if interrupted)")
  R <- terra::rast(lapply(years,\(yy).fetch_year(parm,yy)))
  terra::writeCDF(R,fn,varname=parm,longname=parm,overwrite=TRUE); R}
write_if_missing <- function(x,fn,vn,ln,u){
  if(!file.exists(fn)) terra::writeCDF(x,fn,varname=vn,longname=ln,unit=u,overwrite=TRUE)
  message("✓ ",basename(fn),if(file.exists(fn))" (exists)"else" (written)")}
# 1) Download/load POWER variables
message("=== Downloading (or loading) DNDC Climate Variables ===")
R_TMAX<-.get_or_load_var("T2M_MAX"); R_TMIN<-.get_or_load_var("T2M_MIN"); R_TAVG<-.get_or_load_var("T2M")
R_PREC<-.get_or_load_var("PRECTOTCORR"); R_WIND<-.get_or_load_var("WS2M"); R_ALLSKY<-.get_or_load_var("ALLSKY_SFC_SW_DWN"); R_RH<-.get_or_load_var("RH2M")
# 2) Convert to DNDC units (write only if missing)
message("=== Converting to DNDC Units ===")
write_if_missing(R_TMAX,   file.path(power_dir,"tmax-2005_2024-91.5x101.5x8x29.nc"),"tmax","Daily Max Temperature","C")
write_if_missing(R_TMIN,   file.path(power_dir,"tmin-2005_2024-91.5x101.5x8x29.nc"),"tmin","Daily Min Temperature","C")
write_if_missing(R_TAVG,   file.path(power_dir,"tavg-2005_2024-91.5x101.5x8x29.nc"),"tavg","Daily Avg Temperature","C")
write_if_missing(R_PREC/10,file.path(power_dir,"prec-2005_2024-91.5x101.5x8x29.nc"),"prec","Daily Precipitation","cm")
write_if_missing(R_WIND,   file.path(power_dir,"wind-2005_2024-91.5x101.5x8x29.nc"),"wind","Wind Speed","m/s")
write_if_missing(R_ALLSKY*3.6,file.path(power_dir,"srad-2005_2024-91.5x101.5x8x29.nc"),"srad","Solar Radiation","MJ/m2/day")
write_if_missing(R_RH,     file.path(power_dir,"rhum-2005_2024-91.5x101.5x8x29.nc"),"rhum","Relative Humidity","%")
# 3) Elevation (resample to climate grid; do not redo if exists)
message("=== Processing Elevation ===")
elv_tif <- file.path(raw_dir,"elevation.tif")
if(!file.exists(elv_tif)){
  ref <- terra::rast(file.path(power_dir,"T2M_MAX-2005_2024-91.5x101.5x8x29.nc"))[[1]]
  elv <- geodata::elevation_30s("Myanmar",path=raw_dir)
  terra::resample(elv,ref,"average",filename=elv_tif,overwrite=TRUE); rm(ref,elv)
  message("✓ elevation.tif written")
} else message("✓ elevation.tif exists")
# 4) Soil (do not redo if aggregated exists)
message("=== Processing Soil Data ===")
soil_tif <- file.path(raw_dir,"soil.tif"); soil_agg_tif <- file.path(raw_dir,"soil_agg.tif")
if(!file.exists(soil_agg_tif)){
  elv_ref <- terra::rast(elv_tif)
  if(!file.exists(soil_tif)){
    soil <- geodata::soil_world(c("bdod","clay","sand","silt","soc","phh2o"),c(5,15,30),stat="mean",vsi=TRUE)
    soil <- terra::crop(soil,ext,filename=soil_tif)
  } else soil <- terra::rast(soil_tif)
  getL <- function(r,p){n<-names(r); if(p%in%n) return(r[[p]]); if(paste0(p,"_mean")%in%n) return(r[[paste0(p,"_mean")]]); m<-grep(paste0("^",p),n,value=TRUE); if(length(m)) return(r[[m[1]]]); stop("Missing layer: ",p)}
  soil[["bdod_0-5cm"]]  <- getL(soil,"bdod_0-5cm")/100;   soil[["bdod_5-15cm"]]  <- getL(soil,"bdod_5-15cm")/100;  soil[["bdod_15-30cm"]] <- getL(soil,"bdod_15-30cm")/100
  soil[["clay_0-5cm"]]  <- getL(soil,"clay_0-5cm")/1000;  soil[["sand_0-5cm"]]   <- getL(soil,"sand_0-5cm")/1000;  soil[["silt_0-5cm"]]   <- getL(soil,"silt_0-5cm")/1000
  soil[["soc_0-5cm"]]   <- getL(soil,"soc_0-5cm")/1000;   soil[["pH_0-5cm"]]     <- getL(soil,"phh2o_0-5cm")/10;  soil[["poros_0-5cm"]]  <- 1-(soil[["bdod_0-5cm"]]/2.65)
  cl <- terra::values(soil[["clay_0-5cm"]])*100; sa <- terra::values(soil[["sand_0-5cm"]])*100; si <- terra::values(soil[["silt_0-5cm"]])*100
  tex <- rep(4L,length(cl))
  tex[cl>40&sa>45]<-10L; tex[cl>40&si>40]<-11L; tex[cl>40&sa<=45&si<=40]<-12L
  tex[cl>=27&cl<=40&sa>20&sa<45]<-7L; tex[cl>=27&cl<=40&si>40]<-9L; tex[cl>=27&cl<=40&sa>=20&sa<=45&si<=40]<-8L
  tex[cl<27&cl>=7&sa<50&si>=28]<-5L; tex[cl<27&cl>=7&si>=80]<-6L; tex[cl<27&cl>=7&sa>=23&sa<52&si<50]<-4L
  tex[cl<27&cl>=7&sa>=50]<-3L; tex[cl<7&sa>=85]<-1L; tex[cl<7&sa>=70&sa<85]<-2L
  texture <- soil[["clay_0-5cm"]]; terra::values(texture) <- tex; names(texture) <- "texture_0-5cm"
  soil <- c(soil,texture); rm(cl,sa,si,tex,texture)
  soil2 <- terra::resample(soil,elv_ref,"average",filename=soil_agg_tif,overwrite=TRUE); rm(soil,elv_ref)
  message("✓ soil_agg.tif written")
} else {soil2 <- terra::rast(soil_agg_tif); message("✓ soil_agg.tif exists")}
# 5) Rice mask + cells table (FIXED)
message("=== Creating Rice Cells Table (fixed) ===")
cells_rds <- file.path("data","cells.rds")
if(!file.exists(cells_rds)){
  aoi <- geodata::gadm("Myanmar",level=1,path=raw_dir)
  ref_grid <- terra::rast(file.path(power_dir,"T2M_MAX-2005_2024-91.5x101.5x8x29.nc"))[[1]]
  cell_r <- terra::mask(terra::init(ref_grid,"cell"),aoi,touches=touches_aoi)
  rice5m <- geodata::crop_spam(crop="rice",var="area",raw_dir)|>terra::crop(aoi,mask=TRUE)
  rice_05 <- terra::resample(rice5m[[1]],ref_grid,method="sum")
  cell_r <- terra::mask(cell_r,rice_05>=rice_min_ha,maskvalue=FALSE)
  cells <- data.frame(cell_r,na.rm=TRUE)[,1]; xy <- data.frame(terra::xyFromCell(cell_r,cells)); xy$cell <- cells
  elv <- terra::rast(elv_tif); xy$elevation <- round(elv[cells])
  xy <- cbind(xy,terra::extract(soil2,cells)[,-1,drop=FALSE])
  saveRDS(xy,cells_rds)
  message(sprintf("✓ Created cells table with %d rice cells (min %.0f ha; touches=%s)",nrow(xy),rice_min_ha,touches_aoi))
} else {xy <- readRDS(cells_rds); message(sprintf("✓ Loaded existing cells table with %d rice cells",nrow(xy)))}
message("\n=== SUMMARY ===")
message(sprintf("Climate: 2005-2024 (%d years)",length(years)))
message("Vars: tmax,tmin,tavg,prec,wind,srad,rhum"); message("Grid: 0.5° x 0.5°")
message(sprintf("Rice-cell rule: rice_area >= %.0f ha; touches=%s",rice_min_ha,touches_aoi))
message(sprintf("Rice cells: %d",nrow(xy))); message(sprintf("Saved: %s",cells_rds))
message("\n✓ Step 1 complete. Ready for Step 2 (make_climate_files.R)")
