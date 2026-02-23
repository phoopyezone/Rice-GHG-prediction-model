##### 1_dndc_get_wth_soil_data.R - DNDC VERSION (condensed) #####
# Fix: touches=FALSE + rice_min_ha threshold (not just >0). All 5 steps preserved.
suppressPackageStartupMessages({library(nasapower);library(terra);library(dplyr);library(geodata);library(httr)})
path <- "G:/My Drive/Research/simulation/main_dndc/simulate_dndc"
dir.create(path,FALSE,TRUE); setwd(path)
raw_dir <- file.path("data","raw"); power_dir <- file.path(raw_dir,"weather","power")
dir.create(power_dir,recursive=TRUE,showWarnings=FALSE); dir.create("data",showWarnings=FALSE)
years <- 2005:2024; ext0 <- terra::ext(91.5,101.5,8,29)
rice_min_ha <- 1000; touches_aoi <- FALSE  # key fix
lon_seq <- seq(xmin(ext0),xmax(ext0),by=0.5); lat_seq <- sort(seq(ymin(ext0),ymax(ext0),by=0.5),decreasing=TRUE)
.make_layer <- function(df,val){r<-terra::rast(ncols=length(lon_seq),nrows=length(lat_seq),xmin=xmin(ext0),xmax=xmax(ext0),ymin=ymin(ext0),ymax=ymax(ext0),crs="EPSG:4326"); terra::values(r)<-(df|>dplyr::arrange(dplyr::desc(LAT),LON))[[val]]; r}
.fetch_year <- function(parm,yr){
  d1<-paste0(yr,"-01-01"); d2<-paste0(yr,"-12-31")
  httr::set_config(httr::timeout(300)); on.exit(httr::reset_config(),add=TRUE)
  ts<-4; ms<-2
  lon0s<-seq(xmin(ext0),xmax(ext0),by=ts); lat0s<-seq(ymin(ext0),ymax(ext0),by=ts)
  if(length(lon0s)>1&&(xmax(ext0)-lon0s[length(lon0s)])<ms) lon0s<-lon0s[-length(lon0s)]
  if(length(lat0s)>1&&(ymax(ext0)-lat0s[length(lat0s)])<ms) lat0s<-lat0s[-length(lat0s)]
  out<-list()
  for(lo in lon0s) for(la in lat0s){
    l1<-min(lo+ts,xmax(ext0)); l2<-min(la+ts,ymax(ext0))
    if((l1-lo)<ms||(l2-la)<ms) next
    for(k in 1:5){ok<-tryCatch({out[[length(out)+1]]<-nasapower::get_power(community="ag",lonlat=c(lo,la,l1,l2),pars=parm,dates=c(d1,d2),temporal_api="daily"); Sys.sleep(3); TRUE},error=function(e){msg<-as.character(e); if(grepl("Timeout|timed out",msg,ignore.case=TRUE)){Sys.sleep(30*k);FALSE} else if(grepl("429|rate limit",msg,ignore.case=TRUE)){Sys.sleep(120*k);FALSE} else{Sys.sleep(15);FALSE}}); if(ok) break}
  }
  if(!length(out)) stop("No tiles for ",parm," ",yr)
  dat<-dplyr::bind_rows(out)|>dplyr::rename(DATE=YYYYMMDD); dat$DATE<-as.Date(dat$DATE); dat<-dat[!duplicated(dat[,c("LAT","LON","DATE")]),]
  days<-sort(unique(dat$DATE)); r<-terra::rast(lapply(days,\(dd).make_layer(dat[dat$DATE==dd,c("LAT","LON",parm)],parm))); terra::time(r)<-days; names(r)<-paste0(parm,"_",format(days,"%Y%j")); r}
.get_or_load_var <- function(parm){fn<-file.path(power_dir,sprintf("%s-2005_2024-91.5x101.5x8x29.nc",parm)); if(file.exists(fn)) return(terra::rast(fn)); R<-terra::rast(lapply(years,\(yy).fetch_year(parm,yy))); terra::writeCDF(R,fn,varname=parm,longname=parm,overwrite=TRUE); R}
wim <- function(x,fn,vn,ln,u) if(!file.exists(fn)) terra::writeCDF(x,fn,varname=vn,longname=ln,unit=u,overwrite=TRUE)
# 1) Download/load POWER variables
message("=== 1) POWER climate (load or download) ===")
R_TMAX<-.get_or_load_var("T2M_MAX"); R_TMIN<-.get_or_load_var("T2M_MIN"); R_TAVG<-.get_or_load_var("T2M")
R_PREC<-.get_or_load_var("PRECTOTCORR"); R_WIND<-.get_or_load_var("WS2M"); R_ALLSKY<-.get_or_load_var("ALLSKY_SFC_SW_DWN"); R_RH<-.get_or_load_var("RH2M")
# 2) Convert to DNDC units (skip if exists)
message("=== 2) DNDC-unit netCDFs ===")
wim(R_TMAX,  file.path(power_dir,"tmax-2005_2024-91.5x101.5x8x29.nc"), "tmax","Daily Max Temperature","C")
wim(R_TMIN,  file.path(power_dir,"tmin-2005_2024-91.5x101.5x8x29.nc"), "tmin","Daily Min Temperature","C")
wim(R_TAVG,  file.path(power_dir,"tavg-2005_2024-91.5x101.5x8x29.nc"), "tavg","Daily Avg Temperature","C")
wim(R_PREC/10,file.path(power_dir,"prec-2005_2024-91.5x101.5x8x29.nc"),"prec","Daily Precipitation","cm")
wim(R_WIND,  file.path(power_dir,"wind-2005_2024-91.5x101.5x8x29.nc"), "wind","Wind Speed","m/s")
wim(R_ALLSKY*3.6,file.path(power_dir,"srad-2005_2024-91.5x101.5x8x29.nc"),"srad","Solar Radiation","MJ/m2/day")
wim(R_RH,    file.path(power_dir,"rhum-2005_2024-91.5x101.5x8x29.nc"), "rhum","Relative Humidity","%")
# 3) Elevation (resample to climate grid; write once)
message("=== 3) Elevation ===")
elv_tif <- file.path(raw_dir,"elevation.tif")
if(!file.exists(elv_tif)){ref<-terra::rast(file.path(power_dir,"T2M_MAX-2005_2024-91.5x101.5x8x29.nc"))[[1]]; terra::resample(geodata::elevation_30s("Myanmar",path=raw_dir),ref,"average",filename=elv_tif,overwrite=TRUE)}
# 4) Soil (write aggregated tif once)
message("=== 4) Soil ===")
soil_tif<-file.path(raw_dir,"soil.tif"); soil_agg_tif<-file.path(raw_dir,"soil_agg.tif")
if(!file.exists(soil_agg_tif)){
  if(!file.exists(soil_tif)) soil<-terra::crop(geodata::soil_world(c("bdod","clay","sand","silt","soc","phh2o"),c(5,15,30),stat="mean",vsi=TRUE),ext0,filename=soil_tif) else soil<-terra::rast(soil_tif)
  getL<-function(r,nm){n<-names(r); if(nm%in%n) return(r[[nm]]); if(paste0(nm,"_mean")%in%n) return(r[[paste0(nm,"_mean")]]); m<-grep(paste0("^",nm),n,value=TRUE); if(length(m)) return(r[[m[1]]]); stop("Missing: ",nm)}
  soil[["bdod_0-5cm"]]<-getL(soil,"bdod_0-5cm")/100; soil[["bdod_5-15cm"]]<-getL(soil,"bdod_5-15cm")/100; soil[["bdod_15-30cm"]]<-getL(soil,"bdod_15-30cm")/100
  soil[["clay_0-5cm"]]<-getL(soil,"clay_0-5cm")/1000; soil[["sand_0-5cm"]]<-getL(soil,"sand_0-5cm")/1000; soil[["silt_0-5cm"]]<-getL(soil,"silt_0-5cm")/1000
  soil[["soc_0-5cm"]]<-getL(soil,"soc_0-5cm")/1000; soil[["pH_0-5cm"]]<-getL(soil,"phh2o_0-5cm")/10; soil[["poros_0-5cm"]]<-1-(soil[["bdod_0-5cm"]]/2.65)
  cl<-terra::values(soil[["clay_0-5cm"]])*100; sa<-terra::values(soil[["sand_0-5cm"]])*100; si<-terra::values(soil[["silt_0-5cm"]])*100; tx<-rep(4L,length(cl))
  tx[cl>40&sa>45]<-10L; tx[cl>40&si>40]<-11L; tx[cl>40&sa<=45&si<=40]<-12L; tx[cl>=27&cl<=40&sa>20&sa<45]<-7L; tx[cl>=27&cl<=40&si>40]<-9L; tx[cl>=27&cl<=40&sa>=20&sa<=45&si<=40]<-8L
  tx[cl<27&cl>=7&sa<50&si>=28]<-5L; tx[cl<27&cl>=7&si>=80]<-6L; tx[cl<27&cl>=7&sa>=23&sa<52&si<50]<-4L; tx[cl<27&cl>=7&sa>=50]<-3L; tx[cl<7&sa>=85]<-1L; tx[cl<7&sa>=70&sa<85]<-2L
  tex<-soil[["clay_0-5cm"]]; terra::values(tex)<-tx; names(tex)<-"texture_0-5cm"; soil<-c(soil,tex)
  soil2<-terra::resample(soil,terra::rast(elv_tif),"average",filename=soil_agg_tif,overwrite=TRUE)
} else soil2<-terra::rast(soil_agg_tif)
# 5) Rice mask + cells table (FIXED: touches=FALSE, rice_min_ha threshold)
message("=== 5) Rice cells (fixed) ===")
cells_rds <- file.path("data","cells.rds")
if(!file.exists(cells_rds)){
  aoi<-geodata::gadm("Myanmar",level=1,path=raw_dir)
  ref_grid<-terra::rast(file.path(power_dir,"T2M_MAX-2005_2024-91.5x101.5x8x29.nc"))[[1]]
  cell_r<-terra::mask(terra::init(ref_grid,"cell"),aoi,touches=touches_aoi)
  rice_05<-terra::resample(terra::crop(geodata::crop_spam(crop="rice",var="area",raw_dir)[[1]],aoi,mask=TRUE),ref_grid,method="sum")
  cell_r<-terra::mask(cell_r,rice_05>=rice_min_ha,maskvalue=FALSE)
  cells<-data.frame(cell_r,na.rm=TRUE)[,1]; xy<-data.frame(terra::xyFromCell(cell_r,cells)); xy$cell<-cells
  xy$elevation<-round(terra::rast(elv_tif)[cells]); xy<-cbind(xy,terra::extract(soil2,cells)[,-1,drop=FALSE])
  saveRDS(xy,cells_rds); message(sprintf("✓ Created %d rice cells (>=%.0f ha; touches=%s)",nrow(xy),rice_min_ha,touches_aoi))
} else {xy<-readRDS(cells_rds); message(sprintf("✓ Loaded %d rice cells",nrow(xy)))}
message(sprintf("\n=== SUMMARY: %d years | rice>=%.0f ha | touches=%s | %d cells | saved: %s ===",length(years),rice_min_ha,touches_aoi,nrow(xy),cells_rds))
message("✓ Step 1 complete. Ready for Step 2 (make_climate_files.R)")
