library(metR)
library(data.table)
library(rcdo)

future::plan("multicore", workers = 24)

pattern <- "*/*/*/*/Amon/psl/*/*/*.nc"
experiments <- c("hist-GHG", "hist-totalO3", "hist-stratO3",
                 "historical", "ssp245-stratO3", "ssp245-GHG")

folders <- c("/g/data/oi10/replicas/CMIP6/CMIP/",
             "/g/data/oi10/replicas/CMIP6/DAMIP",
             "/g/data/fs38/publications/CMIP6/CMIP/",
             "/g/data/fs38/publications/CMIP6/DAMIP/")

files <- file.path(folders, pattern) |> 
  Sys.glob()


get_times <- function(files) {
  basename(files) |> 
    gsub("\\.nc", "", x = _) |> 
    strsplit("_") |> 
    lapply(\(x) x[7]) |> 
    unlist() |> 
    strsplit("-") |> 
    transpose()  |> 
    lapply(as.numeric) |> 
    setNames(c("start", "end"))
}



models <- strcapture("([^/]*)/([^/]*)/([^/]*)/([^/]*)/Amon/psl/([^/]*)/v([^/]*)/*", 
                     files, 
                     list(
                       source_id = character(1),
                       model_id = character(1),
                       
                       experiment_id = character(1),
                       
                       member_id = character(1),
                       
                       grid_label = character(1),
                       version_id = numeric(1)
                     ))  |> 
  as.data.table() |> 
  _[, file  := files] |> 
  _[experiment_id %in% experiments]


models[, c("start", "end") := get_times(file)]




era_files <- "/g/data/rt52/era5/single-levels/monthly-averaged/msl/*/*" |> 
  Sys.glob()

era5 <- strcapture("msl_era5_moda_sfc_([^/]*)-([^/]*).nc", 
           basename(era_files), 
           list(start = character(1),
                end = character(1))) |> 
  as.data.table() |> 
  _[, file  := era_files] |> 
  _[, let(source_id = "ERA5", 
          model_id = "ERA5",
          member_id = "r1i1p1f1",
          experiment_id = "historical",
          grid_label = "gr", 
          version_id = 20250213,
          start = as.numeric(substr(start, 1, 6)),
          end = as.numeric(substr(end, 1, 6))
          )] 

models <- rbind(models, era5)


models <- models[end >= 197901] |> 
  _[start <= 201412]

models_latest <- models[, .SD[version_id == max(version_id)], 
                        by = .(source_id, model_id, member_id, experiment_id)]


write_sam <- function(files_in, dir = "sam_indices/cmip6/elio_version") {
  name <- paste0(tools::file_path_sans_ext(basename(files_in[1])), ".csv")
  file_out <- file.path(dir, name)
  
  if (file.exists(file_out)) {
    return(fread(file_out))
  }
  
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  
  # Don't use models with weird grids.
  if (is.null(GlanceNetCDF(files_in[1])$dims$lat)) {
    return(NULL)
  }
  
  if (length(files_in) > 1) {
    Sys.setenv(SKIP_SAME_TIME = 1)
    file_merged <- files_in |> 
      cdo_mergetime() |> 
      cdo_execute(options = "-b f32") 
  } else {
    file_merged <- files_in
  }
  
  subset <- list(lat = list(-40, -65),
                 time = c("1979-01-01", "2014-12-31"))
  
  if ("latitude" %in% names(GlanceNetCDF(files_in[1])$dims)) {
    names(subset)[1] <- "latitude"
  }
  
  sam <- ReadNetCDF(file_merged, 
                    vars = \(x) c(psl = x[x %in% c("psl", "msl")]),
                    subset = subset) |> 
    _[year(time) %between% c(1979, 2014)]
  
  if (nrow(sam) == 0) {
    return(NULL)
  }
  
  sam <- sam |> 
    setnames("latitude", "lat", skip_absent = TRUE) |> 
    _[, .(psl = mean(psl)), by = .(time, lat)] |> 
    _[, psl := (psl - mean(psl))/sd(psl), by = .(lat, month(time))] |> 
    _[order(lat)] |> 
    _[, .(psl = diff(psl)), by = .(time)] 
  
  fwrite(sam, file_out)
  sam
}


sam_data <- models_latest |> 
  split(by = c("source_id", "model_id", "member_id", "experiment_id")) |> 
  furrr::future_map(\(df) df[, write_sam(file),
                             by = .(institute = source_id, 
                                    model = model_id, 
                                    member = member_id, 
                                    forcing = experiment_id)]) |> 
  Filter(f = \(x) ncol(x) == 6) |> 
  rbindlist() |> 
  na.omit() |> 
  setnames("psl", "sam")

fwrite(sam_data, "sam_indices/sam-cmip6.csv")

saveRDS(sam_data, "cmip-damip/shiny/sam-cmip6.Rds")
