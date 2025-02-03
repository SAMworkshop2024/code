library(metR)
library(data.table)
library(rcdo)

pattern <- "*/*/hist*/*/Amon/psl/*/*/*.nc"


folders <- c("/g/data/oi10/replicas/CMIP6/CMIP/",
             "/g/data/oi10/replicas/CMIP6/DAMIP",
             "/g/data/fs38/publications/CMIP6/CMIP/",
             "/g/data/fs38/publications/CMIP6/DAMIP/")

files <- file.path(folders, pattern) |> 
  Sys.glob()

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
  _[experiment_id %in% c("hist-GHG", "hist-stratO3", "historical")]


models_latest <- models[, .SD[version_id == max(version_id)], 
                        by = .(source_id, model_id, member_id, experiment_id)]


write_sam <- function(files_in, dir = "sam_indices/cmip6") {
  name <- paste0(tools::file_path_sans_ext(basename(files_in[1])), ".csv")
  file_out <- file.path(dir, name)
  
  if (file.exists(file_out)) {
    return(file_out)
  }
  
  # Don't use models with weird grids.
  if (is.null(GlanceNetCDF(files_in[1])$dims$lat)) {
    return(NULL)
  }
  
  if (length(files_in) > 1) {
    Sys.setenv(SKIP_SAME_TIME = 1)
    file_merged <- files_in |> 
      cdo_mergetime() |> 
      cdo_execute()
  } else {
    file_merged <- files_in
  }
  
  subset <- list(lat = list(-40, -65),
                 time = c("1979-01-01", "2014-12-31"))
  
  if ("latitude" %in% names(GlanceNetCDF(files_in[1])$dims)) {
    names(subset)[1] <- "latitude"
  }
  
  sam <- ReadNetCDF(file_merged, 
                    vars = "psl", 
                    subset = subset) |> 
    _[year(time) %between% c(1979, 20015)]
  
  if (nrow(sam) == 0) {
    return(NULL)
  }
  
  sam <- sam |> 
    setnames("latitude", "lat", skip_absent = TRUE) |> 
    _[, .(psl = mean(psl)), by = .(lat, time = seasonally(time))] |> 
    _[, psl := (psl - mean(psl))/sd(psl), by = .(lat, season(time))] |> 
    _[order(lat)] |> 
    _[, .(psl = diff(psl)), by = .(time)] 
  
  fwrite(sam, file_out)
  file_out
}

sams <- models_latest[, .(file = write_sam(file)),
                      by = .(source_id, model_id, member_id, experiment_id)]

sam_data <- sams[, fread(file),  by = .(institution = source_id, 
                                        model = model_id, 
                                        member = member_id, 
                                        forcing = experiment_id)] |> 
  setnames("psl", "sam")

fwrite(sam_data, "sam_indices/sam-cmip6.csv")

saveRDS(sam_data, "cmip-damip/shiny/sam-cmip6.Rds")
