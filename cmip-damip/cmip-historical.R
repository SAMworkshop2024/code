
library(data.table)
library(rcdo)

base_folder <- "/g/data/oi10/replicas/CMIP6/CMIP/"

files <- Sys.glob(file.path(base_folder, "*/*/historical/*/Amon/psl/*/*/*.nc"))

models <- strcapture("([^/]*)/([^/]*)/historical/([^/]*)/Amon/psl/([^/]*)/v([^/]*)/*", files, list(
  source_id = character(1),
  model_id = character(1),

  member_id = character(1),

  grid_label = character(1),
  version_id = numeric(1)
))  |> 
    as.data.table() |> 
    _[, file  := files]

models_latest <- models[, .SD[which.max(version_id)], by =  .(source_id, model_id, member_id)]


compute_sam <- function(files_in, dir = "sam_indices/cmip6") {
    name <- basename(files_in[1])
    file_out <- file.path("cmip6-sam", name)

    if (file.exists(file_out)) {
      return(file_out)
    }

    if (length(files_in) > 1) {
      file_in <- tempfile()
      on.exit(unlink(file_in))
      system(paste0("cdo mergetime ", paste(files_in, collapse = " "), " ", file_in))
    } else {
      file_in  <- files_in[1]
    }
    
    dir.create(dir, showWarnings = FALSE)

    system(paste0("bash compute-sam.sh ", file_in, " ", file_out))
    return(file_out)
}

sams <- models[, .(file = compute_sam(file)), by = .(source_id, model_id, member_id)]
