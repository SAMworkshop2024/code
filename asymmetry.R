library(ggplot2)
library(metR)
library(data.table)


# Function to compute p-values
pvaluate <- function (estimate, std.error, df)  {
  2 * stats::pt(abs(estimate)/std.error, df, 
                lower.tail = FALSE)
}

# What happens in the models 


index_file <-  "sam_indices/SAM_GW_1m_1979-2023.nc"

climatology <- c(1979, 2011)
years <- 1979:2023

files <- list.files(paste0("/g/data/rt52/era5/single-levels/monthly-averaged/msl/", years, "/"), 
                    full.names = TRUE)


# Remap to lower resolution 
paste("cdo -L genbil,r144x72", files[1], "remapweights.nc") |> 
  system()

files_out <- file.path("era5temp/", basename(files))
dir.create("era5temp",  showWarnings = FALSE)
# I cannot remap and merge everthing for some reason. 
# Maybe due to limits in the number of arguments. 
for (year in years) {
  paste0("cdo -L -b F32 -O -mergetime -apply,\"-remap,r144x72,remapweights.nc\" [ /g/data/rt52/era5/single-levels/monthly-averaged/msl/",
         year, "/*.nc ] era5temp/era.", year, ".nc") |> 
    system()
  
}

system("cdo -L -O mergetime -apply,\"-sellonlatbox,0.,360.,-90.,0.\" [ era5temp/era.????.nc ] era.nc")
unlink("era5temp", recursive = TRUE)

mslp <- ReadNetCDF("era.nc", subset = list(lat = c(-90, -20))) |> 
  _[, msl := msl - mean(msl[year(time) %between% climatology]), by = .(lon, lat, month(time))] 

sam_index <- mslp[lat %~% c(-40, -65)] |>   # sea level pressure from ERA5 at -65 and -40
  # Compute mean by latitude
  _[, .(msl = mean(msl)), by = .(time, lat)] |>
  # Remove the monthly climatological mean
  _[, msl := msl - mean(msl[year(time) %between% climatology]), by = .(lat, month(time))] |> 
  # Standardise by the monthly climatological standard deviation
  _[, msl := msl/sd(msl[year(time) %between% climatology]), by = .(lat, month(time))] |>
  # Reshape so that each latitude is its own column
  tidyfast::dt_pivot_wider(names_from = lat, values_from = msl) |> 
  setnames(c("time", "-65", "-40")) |> 
  # Sam is the difference of the latitudes
  _[, .(sam = `-40` - `-65` ), by = time] 


sam_map <- mslp |> 
  _[, msl := msl - mean(msl[year(time) %between% climatology]), by = .(lat, lon, season(time))] |> 
  _[sam_index, on = 'time'] |> 
  _[, FitLm(msl, sam), by = .(lat, lon, month(time))] |> 
  _[term == "sam"]


variances <- sam_map |> 
  na.omit() |> 
  _[, estimate_a := estimate - mean(estimate), by = .(month, lat)] |>
  _[, .(asym = modi::weighted.var(estimate_a, sqrt(cos(lat*pi/180)))/modi::weighted.var(estimate, sqrt(cos(lat*pi/180))),
        sym = modi::weighted.var(estimate - estimate_a, sqrt(cos(lat*pi/180)))/modi::weighted.var(estimate, sqrt(cos(lat*pi/180)))),
    by = month] |> 
  tidyfast::dt_pivot_longer(cols = asym:sym) 

variances |> 
  ggplot(aes(month, value)) +
  geom_col(aes(fill = name), position = "fill")  +
  scale_x_continuous(NULL, labels = month.abb, breaks = 1:12, expand = c(0, 0)) +
  scale_y_continuous("Explained variance", limits = c(0, 1), labels = scales::percent_format(),
                     expand = c(0, 0)) +
  scale_fill_brewer("Part", labels = c(asym = "Asymmetric", 
                                       sym = "Symmetric"),
                    guide = "none",
                    palette = "Set1") + 
  geom_text(data = \(x) x[month == 3 & name == "sym"], 
            label = "Asymmetric", angle = 90, color = "white", size = 6, hjust = -0.2) +
  
  geom_text(data = \(x) x[month == 3 & name == "sym"], 
            label = "Symmetric", angle = 90, color = "white", size = 6, hjust = 1.2) +
  
  theme_minimal(base_size = 16)

ggdatasaver::save_plot_data(last_plot(), name = "asymmetry")
