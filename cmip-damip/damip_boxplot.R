library(metR)
library(ggplot2)
library(data.table)


# historical --------------------------------------------------------------

files_historical <- list.files("sam_indices/cmip6", full.names = TRUE)


historical <- basename(files_historical) |>
  strsplit(split = "_") |>
  transpose() |>
  as.data.table() |>
  _[, -c(1, 2, 4, 6, 7)] |>
  setnames(c("model", "ensemble")) |>
  _[, file := files_historical]


read_historical <- function(file) {
  f <- ncdf4::nc_open(file)
  institute <- ncdf4::ncatt_get(f, 0, "institution_id")$value

  ReadNetCDF(f, vars = "sam") |>
    _[year(time) %between% c(1979, 2014)] |>
    _[season(time) %in% c("DJF", "JJA")] |>
    _[, .(sam = mean(sam)), by = .(time = seasonally(time))] |>
    _[, institute := institute][]
}

sam_historical <- historical[, read_historical(file), by = .(model, ensemble)] |>
  _[, forcing := "historical"]



# damip -------------------------------------------------------------------

damip_models <- c(
  "ACCESS-CM2", "ACCESS-ESM1-5", "BCC-CSM2-MR", "CESM2",
  "CMCC-CM2-SR5", "CNRM-CM6-1", "CanESM5", "E3SM-1-0",
  "E3SM-1-1", "E3SM-2-0", "EC-Earth3", "FGOALS-g3", "GFDL-CM4",
  "GFDL-ESM4", "GISS-E2-1-G", "GISS-E2-2-G", "HadGEM3-GC31-LL",
  "IPSL-CM6A-LR", "MIROC-ES2L", "MIROC6", "MPI-ESM1-2-LR",
  "MRI-ESM2-0", "NorESM2-LM", "UKESM1-0-LL"
)

sam_historical <- sam_historical[model %in% damip_models]

file_damip <- "sam_indices/sam_psl.DAMIP.1850-2020.nc"

forcing_names <- c("hist-GHG", "hist-stratO3")

institute_names <- c(
  "BCC", "CAS", "CCCma", "CNRM-CERFACS", "E3SM-Project",
  "EC-Earth-Consortium", "IPSL", "MIROC", "MOHC", "MPI-M",
  "MRI", "NASA-GISS", "NCAR", "NCC", "NOAA-GFDL"
)

# There are some models that are all NAs
# Probably institute-model-forcing-ensemble crossings that didn't exist in the
# original data
frac_valid <- ReadNetCDF(file_damip, vars = "sam") |>
  _[, .(frac_valid = mean(is.finite(sam))), by = .(institute, model, forcing, ensemble)]

sam_damip <- ReadNetCDF(file_damip, vars = "sam") |>
  _[frac_valid[frac_valid != 0], on = .NATURAL] |> # keep only valids
  _[, institute := factor(institute_names[institute])] |>
  _[, forcing := factor(forcing_names[forcing])] |>
  _[, time := lubridate::make_date(year, month)] |>
  _[year %between% c(1979, 2014)] |>
  _[season(time) %in% c("DJF", "JJA")] |>
  _[, .(sam = mean(sam, na.rm = TRUE)),
    by = .(institute, model, ensemble, forcing, time = seasonally(time))] |>
  _[, ensemble := interaction(institute, model, forcing)] |>
  _[, model := interaction(institute, model)]

sam <- rbind(sam_damip, sam_historical)

complete_institute <- sam |>
  dcast(institute ~ forcing, value.var = "sam") |>
  melt(id.vars = "institute") |>
  _[, .SD[all(value != 0)], by = .(institute)] |>
  _[, unique(institute)]

sam <- sam[institute %in% as.character(complete_institute)]

# Function to compute estimates and standard errors of the slopes
# of segmented regression.
fit_segmented <- function(formula, data = NULL) {
  mod <- lm(formula, data = data)
  coef <- cumsum(coef(mod)[-1])
  terms <- names(coef)

  list(
    period = terms,
    estimate = coef,
    se = c(
      sqrt(vcov(mod)[2:2, 2:2]),
      sqrt(sum(vcov(mod)[2:3, 2:3]))
    ),
    df = rep(mod$df.residual, 2)
  )
}

cut_year <- 1999

mmm_trend <- sam |>
  _[, .(sam = mean(sam, na.rm = TRUE)), by = .(model, forcing, time = seasonally(time))] |>
  _[, .(sam = mean(sam, na.rm = TRUE)), by = .(forcing, time = seasonally(time))] |>
  _[, let(
    `1979--1999` = year(time),
    `2000--2014` = pmax(year(time) - cut_year, 0),
    year = year(time)
  )] |>
  na.omit() |>
  _[, fit_segmented(sam ~ `1979--1999` + `2000--2014`), by = .(forcing, season(time))] |>
  _[, expand := qt(1 - 0.025, df)]

model_mean <- sam |>
  _[, .(sam = mean(sam, na.rm = TRUE)), by = .(model, forcing, time = seasonally(time))] |>
  _[, let(
    `1979--1999` = year(time),
    `2000--2014` = pmax(year(time) - cut_year, 0),
    year = year(time)
  )] |>
  na.omit() |>
  _[, fit_segmented(sam ~ `1979--1999` + `2000--2014`), by = .(model, forcing, season(time))] |>
  _[, expand := qt(1 - 0.025, df)]

fwrite(model_mean, "cmip-damip/plot-data/model-mean_trend.csv")

fwrite(mmm_trend, "cmip-damip/plot-data/mmm_trend.csv")


# plot --------------------------------------------------------------------

model_mean |>
  ggplot(aes(forcing, estimate * 10)) +
  geom_tile(
    data = mmm_trend, aes(
      height = se * expand * 2 * 10,
      fill = period
    ),
    width = 1, position = position_dodge2()
  ) +
  geom_tile(
    data = mmm_trend,
    width = 1, height = 0, color = "gray90", linewidth = 2,
    position = position_dodge2()
  ) +
  geom_jitter(aes(fill = period),
    shape = 21, size = 3,
    position = position_jitterdodge(dodge.width = 1, jitter.width = 0.1, seed = 42),
    show.legend = FALSE
  ) +
  geom_hline(yintercept = 0, linewidth = 0.5) +
  scale_color_brewer(
    palette = "Dark2", aesthetics = c("fill", "color"),
    labels = c(
      "1979–1999",
      "2000–2014"
    )
  ) +
  labs(
    x = NULL,
    y = "Trend",
    fill = NULL
  ) +
  guides(fill = guide_legend(position = "inside", direction = "horizontal")) +
  facet_wrap(~season) +
  theme_minimal(base_size = 16) +
  theme(
    legend.position.inside = c(0.7, 0.1),
    panel.background = element_rect(color = NA, fill = "#FAFAFA")
  )

# Save all the data for the plot
ggdatasaver::save_plot_data(last_plot(), name = "cmip-boxplot", dir = "cmip-damip/plot-data")
