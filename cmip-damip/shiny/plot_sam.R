
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


seasonally <- function (x) {
  if (is.character(x)) 
    x <- as.Date(x)
  times <- unique(x)
  m <- data.table::month(times)
  times_org <- times
  lubridate::year(times[m == 12]) <- data.table::year(times[m == 12]) + 1
  s <- season(times)
  lubridate::day(times) <- 15
  lubridate::month(times) <- (as.numeric(s) - 1) * 3 + 1
  as.Date(times[match(x, times_org)])
}

season <- function (x, lang = c("en", "es")) {
  x <- lubridate::month(x)
  if (lang[1] == "en") {
    djf <- "DJF"
  }
  else {
    djf <- "DEF"
  }
  jja <- "JJA"
  mam <- "MAM"
  son <- "SON"
  seasons <- c(djf, djf, rep(c(mam, jja, son), each = 3), 
               djf)
  return(factor(seasons[x], levels = c(djf, mam, jja, son)))
}


plot_sam <- function(data,
                     seasons = c("DJF", "MAM", "JJA", "SON"),
                     cut_year = 1999,
                     common_models = TRUE) {
  
  
  data <- data[season(time) %in% seasons]
  
  if (common_models) {
    models <- data |> 
      dcast(model ~ forcing, fill = NA) |> 
      na.omit() |> 
      _[, model]
    
    data <- data[model %in% models]
  }
  
  mmm_trend <- data |>
    _[, .(sam = mean(sam, na.rm = TRUE)), by = .(model, forcing, time = seasonally(time))] |>
    _[, .(sam = mean(sam, na.rm = TRUE)), by = .(forcing, time = seasonally(time))] |>
    _[, let(
      first_period = year(time),
      second_period = pmax(year(time) - cut_year, 0),
      year = year(time)
    )] |>
    na.omit() |>
    _[, fit_segmented(sam ~ first_period + second_period), by = .(forcing, season(time))] |>
    _[, expand := qt(1 - 0.025, df)]
  
  
  model_mean <- data |>
    _[, .(sam = mean(sam, na.rm = TRUE)), by = .(model, forcing, time = seasonally(time))] |>
    _[, let(
      first_period = year(time),
      second_period = pmax(year(time) - cut_year, 0),
      year = year(time)
    )] |>
    na.omit() |>
    _[, fit_segmented(sam ~ first_period + second_period), by = .(model, forcing, season(time))] |>
    _[, expand := qt(1 - 0.025, df)]
  
  
  trend_labels <- c(first_period = paste0("1979--", cut_year),
                    second_period = paste0(cut_year + 1, "--2014"))
  
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
      labels = trend_labels
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
}
