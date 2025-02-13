
# Function to compute estimates and standard errors of the slopes
# of segmented regression.
# fit_segmented <- function(formula, data = NULL) {
#   mod <- lm(formula, data = data)
#   coef <- cumsum(coef(mod)[-1])
#   terms <- names(coef)
#   
#   list(
#     period = terms,
#     estimate = coef,
#     se = c(
#       sqrt(vcov(mod)[2:2, 2:2]),
#       sqrt(sum(vcov(mod)[2:3, 2:3]))
#     ),
#     df = rep(mod$df.residual, 2)
#   )
# }
# 
fit_segmented_ <- function(time, sam, cut_year) {
  data <- data.table(time = time, sam = sam) |> 
    copy() |> 
    _[, let(
      first_period = year(time),
      second_period = pmax(year(time) - cut_year, 0)
    )] |>
    na.omit()
  
  mod <- lm(sam ~ first_period + second_period, data = data)
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

fit_segmented <- memoise::memoise(fit_segmented_)

fir_normal_impl <- function(time, sam) {
  data <- data.table(time = year(time), sam = sam)
  mod <- lm(sam ~ time, data = data)
  coef <- coef(mod)[-1]
  
  list(
    estimate = coef,
    se = summary(mod)$coefficients[2, 2],
    df = mod$df.residual
  )
}

fit_normal_ <- function(time, sam, cut_year) {
  data <- data.table(time = time, sam = sam) |> 
    copy() |> 
    _[, let(period = ifelse(year(time) < cut_year, "first_period", "second_period"))] |>
    na.omit() |> 
    _[, fir_normal_impl(time, sam), by = period]
  
}

fit_normal <- memoise::memoise(fit_normal_)

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

select_models_ <- function(data, common_models = TRUE) {
  if (common_models) {
    models <- data |> 
      dcast(model ~ forcing, fill = NA) |> 
      na.omit() |> 
      _[, model]
    
    data <- data[model %in% models]
  }
  
  return(data)
}

select_models <- memoise::memoise(select_models_)

select_season <- function(data, seasons = c("DJF", "MAM", "JJA", "SON")) {
  data[season %in% seasons]
}


plot_sam_ <- function(data,
                      obs,
                      cut_year = 1999,
                      continuous = TRUE) {
  
  if (continuous) {
    fit <- fit_segmented
    trends <- geom_smooth(se = FALSE, linewidth = 0.8, 
                          method = "lm", formula = y ~ x + pmax(x - cut_year, 0))
  } else {
    fit <- fit_normal
    trends <- geom_smooth(aes(group = interaction(forcing, year(time) < cut_year)), 
                          se = FALSE, method = "lm", formula = y ~ x) 
  }
  # browser()
  model_mean <- data |>
    _[, .(sam = mean(sam, na.rm = TRUE)), by = .(model, forcing, time = seasonally(time))]
  
  mmm_mean <- model_mean |>
    _[, .(sam = mean(sam, na.rm = TRUE)), by = .(forcing, time = seasonally(time))] 
  
  mmm_trend <- mmm_mean |> 
    _[, fit(time, sam, cut_year), by = .(forcing, season(time))] |>
    _[, expand := qt(1 - 0.025, df)]
  
  obs_trend <- obs |> 
    _[, fit(time, sam, cut_year), by = .(forcing, season(time))] |> 
    _[, expand := qt(1 - 0.025, df)]
  
  model_mean <- model_mean |>
    _[, fit(time, sam, cut_year), by = .(forcing, model, season(time))] |>
    _[, expand := qt(1 - 0.025, df)]
  
  
  trend_labels <- c(first_period = paste0("1979--", cut_year),
                    second_period = paste0(cut_year + 1, "--2014"))
  
  boxplot <- model_mean |>
    ggplot(aes(forcing, estimate * 10)) +
    geom_tile(
      data = mmm_trend, aes(
        height = se * expand * 2 * 10,
        fill = period
      ),
      width = 1, position = position_dodge2(width = 0.75)
    ) +
    geom_tile(
      data = mmm_trend,
      width = 1, height = 0, color = "gray90", linewidth = 2,
      position = position_dodge2(width = 0.75)
    ) +
    geom_jitter(aes(fill = period),
                shape = 21, size = 3,
                position = position_jitterdodge(dodge.width = 1, jitter.width = 0.1, seed = 42),
                show.legend = FALSE
    ) +
    geom_errorbar(data = obs_trend, 
                  position = position_dodge(width = 0.75),
                  aes(color = period,
                      ymax = (estimate + se * expand) * 10,
                      ymin = (estimate - se * expand) * 10), 
                  width = 0.25) +
    geom_point(data = obs_trend, shape = 23,
               position = position_dodge(width = 0.75),
               size = 4, aes(fill = period)) +
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
  
  # browser()
  
  mmm_mean <- rbind(obs[, .(time, sam, forcing = "ERA5")], mmm_mean)
  
  lineplot <- mmm_mean |> 
    ggplot(aes(year(time), sam, color = forcing)) +
    geom_line(alpha = 0.6) +
    trends + 
    facet_wrap(~season(time)) +
    scale_color_brewer(NULL, palette = "Set1") +
    scale_x_continuous(NULL) +
    scale_y_continuous(NULL) +
    guides(color = guide_legend(position = "bottom")) +
    theme_minimal(base_size = 16) +
    theme(
      legend.position.inside = c(0.7, 0.1),
      panel.background = element_rect(color = NA, fill = "#FAFAFA")
    )
  
  
  list(boxplot = boxplot, 
       lineplot = lineplot)
  
}
plot_sam <- memoise::memoise(plot_sam_)
