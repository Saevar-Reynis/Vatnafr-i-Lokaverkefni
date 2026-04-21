data_path <- "lamah_ice/D_gauges/2_timeseries/daily/ID_37.csv"
start_date <- as.Date("1993-10-01")
end_date <- as.Date("2023-09-30")
bootstrap_reps <- 3000
set.seed(37)

df <- read.csv(data_path, sep = ";", stringsAsFactors = FALSE)
names(df) <- c("year", "month", "day", "qobs", "qc_flag")
df$qobs <- as.numeric(df$qobs)
df$date <- as.Date(sprintf("%04d-%02d-%02d", df$year, df$month, df$day))

df <- subset(df, date >= start_date & date <= end_date & !is.na(qobs))
df$water_year <- ifelse(df$month >= 10, df$year + 1, df$year)

annual_peaks <- do.call(
  rbind,
  lapply(split(df, df$water_year), function(d) {
    d[which.max(d$qobs), c("water_year", "date", "month", "qobs")]
  })
)

annual_peaks <- annual_peaks[order(annual_peaks$water_year), ]
row.names(annual_peaks) <- NULL

month_counts <- table(factor(annual_peaks$month, levels = 1:12))

x <- annual_peaks$qobs
n <- length(x)
xs <- sort(x)
gringorten_p <- ((1:n) - 0.44) / (n + 0.12)

negloglik_gumbel <- function(par, x) {
  mu <- par[1]
  beta <- par[2]

  if (!is.finite(beta) || beta <= 0) {
    return(1e100)
  }

  z <- (x - mu) / beta
  length(x) * log(beta) + sum(z + exp(-z))
}

fit_gumbel <- function(x) {
  m <- mean(x)
  s <- sd(x)
  beta0 <- s * sqrt(6) / pi
  mu0 <- m - 0.5772156649 * beta0

  fit <- optim(
    c(mu0, beta0),
    negloglik_gumbel,
    x = x,
    method = "L-BFGS-B",
    lower = c(-Inf, 1e-8)
  )

  fit$par
}

qgumbel <- function(p, par) {
  mu <- par[1]
  beta <- par[2]
  mu - beta * log(-log(p))
}

negloglik_ln3 <- function(par, x) {
  theta <- par[1]

  if (theta >= min(x)) {
    return(1e100)
  }

  y <- log(x - theta)
  mu <- mean(y)
  sigma <- sqrt(mean((y - mu)^2))

  if (!is.finite(sigma) || sigma <= 0) {
    return(1e100)
  }

  -sum(dnorm(y, mean = mu, sd = sigma, log = TRUE) - log(x - theta))
}

fit_ln3 <- function(x) {
  upper <- min(x) - 1e-6
  lower <- upper - 10 * sd(x)

  fit <- optim(
    (lower + upper) / 2,
    negloglik_ln3,
    x = x,
    method = "L-BFGS-B",
    lower = lower,
    upper = upper
  )

  theta <- fit$par[1]
  y <- log(x - theta)
  mu <- mean(y)
  sigma <- sqrt(mean((y - mu)^2))

  c(theta = theta, meanlog = mu, sdlog = sigma)
}

qln3 <- function(p, par) {
  par[1] + qlnorm(p, meanlog = par[2], sdlog = par[3])
}

dpearson3 <- function(y, mu, sigma, skew, log = FALSE) {
  if (abs(skew) < 1e-8) {
    return(dnorm(y, mean = mu, sd = sigma, log = log))
  }

  shape <- 4 / (skew^2)
  scale <- sigma / sqrt(shape)

  if (skew > 0) {
    loc <- mu - shape * scale
    dens <- dgamma(y - loc, shape = shape, scale = scale, log = log)

    if (log) {
      dens[y <= loc] <- -Inf
    } else {
      dens[y <= loc] <- 0
    }
  } else {
    loc <- mu + shape * scale
    dens <- dgamma(loc - y, shape = shape, scale = scale, log = log)

    if (log) {
      dens[y >= loc] <- -Inf
    } else {
      dens[y >= loc] <- 0
    }
  }

  dens
}

qpearson3 <- function(p, mu, sigma, skew) {
  if (abs(skew) < 1e-8) {
    return(qnorm(p, mean = mu, sd = sigma))
  }

  shape <- 4 / (skew^2)
  scale <- sigma / sqrt(shape)

  if (skew > 0) {
    loc <- mu - shape * scale
    loc + qgamma(p, shape = shape, scale = scale)
  } else {
    loc <- mu + shape * scale
    loc - qgamma(1 - p, shape = shape, scale = scale)
  }
}

negloglik_lp3 <- function(par, x) {
  mu <- par[1]
  sigma <- par[2]
  skew <- par[3]

  if (!is.finite(sigma) || sigma <= 0) {
    return(1e100)
  }

  y <- log10(x)
  ll <- dpearson3(y, mu, sigma, skew, log = TRUE) - log(x * log(10))

  if (any(!is.finite(ll))) {
    return(1e100)
  }

  -sum(ll)
}

fit_lp3 <- function(x) {
  y <- log10(x)
  mu0 <- mean(y)
  sigma0 <- sd(y)
  skew0 <- mean((y - mu0)^3) / (sigma0^3)

  fit <- optim(
    c(mu0, sigma0, skew0),
    negloglik_lp3,
    x = x,
    method = "L-BFGS-B",
    lower = c(-Inf, 1e-8, -10),
    upper = c(Inf, Inf, 10)
  )

  setNames(fit$par, c("meanlog10", "sdlog10", "skewlog10"))
}

qlp3 <- function(p, par) {
  10^(qpearson3(p, par[1], par[2], par[3]))
}

fit_g <- fit_gumbel(x)
fit_ln <- fit_ln3(x)
fit_lp <- fit_lp3(x)

fit_rmse <- c(
  Gumbel = sqrt(mean((xs - qgumbel(gringorten_p, fit_g))^2)),
  LogNormal3 = sqrt(mean((xs - qln3(gringorten_p, fit_ln))^2)),
  LogPearson3 = sqrt(mean((xs - qlp3(gringorten_p, fit_lp))^2))
)

best_fit <- names(which.min(fit_rmse))

quantile_function <- switch(
  best_fit,
  Gumbel = function(prob) qgumbel(prob, fit_g),
  LogNormal3 = function(prob) qln3(prob, fit_ln),
  LogPearson3 = function(prob) qlp3(prob, fit_lp)
)

return_periods <- c(10, 50, 100)
non_exceedance_probs <- 1 - 1 / return_periods
design_flows <- sapply(non_exceedance_probs, quantile_function)
names(design_flows) <- paste0("Q", return_periods)

bootstrap_quantiles <- matrix(NA_real_, nrow = bootstrap_reps, ncol = length(return_periods))
colnames(bootstrap_quantiles) <- paste0("Q", return_periods)

for (i in seq_len(bootstrap_reps)) {
  xb <- sample(x, replace = TRUE)

  bootstrap_quantiles[i, ] <- tryCatch(
    {
      if (best_fit == "Gumbel") {
        fitted <- fit_gumbel(xb)
        sapply(non_exceedance_probs, function(prob) qgumbel(prob, fitted))
      } else if (best_fit == "LogNormal3") {
        fitted <- fit_ln3(xb)
        sapply(non_exceedance_probs, function(prob) qln3(prob, fitted))
      } else {
        fitted <- fit_lp3(xb)
        sapply(non_exceedance_probs, function(prob) qlp3(prob, fitted))
      }
    },
    error = function(e) rep(NA_real_, length(return_periods))
  )
}

confidence_intervals <- apply(
  bootstrap_quantiles,
  2,
  function(v) quantile(v[is.finite(v)], probs = c(0.05, 0.95), na.rm = TRUE)
)

png("flood_frequency_fit_id37.png", width = 1200, height = 700, res = 150)
plot(
  xs,
  gringorten_p,
  pch = 19,
  col = "black",
  xlab = "Hæsta rennslisgildi (m3/s)",
  ylab = "Gringorten óyfirstignalíkur",
  main = "Flóðagreining í Hvítá",
  xlim = range(c(xs, design_flows)),
  ylim = c(0, 1)
)

p_seq <- seq(0.01, 0.99, length.out = 300)
lines(qgumbel(p_seq, fit_g), p_seq, col = "firebrick", lwd = 2)
lines(qln3(p_seq, fit_ln), p_seq, col = "steelblue", lwd = 2)
lines(qlp3(p_seq, fit_lp), p_seq, col = "darkgreen", lwd = 2)

legend(
  "bottomright",
  legend = c("Árlegir toppar", "Gumbel", "Log Normal 3", "Log Pearson 3"),
  col = c("black", "firebrick", "steelblue", "darkgreen"),
  pch = c(19, NA, NA, NA),
  lty = c(NA, 1, 1, 1),
  lwd = c(NA, 2, 2, 2),
  bty = "n"
)
dev.off()

cat("\nID 37 - Hvita, Kljafoss\n")
cat("Timabil:", format(start_date), "til", format(end_date), "\n")
cat("Fjoldi annual peaks:", n, "\n\n")

cat("Annual peaks eftir vatnaar:\n")
print(annual_peaks)

cat("\nFjoldi annual peaks i hverjum manudi:\n")
print(month_counts)

cat("\nRMSE fyrir hverja dreifingu:\n")
print(round(fit_rmse, 3))

cat("\nBest fit midad vid Gringorten plotting positions:\n")
cat(best_fit, "\n")

cat("\nQ10, Q50 og Q100 (m3/s):\n")
print(round(design_flows, 2))

cat("\n90% confidence interval (m3/s):\n")
print(round(confidence_intervals, 2))

cat("\nMyndir vistadar sem:\n")
cat("- flood_frequency_fit_id37.png\n")
