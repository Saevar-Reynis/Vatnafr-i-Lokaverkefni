data_path <- "lamah_ice/D_gauges/2_timeseries/daily/ID_37.csv"
start_date <- as.Date("1993-10-01")
end_date <- as.Date("2023-09-30")

df <- read.csv(data_path, sep = ";", stringsAsFactors = FALSE)
names(df) <- c("year", "month", "day", "qobs", "qc_flag")
df$qobs <- as.numeric(df$qobs)
df$date <- as.Date(sprintf("%04d-%02d-%02d", df$year, df$month, df$day))

df <- subset(df, date >= start_date & date <= end_date & !is.na(qobs))
df$analysis_year <- ifelse(df$month >= 10, df$year + 1, df$year)

mk_s <- function(x) {
  n <- length(x)
  s <- 0

  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      s <- s + sign(x[j] - x[i])
    }
  }

  s
}

mk_var_ties <- function(x) {
  n <- length(x)
  ties <- table(x)
  tie_term <- sum(ties * (ties - 1) * (2 * ties + 5))
  (n * (n - 1) * (2 * n + 5) - tie_term) / 18
}

theil_sen <- function(t, x) {
  slopes <- numeric(0)
  n <- length(x)

  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      slopes <- c(slopes, (x[j] - x[i]) / (t[j] - t[i]))
    }
  }

  median(slopes)
}

modified_mk <- function(t, x, alpha = 0.05) {
  n <- length(x)
  slope <- theil_sen(t, x)
  detrended <- x - slope * t
  ranked <- rank(detrended, ties.method = "average")
  lags <- 1:(n - 1)

  acf_vals <- sapply(lags, function(lag) {
    if ((n - lag) < 2) {
      return(0)
    }

    out <- suppressWarnings(cor(ranked[1:(n - lag)], ranked[(lag + 1):n]))
    ifelse(is.na(out), 0, out)
  })

  sig_limit <- qnorm(1 - alpha / 2) / sqrt(n)
  sig_acf <- ifelse(abs(acf_vals) > sig_limit, acf_vals, 0)
  weights <- (n - lags) * (n - lags - 1) * (n - lags - 2)
  correction <- 1 + (2 / (n * (n - 1) * (n - 2))) * sum(weights * sig_acf)
  correction <- max(1, correction)

  s <- mk_s(x)
  var_s <- mk_var_ties(x) * correction

  z <- if (s > 0) {
    (s - 1) / sqrt(var_s)
  } else if (s < 0) {
    (s + 1) / sqrt(var_s)
  } else {
    0
  }

  p_value <- 2 * (1 - pnorm(abs(z)))
  tau <- s / (0.5 * n * (n - 1))

  data.frame(
    n = n,
    slope = slope,
    tau = tau,
    z = z,
    p_value = p_value,
    correction = correction
  )
}

annual <- aggregate(qobs ~ analysis_year, data = df, FUN = mean)
annual_result <- modified_mk(annual$analysis_year, annual$qobs)

monthly <- aggregate(qobs ~ analysis_year + month, data = df, FUN = mean)
monthly_results <- do.call(
  rbind,
  lapply(split(monthly, monthly$month), function(d) {
    out <- modified_mk(d$analysis_year, d$qobs)
    out$month <- unique(d$month)
    out
  })
)

monthly_results <- monthly_results[, c("month", "n", "slope", "tau", "z", "p_value", "correction")]
monthly_results <- monthly_results[order(monthly_results$month), ]
significant_months <- subset(monthly_results, p_value < 0.05)

annual_fit <- median(annual$qobs) + annual_result$slope[1] * (annual$analysis_year - median(annual$analysis_year))

png("trend_analysis_annual_id37.png", width = 1200, height = 700, res = 150)
plot(
  annual$analysis_year,
  annual$qobs,
  pch = 19,
  col = "black",
  xlab = "Ár",
  ylab = "Ársmeðalrennsli (m3/s)",
  main = "Leitnigreining í Hvítá"
)
lines(annual$analysis_year, annual_fit, col = "firebrick", lwd = 2)
legend(
  "topright",
  legend = c("Ársmeðalrennsli", "Theil-Sen trend"),
  col = c("black", "firebrick"),
  pch = c(19, NA),
  lty = c(NA, 1),
  lwd = c(NA, 2),
  bty = "n"
)
dev.off()

cat("\nID 37 - Hvita, Kljafoss\n")
cat("Timabil:", format(start_date), "til", format(end_date), "\n\n")

cat("Arsgrundur (medalrennsli per vatnaar):\n")
print(round(annual_result, 6))

cat("\nManadargrundur (medalrennsli fyrir hvern manud yfir 30 ar):\n")
print(round(monthly_results, 6))

cat("\nMarktaekir manudir vid p < 0.05:\n")
if (nrow(significant_months) == 0) {
  cat("Enginn manudur var marktakur.\n")
} else {
  print(round(significant_months, 6))
}

cat("\nMynd vistud sem:\n")
cat("- trend_analysis_annual_id37.png\n")
