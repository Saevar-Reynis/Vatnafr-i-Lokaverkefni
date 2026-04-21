q_path <- "lamah_ice/D_gauges/2_timeseries/daily/ID_37.csv"
m_path <- "lamah_ice/A_basins_total_upstrm/2_timeseries/daily/meteorological_data/ID_37.csv"
hydro_path <- "lamah_ice/D_gauges/1_attributes/hydro_indices_1981_2018.csv"

# Valinn atburður:
# Flóðlotan í lok febrúar / byrjun mars 2013.
# Hún er notuð hér vegna þess að hún sést skýrt í gögnunum við ID 37
# og tengist opinberum flóðaviðvörunum á vatnasviði Hvítár og Ölfusár.
plot_start <- as.Date("2013-02-20")
plot_end <- as.Date("2013-03-10")
event_start <- as.Date("2013-02-25")
peak_date <- as.Date("2013-03-01")
main_rain_end <- as.Date("2013-03-02")

q <- read.csv(q_path, sep = ";", stringsAsFactors = FALSE)
names(q) <- c("year", "month", "day", "qobs", "qc_flag")
q$qobs <- as.numeric(q$qobs)
q$date <- as.Date(sprintf("%04d-%02d-%02d", q$year, q$month, q$day))

m <- read.csv(m_path, sep = ";", stringsAsFactors = FALSE)
m$date <- as.Date(sprintf("%04d-%02d-%02d", m$YYYY, m$MM, m$DD))

hydro <- read.csv(hydro_path, sep = ";", stringsAsFactors = FALSE)
hydro_37 <- subset(hydro, id == 37)

df <- merge(
  q[, c("date", "qobs")],
  m[, c("date", "prec", "X2m_temp_mean")],
  by = "date",
  all.x = TRUE
)

event_df <- subset(df, date >= plot_start & date <= plot_end)

pre_event_flow <- subset(df, date == (event_start - 1))$qobs[1]
peak_flow <- subset(df, date == peak_date)$qobs[1]
peak_temp <- subset(df, date == peak_date)$X2m_temp_mean[1]
peak_prec <- subset(df, date == peak_date)$prec[1]

return_date <- subset(df, date > peak_date & qobs <= pre_event_flow)$date[1]

time_to_peak_days <- as.numeric(peak_date - event_start)
excess_rain_release_days <- as.numeric(return_date - main_rain_end)
recession_days <- as.numeric(return_date - peak_date)

plot_labels <- list(
  title = "Liður 8: Rennslisatburður í Hvítá við Kljáfoss",
  flow_rise = "Upphaf rennslisaukningar",
  peak = "Qpeak",
  rain_end = "Lok meginúrkoma",
  return_baseflow = "Til baka á grunnstig",
  precipitation = "Úrkoma",
  date = "Dagsetning",
  temperature = "2 m meðalhiti"
)

date_axis_at <- seq(plot_start, plot_end, by = "day")
date_axis_labels <- format(date_axis_at, "%d.%m")

png_args <- list(
  filename = "event_analysis_id37_2013.png",
  width = 1200,
  height = 1000,
  res = 150,
  family = "sans"
)

do.call(png, png_args)
par(mfrow = c(3, 1), mar = c(6, 5, 3, 2), oma = c(2, 0, 0, 0))

plot(
  event_df$date,
  event_df$qobs,
  type = "o",
  pch = 19,
  col = "black",
  xaxt = "n",
  xlab = "",
  ylab = "Q (m3/s)",
  main = plot_labels$title
)
axis.Date(1, at = date_axis_at, labels = date_axis_labels, las = 2, cex.axis = 0.75)
abline(v = event_start, col = "darkorange", lty = 2, lwd = 2)
abline(v = peak_date, col = "firebrick", lty = 2, lwd = 2)
abline(v = main_rain_end, col = "steelblue", lty = 2, lwd = 2)
abline(v = return_date, col = "darkgreen", lty = 2, lwd = 2)
abline(h = pre_event_flow, col = "grey60", lty = 3)
legend(
  "topright",
  legend = c(
    "Q",
    plot_labels$flow_rise,
    plot_labels$peak,
    plot_labels$rain_end,
    plot_labels$return_baseflow
  ),
  col = c("black", "darkorange", "firebrick", "steelblue", "darkgreen"),
  lty = c(1, 2, 2, 2, 2),
  pch = c(19, NA, NA, NA, NA),
  lwd = c(1, 2, 2, 2, 2),
  bty = "n"
)

barplot(
  event_df$prec,
  names.arg = format(event_df$date, "%m-%d"),
  col = "steelblue",
  border = NA,
  ylab = "P (mm/dag)",
  main = plot_labels$precipitation
)

plot(
  event_df$date,
  event_df$X2m_temp_mean,
  type = "o",
  pch = 19,
  col = "firebrick",
  xaxt = "n",
  xlab = plot_labels$date,
  ylab = "T (degC)",
  main = plot_labels$temperature
)
axis.Date(1, at = date_axis_at, labels = date_axis_labels, las = 2, cex.axis = 0.75)
abline(h = 0, col = "grey60", lty = 3)
abline(v = event_start, col = "darkorange", lty = 2, lwd = 2)
abline(v = peak_date, col = "firebrick", lty = 2, lwd = 2)
abline(v = main_rain_end, col = "steelblue", lty = 2, lwd = 2)
abline(v = return_date, col = "darkgreen", lty = 2, lwd = 2)

dev.off()

cat("\nID 37 - Hvita, Kljafoss\n")
cat("Valinn atburdur:", format(event_start), "til", format(return_date), "\n")
cat("Plottad timabil:", format(plot_start), "til", format(plot_end), "\n\n")

cat("Atburdaskilgreining:\n")
cat("- Rennslisaukning byrjar:", format(event_start), "\n")
cat("- Qpeak:", format(peak_date), "(", round(peak_flow, 2), "m3/s )\n")
cat("- Lok meginurkoma:", format(main_rain_end), "\n")
cat("- Rennsli aftur likt og fyrir atburd:", format(return_date), "\n")
cat("- Rennsli fyrir atburd:", round(pre_event_flow, 2), "m3/s\n\n")

cat("Timametrar:\n")
cat("- Time-to-peak:", time_to_peak_days, "dagar\n")
cat("- Excess rain release time:", excess_rain_release_days, "dagar\n")
cat("- Recession time:", recession_days, "dagar\n\n")

cat("Veður og vatnafar vid topp:\n")
cat("- Urkoma a toppdegi:", round(peak_prec, 2), "mm/dag\n")
cat("- Medalhiti a toppdegi:", round(peak_temp, 2), "degC\n\n")

cat("Tenging vid hydro-indices fyrir ID 37 (1981-2018):\n")
cat("- baseflow_index_ladson:", round(hydro_37$baseflow_index_ladson, 3), "\n")
cat("- slope_fdc:", round(hydro_37$slope_fdc, 3), "\n\n")

cat("Mynd vistud sem:\n")
cat("- event_analysis_id37_2013.png\n")
