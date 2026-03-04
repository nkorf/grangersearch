#!/usr/bin/env Rscript
# Generate CRAN download statistics plot for README (base R graphics)

if (!requireNamespace("cranlogs", quietly = TRUE)) {
  install.packages("cranlogs", repos = "https://cloud.r-project.org")
}

library(cranlogs)

# Get daily downloads since first available date
downloads <- cran_downloads(
  packages = "grangersearch",
  from     = "2025-01-01",
  to       = Sys.Date()
)

# Remove days with zero downloads before the package was on CRAN
first_nonzero <- min(which(downloads$count > 0))
downloads <- downloads[first_nonzero:nrow(downloads), ]

# Compute cumulative downloads
downloads$cumulative <- cumsum(downloads$count)
total <- sum(downloads$count)

# Output file
outpath <- file.path("man", "figures", "cran-downloads.png")
png(outpath, width = 1200, height = 525, res = 150, bg = "white")

par(mar = c(4, 4.5, 3, 1))

# Cumulative area plot (extra y headroom for the total label)
plot(downloads$date, downloads$cumulative,
     type = "n",
     xlab = "", ylab = "Cumulative Downloads",
     main = "grangersearch CRAN Downloads",
     ylim = c(0, total * 1.15),
     las  = 1,
     bty  = "l",
     xaxt = "n")

# Shaded area
polygon(
  c(downloads$date, rev(downloads$date)),
  c(downloads$cumulative, rep(0, nrow(downloads))),
  col = adjustcolor("#2c7bb6", alpha.f = 0.25),
  border = NA
)

# Line
lines(downloads$date, downloads$cumulative, col = "#2c7bb6", lwd = 2)

# X-axis with month labels
month_seq <- seq(
  from = as.Date(format(min(downloads$date), "%Y-%m-01")),
  to   = max(downloads$date),
  by   = "month"
)
axis(1, at = month_seq, labels = format(month_seq, "%b %Y"))

# Total marker at end point
last_date <- downloads$date[nrow(downloads)]
last_cum  <- downloads$cumulative[nrow(downloads)]
points(last_date, last_cum, pch = 19, col = "#d7191c", cex = 1.3)
text(last_date, last_cum,
     labels = format(total, big.mark = ","),
     col = "#d7191c", font = 2, cex = 1.1,
     pos = 3, offset = 0.6)

# Subtitle and caption
mtext(paste0(min(downloads$date), " to ", max(downloads$date)),
      side = 3, line = 0, cex = 0.85, col = "grey40")
mtext("Source: cranlogs", side = 1, line = 2.8, adj = 1, cex = 0.7, col = "grey50")

dev.off()
message("Saved: ", outpath)
