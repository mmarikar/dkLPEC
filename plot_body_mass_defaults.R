#------------------------------------------------------------------------------
# plot_body_mass_defaults.R
#
# This file contains an R script for plotting default body mass functions for
# various animal species.
#
# Author: Dustin Kapraun, U.S. EPA, September 2020.
#------------------------------------------------------------------------------

# Save plots to files in the output directory.
write_output = TRUE

# Set working directory to the directory containing this file.
script.dir = dirname(sys.frame(1)$ofile)
setwd(script.dir)

# Load simulation functions.
source("sim_functions.R")

#
# Plot rat body mass functions.
#

# Define simulation parameters.
p = c(A_m_init = 100)     # Initial amount (mg) in mother.
p = init_rat_p(p)

# Define times for simulation.
times = seq(from=-20, to=100, by=0.1)

# Run simulation.
out <- rat_sim(times=times, p=p)

# Find row index of the output data frame that corresponds to the instant of
# birth.
birth_idx = match(p[["t_gest"]], out[, 1])

# Plot maternal and infant mass vs. time on the same axes.
ymax = max(c(out[, "M_mf"], out[, "M_i"]))
if (write_output) {
    file_name = "./Output/rat_body_mass.tif"
    tiff(file_name, units="in", width=6, height=5, res=300)
}
plot(out[, 1], out[, "M_mf"], type='l', lty=1, col=vcol[1], lwd=2,
     ylim=c(0, ymax),
     xlab="Time Since Conception (d)",
     ylab="Mass (kg)")
lines(out[birth_idx:nrow(out), 1], out[birth_idx:nrow(out), "M_i"], lty=2,
      col=vcol[2], lwd=2)
lines(c(0, 0), c(-20, 1.5*ymax), lty=3, col="black", lwd=2)
lines(c(p[["t_gest"]], p[["t_gest"]]), c(-20, 1.5*ymax), lty=3, col="black",
      lwd=2)
legend("bottomright", c("mother", "infant"), lty=c(1, 2), col=vcol[1:2],
       lwd=c(2, 2))
if (write_output) {
    dev.off()
}


# # Plot maternal and infant mass vs. time on the same axes.
# ymax = max(c(out[, "M_mf"], out[, "M_i"]))
# if (write_output) {
#   file_name = "./Output/rat_body_mass_default.tif"
#   tiff(file_name, units="in", width=6, height=5, res=300)
# }
# plot(out[, 1], out[, "M_mf"], type='l', lty=1, col=vcol[1], lwd=2,
#      ylim=c(0, ymax),
#      xlab="Time Since Conception (d)",
#      ylab="Mass (kg)")
# points(c(p["t_m_1"], p["t_m_2"]), c(p["M_m_1"], p["M_m_2"]))
# lines(out[birth_idx:nrow(out), 1], out[birth_idx:nrow(out), "M_i"], lty=2,
#       col=vcol[2], lwd=2)
# lines(c(0, 0), c(-20, 1.5*ymax), lty=3, col="black", lwd=2)
# lines(c(p[["t_gest"]], p[["t_gest"]]), c(-20, 1.5*ymax), lty=3, col="black",
#       lwd=2)
# legend("bottomright", c("mother", "infant"), lty=c(1, 2), col=vcol[1:2],
#        lwd=c(2, 2))
# if (write_output) {
#   dev.off()
# }


#
# Plot human body mass functions.
#

# Define simulation parameters.
p = c(A_m_init = 100)     # Initial amount (mg) in mother.
p = init_human_p(p)

# Define times for simulation.
t_start = -24.25 * 365
t_end = 78.75 * 365
times = seq(from=t_start, to=t_end, by=0.25)

# Run simulation.
out <- human_sim(times=times, p=p)

# Find row index of the output data frame that corresponds to the instant of
# birth.
birth_idx = match(p[["t_gest"]], out[, 1])

# Plot maternal and infant mass vs. time on the same axes.
ymax = max(c(out[, "M_mf"], out[, "M_i"]))
if (write_output) {
    file_name = "./Output/human_body_mass.tif"
    tiff(file_name, units="in", width=6, height=5, res=300)
}
plot(out[, 1] / 365, out[, "M_mf"], type='l', lty=1, col=vcol[1], lwd=2,
     ylim=c(0, ymax),
     xlab="Time Since Conception (y)",
     ylab="Mass (kg)")
lines(out[birth_idx:nrow(out), 1] / 365, out[birth_idx:nrow(out), "M_i"],
      lty=2, col=vcol[2], lwd=2)
lines(c(0, 0), c(-20, 1.5*ymax), lty=3, col="black", lwd=2)
lines(c(p[["t_gest"]], p[["t_gest"]]) / 365, c(-20, 1.5*ymax), lty=3,
      col="black", lwd=2)
legend("bottomright", c("mother", "infant"), lty=c(1, 2), col=vcol[1:2],
       lwd=c(2, 2))
if (write_output) {
    dev.off()
}