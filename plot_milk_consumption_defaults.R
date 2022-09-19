#------------------------------------------------------------------------------
# plot_milk_consumption_defaults.R
#
# This file contains an R script for plotting default milk consumption rates
# for various animal species.
#
# Author: Dustin Kapraun, U.S. EPA, September 2020.
#------------------------------------------------------------------------------

# Save plots to files in the output directory.
write_output = TRUE

# Set working directory to the directory containing this file.
#script.dir = dirname(sys.frame(1)$ofile)
#setwd(script.dir)

# Load simulation functions.
source("sim_functions.R")

#
# Plot rat milk consumption rate.
#

# Define simulation parameters.
p = c(A_m_init = 100)     # Initial amount (mg) in mother.
p = init_rat_p(p)

# Define times for simulation.
times = seq(from=0, to=50, by=0.1)

# Run simulation.
out <- rat_sim(times=times, p=p)

# Find row index of the output data frame that corresponds to the instant of
# birth.
birth_idx = match(p[["t_gest"]], out[, 1])

# Plot milk consumption rate vs. time.
ymax = max(out[, "R_milk"])
if (write_output) {
    file_name = "./Output/rat_milk_consumption.tif"
    tiff(file_name, units="in", width=6, height=5, res=300)
}
plot(out[birth_idx:nrow(out), 1] - p[["t_gest"]],
     out[birth_idx:nrow(out), "R_milk"],
     type='l', lty=1, col=vcol[1], lwd=2,
     xlab="Time Since Birth (d)",
     ylab="Consumption Rate (kg/d)",
     # xlim=c(20, max(times)),
     ylim=c(0, ymax))
# lines(c(p[["t_gest"]], p[["t_gest"]]), c(-20, 1.5*ymax), lty=3,
#       col="black", lwd=2)
if (write_output) {
    dev.off()
}


#
# Plot human milk consumption rate.
#

# Define simulation parameters.
p = c(A_m_init = 100)     # Initial amount (mg) in mother.
p = init_human_p(p)

# Define times for simulation.
t_start = 0
# t_end = 1.75 * 365
t_end = 2.0 * 365
times = seq(from=t_start, to=t_end, by=0.25)

# Run simulation.
out <- human_sim(times=times, p=p)

# Find row index of the output data frame that corresponds to the instant of
# birth.
birth_idx = match(p[["t_gest"]], out[, 1])

# Plot milk consumption rate vs. time.
ymax = max(out[, "R_milk"])
if (write_output) {
    file_name = "./Output/human_milk_consumption.tif"
    tiff(file_name, units="in", width=6, height=5, res=300)
}
plot(out[birth_idx:nrow(out), 1] - p[["t_gest"]],
     out[birth_idx:nrow(out), "R_milk"],
     type='l', lty=1, col=vcol[1], lwd=2,
     xlab="Time Since Birth (d)",
     ylab="Consumption Rate (kg/d)",
     # xlim=c(260, max(times)),
     ylim=c(0, ymax))
# lines(c(p[["t_gest"]], p[["t_gest"]]), c(-20, 1.5*ymax), lty=3,
#       col="black", lwd=2)
if (write_output) {
    dev.off()
}