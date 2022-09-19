#------------------------------------------------------------------------------
# build_lipophilic_tk.R
#
# This R source code file reads the file "lipophilic_tk.model" and creates the
# files "lipophilic_tk_model.c" (a C implementation of the model),
# "lipophilic_tk_model.o", (a C object file), "lipophilic_tk_model_inits.R" (an
# R source code file that defines initialization functions for the model), and
# a compiled version of the model called either "lipophilic_tk_model.dll" (on a
# Windows system) or "lipophilic_tk_model.so" (on a Unix system) that can be
# used by the R package "deSolve".
#
# Author: Dustin Kapraun, U.S. EPA, January 2022.
#------------------------------------------------------------------------------

# Set working directory to the directory containing this file.
#script.dir = dirname(sys.frame(1)$ofile)
#setwd(script.dir)

# Load relevant functions.
source("RMCSim.R")

# Load the lipophilic TK model.
compile_model("lipophilic_tk")