## ---- include = FALSE------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup-----------------------------------------------------------------
library(simpleEpiModel)
library(ggplot2)

## --------------------------------------------------------------------------
# Run the SEIR function
# This returns a dataframe, which gets stored in SEIR_df
# Each row in the dataframe corresponds to a single time step (e.g., day)
# Each column contains the value for that model parameter at that step
#
SEIR_flu <- SEIR(infection="influenza", N=10^6, duration=365, number_infectious=1, percent_susceptible=100)

# You can easily visualize the response of individual model parameters over hte course of the run using ggplot
#
plot_SEIR_flu <- ggplot(SEIR_flu) + geom_line(aes(x=t, y=S), color="green") + geom_line(aes(x=t, y=E), color="orange") + geom_line(aes(x=t, y=I), color="red") + geom_line(aes(x=t, y=R), color="blue")

# And visualize it
plot_SEIR_flu


