
#' Some built-in data sets for common infections
simple_infections <- data.frame(
  "infection_name" = c("influenza", "COVID"),
  "latent_period" = c(2,3), #D'
  "infectious_period" = c(2,10), #D
  "subclinical_period" = c(0,5), #w
  "basic_reproduction_number" = c(2,3.2), #R0
  "time_units" = c("days", "days"),
  "immune_period" = c(30, 30)
)

#' Generate an SIR model for an infection based on input parameters
#' @param infection The name of the infection to model c("influenza", "COVID").
#' @param duration The length of time to model, in days (optional; 120).
#' @param time_step The time between each model step, in days (optional; 1).
#' @param N The total population size (optional; 10^6).
#' @param percent_susceptible The initial susceptible population, as a percentage of N (optional; 100).
#' @param number_infectious The initial number of infectious individuals (optional; 1).
#' @param infection_data A list of key-value pairs with basic infection data. See simple_infections.
SIR <- function(infection, duration, time_step, N, percent_susceptible, number_infectious, infection_data) {

  model = "SIR"

  if (missing(infection)) infection = "influenza"
  if (missing(duration)) duration = 120
  if (missing(time_step)) time_step = 1
  if (missing(N)) N = 10^6
  if (missing(percent_susceptible)) percent_susceptible = 100
  if (missing(number_infectious)) number_infectious = 1

  if (missing(infection_data)) {
    if (infection %in% simple_infections$infection_name) {
      infection_data <- as.list(subset(simple_infections, infection_name==infection))
    } else {
      infection_data <- as.list(subset(simple_infections, infection_name=="influenza"))
    }
  }

  #
  # Need to check user-provided infection_data for required values
  #

  S = (percent_susceptible/100) * N
  I = number_infectious
  R = 0
  beta = infection_data$basic_reproduction_number/(N*infection_data$infectious_period)
  f = 1/infection_data$latent_period
  r = 1/infection_data$infectious_period
  t = 0

  #  print(paste0("beta = ", beta, sep=""))
  #  print(paste0("r = ", r, sep=""))
  #  print(paste0("f = ", f, sep=""))

  df <- data.frame(
    "t" = c(t),
    "t" = c(t),
    "S" = c(S),
    "I" = c(I),
    "R" = c(R)
  )
  steps <- seq(1, duration, by=time_step)
  for (step in steps) {
    St = S - beta*I*S
    It = I + f*S - r*I
    Rt = R + r*I

    df <- rbind(df, c(step, St, It, Rt))

    S = St
    I = It
    R = Rt
  }

  return(df);
}



#' Generate an SIRS model for an infection based on input parameters
#' @param infection The name of the infection to model c("influenza", "COVID").
#' @param duration The length of time to model, in days (optional; 120).
#' @param time_step The time between each model step, in days (optional; 1).
#' @param N The total population size (optional; 10^6).
#' @param percent_susceptible The initial susceptible population, as a percentage of N (optional; 100).
#' @param number_infectious The initial number of infectious individuals (optional; 1).
#' @param infection_data A list of key-value pairs with basic infection data. See simple_infections.
SIRS <- function(infection, duration, time_step, N, percent_susceptible, number_infectious, infection_data) {

  model = "SIRS"

  if (missing(infection)) infection = "influenza"
  if (missing(duration)) duration = 120
  if (missing(time_step)) time_step = 1
  if (missing(N)) N = 10^6
  if (missing(percent_susceptible)) percent_susceptible = 100
  if (missing(number_infectious)) number_infectious = 1

  if (missing(infection_data)) {
    if (infection %in% simple_infections$infection_name) {
      infection_data <- as.list(subset(simple_infections, infection_name==infection))
    } else {
      infection_data <- as.list(subset(simple_infections, infection_name=="influenza"))
    }
  }

  #
  # Need to check user-provided infection_data for required values
  #

  S = (percent_susceptible/100) * N
  I = number_infectious
  R = 0
  beta = infection_data$basic_reproduction_number/(N*infection_data$infectious_period)
  f = 1/infection_data$latent_period
  r = 1/infection_data$infectious_period
  gamma = 1/infection_data$immune_period
  t = 0

  #  print(paste0("beta = ", beta, sep=""))
  #  print(paste0("r = ", r, sep=""))
  #  print(paste0("f = ", f, sep=""))

  df <- data.frame(
    "t" = c(t),
    "S" = c(S),
    "I" = c(I),
    "R" = c(R)
  )
  steps <- seq(1, duration, by=time_step)
  for (step in steps) {
    St = S - beta*I*S + gamma*R
    It = I + f*S - r*I
    Rt = R + r*I - gamma*R

    df <- rbind(df, c(step, St, It, Rt))

    S = St
    I = It
    R = Rt
  }

  return(df);
}



#' Generate an SEIR model for an infection based on input parameters
#' @param infection The name of the infection to model c("influenza", "COVID").
#' @param duration The length of time to model, in days (optional; 120).
#' @param time_step The time between each model step, in days (optional; 1).
#' @param N The total population size (optional; 10^6).
#' @param percent_susceptible The initial susceptible population, as a percentage of N (optional; 100).
#' @param number_infectious The initial number of infectious individuals (optional; 1).
#' @param infection_data A list of key-value pairs with basic infection data. See simple_infections.
SEIR <- function(infection, duration, time_step, N, percent_susceptible, number_infectious, infection_data) {

  model = "SEIR"

  if (missing(infection)) infection = "influenza"
  if (missing(duration)) duration = 120
  if (missing(time_step)) time_step = 1
  if (missing(N)) N = 10^6
  if (missing(percent_susceptible)) percent_susceptible = 100
  if (missing(number_infectious)) number_infectious = 1

  if (missing(infection_data)) {
    if (infection %in% simple_infections$infection_name) {
      infection_data <- as.list(subset(simple_infections, infection_name==infection))
    } else {
      infection_data <- as.list(subset(simple_infections, infection_name=="influenza"))
    }
  }

  #
  # Need to check user-provided infection_data for required values
  #

  S = (percent_susceptible/100) * N
  E = 0
  I = number_infectious
  R = 0
  beta = infection_data$basic_reproduction_number/(N*infection_data$infectious_period)
  f = 1/infection_data$latent_period
  r = 1/infection_data$infectious_period
  t = 0

#  print(paste0("beta = ", beta, sep=""))
#  print(paste0("r = ", r, sep=""))
#  print(paste0("f = ", f, sep=""))

  df <- data.frame(
    "t" = c(t),
    "S" = c(S),
    "E" = c(E),
    "I" = c(I),
    "R" = c(R)
  )
  steps <- seq(1, duration, by=time_step)
  for (step in steps) {
    St = S - beta*I*S
    Et = E + beta*I*S - f*E
    It = I + f*E - r*I
    Rt = R + r*I

    df <- rbind(df, c(step, St, Et, It, Rt))

    S = St
    E = Et
    I = It
    R = Rt
  }

  return(df);
}


#' Generate an SEIRS model for an infection based on input parameters
#' @param infection The name of the infection to model c("influenza", "COVID").
#' @param duration The length of time to model, in days (optional; 120).
#' @param time_step The time between each model step, in days (optional; 1).
#' @param N The total population size (optional; 10^6).
#' @param percent_susceptible The initial susceptible population, as a percentage of N (optional; 100).
#' @param number_infectious The initial number of infectious individuals (optional; 1).
#' @param infection_data A list of key-value pairs with basic infection data. See simple_infections.
SEIRS <- function(infection, duration, time_step, N, percent_susceptible, number_infectious, infection_data) {

  model = "SEIRS"

  if (missing(infection)) infection = "influenza"
  if (missing(duration)) duration = 120
  if (missing(time_step)) time_step = 1
  if (missing(N)) N = 10^6
  if (missing(percent_susceptible)) percent_susceptible = 100
  if (missing(number_infectious)) number_infectious = 1

  if (missing(infection_data)) {
    if (infection %in% simple_infections$infection_name) {
      infection_data <- as.list(subset(simple_infections, infection_name==infection))
    } else {
      infection_data <- as.list(subset(simple_infections, infection_name=="influenza"))
    }
  }

  #
  # Need to check user-provided infection_data for required values
  #

  S = (percent_susceptible/100) * N
  E = 0
  I = number_infectious
  R = 0
  beta = infection_data$basic_reproduction_number/(N*infection_data$infectious_period)
  f = 1/infection_data$latent_period
  r = 1/infection_data$infectious_period
  gamma = 1/infection_data$immune_period
  t = 0

#    print(paste0("beta = ", beta, sep=""))
#    print(paste0("r = ", r, sep=""))
#    print(paste0("f = ", f, sep=""))
#    print(paste0("gamma = ", gamma, sep=""))

  df <- data.frame(
    "t" = c(t),
    "S" = c(S),
    "E" = c(E),
    "I" = c(I),
    "R" = c(R)
  )
  steps <- seq(1, duration, by=time_step)
  for (step in steps) {
    St = S - beta*I*S + gamma*R
    Et = E + beta*I*S - f*E
    It = I + f*E - r*I
    Rt = R + r*I - gamma*R

    df <- rbind(df, c(step, St, Et, It, Rt))

    S = St
    E = Et
    I = It
    R = Rt
  }

  return(df);
}


globalVariables(c("simple_infections", names(simple_infections)))

