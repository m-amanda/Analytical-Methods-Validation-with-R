# Load library
library(writexl) # Used to export data
library(reshape2) # Used to melt data.frames

set.seed(123) # Set seed for reproducibility

# Calibrations
# Set calibration values
concentrations <- rep(seq(3, 15, by = 3), times = 3)

# Create a data frame with concentrations and SA and ES responses
calibrations <- data.frame(
  concentrations = concentrations,
  SA_responses = 2.03 * concentrations + 0.96 + rnorm(15, 0, 0.06) + 11, # Regression function to define responses with slope, intercept, error, and a constant initial value simulating an increase in response due to sample analyte.
  ES_responses = 3.05 * concentrations + 1.96 + rnorm(15, 0, 0.03) # Regression function to define responses with slope, intercept, and error.
)

YC_volume <- seq(100, 500, by = 100) # Set Youden calibration volumes

# Create a data frame for Youden calibration
youden <- data.frame(
  YC_volume = YC_volume,
  YC_responses = 1.60 * YC_volume + 2.05 # Regression with slope and intercept to define Youden responses
)

# ES accuracy
# Generate control points
control_points <- rep(c(0, 3, 9, 15), each = 3)

# Function to calculate responses with varied deviations and a constant to simule sample analyte response
generate_responses <- function(control_points, sd) {
  3.05 * control_points + 1.96 + rnorm(length(control_points), 0, sd) + 11
}

# Define standard deviations for each day
sds <- c(0.07, 0.08, 0.1, 0.09, 0.05)

# Generate responses for days 1 to 5
days_responses <- sapply(sds, function(sd) generate_responses(control_points, sd))

ES_accuracy <- data.frame(control_points, days_responses)
colnames(ES_accuracy) <- c("control_points", "1", "2", "3", "4", "5")

# Melt the data frame
ES_accuracy <- melt(ES_accuracy, id = "control_points")
colnames(ES_accuracy) <- c("concentrations", "days", "response")

# SA accuracy
# Control points for concentration (sequence from 0 to 15 with increment of 3)
concentrations <- rep(seq(0, 15, by = 3), times = 3)

# Function to generate responses based on concentration and deviation, and response increment from fortification
generate_responses <- function(concentrations, sd, Qr) {
  2.03 * concentrations + 0.96 + rnorm(length(concentrations), 0, sd) + (Qr + 11)
}

# Define standard deviations for each day
sds <- c(0.07, 0.05, 0.09, 0.03, 0.05)

# Generate responses for the 5 days
days_responses_non_fortified <- sapply(sds, function(sd) generate_responses(concentrations, sd, 0))
days_responses_Q1_fortified <- sapply(sds, function(sd) generate_responses(concentrations, sd, 7.05))
days_responses_Q2_fortified <- sapply(sds, function(sd) generate_responses(concentrations, sd, 19.23))
days_responses_Q3_fortified <- sapply(sds, function(sd) generate_responses(concentrations, sd, 31.41))

# Function to process data
process_data <- function(concentrations, responses) {
  df <- data.frame(concentrations, responses)
  colnames(df) <- c("control_points", as.character(1:5))
  melted_df <- melt(df, id = "control_points")
  colnames(melted_df) <- c("concentrations", "days", "response")
  return(melted_df)
}

# Apply the function to each data set
SA_accuracy_non_fortified <- process_data(concentrations, days_responses_non_fortified)
SA_accuracy_Q1_fortified <- process_data(concentrations, days_responses_Q1_fortified)
SA_accuracy_Q2_fortified <- process_data(concentrations, days_responses_Q2_fortified)
SA_accuracy_Q3_fortified <- process_data(concentrations, days_responses_Q3_fortified)

# Define path to export data
setwd("~/Dropbox/AMV-app")

# Export data to .xlsx documents
write_xlsx(calibrations, path = "calibrations.xlsx")
write_xlsx(youden, path = "youden.xlsx")
write_xlsx(ES_accuracy, path = "ES_accuracy.xlsx")
write_xlsx(SA_accuracy_non_fortified, path = "SA_accuracy_non_fortified.xlsx")
write_xlsx(SA_accuracy_Q1_fortified, path = "SA_accuracy_Q1_fortified.xlsx")
write_xlsx(SA_accuracy_Q2_fortified, path = "SA_accuracy_Q2_fortified.xlsx")
write_xlsx(SA_accuracy_Q3_fortified, path = "SA_accuracy_Q3_fortified.xlsx")
