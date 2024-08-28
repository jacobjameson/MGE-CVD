################################################################################
# AUTHOR:             Jacob Jameson
#
# DESCRIPTION:        This script loads all necessary packages for the project 
#                     and runs all other scripts necessary for the project.
#
# R version:           4.3.3 (2024-8-28)
################################################################################

rm(list = ls()) # Clear workspace

# Start logging ---------------------------------------------------------------
log_file <- "outputs/run_all_log.txt"
sink(log_file) # Start diverting R output to log file

cat("Log started: ", date(), "\n\n")

library(lmtest)
library(sandwich)
library(broom)
library(tidyverse)
library(sjstats)
library(stargazer)
library(marginaleffects)
library(gtsummary)
library(ggeffects)
library(ggsci)
library(haven)
library(labelled)
library(scales)
library(knitr)
library(kableExtra)
library(survey)
library(margins)
library(gtsummary)
library(ggthemes)


# Run all scripts -------------------------------------------------------------
## Generate data
options(warn = -1) # Suppress warnings

source("src/prepare data/Construct Analytical dataset.R",
       echo = TRUE, max.deparse.length = 1000)

## Run models
source("src/gen tables.R", echo = TRUE, max.deparse.length = 1000)

## Generate figures
source("src/gen figures.R", echo = TRUE,max.deparse.length = 1000)

# End logging -----------------------------------------------------------------
cat("\nAll scripts have been executed successfully. 
    Outputs are available in the outputs directory.\n")

cat("Log ended: ", date(), "\n")

# Stop logging
sink(type = "message") # Stop diverting messages to log file
sink() # Stop diverting output to log file
################################################################################