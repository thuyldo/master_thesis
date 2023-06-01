# packages

library(survival)
library(PBIR)
library(purrr)
library(ggplot2)
library(tidyr)
library(latex2exp)
library(dplyr)
library(plyr)
library(haven)
library(epitools)
library(survminer)
library(gridExtra)

# load data

demo <- read_sas("Datasphère/zip/AllProvidedFiles_265/demo.sas7bdat", NULL)
aeendpt <- read_sas("Datasphère/zip/AllProvidedFiles_265/aeendpt.sas7bdat", NULL)
asendpt <- read_sas("Datasphère/zip/AllProvidedFiles_265/asendpt.sas7bdat",  NULL)
medhist <- read_sas("Datasphère/zip/AllProvidedFiles_265/medhist.sas7bdat", NULL)
respeval <- read_sas("Datasphère/zip/AllProvidedFiles_265/respeval.sas7bdat", NULL)
radiotx <- read_sas("Datasphère/zip/AllProvidedFiles_265/radiotx.sas7bdat", NULL)
chemotx <- read_sas("Datasphère/zip/AllProvidedFiles_265/chemotx.sas7bdat", NULL)
surghist <- read_sas("Datasphère/zip/AllProvidedFiles_265/surghist.sas7bdat", NULL)
