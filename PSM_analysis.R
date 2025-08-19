## ---- 0) Install & Load packages ----
pkgs <- c("readxl", "dplyr", "MatchIt", "tableone")
to_install <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(to_install) > 0) install.packages(to_install)

library(readxl)
library(dplyr)
library(MatchIt)
library(tableone)

## ---- 1) Read data ----
my_data <- read_excel("D:/my_data.xlsx")
df <- as.data.frame(my_data)
print(names(df))
print(head(df))

## ---- 2) Define treatment variable ----
df <- df %>%
  mutate(
    treatment = case_when(
      treatment %in% c(1, "1", "treated", "Treatment") ~ 1,
      TRUE ~ 0
    ),
    treatment = factor(treatment, levels = c(0, 1),
                       labels = c("Control", "Treatment"))
  )

## ---- 3) Define covariates ----
catVars <- c("Age", "Sex", "Smoking", "Alcohol", "Hypertension",
             "Diabetes_mellitus", "BMI", "Tumor_length",
             "Tumor_location", "TNM")

vars_in_data <- intersect(catVars, names(df))
validCatVars <- vars_in_data[sapply(df[vars_in_data], function(x) {
  length(unique(x[!is.na(x)])) > 1
})]
df[validCatVars] <- lapply(df[validCatVars], as.factor)

## ---- 4) Pre-matching balance ----
cat("\n===== Pre-matching sample sizes =====\n")
print(table(df$treatment))

tab1_pre <- CreateTableOne(
  vars       = vars_in_data,
  strata     = "treatment",
  data       = df,
  factorVars = intersect(validCatVars, vars_in_data),
  test       = FALSE
)
cat("\n===== Table 1 (Pre-matching) with SMD =====\n")
print(tab1_pre, showAllLevels = TRUE, quote = TRUE, noSpaces = TRUE, smd = TRUE)

tab1_pre_df <- print(tab1_pre, showAllLevels = TRUE, quote = TRUE, noSpaces = TRUE,
                     smd = TRUE, printToggle = FALSE)
write.csv(as.data.frame(tab1_pre_df), "table1_pre_matching_with_SMD.csv", row.names = FALSE)

## ---- 5) Propensity Score Matching ----
m.formula <- as.formula(paste("treatment ~", paste(vars_in_data, collapse = " + ")))

m.out <- matchit(
  formula = m.formula,
  data    = df,
  method  = "nearest",
  caliper = 0.1,
  ratio   = 1
)

m.data <- match.data(m.out)

## ---- 6) Post-matching balance ----
cat("\n===== Post-matching sample sizes =====\n")
print(table(m.data$treatment))

tab1_post <- CreateTableOne(
  vars       = vars_in_data,
  strata     = "treatment",
  data       = m.data,
  factorVars = intersect(validCatVars, vars_in_data),
  test       = FALSE
)
cat("\n===== Table 1 (Post-matching) with SMD =====\n")
print(tab1_post, showAllLevels = TRUE, quote = TRUE, noSpaces = TRUE, smd = TRUE)

tab1_post_df <- print(tab1_post, showAllLevels = TRUE, quote = TRUE, noSpaces = TRUE,
                      smd = TRUE, printToggle = FALSE)
write.csv(as.data.frame(tab1_post_df), "table1_post_matching_with_SMD.csv", row.names = FALSE)

## ---- 7) Save matched dataset ----
write.csv(m.data, "match.csv", row.names = FALSE)

cat("\nOutputs saved:\n",
    "- table1_pre_matching_with_SMD.csv\n",
    "- table1_post_matching_with_SMD.csv\n",
    "- match.csv (matched dataset)\n")