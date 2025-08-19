###########################################
# Statistical Analysis Script
# Purpose: Univariate and Multivariate Cox regression with FDR adjustment
#          + Proportional hazards assumption test
###########################################
#-----------------------------------------
# 0. Install required packages (run once)
#-----------------------------------------
# Uncomment the following lines if the packages are not installed:
install.packages("survival")
install.packages("dplyr")
install.packages("readxl")
install.packages("forestmodel")

#-----------------------------------------
# Load required packages
#-----------------------------------------
# Load required packages
library(survival)      # For Cox models
library(forestmodel)   # For forest plots
library(dplyr)         # For data manipulation
library(readxl)        # For reading Excel files

#-----------------------------------------
# 1. Data Import
#-----------------------------------------
# Replace with your actual data file name
my_data <- read_excel("my_data.xlsx")

# Ensure time and status are numeric
my_data$time <- as.numeric(my_data$time)
my_data$status <- as.numeric(my_data$status)

#-----------------------------------------
# 2. Define categorical variables
#-----------------------------------------
factor_vars <- c("Age", "Sex", "Smoking", "Alcohol", "Hypertension",
                 "Diabetes mellitus", "BMI", "T stage", "N stage",  
                 "Tumor location", "Tumor length(cm)", "Induction therapy", 
                 "Lymphopenia", "Clinical response")

# Convert to factor and set reference levels
my_data$`Tumor length(cm)` <- factor(my_data$`Tumor length(cm)`)
my_data$`Tumor length(cm)` <- relevel(my_data$`Tumor length(cm)`, ref = "≤5")

my_data$Age <- factor(my_data$Age)
my_data$Age <- relevel(my_data$Age, ref = "<60")

my_data$BMI <- factor(my_data$BMI)
my_data$BMI <- relevel(my_data$BMI, ref = "18.5-24.9")

my_data$`Tumor location` <- factor(my_data$`Tumor location`)
my_data$`Tumor location` <- relevel(my_data$`Tumor location`, ref = "Upper")

# Convert other variables to factor if not already
for (var in factor_vars) {
  if (!is.factor(my_data[[var]])) {
    my_data[[var]] <- factor(my_data[[var]])
  }
}

# Function to safely handle variable names with spaces/special characters
safe_var <- function(x) {
  if (grepl("[^a-zA-Z0-9_]", x)) paste0("`", x, "`") else x
}

#-----------------------------------------
# 3. Univariate Cox regression
#-----------------------------------------
univ_results <- lapply(factor_vars, function(var) {
  var_safe <- safe_var(var)
  f <- as.formula(paste("Surv(time, status) ~", var_safe))
  model <- coxph(f, data = my_data)
  s <- summary(model)
  data.frame(
    Variable = var,
    HR = round(s$coefficients[,"exp(coef)"], 2),
    Lower_CI = round(s$conf.int[,"lower .95"], 2),
    Upper_CI = round(s$conf.int[,"upper .95"], 2),
    P_value = round(s$coefficients[,"Pr(>|z|)"], 4),
    row.names = NULL
  )
})

# Combine results
univ_df <- do.call(rbind, univ_results)

# Adjust p-values with Benjamini-Hochberg (FDR)
univ_df$Adj_P_value <- p.adjust(univ_df$P_value, method = "BH")

# Print univariate results
print(univ_df)

#-----------------------------------------
# 4. Multivariate Cox regression
#-----------------------------------------
# Select variables with adjusted p < 0.2 (threshold can be changed)
selected_vars <- univ_df %>% filter(Adj_P_value < 0.2) %>% pull(Variable)
selected_vars_safe <- sapply(selected_vars, safe_var)

if (length(selected_vars_safe) > 0) {
  multiv_formula <- as.formula(paste("Surv(time, status) ~", paste(selected_vars_safe, collapse = " + ")))
  multiv_model <- coxph(multiv_formula, data = my_data)
  
  # Plot forest model
  forest_model(multiv_model)
  
  # Test proportional hazards assumption
  ph_test <- cox.zph(multiv_model)
  print(ph_test)
  
  # Plot Schoenfeld residuals
  plot(ph_test)
  
} else {
  cat("No variables with adjusted p < 0.2. Multivariate model not performed.\n")
}