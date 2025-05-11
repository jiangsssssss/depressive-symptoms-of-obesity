# Packages
library(foreign)
library(tidyverse)
library(tidyselect)
library(dplyr)
library(broom)
library(devtools)
library(dietaryindex)
library(caret)
library(VIM)
library(Boruta)
library(glmnet)
library(smotefamily)
library(pROC)
library(PRROC)
library(rms)
library(pmsampsize)
library(randomForest)
library(purrr)
library(boot)
library(DescTools)
library(dcurves)
library(e1071)    
library(ggplot2)  
library(rmda) 
library(xgboost)
library(ggplot2)
library(patchwork)
library(shapviz)
library(fastshap)

# Participants
## Import data(exclude dietary data)
DATA = read.csv("D:/Clinical prediction model/01data/clean/DATA_without_DII.csv")
## Import 15-16 dietary data (dietaryindex comes with 17-18 dietary data)
setwd("D:/Clinical prediction model/01data/raw/I/")
load("NHANES_20152016.RDA")
NHANES_20152016 <- NHANES_20152016 
## Calculate Dietary Inflammation Index
data("NHANES_20152016")
DII_I <- DII_NHANES_FPED(FPED_PATH = NHANES_20152016$FPED, NUTRIENT_PATH = NHANES_20152016$NUTRIENT, DEMO_PATH = NHANES_20152016$DEMO, FPED_PATH2 = NHANES_20152016$FPED2, NUTRIENT_PATH2 = NHANES_20152016$NUTRIENT2)
data("NHANES_20172018")
DII_J <- DII_NHANES_FPED(FPED_PATH = NHANES_20172018$FPED, NUTRIENT_PATH = NHANES_20172018$NUTRIENT, DEMO_PATH = NHANES_20172018$DEMO, FPED_PATH2 = NHANES_20172018$FPED2, NUTRIENT_PATH2 = NHANES_20172018$NUTRIENT2)
## Merge data and include samples with BMI >30
DII_Total <- bind_rows(DII_I, DII_J)
DATA <- DATA %>%
  left_join(DII_Total, by = "SEQN")
DATA_OBS <- subset(DATA, BMXBMI > 30)

# Predictors and Outcome
## Alcohol_intake
DATA_OBS$`Alcohol_intake` <- NA
DATA_OBS$`Alcohol_intake`[DATA_OBS$ALQ111 == 2] <- "No"
DATA_OBS$`Alcohol_intake`[DATA_OBS$ALQ111 == 1 & DATA_OBS$ALQ121 == 0] <- "Past"
DATA_OBS$`Alcohol_intake`[DATA_OBS$ALQ121 > 0 & DATA_OBS$ALQ121 <= 10 & DATA_OBS$ALQ142 == 0 &
                            ((DATA_OBS$RIAGENDR == 1 & DATA_OBS$ALQ130 <= 2) | (DATA_OBS$RIAGENDR == 2 & DATA_OBS$ALQ130 <= 1))] <- "Normal"
DATA_OBS$`Alcohol_intake`[DATA_OBS$ALQ142 > 0 & DATA_OBS$ALQ142 <= 10] <- "Over"
DATA_OBS$`Alcohol_intake`[DATA_OBS$ALQ110 == 2] <- "No"
DATA_OBS$`Alcohol_intake`[(DATA_OBS$ALQ110 == 1 | DATA_OBS$ALQ101 == 1) & DATA_OBS$ALQ120Q == 0] <- "Past"
DATA_OBS$`Alcohol_intake`[DATA_OBS$ALQ120Q > 0 & DATA_OBS$ALQ120Q <= 365 & DATA_OBS$ALQ141Q == 0 &
                            ((DATA_OBS$RIAGENDR == 1 & DATA_OBS$ALQ130 <= 2) | (DATA_OBS$RIAGENDR == 2 & DATA_OBS$ALQ130 <= 1))] <- "Normal"
DATA_OBS$`Alcohol_intake`[(DATA_OBS$RIAGENDR == 1 & DATA_OBS$ALQ130 > 2) |
                            (DATA_OBS$RIAGENDR == 2 & DATA_OBS$ALQ130 > 1)] <- "Over"
DATA_OBS$`Alcohol_intake`[DATA_OBS$ALQ141Q > 0 & DATA_OBS$ALQ141Q <= 365] <- "Over"
## Smoking
DATA_OBS <- DATA_OBS %>%
  mutate(Smoking = case_when(SMQ020 == 2 ~ "Never", SMQ020 == 1 & SMQ040 == 3 ~ "Past", SMQ040 %in% c(1, 2) ~ "Now", TRUE ~ NA_character_))
## Moderate recreational activities
DATA_OBS <- DATA_OBS %>%
  mutate(MRA = case_when(PAQ665 == 1 ~ "Yes", PAQ665 == 2 ~ "No", PAQ665 %in% c(7, 9) | is.na(PAQ665) ~ NA_character_, TRUE ~ NA_character_))
## Sedentary behavior
DATA_OBS <- DATA_OBS %>%
  mutate(`SB` = case_when(PAD680 != 7777 & PAD680 != 9999 & PAD680 >= 480 ~ "Yes", PAD680 != 7777 & PAD680 != 9999 & PAD680 < 480 ~ "No", TRUE ~ NA_character_))
## Health_condition
DATA_OBS <- DATA_OBS %>%
  mutate(`Health_condition` = case_when(HSD010 %in% c(1, 2, 3) ~ "Yes", HSD010 %in% c(4, 5) ~ "No", TRUE ~ NA_character_))
## Trouble_sleeping
DATA_OBS <- DATA_OBS %>%
  mutate(Trouble_sleeping = case_when(SLQ050 == 1 ~ "Yes", SLQ050 == 2 ~ "No", SLQ050 %in% c(7, 9) | is.na(SLQ050) ~ NA_character_, TRUE ~ NA_character_))
## Seen_mental_health_professional
DATA_OBS <- DATA_OBS %>%
  mutate(Seen_mental_health_professional = case_when(HUQ090 == 1 ~ "Yes", HUQ090 == 2 ~ "No", HUQ090 %in% c(7, 9) | is.na(SLQ050) ~ NA_character_, TRUE ~ NA_character_))
## Hypertension
DATA_OBS <- DATA_OBS %>%
  mutate(Sys_Avg = (BPXSY1 + BPXSY2 + BPXSY3) / 3, Dia_Avg = (BPXDI1 + BPXDI2 + BPXDI3) / 3, Hypertension = if_else(BPQ020 == 1 | Sys_Avg > 140 | Dia_Avg > 90, "Yes", "No"))
## Diabetes
DATA_OBS$Diabetes <- apply(DATA_OBS, 1, function(row) {
  if (!is.na(row["LBXGLU"]) && row["LBXGLU"] >= 126) {
    return("Yes")
  } else if (!is.na(row["LBXGH"]) && row["LBXGH"] >= 6.5) {
    return("Yes")
  } else if (!is.na(row["LBXGLT"]) && row["LBXGLT"] >= 198) {
    return("Yes")
  } else if (!is.na(row["DIQ010"]) && row["DIQ010"] == 1) {
    return("Yes")
  } else if (!is.na(row["DIQ050"]) && row["DIQ050"] == 1) {
    return("Yes")
  } else if (!is.na(row["DIQ070"]) && row["DIQ070"] == 1) {
    return("Yes")
  } else if (is.na(row["LBXGLU"]) && is.na(row["LBXGH"]) && is.na(row["LBXGLT"]) &&
             (is.na(row["DIQ010"]) || row["DIQ010"] %in% c(7, 9)) &&
             (is.na(row["DIQ050"]) || row["DIQ050"] %in% c(7, 9)) &&
             (is.na(row["DIQ070"]) || row["DIQ070"] %in% c(7, 9))) {
    return(NA)
  } else {
    return("No")
  }
})
## Coronary heart disease
DATA_OBS <- DATA_OBS %>%
  mutate(CHD = if_else(MCQ160C %in% c(7, 9), NA_real_, MCQ160C))
DATA_OBS <- DATA_OBS %>%
  mutate(CHD = case_when(CHD == 2 ~ "No", CHD == 1 ~ "Yes", TRUE ~ as.character(CHD)))
## Congestive heart failure
DATA_OBS <- DATA_OBS %>%
  mutate(CHF = if_else(MCQ160B %in% c(7, 9), NA_real_, MCQ160B))
DATA_OBS <- DATA_OBS %>%
  mutate(CHF = case_when(CHF == 2 ~ "No", CHF == 1 ~ "Yes", TRUE ~ as.character(CHD)))
## Angina
DATA_OBS <- DATA_OBS %>%
  mutate(Angina = if_else(MCQ160D %in% c(7, 9), NA_real_, MCQ160D))
DATA_OBS <- DATA_OBS %>%
  mutate(Angina = case_when(Angina == 2 ~ "No", Angina == 1 ~ "Yes", TRUE ~ as.character(Angina)))
## Heart attack
DATA_OBS <- DATA_OBS %>%
  mutate(Heartattack = if_else(MCQ160E %in% c(7, 9), NA_real_, MCQ160E))
DATA_OBS <- DATA_OBS %>%
  mutate(Heartattack = case_when(Heartattack == 2 ~ "No", Heartattack == 1 ~ "Yes", TRUE ~ as.character(Heartattack)))
## Stroke
DATA_OBS <- DATA_OBS %>%
  mutate(Stroke = if_else(MCQ160F %in% c(7, 9), NA_real_, MCQ160F))
DATA_OBS <- DATA_OBS %>%
  mutate(Stroke = case_when(Stroke == 2 ~ "No", Stroke == 1 ~ "Yes", TRUE ~ as.character(Stroke)))
## Asthma
DATA_OBS <- DATA_OBS %>%
  mutate(Asthma = if_else(MCQ010 %in% c(7, 9), NA_real_, MCQ010))
DATA_OBS <- DATA_OBS %>%
  mutate(Asthma = case_when(Asthma == 2 ~ "No", Asthma == 1 ~ "Yes", TRUE ~ as.character(Asthma)))
# Arthritis
DATA_OBS <- DATA_OBS %>%
  mutate(Arthritis = if_else(MCQ160A %in% c(7, 9), NA_real_, MCQ160A))
DATA_OBS <- DATA_OBS %>%
  mutate(Arthritis = case_when(Arthritis == 2 ~ "No", Arthritis == 1 ~ "Yes", TRUE ~ as.character(Arthritis)))
## Obstructive sleep apnea
DATA_OBS <- DATA_OBS %>%
  mutate(OSA = case_when(
    SLQ030 %in% c(2, 3) | SLQ040 %in% c(2, 3) | SLQ120 == 4 ~ "Yes",
    (SLQ120 %in% c(7, 9) | is.na(SLQ120)) & (SLQ040 %in% c(7, 9) | is.na(SLQ040)) & (SLQ030 %in% c(7, 9) | is.na(SLQ030)) ~ NA_character_, TRUE ~ "No"
  ))
## Gender
DATA_OBS <- DATA_OBS %>%
  mutate(Gender = case_when(RIAGENDR == 1 ~ "Male", RIAGENDR == 2 ~ "Female", TRUE ~ NA_character_))
## Race
DATA_OBS <- DATA_OBS %>%
  mutate(Race = case_when(
    RIDRETH1 == 1 ~ "Mexican_American", RIDRETH1 == 2 ~ "Other_Hispanic", RIDRETH1 == 3 ~ "Non-Hispanic_White", RIDRETH1 == 4 ~ "Non-Hispanic_Black", RIDRETH1 == 5 ~ "Other_Race", TRUE ~ NA_character_
  ))
## Education
DATA_OBS$Education <- NA
DATA_OBS$Education[DATA_OBS$DMDEDUC2 %in% c(1, 2)] <- "less than high school"
DATA_OBS$Education[DATA_OBS$DMDEDUC2 == 3] <- "completed high school"
DATA_OBS$Education[DATA_OBS$DMDEDUC2 %in% c(4, 5)] <- "more than high school"
## Marital status
DATA_OBS$Marital_status <- NA
DATA_OBS$Marital_status[DATA_OBS$DMDMARTL %in% c(1, 6)] <- "Married/Living_with_partner"
DATA_OBS$Marital_status[DATA_OBS$DMDMARTL == 5] <- "Never"
DATA_OBS$Marital_status[DATA_OBS$DMDMARTL %in% c(2, 3, 4)] <- "Widowed/Divorced/Separated"
## Poverty-to-income ratio
DATA_OBS$PIR <- NA
DATA_OBS$PIR[DATA_OBS$INDFMPIR < 1.3] <- "<1.3"
DATA_OBS$PIR[DATA_OBS$INDFMPIR >= 1.3 & DATA_OBS$INDFMPIR <= 3.5] <- "1.3-3.5"
DATA_OBS$PIR[DATA_OBS$INDFMPIR > 3.5] <- ">3.5"
## Health_insurance
DATA_OBS <- DATA_OBS %>%
  mutate(Health_insurance = if_else(HIQ011 %in% c(7, 9), NA_real_, HIQ011))
DATA_OBS <- DATA_OBS %>%
  mutate(Health_insurance = case_when(
    Health_insurance == 2 ~ "No",
    Health_insurance == 1 ~ "Yes",
    TRUE ~ as.character(Health_insurance)
  ))
names(DATA_OBS)[names(DATA_OBS) == "RIDAGEYR"] <- "Age"
names(DATA_OBS)[names(DATA_OBS) == "BMXBMI"] <- "BMI"
colnames(DATA_OBS)[colnames(DATA_OBS) == "LBXTC"] <- "TC"
names(DATA_OBS)[names(DATA_OBS) == "LBDHDD"] <- "HDL"
names(DATA_OBS)[names(DATA_OBS) == "LBXHSCRP"] <- "HS_CRP"
names(DATA_OBS)[names(DATA_OBS) == "DII_ALL"] <- "DII"
## Depression
sum_vars <- rowSums(DATA_OBS[, c("DPQ010", "DPQ020", "DPQ030", "DPQ040", "DPQ050", "DPQ060", "DPQ070", "DPQ080", "DPQ090")], na.rm = TRUE)
has_7_or_9 <- apply(DATA_OBS[, c("DPQ010", "DPQ020", "DPQ030", "DPQ040", "DPQ050", "DPQ060", "DPQ070", "DPQ080", "DPQ090")], 1, function(x) any(x %in% c(7, 9)))
has_na <- apply(DATA_OBS[, c("DPQ010", "DPQ020", "DPQ030", "DPQ040", "DPQ050", "DPQ060", "DPQ070", "DPQ080", "DPQ090")], 1, function(x) any(is.na(x)))
DATA_OBS$Depression <- ifelse(has_na | has_7_or_9, NA, ifelse(sum_vars >= 10, "Yes", "No"))

# Exclusion
Data_OBS <- DATA_OBS[, c(
  "SEQN", "Gender", "Age", "BMI", "Race", "DII", "Alcohol_intake",
  "Smoking", "Hypertension", "Diabetes", "CHD", "CHF", "Angina",
  "Heartattack", "Stroke", "Asthma", "Arthritis", "OSA","Education",
  "Marital_status", "PIR", "Health_insurance", "MRA", "SB", "TC", "HDL","Trouble_sleeping","Health_condition","Seen_mental_health_professional","Depression"
)]
Data_OBS <- Data_OBS[!is.na(Data_OBS$Depression), ]
Data_OBS <- Data_OBS[Data_OBS$Age >= 18, ]
Data_OBS <- Data_OBS %>% dplyr::select(-SEQN)
missing_values <- sapply(Data_OBS, function(x) sum(is.na(x)))
missing_values_df <- data.frame(Variable = names(missing_values), MissingValues = missing_values)
print(missing_values_df)
# Interpolation of missing data
factor_cols <- c( "Gender","Race", "Alcohol_intake",
                  "Smoking", "Hypertension", "Diabetes", "CHD", "CHF", "Angina",
                  "Heartattack", "Stroke", "Asthma", "Arthritis", "OSA","Education",
                  "Marital_status", "PIR", "Health_insurance", "MRA", "SB", "Trouble_sleeping","Health_condition",
                  "Seen_mental_health_professional","Depression")
Data_OBS[factor_cols] <- lapply(Data_OBS[factor_cols], as.factor)
numeric_cols <- c("Age", "BMI", "DII", "TC", "HDL")
vars_to_impute <- names(Data_OBS)[names(Data_OBS) != "Depression"]
data_OBS <- kNN(Data_OBS, variable = vars_to_impute, k = 5)
summary(data_OBS) 
data_OBS <- data_OBS[, -c(30:57)]
str(data_OBS)

# Sample description
data_OBS <- data_OBS %>%
  mutate(across(c(
    Gender, Race, Alcohol_intake, Hypertension, Diabetes, CHD,
    CHF, Angina, Heartattack, Stroke, OSA, Asthma, Trouble_sleeping, Education,
    Marital_status, PIR, Health_insurance, MRA, SB, Smoking, Arthritis,
    Health_condition, Seen_mental_health_professional, Depression
  ), as.factor))
categorical_vars <- c("Gender", "Race", "Alcohol_intake", "Smoking", "Hypertension", "Diabetes", "CHD", "OSA","CHF", "Angina", "Heartattack", "Stroke", "Asthma", "Arthritis","Health_condition", "Trouble_sleeping","Education", "Marital_status", "PIR", "Health_insurance", "Seen_mental_health_professional","MRA", "SB")
continuous_vars <- c("Age", "BMI", "DII", "TC", "HDL")
calc_mean_ci <- function(data, variable) {
  est <- mean(data, na.rm = TRUE)
  se <- sd(data, na.rm = TRUE) / sqrt(length(data))
  lower <- est - 1.96 * se
  upper <- est + 1.96 * se
  return(c(round(est, 3), round(lower, 3), round(upper, 3)))
}
calc_percent_and_ci <- function(x, n) {
  percent <- round(x / n * 100, 3)
  lower <- round(percent - 1.96 * sqrt((x / n) * (1 - x / n) / n) * 100, 3)
  upper <- round(percent + 1.96 * sqrt((x / n) * (1 - x / n) / n) * 100, 3)
  return(c(percent, lower, upper))
}
Sample_description_results <- data.frame(
  Variable = character(),
  No_Percent = numeric(), No_CI_Lower = numeric(), No_CI_Upper = numeric(),
  Yes_Percent = numeric(), Yes_CI_Lower = numeric(), Yes_CI_Upper = numeric(),
  P_Value = numeric(), stringsAsFactors = FALSE
)
for (var in continuous_vars) {
  no_data <- data_OBS[data_OBS$Depression == "No", var]
  yes_data <- data_OBS[data_OBS$Depression == "Yes", var]
  no_mean_ci <- calc_mean_ci(no_data, var)
  yes_mean_ci <- calc_mean_ci(yes_data, var)
  t_test <- t.test(no_data, yes_data, var.equal = TRUE)
  Sample_description_results <- rbind(Sample_description_results, data.frame(
    Variable = var,
    No_Percent = no_mean_ci[1], No_CI_Lower = no_mean_ci[2], No_CI_Upper = no_mean_ci[3],
    Yes_Percent = yes_mean_ci[1], Yes_CI_Lower = yes_mean_ci[2], Yes_CI_Upper = yes_mean_ci[3],
    P_Value = round(t_test$p.value, 4)
  ))
}
for (var in categorical_vars) {
  table_result <- table(data_OBS$Depression, data_OBS[, var])
  chi_result <- chisq.test(table_result)
  levels <- levels(data_OBS[, var])
  for (level in levels) {
    no_count <- sum(data_OBS[data_OBS$Depression == "No", var] == level, na.rm = TRUE)
    yes_count <- sum(data_OBS[data_OBS$Depression == "Yes", var] == level, na.rm = TRUE)
    no_total <- sum(data_OBS$Depression == "No", na.rm = TRUE)
    yes_total <- sum(data_OBS$Depression == "Yes", na.rm = TRUE)
    no_percent_ci <- calc_percent_and_ci(no_count, no_total)
    yes_percent_ci <- calc_percent_and_ci(yes_count, yes_total)
    Sample_description_results <- rbind(Sample_description_results, data.frame(
      Variable = paste(var, level, sep = ": "),
      No_Percent = no_percent_ci[1], No_CI_Lower = no_percent_ci[2], No_CI_Upper = no_percent_ci[3],
      Yes_Percent = yes_percent_ci[1], Yes_CI_Lower = yes_percent_ci[2], Yes_CI_Upper = yes_percent_ci[3],
      P_Value = round(chi_result$p.value, 4)
    ))
  }
}
print(Sample_description_results)


#Feature selection
##Boruta
data_OBS_log <- data_OBS[, c("Gender", "Alcohol_intake", "Smoking", "Diabetes", "Hypertension","CHF", "CHD",
                             "Angina", "Stroke", "Asthma", "Arthritis", 
                             "Education", "Marital_status", "PIR", "MRA", "OSA", "Health_condition","Trouble_sleeping","Seen_mental_health_professional",
                             "BMI","DII", "Depression")]
data_OBS_log <- data_OBS_log %>%
  mutate(across(c(
    Gender, Alcohol_intake, Hypertension, Diabetes, CHD,
    CHF, Angina, Stroke, OSA, Asthma, Trouble_sleeping, Education,
    Marital_status, PIR, MRA, Smoking, Arthritis,
    Health_condition, Seen_mental_health_professional, Depression
  ), as.factor))
str(data_OBS_log)
set.seed(123)
trains <- createDataPartition(y = data_OBS_log$Depression, p = 0.7, list = F)
data_OBS_log_train <- data_OBS_log[trains, ]
data_OBS_log_test <- data_OBS_log[-trains, ]
table(data_OBS_log_train$Depression)
table(data_OBS_log_test$Depression)
set.seed(123)
feature.selection <- Boruta(Depression ~ ., data = data_OBS_log_train, doTrace = 1,maxRuns = 100)
table(feature.selection$finalDecision)
par(mar = c(7.5, 4, 0.5, 2) + 0.1)
plot(feature.selection, cex.axis = 0.50, las = 3, xlab = "")
Feature2 <- getSelectedAttributes(feature.selection)
Feature2


#Modeling
## Data preparation
data_simplify <- data_OBS_log[, c("DII","Health_condition",
                                  "Seen_mental_health_professional", "Trouble_sleeping",
                                  "Marital_status", "Gender","Alcohol_intake","Smoking",
                                  "Depression")]
vars <- c("DII","Depression")
dummy_matrix <- model.matrix(~ -1 + Alcohol_intake + Gender + Smoking + 
                               Marital_status + Health_condition + 
                               Trouble_sleeping + Seen_mental_health_professional, 
                             data = data_simplify)
data_matrix <- cbind(dummy_matrix, data_simplify[, vars])
data_simplify <- as.data.frame(data_matrix)
data_simplify$Depression <- ifelse(data_simplify$Depression == "Yes", 1, 0)
names(data_simplify)[names(data_simplify) == "Marital_statusNever"] <- "Marital_statusNever_married"
names(data_simplify)[names(data_simplify) == "Marital_statusWidowed/Divorced/Separated"] <- "Marital_status_Widowed_Divorced_Separated"
data_simplify <- data_simplify %>% mutate(across(c(Depression), as.factor))
data_simplify$Alcohol_intakeNo <- NULL
set.seed(123)
trains <- createDataPartition(y = data_simplify$Depression, p = 0.7, list = F)
data_simplify_train <- data_simplify[trains, ]
data_simplify_test <- data_simplify[-trains, ]
table(data_simplify_train$Depression)
table(data_simplify_test$Depression)
str(data_simplify_train)
##Sample size calculation
pmsampsize(type = "b",csrsquared = 0.072,parameters = 12,prevalence = 0.1,seed = 123)

##Random Forests
###Imbalance data
set.seed(123)
param_grid <- expand.grid(mtry = 3:5,ntree = seq(200, 500, by = 50))
folds <- createFolds(data_simplify_train$Depression,k = 10,list = TRUE,returnTrain = FALSE)
results <- data.frame()
for (i in 1:nrow(param_grid)) {
  current_mtry <- param_grid$mtry[i]
  current_ntree <- param_grid$ntree[i]
  
  fold_aucs <- numeric(10)  
  
  for (fold in 1:10) {
    val_indices <- folds[[fold]]
    train_data <- data_simplify_train[-val_indices, ]
    val_data <- data_simplify_train[val_indices, ]
    
    
    model <- randomForest(
      Depression ~ .,
      data = train_data,
      mtry = current_mtry,
      ntree = current_ntree
    )
    
    
    pred_probs <- predict(model, val_data, type = "prob")[, 2]  
    
    
    roc_obj <- roc(
      response = val_data$Depression,
      predictor = pred_probs,
      levels = levels(val_data$Depression)  
    )
    fold_aucs[fold] <- auc(roc_obj)
  }
  
  
  mean_auc <- mean(fold_aucs)
  
  
  results <- rbind(results, data.frame(
    mtry = current_mtry,
    ntree = current_ntree,
    auc = mean_auc
  ))
  
  cat(sprintf("mtry=%d, ntree=%d | Mean AUC=%.4f\n", current_mtry, current_ntree, mean_auc))
}
results_sorted <- results[order(-results$auc), ]
best_params_RF <- results_sorted[1, ]
print(best_params_RF)
set.seed(123)
rf.OBS.2 <- randomForest(Depression ~ .,data = data_simplify_train,mtry = best_params_RF$mtry,ntree = best_params_RF$ntree)
print(rf.OBS.2)
rf.OBS.trainprob <- predict(rf.OBS.2, newdata = data_simplify_train, type = "prob")[, 2]
rf_train_roc <- roc(response = data_simplify_train$Depression, predictor = rf.OBS.trainprob)
rf.OBS.testprob <- predict(rf.OBS.2, newdata = data_simplify_test, type = "prob")[, 2]
rf_test_roc <- roc(response = data_simplify_test$Depression, predictor = rf.OBS.testprob)
###Balance data
set.seed(123)
param_grid <- expand.grid(mtry = 3:5, ntree = seq(200, 500, by = 50))
folds <- createFolds(data_simplify_train$Depression, k = 10, list = TRUE, returnTrain = FALSE)
results <- data.frame()
for (i in 1:nrow(param_grid)) {
  current_mtry <- param_grid$mtry[i]
  current_ntree <- param_grid$ntree[i]
  fold_aucs <- numeric(10)
  for (fold in 1:10) {
    val_indices <- folds[[fold]]
    raw_train <- data_simplify_train[-val_indices, ]
    val_data <- data_simplify_train[val_indices, ]
    smote_data <- SMOTE(
      X = raw_train[, -which(names(raw_train) == "Depression")],
      target = raw_train$Depression,
      K = 5,
      dup_size = 3
    )
    balanced_train <- smote_data$data
    colnames(balanced_train)[ncol(balanced_train)] <- "Depression"  
    balanced_train <- balanced_train %>% mutate(across(c(Depression), as.factor))
    model <- randomForest(
      Depression ~ .,
      data = balanced_train,
      mtry = current_mtry,
      ntree = current_ntree
    )
    pred_probs <- predict(model, val_data, type = "prob")[, 2]
    roc_obj <- roc(
      response = val_data$Depression,
      predictor = pred_probs,
      levels = levels(val_data$Depression)
    )
    fold_aucs[fold] <- auc(roc_obj)
  }
  mean_auc <- mean(fold_aucs)
  
  results <- rbind(results, data.frame(
    mtry = current_mtry,
    ntree = current_ntree,
    auc = mean_auc
  ))
  cat(sprintf("mtry=%d, ntree=%d | Mean AUC=%.4f\n", current_mtry, current_ntree, mean_auc))}
results_sorted <- results[order(-results$auc), ]
best_params_RF_balanced <- results_sorted[1, ]
print(best_params_RF_balanced)
final_smote <- SMOTE(X = data_simplify_train[, -which(names(data_simplify_train) == "Depression")],
  target = data_simplify_train$Depression,K = 5,dup_size = 3)
balanced_full_train <- final_smote$data
colnames(balanced_full_train)[ncol(balanced_full_train)] <- "Depression"
balanced_full_train <- balanced_full_train %>% mutate(across(c(Depression), as.factor))
set.seed(123)
rf.OBS.2_balanced <- randomForest(
  Depression ~ .,
  data = balanced_full_train,
  mtry = best_params_RF_balanced$mtry,
  ntree = best_params_RF_balanced$ntree)
rf.OBS_balanced_trainprob <- predict(rf.OBS.2_balanced, newdata = balanced_full_train, type = "prob")[, 2]
rf_train_balanced_roc <- roc(response = balanced_full_train$Depression, predictor = rf.OBS_balanced_trainprob)
rf.OBS_balanced_testprob <- predict(rf.OBS.2_balanced, newdata = data_simplify_test, type = "prob")[, 2]
rf_test_balanced_roc <- roc(response = data_simplify_test$Depression, predictor = rf.OBS_balanced_testprob)

##SVM
###Imbalance data
continuous_vars <- c("DII")
categorical_vars <- c("Alcohol_intakeNormal","Alcohol_intakeOver","Alcohol_intakePast","GenderMale","SmokingNow","SmokingPast",
                      "Marital_statusNever_married","Marital_status_Widowed_Divorced_Separated",
                      "Health_conditionYes","Trouble_sleepingYes","Seen_mental_health_professionalYes")
preProc <- preProcess(data_simplify_train[, continuous_vars, drop = FALSE], method = c("center", "scale"))
x_train_cont_scaled <- predict(preProc, data_simplify_train[, continuous_vars, drop = FALSE])
x_test_cont_scaled <- predict(preProc, data_simplify_test[, continuous_vars, drop = FALSE])
x_train_scaled <- cbind(x_train_cont_scaled, data_simplify_train[, categorical_vars])
x_test_scaled <- cbind(x_test_cont_scaled, data_simplify_test[, categorical_vars])
x_train <- as.matrix(x_train_scaled)
x_test <- as.matrix(x_test_scaled)
y_train <- data_simplify_train$Depression
y_test <- data_simplify_test$Depression
class_weights <- c("1" = 1, "0" = 8) 
cost_values <- 10^seq(-1, 2, length = 5)  
gamma_values <- 10^seq(-3, 1, length = 5) 
param_grid <- expand.grid(cost = cost_values, gamma = gamma_values)
set.seed(123)
folds <- createFolds(data_simplify_train$Depression, k = 10, list = TRUE, returnTrain = FALSE)
results <- data.frame()
for (i in 1:nrow(param_grid)) {
  current_cost <- param_grid$cost[i]
  current_gamma <- param_grid$gamma[i]
  
  fold_aucs <- numeric(10) 
  for (fold in 1:10) {
    val_indices <- folds[[fold]]
    train_indices <- setdiff(1:nrow(data_simplify_train), val_indices)
    train_data <- data_simplify_train[train_indices, ]
    val_data <- data_simplify_train[val_indices, ]
    preProc <- preProcess(train_data[, continuous_vars, drop = FALSE], method = c("center", "scale"))
    xtrain_cont_scaled <- predict(preProc, train_data[, continuous_vars, drop = FALSE])
    x_val_cont_scaled <- predict(preProc, val_data[, continuous_vars, drop = FALSE])
    xtrain_scaled <- cbind(xtrain_cont_scaled, train_data[, categorical_vars])
    x_val_scaled <- cbind(x_val_cont_scaled, val_data[, categorical_vars])
    xtrain_mat <- as.matrix(xtrain_scaled)
    x_val_mat <- as.matrix(x_val_scaled)
    ytrain_fold <- as.factor(train_data$Depression)
    y_val_fold <- as.factor(val_data$Depression)
    svm_model <- svm(
      x = xtrain_mat,
      y = ytrain_fold,
      kernel = "radial",
      probability = TRUE,
      class.weights = class_weights,  
      cost = current_cost,
      gamma = current_gamma
    )
    val_probs <- attr(predict(svm_model, x_val_mat, probability = TRUE), "prob")[, "1"]  
    roc_obj <- roc(response = y_val_fold, predictor = val_probs, levels = levels(y_val_fold))
    fold_aucs[fold] <- auc(roc_obj)
  }
  mean_auc <- mean(fold_aucs)
  results <- rbind(results, data.frame(
    cost = current_cost,
    gamma = current_gamma,
    auc = mean_auc
  ))
  
  cat(sprintf("cost=%.3f, gamma=%.3f | Mean AUC=%.4f\n", current_cost, current_gamma, mean_auc))
}
results_sorted <- results[order(-results$auc), ]
best_params_svm <- results_sorted[1, ]
print(best_params_svm)
set.seed(123)
svm_model <- svm(
  x = x_train_scaled,
  y = as.factor(y_train),
  kernel = "radial",
  probability = TRUE,
  class.weights = class_weights,
  cost = best_params_svm$cost,
  gamma = best_params_svm$gamma
)
svm_train_probs <- attr(predict(svm_model, x_train, probability = TRUE), "prob")[, 2]
svm_test_probs <- attr(predict(svm_model, x_test, probability = TRUE), "prob")[, 2]
svm_train_roc <- roc(y_train, svm_train_probs)
svm_test_roc <- roc(y_test, svm_test_probs)
###Balance data
continuous_vars <- c("DII")
categorical_vars <- c("Alcohol_intakeNormal","Alcohol_intakeOver","Alcohol_intakePast","GenderMale","SmokingNow","SmokingPast",
                      "Marital_statusNever_married","Marital_status_Widowed_Divorced_Separated",
                      "Health_conditionYes","Trouble_sleepingYes","Seen_mental_health_professionalYes")
preProc <- preProcess(data_simplify_train[, continuous_vars, drop = FALSE], method = c("center", "scale"))
x_train_cont_scaled <- predict(preProc, data_simplify_train[, continuous_vars, drop = FALSE])
x_test_cont_scaled <- predict(preProc, data_simplify_test[, continuous_vars, drop = FALSE])
x_train_scaled <- cbind(x_train_cont_scaled, data_simplify_train[, categorical_vars])
x_test_scaled <- cbind(x_test_cont_scaled, data_simplify_test[, categorical_vars])
x_train <- as.matrix(x_train_scaled)
x_test <- as.matrix(x_test_scaled)
y_train <- data_simplify_train$Depression
y_test <- data_simplify_test$Depression
class_weights <- c("1" = 1, "0" = 8) 
cost_values <- 10^seq(-1, 2, length = 5)  
gamma_values <- 10^seq(-3, 1, length = 5) 
param_grid <- expand.grid(cost = cost_values, gamma = gamma_values)
set.seed(123)
folds <- createFolds(data_simplify_train$Depression, k = 10, list = TRUE, returnTrain = FALSE)
results <- data.frame()

for (i in 1:nrow(param_grid)) {
  current_cost <- param_grid$cost[i]
  current_gamma <- param_grid$gamma[i]
  
  fold_aucs <- numeric(10)
  
  for (fold in 1:10) {
    val_indices <- folds[[fold]]
    train_indices <- setdiff(1:nrow(data_simplify_train), val_indices)
    raw_train <- data_simplify_train[train_indices, ]
    val_data <- data_simplify_train[val_indices, ]
    preProc <- preProcess(raw_train[, continuous_vars, drop = FALSE], method = c("center", "scale"))
    xtrain_cont_scaled <- predict(preProc, raw_train[, continuous_vars, drop = FALSE])
    x_val_cont_scaled <- predict(preProc, val_data[, continuous_vars, drop = FALSE])
    train_for_smote <- cbind(xtrain_cont_scaled, raw_train[, categorical_vars])
    train_for_smote$Depression <- raw_train$Depression
    smote_data <- SMOTE(
      X = train_for_smote[, -which(names(train_for_smote) == "Depression")],
      target = train_for_smote$Depression,
      K = 5,
      dup_size = 3
    )
    balanced_train <- smote_data$data 
    colnames(balanced_train)[ncol(balanced_train)] <- "Depression" 
    balanced_train <- balanced_train %>% mutate(Depression = as.factor(Depression))  
    xtrain_smote <- as.matrix(balanced_train[, -which(names(balanced_train) == "Depression")])
    ytrain_smote <- balanced_train$Depression
    x_val_scaled <- cbind(x_val_cont_scaled, val_data[, categorical_vars])
    x_val_mat <- as.matrix(x_val_scaled)
    y_val_fold <- as.factor(val_data$Depression)
    svm_model <- svm(
      x = xtrain_smote,
      y = ytrain_smote,
      kernel = "radial",
      probability = TRUE,
      class.weights = class_weights,
      cost = current_cost,
      gamma = current_gamma
    )
    val_probs <- attr(predict(svm_model, x_val_mat, probability = TRUE), "prob")[, "1"] 
    roc_obj <- roc(response = y_val_fold, predictor = val_probs, levels = levels(y_val_fold))
    fold_aucs[fold] <- auc(roc_obj)
  }
  mean_auc <- mean(fold_aucs)
  results <- rbind(results, data.frame(
    cost = current_cost,
    gamma = current_gamma,
    auc = mean_auc
  ))
  
  cat(sprintf("cost=%.3f, gamma=%.3f | Mean AUC=%.4f\n", current_cost, current_gamma, mean_auc))
}
results_sorted <- results[order(-results$auc), ]
best_params_svm_balanced <- results_sorted[1, ]
print(best_params_svm_balanced)
full_train_for_smote <- cbind(x_train_cont_scaled, data_simplify_train[, categorical_vars])
full_train_for_smote$Depression <- data_simplify_train$Depression
final_smote <- SMOTE(
  X = full_train_for_smote[, -which(names(full_train_for_smote) == "Depression")],
  target = full_train_for_smote$Depression,
  K = 5,
  dup_size = 3
)
balanced_full_train <- final_smote$data 
colnames(balanced_full_train)[ncol(balanced_full_train)] <- "Depression"  
balanced_full_train <- balanced_full_train %>% mutate(Depression = as.factor(Depression))  

x_train_final <- as.matrix(balanced_full_train[, -which(names(balanced_full_train) == "Depression")])
y_train_final <- balanced_full_train$Depression

set.seed(123)
svm_model_balanced <- svm(
  x = x_train_final,
  y = y_train_final,
  kernel = "radial",
  probability = TRUE,
  class.weights = class_weights,
  cost = best_params_svm_balanced$cost,
  gamma = best_params_svm_balanced$gamma
)
svm_train_probs_balanced <- attr(predict(svm_model_balanced, x_train_final, probability = TRUE), "prob")[, 2]
svm_test_probs_balanced <- attr(predict(svm_model_balanced, x_test, probability = TRUE), "prob")[, 2]
svm_train_roc_balanced <- roc(y_train_final, svm_train_probs_balanced)
svm_test_roc_balanced <- roc(y_test, svm_test_probs_balanced)

##KNN
###Imbalance data
set.seed(123)
print(table(data_simplify_train$Depression))
ctrl <- trainControl(method = "cv",number = 10,classProbs = TRUE,summaryFunction = twoClassSummary  )
k_grid <- expand.grid(k = seq(1, 70, by = 2))
data_simplify_train$Depression <- factor(data_simplify_train$Depression,levels = c("0", "1"),labels = c("No", "Yes"))
data_simplify_test$Depression <- factor(data_simplify_test$Depression,levels = c("0", "1"),labels = c("No", "Yes"))
set.seed(123)
knn_model <- train(Depression ~ ., data = data_simplify_train,method = "knn",trControl = ctrl,                 
                   preProcess = c("center", "scale"),tuneGrid = k_grid,metric = "ROC")
print(knn_model$bestTune)
print(knn_model)
knn_train_prob <- predict(knn_model, newdata = data_simplify_train, "prob")[, 2]
knn_test_prob <- predict(knn_model, newdata = data_simplify_test, "prob")[, 2]
knn_train_roc <- roc(response = data_simplify_train$Depression, predictor = knn_train_prob)
knn_test_roc <- roc(response = data_simplify_test$Depression, predictor = knn_test_prob)
###Balance data
categorical_vars <- c("Alcohol_intakeNormal","Alcohol_intakeOver","Alcohol_intakePast","GenderMale","SmokingNow","SmokingPast",
                      "Marital_statusNever_married","Marital_status_Widowed_Divorced_Separated",
                      "Health_conditionYes","Trouble_sleepingYes","Seen_mental_health_professionalYes")
data_simplify_train$Depression <- factor(data_simplify_train$Depression,levels = c("No", "Yes"), labels = c("0", "1"))  
data_simplify_test$Depression <- factor(data_simplify_test$Depression,levels = c("No", "Yes"),labels = c("0", "1"))
k_grid <- expand.grid(k = seq(1, 70, by = 2))
set.seed(123)
folds <- createFolds(data_simplify_train$Depression, k = 10, list = TRUE)
results <- data.frame(k = integer(), auc = numeric())
for (k_val in k_grid$k) {
  fold_aucs <- numeric(10)
  for (fold in 1:10) {
    val_indices <- folds[[fold]]
    raw_train <- data_simplify_train[-val_indices, ]
    val_data <- data_simplify_train[val_indices, ]
    preProc <- preProcess(raw_train, method = c("center", "scale"))
    train_processed <- predict(preProc, raw_train)
    smote_data <- SMOTE(X = train_processed[, -which(names(train_processed) == "Depression")],
                        target = train_processed$Depression,K = 5,dup_size = 3)
    balanced_train <- smote_data$data 
    colnames(balanced_train)[ncol(balanced_train)] <- "Depression" 
    balanced_train <- balanced_train %>% 
      mutate(Depression = factor(Depression, levels = c("0", "1"))) %>%
      mutate(across(all_of(categorical_vars), ~ round(.x))) 
    val_processed <- predict(preProc, val_data)
    knn_model <- knn3(Depression ~ .,data = balanced_train,k = k_val)
    val_probs <- predict(knn_model, val_processed, type = "prob")[, "1"]
    roc_obj <- roc(response = val_data$Depression, predictor = val_probs)
    fold_aucs[fold] <- auc(roc_obj)}
  results <- rbind(results, data.frame(k = k_val, auc = mean(fold_aucs)))
  cat(sprintf("k=%d | Mean AUC=%.4f\n", k_val, mean(fold_aucs)))}
best_k <- results[which.max(results$auc), "k"]
preProc_full <- preProcess(data_simplify_train, method = c("center", "scale"))
train_full_processed <- predict(preProc_full, data_simplify_train)
final_smote <- SMOTE(
  X = train_full_processed[, -which(names(train_full_processed) == "Depression")],
  target = train_full_processed$Depression,K = 5,dup_size = 3)
balanced_full <- final_smote$data
colnames(balanced_full)[ncol(balanced_full)] <- "Depression"  
balanced_full <- balanced_full %>% 
  mutate(Depression = factor(Depression, levels = c("0", "1"))) %>%
  mutate(across(all_of(categorical_vars), ~ round(.x)))
knn_final <- knn3(Depression ~ .,data = balanced_full,k = best_k)
test_processed <- predict(preProc_full, data_simplify_test)
knn_train_prob_balanced <- predict(knn_final, train_full_processed, type = "prob")[, "1"]
knn_test_prob_balanced <- predict(knn_final, test_processed, type = "prob")[, "1"]
knn_train_roc_balanced <- roc(response = data_simplify_train$Depression, predictor = knn_train_prob)
knn_test_roc_balanced <- roc(response = data_simplify_test$Depression, predictor = knn_test_prob)




##XGBoost
###Imbalance data
data_simplify_test$Depression <- ifelse(data_simplify_test$Depression == "1", 1, 0)
data_simplify_train$Depression <- ifelse(data_simplify_train$Depression == "1", 1, 0)
dtrain <- xgb.DMatrix(
  data = as.matrix(data_simplify_train[, -which(colnames(data_simplify_train) == "Depression")]),
  label = data_simplify_train$Depression
)
dtest <- xgb.DMatrix(
  data = as.matrix(data_simplify_test[, -which(colnames(data_simplify_test) == "Depression")]),
  label = data_simplify_test$Depression
)
scale_pos_weight <- sum(data_simplify_train$Depression == 0) / sum(data_simplify_train$Depression == 1)
max_depth_values <- c(3, 5, 7)
eta_values <- c(0.01, 0.1, 0.3)
best_auc <- 0
for (max_depth in max_depth_values) {
  for (eta in eta_values) {params <- list(objective = "binary:logistic",eval_metric = "auc",max_depth = max_depth,
                                          eta = eta,subsample = 0.8,colsample_bytree = 0.8)
  set.seed(123)
  cv <- xgb.cv(params, dtrain, nrounds = 500, nfold = 10,early_stopping_rounds = 20, 
               scale_pos_weight = scale_pos_weight,verbose = FALSE)
  current_auc <- max(cv$evaluation_log$test_auc_mean)
  if (current_auc > best_auc) {
    best_auc <- current_auc
    best_params <- params
    best_nrounds <- cv$best_iteration}}}
set.seed(123)
xgb_model <- xgb.train(best_params, dtrain, best_nrounds, print_every_n = 10) 
xgb_testpredprob <- predict(xgb_model, dtest)
xgb_trainpredprob <- predict(xgb_model, dtrain)
xgb_train_roc <- roc(data_simplify_train$Depression, xgb_trainpredprob)
xgb_test_roc <- roc(data_simplify_test$Depression, xgb_testpredprob)
##Balanced data
max_depth_values <- c(3, 5, 7)
eta_values <- c(0.01, 0.1, 0.3)
best_auc <- 0
best_params <- list()
best_nrounds <- 0
set.seed(123)
folds <- createFolds(data_simplify_train$Depression, k = 10, list = TRUE)
for (max_depth in max_depth_values) {
  for (eta in eta_values) {
    fold_aucs <- numeric(10)
    for (fold in 1:10) {
      val_indices <- folds[[fold]]
      raw_train <- data_simplify_train[-val_indices, ]
      val_data <- data_simplify_train[val_indices, ]
      smote_data <- SMOTE(
        X = raw_train[, -which(names(raw_train) == "Depression")],
        target = raw_train$Depression,
        K = 5,
        dup_size = 3
      )
      balanced_train <- smote_data$data       
      colnames(balanced_train)[ncol(balanced_train)] <- "Depression"  
      balanced_train <- balanced_train %>%
        mutate(Depression = as.numeric(Depression)) 
      dtrain_smote <- xgb.DMatrix(
        data = as.matrix(balanced_train[, -which(names(balanced_train) == "Depression")]),
        label = balanced_train$Depression
      )
      
      dval <- xgb.DMatrix(
        data = as.matrix(val_data[, -which(names(val_data) == "Depression")]),
        label = val_data$Depression
      )
      params <- list(
        objective = "binary:logistic",
        eval_metric = "auc",
        max_depth = max_depth,
        eta = eta,
        subsample = 0.8,
        colsample_bytree = 0.8
      )
      set.seed(123)
      model <- xgb.train(
        params,
        dtrain_smote,
        nrounds = 500,
        early_stopping_rounds = 20,
        watchlist = list(val = dval),
        verbose = 0
      )
      fold_aucs[fold] <- model$best_score
    }
    current_auc <- mean(fold_aucs)
    if (current_auc > best_auc) {
      best_auc <- current_auc
      best_params_balanced <- params
      best_nrounds_balanced <- which.max(model$evaluation_log$val_auc)
    }
  }
}
final_smote <- SMOTE(
  X = data_simplify_train[, -which(names(data_simplify_train) == "Depression")],
  target = data_simplify_train$Depression,
  K = 5,
  dup_size = 3
)

balanced_full_train <- final_smote$data       
colnames(balanced_full_train)[ncol(balanced_full_train)] <- "Depression"  
balanced_full_train <- balanced_full_train %>%
  mutate(Depression = as.numeric(Depression)) 

dtrain_final <- xgb.DMatrix(
  data = as.matrix(balanced_full_train[, -which(names(balanced_full_train) == "Depression")]),
  label = balanced_full_train$Depression
)
dtest <- xgb.DMatrix(
  data = as.matrix(data_simplify_test[, -which(names(data_simplify_test) == "Depression")]),
  label = data_simplify_test$Depression
)

set.seed(123)
xgb_model_balanced <- xgb.train(
  best_params_balanced,
  dtrain_final,
  nrounds = best_nrounds_balanced,
  print_every_n = 10
)
xgb_testpredprob_balanced <- predict(xgb_model_balanced, dtest)
xgb_trainpredprob_balanced <- predict(xgb_model_balanced, dtrain_final)
xgb_train_roc_balanced <- roc(balanced_full_train$Depression, xgb_trainpredprob_balanced)
xgb_test_roc_balanced <- roc(data_simplify_test$Depression, xgb_testpredprob_balanced)


##logistic regression model
###Imbalance data
data_simplify_test$Depression <- ifelse(data_simplify_test$Depression == 1, "Yes", "No")
data_simplify_train$Depression <- ifelse(data_simplify_train$Depression == 1, "Yes", "No")
set.seed(123) 
ctrl <- trainControl(method = "cv",number = 10,savePredictions = "final",classProbs = TRUE,summaryFunction = twoClassSummary)
model <- train(Depression ~ .,data = data_simplify_train,method = "glm",family = binomial,
               trControl = ctrl,metric = "ROC")
cv_results <- model$results
print(cv_results)
lr_prob_train <- predict(model, data_simplify_train, type = "prob")[ ,2]
lr_prob_test <- predict(model, data_simplify_test, type = "prob")[ ,2]
lr_train_roc <- roc(response = data_simplify_train$Depression, predictor = lr_prob_train)
lr_test_roc <- roc(response = data_simplify_test$Depression, predictor = lr_prob_test)
###Balance data
#upsampling
data_simplify_train$Depression <- ifelse(data_simplify_train$Depression == "Yes", 1, 0)
set.seed(123)
smote_result <- SMOTE(X = data_simplify_train[, -13],target = data_simplify_train$Depression,
                      K = 5,dup_size = 3)
balanced_data <- smote_result$data
colnames(balanced_data)[13] <- "Depression"
balanced_data$Depression <- ifelse(balanced_data$Depression == 1, "Yes", "No")
set.seed(123) 
ctrl <- trainControl(method = "cv",number = 10,savePredictions = "final",classProbs = TRUE,summaryFunction = twoClassSummary)
model_balanced <- train(Depression ~ .,data = balanced_data,method = "glm",family = binomial, 
                        trControl = ctrl,metric = "ROC")
cv_results_balanced <- model_balanced$results
print(cv_results_balanced)
lr_prob_train_balanced <- predict(model_balanced, balanced_data, type = "prob")[ ,2]
lr_prob_test_balanced <- predict(model_balanced, data_simplify_test, type = "prob")[ ,2]
lr_train_roc_balanced <- roc(response = balanced_data$Depression, predictor = lr_prob_train_balanced)
lr_test_roc_balanced <- roc(response = data_simplify_test$Depression, predictor = lr_prob_test_balanced)



##ROC curves of 10 models
train_rocs <- list(
  RF = rf_train_roc,
  XGBoost = xgb_train_roc,
  SVM = svm_train_roc,
  KNN = knn_train_roc,
  LR = lr_train_roc,
  RF_SM = rf_train_balanced_roc,
  XGBoost_SM = xgb_train_roc_balanced,
  SVM_SM = svm_train_roc_balanced,
  KNN_SM = knn_train_roc_balanced,
  LR_SM = lr_train_roc_balanced
)
test_rocs <- list(
  RF = rf_test_roc,
  XGBoost = xgb_test_roc,
  SVM = svm_test_roc,
  KNN = knn_test_roc,
  LR = lr_test_roc,
  RF_SM = rf_test_balanced_roc,
  XGBoost_SM = xgb_test_roc_balanced,
  SVM_SM = svm_test_roc_balanced,
  KNN_SM = knn_test_roc_balanced,
  LR_SM = lr_test_roc_balanced
)
model_colors <- c(
  RF          = "#4E79A7",  
  XGBoost     = "#F28E2B",  
  SVM         = "#E15759",  
  KNN         = "#76B7B2",  
  LR          = "#59A14F",  
  RF_SM       = "#EDC948",  
  XGBoost_SM  = "#B07AA1",  
  SVM_SM      = "#FF9DA7",  
  KNN_SM      = "#9C755F",  
  LR_SM       = "#BAB0AC"   
)
create_roc_plot <- function(roc_list) {
  roc_data <- imap_dfr(roc_list, function(roc_obj, model_name) {
    tibble(
      Model = model_name,
      FPR = 1 - roc_obj$specificities,
      TPR = roc_obj$sensitivities,
      AUC = auc(roc_obj) %>% round(3)
    )
  })
  auc_labels <- roc_data %>%
    distinct(Model, AUC) %>%  
    mutate(Label = paste0(Model, " = ", AUC))
  
  
  label_y <- seq(0.4, 0.002, length.out = 10)
  
  ggplot(roc_data, aes(x = FPR, y = TPR, color = Model)) +
    geom_line(size = 0.8) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +
    scale_color_manual(values = model_colors) +
    coord_equal() +
    theme_bw() +
    labs(x = "1 - Specificity",
         y = "Sensitivity") +
    theme(legend.position = "none",
          plot.tag = element_text(size = 12, face = "bold", colour = "black"),  
          plot.tag.position = c(0, 1.05),
          plot.title = element_text(hjust = 0.5)) +
    annotate("text", 
             x = rep(0.6, 10), y = label_y,
             label = auc_labels$Label,
             color = model_colors[auc_labels$Model],
             hjust = 0,
             size = 3,
             fontface = "bold")
}
ROC_train_plot <- create_roc_plot(train_rocs)
ROC_test_plot <- create_roc_plot(test_rocs)
print(ROC_train_plot)
print(ROC_test_plot)




##LR
###calibration curve
data_simplify_test <- data_simplify_test %>% mutate(across(c(Depression), as.factor))
calib_data_test <- data.frame(Observed = as.numeric(data_simplify_test$Depression)-1,
                              Predicted = lr_prob_test)
calib_test <- val.prob(calib_data_test$Predicted, calib_data_test$Observed)
###Decision curve
LR <- lr_prob_test
test_dca_result <- dca(Depression ~ LR, data = data_simplify_test)
DCA_plot <- plot(test_dca_result, colorize = TRUE) + theme_minimal() +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5))
print(DCA_plot)


###line graph
data_simplify <- data_OBS_log[, c("DII","Health_condition",
                                  "Seen_mental_health_professional", "Trouble_sleeping",
                                  "Marital_status", "Gender","Alcohol_intake","Smoking",
                                  "Depression")]
set.seed(123)
trains <- createDataPartition(y = data_simplify$Depression, p = 0.7, list = F)
data_simplify_train <- data_simplify[trains, ]
data_simplify_test <- data_simplify[-trains, ]
ddist <- datadist(data_simplify)
options(datadist = "ddist")
model_lrm <- lrm(Depression ~ Trouble_sleeping + Health_condition + Seen_mental_health_professional + 
                   DII + Smoking + Marital_status + Gender + Alcohol_intake, 
                 data = data_simplify_train,
                 x = TRUE, y = TRUE 
)
nom <- nomogram(model_lrm, 
                fun = plogis,
                funlabel = "Depression Risk")
plot(nom,
     xfrac = 0.25,    
     cex.var = 0.8,   
     cex.axis = 0.7,  
     lmgp = 0.1       
)

str(data_simplify_train)

##SHAP interpretation
pfun <- function(object, newdata) {predict(object, newdata = newdata, type = "fitted")}
X <- data_simplify_train[, names(data_simplify_train) != "Depression"]
shap <- explain(model_lrm,X = X,nsim = 50,pred_wrapper = pfun)
shap_plots <- shapviz(shap, X = X)
##Global explanation
sv_importance(shap_plots, kind = "beeswarm") + theme_bw() + theme(
  axis.title.y = element_text(face = "bold", size = 12),
  axis.text.y = element_text(face = "bold", size = 10))
sv_importance(shap_plots) + theme_bw() + theme(
  axis.title.y = element_text(face = "bold", size = 12),
  axis.text.y = element_text(face = "bold", size = 10))
sv_dependence(shap_plots, v = "DII")+theme_bw()

##Individual Force Plot
positive_samples <- which(data_simplify_train$Depression == "Yes")  
negative_samples <- which(data_simplify_train$Depression == "No")  
set.seed(22)  
selected_id <- sample(positive_samples, 1)
sv_force(shap_plots, row_id = selected_id) + theme_bw()
set.seed(123)  
selected_id <- sample(negative_samples, 1)
sv_force(shap_plots, row_id = selected_id) + theme_bw()
