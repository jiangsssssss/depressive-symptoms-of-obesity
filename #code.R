#Packages
library(foreign)
library(tidyverse)
library(tidyselect)
library(dplyr)
library(broom)
library(dietaryindex)
library(xgboost)
library(caret)
library(pROC)
library(viridis)
library(shapviz)
library(ggforce)
library(rms)
library(pmsampsize)
library(randomForest)
library(purrr)
library(boot)
library(DescTools)
library(dcurves)
library(treeshap)

#Participants
##Import data(exclude dietary data)
DATA=read.csv("D:/Clinical prediction model/01data/clean/DATA.csv")
##Import 15-16 dietary data (dietaryindex comes with 17-18 dietary data)
setwd("D:/Clinical prediction model/01data/raw/I/")
load("NHANES_20152016.RDA")
##Calculate Dietary Inflammation Index
data("NHANES_20152016")
DII_I=DII_NHANES_FPED(FPED_PATH=NHANES_20152016$FPED, NUTRIENT_PATH=NHANES_20152016$NUTRIENT, DEMO_PATH=NHANES_20152016$DEMO, FPED_PATH2=NHANES_20152016$FPED2, NUTRIENT_PATH2=NHANES_20152016$NUTRIENT2)
data("NHANES_20172018")
DII_J=DII_NHANES_FPED(FPED_PATH=NHANES_20172018$FPED, NUTRIENT_PATH=NHANES_20172018$NUTRIENT, DEMO_PATH=NHANES_20172018$DEMO, FPED_PATH2=NHANES_20172018$FPED2, NUTRIENT_PATH2=NHANES_20172018$NUTRIENT2)
##Merge data and include samples with BMI >30
DII_Total=bind_rows(DII_I,DII_J)
DATA <- DATA %>%
  left_join(DII_Total, by = "SEQN")
DATA_OBS <- subset(DATA, BMXBMI > 30)

#Predictors and Outcome
##Alcohol_intake
DATA_OBS$`Alcohol_intake`<- NA
DATA_OBS$`Alcohol_intake`[DATA_OBS$ALQ110 == 2] <- "No"
DATA_OBS$`Alcohol_intake`[(DATA_OBS$ALQ110 == 1 | DATA_OBS$ALQ101 == 1) & DATA_OBS$ALQ120Q == 0] <- "Past"
DATA_OBS$`Alcohol_intake`[DATA_OBS$ALQ120Q > 0 & DATA_OBS$ALQ120Q <= 365 & DATA_OBS$ALQ141Q == 0 & 
    ((DATA_OBS$RIAGENDR == 1 & DATA_OBS$ALQ130 <= 2) |(DATA_OBS$RIAGENDR == 2 & DATA_OBS$ALQ130 <= 1))] <- "Normal"
DATA_OBS$`Alcohol_intake`[(DATA_OBS$RIAGENDR == 1 & DATA_OBS$ALQ130 > 2) |
    (DATA_OBS$RIAGENDR == 2 & DATA_OBS$ALQ130 > 1)] <- "Over"
DATA_OBS$`Alcohol_intake`[DATA_OBS$ALQ141Q > 0 & DATA_OBS$ALQ141Q <= 365] <- "Over"
##Smoking
DATA_OBS <- DATA_OBS %>%
  mutate(Smoking = case_when(SMQ020 == 2 ~ "Never",SMQ020 == 1 & SMQ040 == 3 ~ "Past",SMQ040 %in% c(1, 2) ~ "Now",TRUE ~ NA_character_  ))
##Moderate recreational activities
DATA_OBS <- DATA_OBS %>%
  mutate(MRA = case_when(PAQ665 == 1 ~ "Yes",PAQ665 == 2 ~ "No",PAQ665 %in% c(7, 9) | is.na(PAQ665) ~ NA_character_,TRUE ~ NA_character_  ))
##Sedentary behavior
DATA_OBS <- DATA_OBS %>%
  mutate(`SB` = case_when(PAD680 != 7777 & PAD680 != 9999 & PAD680 >= 480 ~ "Yes", PAD680 != 7777 & PAD680 != 9999 & PAD680 < 480 ~ "No",TRUE ~ NA_character_))
##Health_condition
DATA_OBS <- DATA_OBS %>%
  mutate(`Health_condition` = case_when(HSD010 %in% c(1, 2, 3) ~ "Yes",HSD010 %in% c(4, 5) ~ "No",TRUE ~ NA_character_))
##Trouble_sleeping
DATA_OBS <- DATA_OBS %>%
  mutate(Trouble_sleeping = case_when(SLQ050 == 1 ~ "Yes",SLQ050 == 2 ~ "No",SLQ050 %in% c(7, 9) | is.na(SLQ050) ~ NA_character_,TRUE ~ NA_character_  ))
##Seen_mental_health_professional
DATA_OBS <- DATA_OBS %>%
  mutate(Seen_mental_health_professional = case_when(HUQ090 == 1 ~ "Yes",HUQ090 == 2 ~ "No",HUQ090 %in% c(7, 9) | is.na(SLQ050) ~ NA_character_,TRUE ~ NA_character_  ))
##Hypertension
DATA_OBS <- DATA_OBS %>%
  mutate(Sys_Avg = (BPXSY1 + BPXSY2 + BPXSY3) / 3,Dia_Avg = (BPXDI1 + BPXDI2 + BPXDI3) / 3,Hypertension = if_else(BPQ020 == 1 | Sys_Avg > 140 | Dia_Avg > 90, "Yes","No")) 
##Diabetes
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
##Coronary heart disease
DATA_OBS <- DATA_OBS %>%
  mutate(CHD = if_else(MCQ160C %in% c(7, 9), NA_real_, MCQ160C))
DATA_OBS <- DATA_OBS %>%
  mutate(CHD = case_when(CHD == 2 ~ "No",CHD == 1 ~ "Yes",TRUE ~ as.character(CHD)))
##Congestive heart failure
DATA_OBS <- DATA_OBS %>%
  mutate(CHF = if_else(MCQ160B %in% c(7, 9), NA_real_, MCQ160B))
DATA_OBS <- DATA_OBS %>%
  mutate(CHF = case_when(CHF == 2 ~ "No",CHF == 1 ~ "Yes",TRUE ~ as.character(CHD)))
##Angina
DATA_OBS <- DATA_OBS %>%
  mutate(Angina = if_else(MCQ160D %in% c(7, 9), NA_real_, MCQ160D))
DATA_OBS <- DATA_OBS %>%
  mutate(Angina = case_when(Angina == 2 ~ "No",Angina == 1 ~ "Yes",TRUE ~ as.character(Angina)  ))
##Heart attack
DATA_OBS <- DATA_OBS %>%
  mutate(Heartattack = if_else(MCQ160E %in% c(7, 9), NA_real_, MCQ160E))
DATA_OBS <- DATA_OBS %>%
  mutate(Heartattack = case_when(Heartattack == 2 ~ "No",Heartattack == 1 ~ "Yes", TRUE ~ as.character(Heartattack)))
##Stroke
DATA_OBS <- DATA_OBS %>%
  mutate(Stroke = if_else(MCQ160F %in% c(7, 9), NA_real_, MCQ160F))
DATA_OBS <- DATA_OBS %>%
  mutate(Stroke = case_when(Stroke == 2 ~ "No",Stroke == 1 ~ "Yes",TRUE ~ as.character(Stroke)))
##Asthma
DATA_OBS <- DATA_OBS %>%
  mutate(Asthma = if_else(MCQ010 %in% c(7, 9), NA_real_, MCQ010))
DATA_OBS <- DATA_OBS %>%
  mutate(Asthma = case_when(Asthma == 2 ~ "No",Asthma == 1 ~ "Yes",TRUE ~ as.character(Asthma)))
#Arthritis
DATA_OBS <- DATA_OBS %>%
  mutate(Arthritis = if_else(MCQ160A %in% c(7, 9), NA_real_, MCQ160A))
DATA_OBS <- DATA_OBS %>%
  mutate(Arthritis = case_when(Arthritis == 2 ~ "No",Arthritis == 1 ~ "Yes",TRUE ~ as.character(Arthritis)))
##Obstructive sleep apnea
DATA_OBS <- DATA_OBS %>%
  mutate(OSA = case_when(SLQ030 %in% c(2, 3) |SLQ040 %in% c(2, 3) |SLQ120 == 4 ~ "Yes",
  (SLQ120 %in% c(7, 9) | is.na(SLQ120)) &(SLQ040 %in% c(7, 9) | is.na(SLQ040)) & (SLQ030 %in% c(7, 9) | is.na(SLQ030)) ~ NA_character_,TRUE ~ "No"))
##Gender
DATA_OBS <- DATA_OBS %>%
  mutate(Gender = case_when(RIAGENDR == 1 ~ "Male",RIAGENDR == 2 ~ "Female",TRUE ~ NA_character_ ))
##Race
DATA_OBS <- DATA_OBS %>%
  mutate(Race = case_when(
    RIDRETH1 == 1 ~ "Mexican American",RIDRETH1 == 2 ~ "Other Hispanic",RIDRETH1 == 3 ~ "Non-Hispanic White",RIDRETH1 == 4 ~ "Non-Hispanic Black",RIDRETH1 == 5 ~ "Other Race",TRUE ~ NA_character_  ))
##Education
DATA_OBS$Education <- NA
DATA_OBS$Education[DATA_OBS$DMDEDUC2 %in% c(1, 2)] <- "less than high school"
DATA_OBS$Education[DATA_OBS$DMDEDUC2 == 3] <- "completed high school"
DATA_OBS$Education[DATA_OBS$DMDEDUC2 %in% c(4, 5)] <- "more than high school"
##Marital status
DATA_OBS$Marital_status <- NA
DATA_OBS$Marital_status[DATA_OBS$DMDMARTL %in% c(1, 6)] <- "Married / Living with partner"
DATA_OBS$Marital_status[DATA_OBS$DMDMARTL == 5] <- "Never married"
DATA_OBS$Marital_status[DATA_OBS$DMDMARTL %in% c(2, 3, 4)] <- "Widowed / Divorced/Separated"
##Poverty-to-income ratio 
DATA_OBS$PIR <- NA
DATA_OBS$PIR[DATA_OBS$INDFMPIR < 1.3] <- "<1.3"
DATA_OBS$PIR[DATA_OBS$INDFMPIR >= 1.3 & DATA_OBS$INDFMPIR <= 3.5] <- "1.3-3.5"
DATA_OBS$PIR[DATA_OBS$INDFMPIR > 3.5] <- ">3.5"
##Health_insurance
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
##Depression
sum_vars <- rowSums(DATA_OBS[, c("DPQ010", "DPQ020", "DPQ030", "DPQ040", "DPQ050", "DPQ060", "DPQ070", "DPQ080", "DPQ090")], na.rm = TRUE)
has_7_or_9 <- apply(DATA_OBS[, c("DPQ010", "DPQ020", "DPQ030", "DPQ040", "DPQ050", "DPQ060", "DPQ070", "DPQ080", "DPQ090")], 1, function(x) any(x %in% c(7, 9)))
has_na <- apply(DATA_OBS[, c("DPQ010", "DPQ020", "DPQ030", "DPQ040", "DPQ050", "DPQ060", "DPQ070", "DPQ080", "DPQ090")], 1, function(x) any(is.na(x)))
DATA_OBS$Depression <- ifelse(has_na | has_7_or_9, NA, ifelse(sum_vars >= 10, "Yes", "No"))

#Exclusion
data_OBS <- DATA_OBS[, c("SEQN", "Gender", "Age", "BMI", "Race","DII", "Alcohol_intake", 
                         "Smoking", "Hypertension", "Diabetes","CHD", "CHF", "Angina", 
                         "Heartattack", "Stroke", "OSA", "Asthma", "Arthritis", "Education", 
                         "Marital_status", "PIR", "Health_insurance", "MRA", "SB", 
                         "Trouble_sleeping", "Health_condition", "Seen_mental_health_professional",
                         "HS_CRP", "TC","HDL","Depression")]
missing_values <- sapply(data_OBS, function(x) sum(is.na(x)))
missing_values_df <- data.frame(Variable = names(missing_values), MissingValues = missing_values)
print(missing_values_df)
data_OBS <- data_OBS[data_OBS$Age >= 18, ]
data_OBS <- na.omit(data_OBS)
data_OBS <- data_OBS %>% select(-SEQN)

#Sample description
data_OBS <- data_OBS %>%
  mutate(across(c(Gender, Race, Alcohol_intake, Hypertension, Diabetes, CHD, 
                  CHF,  Angina, Heartattack, Stroke, OSA, Asthma, Trouble_sleeping, Education, 
                  Marital_status, PIR, Health_insurance, MRA, SB, Smoking, Arthritis, 
                  Health_condition, Seen_mental_health_professional,Depression), as.factor))
categorical_vars <- c("Gender", "Race", "Alcohol_intake", "Smoking", "Hypertension", "Diabetes", "CHD", "CHF", "Angina", "Heartattack", "Stroke", "OSA", "Asthma", "Arthritis", "Education", "Marital_status", "PIR", "Health_insurance", "MRA", "SB", "Trouble_sleeping", "Seen_mental_health_professional", "Health_condition")
continuous_vars <- c("Age", "BMI", "DII", "TC", "HDL", "HS_CRP")
calc_mean_ci <- function(data, variable) {est <- mean(data, na.rm = TRUE)
  se <- sd(data, na.rm = TRUE) / sqrt(length(data))
  lower <- est - 1.96 * se
  upper <- est + 1.96 * se
  return(c(round(est, 3), round(lower, 3), round(upper, 3)))}
calc_percent_and_ci <- function(x, n) {
  percent <- round(x / n * 100, 3)  
  lower <- round(percent - 1.96 * sqrt((x / n) * (1 - x / n) / n) * 100, 3)  
  upper <- round(percent + 1.96 * sqrt((x / n) * (1 - x / n) / n) * 100, 3)  
  return(c(percent, lower, upper))}
Sample_description_results <- data.frame(Variable = character(),
  No_Percent = numeric(),No_CI_Lower = numeric(),No_CI_Upper = numeric(),
  Yes_Percent = numeric(),Yes_CI_Lower = numeric(),Yes_CI_Upper = numeric(),
  P_Value = numeric(), stringsAsFactors = FALSE)
for (var in continuous_vars) {
  no_data <- data_OBS[data_OBS$Depression == "No", var]
  yes_data <- data_OBS[data_OBS$Depression == "Yes", var]
  no_mean_ci <- calc_mean_ci(no_data, var)
  yes_mean_ci <- calc_mean_ci(yes_data, var)
  t_test <- t.test(no_data, yes_data, var.equal = TRUE)
  Sample_description_results <- rbind(Sample_description_results, data.frame(Variable = var,
  No_Percent = no_mean_ci[1],No_CI_Lower = no_mean_ci[2],No_CI_Upper = no_mean_ci[3],
  Yes_Percent = yes_mean_ci[1],Yes_CI_Lower = yes_mean_ci[2],Yes_CI_Upper = yes_mean_ci[3],
  P_Value = round(t_test$p.value, 4)))}
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
  No_Percent = no_percent_ci[1],No_CI_Lower = no_percent_ci[2],No_CI_Upper = no_percent_ci[3],
  Yes_Percent = yes_percent_ci[1],Yes_CI_Lower = yes_percent_ci[2],Yes_CI_Upper = yes_percent_ci[3],
  P_Value = round(chi_result$p.value, 4)  ))}}
print(Sample_description_results)

#Feature Selection
##XGBoost
###Data preparation
set.seed(42)
trains <- createDataPartition(y = data_OBS$Depression,p =0.7,list = F)
data_train <- data_OBS[trains, ]
data_test <- data_OBS[-trains, ]
table(data_train$Depression)
table(data_test$Depression)
dvfunc <- dummyVars(~., data = data_train[, 1:29], fullRank = T)
data_trainx <- predict(dvfunc, newdata = data_train[, 1:29]) 
data_trainy <- ifelse(data_train$Depression == "Yes", "1", "0")
data_testx <- predict(dvfunc, newdata = data_test[, 1:29]) 
data_testy <- ifelse(data_test$Depression == "Yes", "1", "0")
dtrain <- xgb.DMatrix(data = data_trainx,label = data_trainy)
dtest <- xgb.DMatrix(data = data_testx,label = data_testy)
###Modeling
params <- list(
  objective ="binary:logistic",
  eval_metric ="auc",
  eta = 0.01,
  max_depth = 4,
  subsample = 0.8,
  colsample_bytree = 0.8)
set.seed(123)
cv <- xgb.cv(params = params,
             data = dtrain,
             nrounds =1000,
             nfold =10,
             early_stopping_rounds =10)
best_iter <- which.max(cv$evaluation_log$test_auc_mean) 
model <- xgb.train(params = params,
                   data = dtrain, nrounds = best_iter)
testpredprob <- predict (model, dtest)
roc_xgb <- roc(data_test$Depression, testpredprob)
auc <- auc (roc_xgb)
ROC_of_XGBoost<- ggroc (roc_xgb,legacy.axes = TRUE,size = 1,color = "#69b3a2")+
  geom_abline (intercept = 0,slope = 1,linetype ="dashed",color = "grey")+
  theme_bw()+annotate(geom = "text",label = paste0("AUC in test set:", round(auc,3)), x = 0.7,y = 0.05)
ROC_of_XGBoost
###Variables importance
importance <- xgb.importance(model = model)
importance <- importance [order(importance$Gain,
                                decreasing = TRUE),][1:10,1:2]
Variables_importance_plot_of_XGBoost <- ggplot(importance,aes(y = Gain,x = reorder(Feature, Gain))) +
  geom_bar(stat = "identity",fill = "#69b3a2", alpha = 0.8, width = 0.6) +scale_fill_viridis() +
  labs(y = "Importance Score",x = "") +coord_flip()+theme_bw()+theme(axis.text.y = element_text(angle = 0, hjust =0.9))
Variables_importance_plot_of_XGBoost
importance$Feature
Feature1 = c("Seen_mental_health_professional", "DII", "Health_condition", 
             "Trouble_sleeping", "BMI", "TC", "Age", "Arthritis", 
             "HS_CRP", "HDL")
##SHAP
shap_long_hd <- shapviz(model, X_pred = data_trainx)
sv_importance(shap_long_hd)+theme_bw()
Feature2 = c("Seen_mental_health_professional", "DII", "Health_condition", 
             "Trouble_sleeping", "BMI", "TC", "Age", "Arthritis", 
             "HS_CRP", "HDL")
##Boruta
library(Boruta)
set.seed(1)
feature.selection <- Boruta(Depression ~ ., data = data_OBS, doTrace = 1)
table(feature.selection$finalDecision)
par(mar=c(7.5, 4, 0.5, 2) + 0.1) 
plot(feature.selection, cex.axis=0.5, las=3, xlab="")
Feature3 <- getSelectedAttributes(feature.selection) 
Feature3
##Intersection
intersection <- Reduce(intersect, list(Feature1, Feature2, Feature3))
print(intersection)
data_feature <- data_OBS[, c(intersection,"Depression")]

#Modeling
##ample size estimation
pmsampsize(type = "b",csrsquared = 0.072,parameters = 7,prevalence = 0.1,seed = 123)
##Data preparation
data_feature$Seen_mental_health_professional <- ifelse(data_feature$Seen_mental_health_professional == "Yes", 1, 0)
data_feature$Health_condition <- ifelse(data_feature$Health_condition == "Yes", 1, 0)
data_feature$Trouble_sleeping <- ifelse(data_feature$Trouble_sleeping == "Yes", 1, 0)
data_feature$Arthritis <- ifelse(data_feature$Arthritis == "Yes", 1, 0)
data_feature$Depression <- ifelse(data_feature$Depression == "Yes", 1, 0)
data_feature <- data_feature %>% mutate(across(c(Depression), as.factor))
str(data_feature)
set.seed(123)
trains <- createDataPartition(y = data_feature$Depression,p =0.7,list = F )
data_feature_train <- data_feature[trains, ]
data_feature_test <- data_feature[-trains, ]
table(data_feature_train$Depression)
table(data_feature_test$Depression)
##Modeling Random Forests
set.seed(123)
rf.OBS <- randomForest(Depression ~ ., data = data_feature_train)
plot(rf.OBS)
which.min(rf.OBS$err.rate[, 1])
###ntree = 131
set.seed(123)
rf.OBS.2 <- randomForest(Depression ~ ., data = data_feature_train, ntree = 131)
rf.OBS.2
###ROC of training set  
rf.OBS.trainprob <- predict(rf.OBS.2,newdata = data_feature_train,type = "prob")
trainroc <- roc(response = data_feature_train$Depression,predictor = rf.OBS.trainprob[ ,2])
ROC_of_Training_set <- ggroc(trainroc,legacy.axes = TRUE,size = 1,color = "black") +
  geom_abline(intercept = 0,slope = 1,linetype = "dashed",color = "grey") + 
  theme_bw() + annotate(geom = "text",label = paste0("AUC in training set: ", round(auc(trainroc), 3)),x = 0.7,y = 0.05)  
ROC_of_Training_set
###ROC of test set 
rf.OBS.testprob <- predict(rf.OBS.2,newdata = data_feature_test,type = "prob")
testroc <- roc(response = data_feature_test$Depression,predictor = rf.OBS.testprob[ ,2])
ROC_of_Test_set <- ggroc(testroc,legacy.axes = TRUE,size = 1,color = "black") +  
  geom_abline(intercept = 0,slope = 1,linetype = "dashed",color = "grey") +  # 参考线
  theme_bw() + annotate(geom = "text",label = paste0("AUC in test set: ", round(auc(testroc), 3)), x = 0.7,y = 0.05)  
ROC_of_Test_set
###Calibration curve of training set
y_train <- as.numeric(data_feature_train$Depression)-1
prob.rf.train <- as.numeric(predict(rf.OBS.2, data_feature_train, type="prob")[,2])
breaks <- seq(0, 1, by=0.1)
ypb.rf.train <- cut(prob.rf.train, breaks=breaks, include.lowest = TRUE, labels=FALSE)
am.rf.train  <- tapply(y_train, ypb.rf.train, mean)
pm.rf.train <- tapply(prob.rf.train, ypb.rf.train, mean)
df.rf.train <- data.frame(pm = pm.rf.train, am = am.rf.train)
ggplot(df.rf.train, aes(x=pm, y=am)) +
  geom_point() +
  geom_line() +
  geom_abline(intercept=0, slope=1, color="black", linetype="dashed")+
  xlab("Mean Predicted Probability") + ylab("True Fraction of Positives") +
  theme(plot.title = element_text(hjust = 0.5))+
  theme_bw() + 
  scale_y_continuous(limits=c(0,1))
###Calibration curve of test set
y_test <- as.numeric(data_feature_test$Depression)-1
prob.rf.test <- as.numeric(predict(rf.OBS.2, data_feature_test, type="prob")[,2])
breaks <- seq(0, 1, by=0.1)
ypb.rf.test <- cut(prob.rf.test, breaks=breaks, include.lowest = TRUE, labels=FALSE)
am.rf.test  <- tapply(y_test, ypb.rf.test, mean)
pm.rf.test <- tapply(prob.rf.test, ypb.rf.test, mean)
df.rf.test <- data.frame(pm = pm.rf.test, am = am.rf.test)
ggplot(df.rf.test, aes(x=pm, y=am)) +
  geom_point() +
  geom_line() +
  geom_abline(intercept=0, slope=1, color="black", linetype="dashed")+
  xlab("Mean Predicted Probability") + ylab("True Fraction of Positives") +
  theme(plot.title = element_text(hjust = 0.5))+
  theme_bw() + 
  scale_y_continuous(limits=c(0,1))
###AUC and brier score
boot_auc <- function(data, indices) {
  d <- data[indices, ]
  pred <- predict(rf.OBS.2, newdata = d, type = "prob")[,2]
  roc_obj <- roc(response = d$Depression, predictor = pred)
  return(auc(roc_obj))
}
set.seed(123)
boot_train <- boot(data_feature_train, boot_auc, R = 1000)
ci_train <- boot.ci(boot_train, type = "perc")
set.seed(123)
boot_test <- boot(data_feature_test, boot_auc, R = 1000)
ci_test <- boot.ci(boot_test, type = "perc")
boot_brier <- function(data, indices) {
  d <- data[indices, ]
  pred <- predict(rf.OBS.2, newdata = d, type = "prob")[,2]
  return(BrierScore(as.numeric(d$Depression)-1, pred, scaled = FALSE))
}
set.seed(123)
boot_brier_train <- boot(data_feature_train, boot_brier, R = 1000)
ci_brier_train <- boot.ci(boot_brier_train, type = "perc")
set.seed(123)
boot_brier_test <- boot(data_feature_test, boot_brier, R = 1000)
ci_brier_test <- boot.ci(boot_brier_test, type = "perc")
cat("Training set AUC: ", auc(trainroc), "\n")
cat("Training set AUC 95% CI: [", ci_train$percent[4], ", ", ci_train$percent[5], "]\n")
cat("Test set AUC: ", auc(testroc), "\n")
cat("Test set AUC 95% CI: [", ci_test$percent[4], ", ", ci_test$percent[5], "]\n")
cat("Training set Brier Score: ", BrierScore(y_train, prob.rf.train, scaled = FALSE), "\n")
cat("Training set Brier Score 95% CI: [", ci_brier_train$percent[4], ", ", ci_brier_train$percent[5], "]\n")
cat("Test set Brier Score: ", BrierScore(y_test, prob.rf.test, scaled = FALSE), "\n")
cat("Test set Brier Score 95% CI: [", ci_brier_test$percent[4], ", ", ci_brier_test$percent[5], "]\n")
###Decision curves analysis
data_feature_train$train_prob <- rf.OBS.trainprob[,2]
data_feature_test$test_prob <- rf.OBS.testprob[,2]
train_dca_result <- dca(Depression ~ train_prob, data = data_feature_train)
test_dca_result <- dca(Depression ~ test_prob, data = data_feature_test)
plot(train_dca_result, colorize = TRUE)
plot(test_dca_result, colorize = TRUE)

#Model interpretation by SHAP
unified_model <- randomForest.unify(rf.OBS.2, data_feature_train)
shaps <- treeshap(unified_model, data_feature_train)
shp <- shapviz(shaps)
sv_importance(shp, kind = "beeswarm")+theme_bw()
sv_dependence(shp, v = "DII",color_var = NULL)+theme_bw()

