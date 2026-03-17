# Diagnostic-Accuracy-meta-analyses
Diagnostic accuracy meta-analysis of miRNA biomarkers in periodontitis

Bivariate Reitsma with full co variance 
library(mada)

############################
# miR-146
############################
miR146_data <- data.frame(
  Study = c("Alminana_Pastor_2023", "Radovic_2018", "Rovas_2022"),
  TP = c(8, 48, 50),
  FN = c(3, 0, 26),
  FP = c(3, 0, 46),
  TN = c(8, 48, 88)
)

miR146_model <- reitsma(
  cbind(TP, FN, FP, TN) ~ 1,
  data = miR146_data,
  method = "reml"
)

summary(miR146_model)
miR146_model$vcov


############################
# miR-155
############################
miR155_data <- data.frame(
  Study = c("Appukuttan_2025", "Baru_2025", "Dailly_2023",
            "Nandipati_2022", "Radovic_2018"),
  TP = c(20, 12, 54, 46, 46),
  FN = c(4, 5, 6, 3, 10),
  FP = c(3, 9, 0, 4, 2),
  TN = c(9, 24, 60, 45, 38)
)

miR155_model <- reitsma(
  cbind(TP, FN, FP, TN) ~ 1,
  data = miR155_data,
  method = "reml"
)

summary(miR155_model)
miR155_model$vcov


############################
# miR-223 (continuity correction applied)
############################
miR223_data <- data.frame(
  Study = c("AbdelKawy_2024", "Bandi_2023a", "Bandi_2023b",
            "Bandi_2024", "Elazazy_2021", "Liu_2022"),
  TP = c(20.5, 37, 15, 40, 35, 19),
  FN = c(8.5, 13, 10, 10, 5, 29),
  FP = c(0.5, 11, 1, 11, 5, 17),
  TN = c(14.5, 39, 24, 39, 15, 32)
)

miR223_model <- reitsma(
  cbind(TP, FN, FP, TN) ~ 1,
  data = miR223_data,
  method = "reml"
)

summary(miR223_model)
miR223_model$vcov



plot(miR146_model, sroclwd = 2)
plot(miR155_model, sroclwd = 2)
plot(miR223_model, sroclwd = 2)


# Covariance matrices
vcov_miR146 <- miR146_model$vcov
vcov_miR155 <- miR155_model$vcov
vcov_miR223 <- miR223_model$vcov

COMBINE HSROC CURVE 
library(mada)

# Example for miR-146
m146 <- reitsma(cbind(TP,FN) ~ cbind(FP,TN), data=miR146, studlab=Study)

# miR-155
m155 <- reitsma(cbind(TP,FN) ~ cbind(FP,TN), data=miR155, studlab=Study)

# miR-223
m223 <- reitsma(cbind(TP,FN) ~ cbind(FP,TN), data=miR223, studlab=Study)


add_hsroc <- function(model, data, col) {
  # Add study points
  points(x = 1 - data$TN / (data$TN + data$FP),
         y = data$TP / (data$TP + data$FN),
         pch=19, col=col)
  
  # Add HSROC summary curve
  sroc_curve <- sroc(model, type="ruttergatsonis")
  lines(sroc_curve$FP, sroc_curve$Sens, col=col, lwd=2)
}




# Your miR-155 data
miR155 <- data.frame(
  Study = c("Appukuttan 2025","Baru 2025","Dailly 2023","Nandipati 2022","Radovic 2018"),
  TP = c(20,12,54,46,46),
  FN = c(4,5,6,3,10),
  FP = c(3,9,0,4,2),
  TN = c(9,24,60,45,38)
)

# Fit bivariate Reitsma model
miR155_model <- reitsma(cbind(TP,FN) ~ cbind(FP,TN), data=miR155, studlab=Study)

# Plot HSROC curve
plot(miR155_model, type="SROC", 
     main="HSROC Curve for miR-155",
     xlab="1 - Specificity", ylab="Sensitivity",
     col.points="darkgreen", cex.points=1.2)

# Add 95% confidence region
add.confidence(miR155_model)





# Example for miR-146
plot(miR146_model, type="ruttergatsonis")  # HSROC curve using Rutter & Gatsonis method
title(main="HSROC Curve for miR-146")





plot(miR155_model, type="ruttergatsonis")
title(main="HSROC Curve for miR-155")




plot(miR223_model, type="ruttergatsonis")
title(main="HSROC Curve for miR-223")







plot(NA, xlim=c(0,1), ylim=c(0,1),
     xlab="1 - Specificity", ylab="Sensitivity",
     main="HSROC Curves for miR-146, miR-155, miR-223")

add_hsroc(m146, miR146, "blue")
add_hsroc(m155, miR155, "darkgreen")
add_hsroc(m223, miR223, "red")

legend("bottomright", legend=c("miR-146","miR-155","miR-223"),
       col=c("blue","darkgreen","red"), lwd=2, pch=19)

       STUDY LEVEL ACCURACY AND PRIORI INFLUNCE EXAMPLE ONE miRNA (biomarker evaluated)

       library(dplyr)

# Input data
miR146_data <- data.frame(
  study = c("Alminana Pastor 2023", "Radovic 2018", "Rovas 2022"),
  TP = c(8, 48, 50),
  FN = c(3, 0, 26),
  FP = c(3, 0, 46),
  TN = c(8, 48, 88)
)

# Calculate study-level accuracy
miR146_table <- miR146_data %>%
  mutate(
    sensitivity = TP/(TP+FN),
    specificity = TN/(TN+FP),
    total_error = FN + FP
  )

# Apply a priori influence rules
miR146_table <- miR146_table %>%
  mutate(
    influential_flag =
      case_when(
        sensitivity == 1 | specificity == 1 ~ "Influential: Perfect accuracy",
        sensitivity < 0.50 | specificity < 0.50 ~ "Influential: Low performance",
        TRUE ~ "No major influence"
      )
  )

# Print final table
print(miR146_table)

cat("\nNote: Influence defined by a priori outlier rules, Baujat/influence plots skipped due to only 3 studies with extreme values.\n")



PPV AND NPV

#---------------------------------------------------
# Function to calculate PPV and NPV with confidence intervals
#---------------------------------------------------
calc_PPV_NPV <- function(sens, sens_lower, sens_upper,
                         spec, spec_lower, spec_upper,
                         prevalence) {
  
  # Point estimates
  PPV <- (sens * prevalence) / (sens * prevalence + (1 - spec) * (1 - prevalence))
  NPV <- (spec * (1 - prevalence)) / ((1 - sens) * prevalence + spec * (1 - prevalence))
  
  # Approximate confidence intervals using delta method
  PPV_Lower <- (sens_lower * prevalence) / (sens_lower * prevalence + (1 - spec_upper) * (1 - prevalence))
  PPV_Upper <- (sens_upper * prevalence) / (sens_upper * prevalence + (1 - spec_lower) * (1 - prevalence))
  
  NPV_Lower <- (spec_lower * (1 - prevalence)) / ((1 - sens_upper) * prevalence + spec_lower * (1 - prevalence))
  NPV_Upper <- (spec_upper * (1 - prevalence)) / ((1 - sens_lower) * prevalence + spec_upper * (1 - prevalence))
  
  # Return results as a data frame
  data.frame(
    Prevalence = prevalence,
    PPV = PPV, PPV_Lower = PPV_Lower, PPV_Upper = PPV_Upper,
    NPV = NPV, NPV_Lower = NPV_Lower, NPV_Upper = NPV_Upper
  )
}

#---------------------------------------------------
# Define prevalence scenarios
#---------------------------------------------------
prevalences <- c(0.05, 0.1, 0.2, 0.5)

#---------------------------------------------------
# miR-146
#---------------------------------------------------
sens <- 0.866; sens_lower <- 0.385; sens_upper <- 0.985
spec <- 0.865; spec_lower <- 0.015; spec_upper <- 0.614

results_miR146 <- do.call(rbind, lapply(prevalences, function(p) {
  calc_PPV_NPV(sens, sens_lower, sens_upper, spec, spec_lower, spec_upper, p)
}))
results_miR146

#---------------------------------------------------
# miR-155
#---------------------------------------------------
sens <- 0.836; sens_lower <- 0.748; sens_upper <- 0.898
spec <- 0.881; spec_lower <- 0.738; spec_upper <- 0.948

results_miR155 <- do.call(rbind, lapply(prevalences, function(p) {
  calc_PPV_NPV(sens, sens_lower, sens_upper, spec, spec_lower, spec_upper, p)
}))
results_miR155

#---------------------------------------------------
# miR-223
#---------------------------------------------------
sens <- 0.699; sens_lower <- 0.548; sens_upper <- 0.817
spec <- 0.773; spec_lower <- 0.626; spec_upper <- 0.876

results_miR223 <- do.call(rbind, lapply(prevalences, function(p) {
  calc_PPV_NPV(sens, sens_lower, sens_upper, spec, spec_lower, spec_upper, p)
}))
results_miR223



Specify by specimen and assay platform 
library(mada)

############################
# miR-146
############################
miR146_data <- data.frame(
  Study = c("Alminana_Pastor_2023","Radovic_2018","Rovas_2022"),
  Specimen = c("GCF","GCF","Gingival_tissue"),
  Platform = c("RT-qPCR","RT-qPCR","RT-qPCR"),
  TP = c(8,48,50),
  FN = c(3,0,26),
  FP = c(3,0,46),
  TN = c(8,48,88)
)

# Stratify by Specimen
miR146_gcf <- subset(miR146_data, Specimen == "GCF")
miR146_tissue <- subset(miR146_data, Specimen == "Gingival_tissue")

# Stratify by Platform (all RT-qPCR, but included for completeness)
miR146_rtpcr <- subset(miR146_data, Platform == "RT-qPCR")

# Fit Reitsma models
miR146_gcf_model <- reitsma(cbind(TP,FN,FP,TN) ~ 1, data=miR146_gcf, method="reml")
miR146_tissue_model <- reitsma(cbind(TP,FN,FP,TN) ~ 1, data=miR146_tissue, method="reml")
miR146_rtpcr_model <- reitsma(cbind(TP,FN,FP,TN) ~ 1, data=miR146_rtpcr, method="reml")

# Summaries
summary(miR146_gcf_model)
summary(miR146_tissue_model)
summary(miR146_rtpcr_model)

############################
# miR-155
############################
miR155_data <- data.frame(
  Study = c("Appukuttan_2025","Baru_2025","Dailly_2023","Nandipati_2022","Radovic_2018"),
  Specimen = c("Gingival_tissue","Gingival_tissue","Blood","Gingival_tissue","GCF"),
  Platform = c("RT-qPCR","RT-qPCR","RT-qPCR","RT-qPCR","RT-qPCR"),
  TP = c(20,12,54,46,46),
  FN = c(4,5,6,3,10),
  FP = c(3,9,0,4,2),
  TN = c(9,24,60,45,38)
)

# Stratify by Specimen
miR155_gcf <- subset(miR155_data, Specimen == "GCF")
miR155_tissue <- subset(miR155_data, Specimen == "Gingival_tissue")
miR155_blood <- subset(miR155_data, Specimen == "Blood")

# Stratify by Platform
miR155_rtpcr <- subset(miR155_data, Platform == "RT-qPCR")

# Fit Reitsma models
miR155_gcf_model <- reitsma(cbind(TP,FN,FP,TN) ~ 1, data=miR155_gcf, method="reml")
miR155_tissue_model <- reitsma(cbind(TP,FN,FP,TN) ~ 1, data=miR155_tissue, method="reml")
miR155_blood_model <- reitsma(cbind(TP,FN,FP,TN) ~ 1, data=miR155_blood, method="reml")
miR155_rtpcr_model <- reitsma(cbind(TP,FN,FP,TN) ~ 1, data=miR155_rtpcr, method="reml")

# Summaries
summary(miR155_gcf_model)
summary(miR155_tissue_model)
summary(miR155_blood_model)
summary(miR155_rtpcr_model)

############################
# miR-223
############################
miR223_data <- data.frame(
  Study = c("AbdelKawy_2024","Bandi_2023a","Bandi_2023b","Bandi_2024","Elazazy_2021","Liu_2022"),
  Specimen = c("GCF","Saliva","Saliva","GCF","GCF","GCF"),
  Platform = c("RT-qPCR","RT-qPCR","RT-qPCR","RT-qPCR","RT-qPCR","RT-qPCR"),
  TP = c(20.5,37,15,40,35,19),
  FN = c(8.5,13,10,10,5,29),
  FP = c(0.5,11,1,11,5,17),
  TN = c(14.5,39,24,39,15,32)
)

# Stratify by Specimen
miR223_gcf <- subset(miR223_data, Specimen == "GCF")
miR223_saliva <- subset(miR223_data, Specimen == "Saliva")

# Stratify by Platform
miR223_rtpcr <- subset(miR223_data, Platform == "RT-qPCR")

# Fit Reitsma models
miR223_gcf_model <- reitsma(cbind(TP,FN,FP,TN) ~ 1, data=miR223_gcf, method="reml")
miR223_saliva_model <- reitsma(cbind(TP,FN,FP,TN) ~ 1, data=miR223_saliva, method="reml")
miR223_rtpcr_model <- reitsma(cbind(TP,FN,FP,TN) ~ 1, data=miR223_rtpcr, method="reml")

# Summaries
summary(miR223_gcf_model)
summary(miR223_saliva_model)
summary(miR223_rtpcr_model)
