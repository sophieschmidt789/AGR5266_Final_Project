##############################
### Neo Compilation Script ###
##############################
library(reshape)
library(reshape2)
library(lme4)
library(zoo)
library(data.table)
library(dplyr)
library(SpATS)
library(BGLR)
library(randomForest)
library(rrBLUP)
library(sommer)
library(tidyr)
library(ggplot2)

setwd("C:/Users/s.schmidt1/OneDrive - University of Florida/Fall 2025/AGR5266c_Field_Plot_Tech/Final Project")
getwd()
df = read.csv("C:/Users/s.schmidt1/OneDrive - University of Florida/Fall 2025/AGR5266c_Field_Plot_Tech/Final Project/2024-25NeoStage2Data_CSV.csv", header = T, as.is = T, fileEncoding ="latin1")
#remove pre-calculated AUDPC column
df$AUDPC <- NULL
unique <- unique(df$Genotype)
##############################################Spatial Adjustment#############################################################################
#column names
all_cols <- colnames(df)
all_cols

#grab just the date columns
date_cols <- grep("^X\\d{1,2}\\.\\d{1,2}\\.\\d{4}$", colnames(df), value = TRUE)
if(length(date_cols) == 0) stop("No date columns detected. Make sure column names look like X1.8.2025 etc.")

#melt data for long format
long <- reshape2::melt(df,
             id.vars = c("ï..SampleNo.","Pos","Rep","Bed","Genotype"),
             measure.vars = date_cols,
             variable.name = "DateString",
             value.name = "Value")


#convert date to correct time
long$DateString <- gsub("^X", "", long$DateString)
long$Date <- as.Date(long$DateString, format = "%m.%d.%Y")
print(unique(long$Date))

#calculate days from first obs for AUDPC
long <- long %>% arrange(Genotype, Date)
start_date <- min(long$Date, na.rm = TRUE)
long$Days <- as.numeric(long$Date - start_date)
str(long)

#spatial coordinats for SPaTS
### we assume that the Bed = column (11 beds per rep)
### we assume that the Pos = row ( beds per rep)
long$Column <- as.numeric(long$Bed)
long$Row <- as.numeric(long$Pos)
long$R <- as.factor(long$Row)
long$C <- as.factor(long$Column)

#run spatial model assuming bed as fixed and row and column as random. 

model.F = SpATS(data = long, response = "Value", genotype = "Genotype", genotype.as.random = FALSE,
                spatial = ~ SAP(Column, Row, nseg = c(10,20), degree = 3, pord = 2), fixed = ~ Bed, random = ~ R + C,
                family = gaussian(), offset = 0, control = list(tolerance = 1e-03))


long$BLUE = model.F$fitted;
temp = aggregate(data = long, BLUE ~ Genotype, FUN = mean);

###Check Residuals
hist(model.F$residuals, main = "Residuals", xlab = "Residual value")
### Residuals v fittted
plot(model.F$fitted, model.F$residuals,
     xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, col = "red")
### Examine Blues - looks good, not all the same
blues <- fitted(model.F, which = "genotype")

out.fit = long[,c("Date", "Rep", "Genotype", "BLUE")]; 
total.fit = c();
total.fit = rbind(total.fit, out.fit);

#Save csv that is organized as Date, Rep, Genotype, and Adjusted BlUE
write.csv(total.fit, "NeoAUDPC_SpATS.csv", row.names = F, quote = F);

###################################Calculate AUDPC on spatially adjusted values#####################################
#prep data frame total.fit
str(total.fit)
total.fit$Rep <- as.factor(total.fit$Rep)
total.fit$Genotype <- as.factor(total.fit$Genotype)
total.fit$Date <- as.factor(total.fit$Date)

#Average BLUE across reps within each date
blues_avg <- total.fit %>%
  group_by(Genotype, Date) %>%
  summarize(MeanBLUE = mean(BLUE, na.rm = TRUE)) %>%
  ungroup()

#Wide format
blues_wide <- blues_avg %>%
  pivot_wider(names_from = Date, values_from = MeanBLUE)


#Caluclate days between observations
date_cols <- sort(unique(blues_avg$Date))
days <- unique(long$Days)

#Function for calculating AUDPC
calc_audpc <- function(values, days) {
  sum((head(values, -1) + tail(values, -1)) / 2 * diff(days), na.rm = TRUE)
}

#Calculate AUDPC using
audpc_results <- blues_wide %>%
  rowwise() %>%
  mutate(AUDPC = calc_audpc(c_across(all_of(as.character(date_cols))), days)) %>%
  ungroup()

write.csv(audpc_results, "NeoAUDPC_correct.csv", row.names = F, quote = F);

# Result: one row per genotype with mean AUDPC and SE
head(audpc_summary)

#Check to include standard error
blues_wide1 <- total.fit %>%
  pivot_wider(names_from = Date, values_from = BLUE)
audpc_results1 <- blues_wide1 %>%
  rowwise() %>%
  mutate(AUDPC = calc_audpc(c_across(all_of(as.character(date_cols))), days)) %>%
  ungroup()
audpc_summary1 <- audpc_results1 %>%
  group_by(Genotype) %>%
  summarize(
    AUDPC_mean = mean(AUDPC, na.rm = TRUE),
    AUDPC_se   = sd(AUDPC, na.rm = TRUE)/sqrt(n()),
    .groups = "drop"
  )

#merge with df with standard error and compare
audpc_se <- merge(audpc_results, audpc_summary1, "Genotype")
write.csv(audpc_se, "NeoAUDPC_correct-SE_average.csv", row.names = F, quote = F);

library(ggplot2)
library(dplyr)

# 1️⃣ Histogram of AUDPC
ggplot(audpc_se, aes(x = AUDPC_mean)) +
  geom_histogram(binwidth = 20, fill = "steelblue", color = "black") +
  theme_minimal() +
  xlab("AUDPC") +
  ylab("Number of Genotypes") +
  ggtitle("Distribution of AUDPC Across Genotypes")

top20 <- audpc_se %>% arrange(AUDPC_mean) %>% slice(1:20)

ggplot(top20, aes(x = reorder(Genotype, AUDPC_mean), y = AUDPC_mean)) +
  geom_point(size = 3, color = "darkgreen") +
  geom_errorbar(aes(ymin = AUDPC_mean - AUDPC_se, ymax = AUDPC_mean + AUDPC_se), width = 0.3) +
  coord_flip() +
  theme_minimal() +
  xlab("Genotype") +
  ylab("AUDPC ± SE") +
  ggtitle("Top 20 Genotypes by AUDPC")

top10 <- audpc_se %>% arrange(AUDPC_mean) %>% slice(1:10)
bottom10 <- audpc_se %>% arrange(desc(AUDPC_mean)) %>% slice(1:10)

# Combine into one dataframe
top_bottom <- bind_rows(top10, bottom10) %>%
  mutate(Group = ifelse(Genotype %in% top10$Genotype, "Top 10", "Bottom 10"))

# Plot
ggplot(top_bottom, aes(x = reorder(Genotype, AUDPC_mean), y = AUDPC_mean, color = Group)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = AUDPC_mean - AUDPC_se, ymax = AUDPC_mean + AUDPC_se), width = 0.3) +
  coord_flip() +
  scale_color_manual(values = c("Top 10" = "darkgreen", "Bottom 10" = "red")) +
  theme_minimal() +
  xlab("Genotype") +
  ylab("AUDPC ± SE") +
  ggtitle("Top 10 and Bottom 10 Genotypes by AUDPC") +
  theme(legend.title = element_blank())

###########################################################################rrBLUP to merge genotype data with phenotype AUDPC#################################################
#SNP Data:
geno = fread("C:/Users/s.schmidt1/OneDrive - University of Florida/Fall 2025/AGR5266c_Field_Plot_Tech/Final Project/GenoDatabase_Fana_iStraw_AgriSeq_Imputed_4.8.25.csv", header = T, sep = ",");
#Map:
map = read.table("C:/Users/s.schmidt1/OneDrive - University of Florida/Fall 2025/AGR5266c_Field_Plot_Tech/Final Project/map_Zhen.txt", header = T, sep = "\t");
#AUDPC Data:
pheno = fread("C:/Users/s.schmidt1/OneDrive - University of Florida/Fall 2025/AGR5266c_Field_Plot_Tech/Final Project/NeoAUDPC_correct-SE_average.csv", header = T, sep = ",");

#Pull out matching genotypes from marker data and my audpc data
common_genos <- intersect(pheno$Genotype, geno$Genotype) 
common_genos #590 common genos - makes sense as these are UFL varieties not UC Davis Varieties
common_pheno <- pheno %>% filter(Genotype %in% common_genos) #590 genotypes
common_marker_geno <- geno %>% filter(Genotype %in% common_genos) #611 must have some duplicates
sum(duplicated(common_marker_geno$Genotype)) #21 duplicates
common_marker_geno <- common_marker_geno[!duplicated(common_marker_geno$Genotype), ] #remove duplicates and keep 1st occurance
common_marker_geno #now 590 obs


#Check if all the markers are in my map file - not sure about this
geno1stcolumn <- common_marker_geno[, -1] #get list of marker names in collumns
all(colnames(geno1stcolumn) %in% map$Probe_ID) #True so yes they are

#make Phenotype Matrix for response vector
y <- common_pheno$AUDPC_mean
names(y) <- common_pheno$Genotype #assign the audpc values to genotype
y #looks good

#make Marker matrix
M <- as.matrix(geno1stcolumn) #don't need genotype column at first
rownames(M) <- common_marker_geno$Genotype #add genotypes back in as rownames
str(M)
dim(M) #looks good
sum_na <- sum(is.na(M)) #no NAs
M <- data.matrix(M) #is numeric

#confirm everything is ordered correctly
M <- M[common_pheno$Genotype, ]
stopifnot(identical(rownames(M), common_pheno$Genotype))

# RR-BLUP: estimate marker effects and predict GEBVs
rr_model <- mixed.solve(y = y, Z = M)

#Results
blups <- rr_model$u #blups = predicted genetic values
mean_effect <- rr_model$beta #mean effects for overall mean
###Variance components
genetic_var <- rr_model$Vu
residual_var <- rr_model$Ve
heritability <- genetic_var / (genetic_var + residual_var)
cat("Heritability (H2):", round(heritability, 3), "\n")

# Combine into a results table
GEBV <- M %*% rr_model$u
blup_results <- data.frame(
  Genotype = rownames(M),
  GEBV = as.numeric(GEBV)
)

#Kinship method
K <- A.mat(M)  # make kinship matrix
rr_modelK <- mixed.solve(y = y, K = K)
blup_resultsK <- data.frame(
  Genotype = names(y),
  GEBV = rr_modelK$u
)

#Merge results
results <- merge(blup_results, pheno, by = "Genotype")
resultsK <- merge(blup_resultsK, pheno, by = "Genotype")

#Check peformance of the model
cor(results$GEBV, results$AUDPC_mean, use = "complete.obs") #0.7689853 High correlation means markers explain lots of phenotypic variation

#Plot
plot(results$AUDPC_mean, results$GEBV,
     xlab = "Observed BLUE AUDPC",
     ylab = "Predicted GEBV",
     main = "Genomic Prediction Accuracy")
abline(0, 1, col = "red")

#Check peformance of the kinship model
cor(results$GEBV, resultsK$AUDPC_mean, use = "complete.obs") #0.7689853 High correlation means markers explain lots of phenotypic variation

#Plot - basically same graph y-value is just on smaller scale
plot(resultsK$AUDPC_mean, resultsK$GEBV,
     xlab = "Observed BLUE AUDPC",
     ylab = "Predicted GEBV",
     main = "Genomic Prediction Accuracy")
abline(0, 0.7689853, col = "red")

#Nice Plot
acc <- cor(results$GEBV, results$AUDPC_mean, use = "complete.obs")
plot(results$AUDPC_mean, results$GEBV,
    xlab = "Observed BLUE AUDPC",
    ylab = "Predicted GEBV",
    main = paste("Genomic Prediction Accuracy (r =", round(acc, 2), ")"),
    pch = 19, col = "blue")
abline(0, 1, col = "red", lwd = 2)

#Rank genotypes: we want genotypes with lowest AUDPC = more resistance
results <- results[order(results$GEBV), ] #order from lowest GEBV to smallest
results$Rank <- rank(results$GEBV, ties.method = "first") #rank values
results

#visualize: top 20 disese resistant genotypes

top20 <- head(results, 20)
ggplot(top20, aes(x = reorder(Genotype, GEBV), y = GEBV)) +
  geom_bar(stat = "identity", fill = "darkgreen") +
  coord_flip() +
  labs(
    title = "Top 20 Most Disease-Resistant Genotypes (Lowest GEBV)",
    x = "Genotype",
    y = "Predicted GEBV (AUDPC)"
  ) +
  theme_minimal(base_size = 14)

#############################################Statistical method##################################################
results$upper <- results$AUDPC_mean + 1.96*results$AUDPC_se
results$lower <- results$AUDPC_mean - 1.96*results$AUDPC_se

ggplot(results, aes(x = AUDPC_mean, y = GEBV)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", color = "red") +
  xlab("Observed AUDPC") +
  ylab("Predicted Resistance (flipped GEBV)") +
  ggtitle("Prediction of AUDPC using SNP markers")
library(ggplot2)
library(ggpubr)  # for easy stat_cor

# Compute correlation
cor_val <- cor(results$GEBV, results$AUDPC_mean, use = "complete.obs")
cor_text <- paste0("r = ", round(cor_val, 2))

# Plot with correlation annotation
ggplot(results, aes(x = AUDPC_mean, y = GEBV)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", color = "red") +
  annotate("text", x = max(results$AUDPC_mean, na.rm = TRUE)*0.7,
           y = max(-results$GEBV, na.rm = TRUE)*0.9,
           label = cor_text, size = 5, color = "black") +
  xlab("Observed AUDPC") +
  ylab("Predicted Resistance (flipped GEBV)") +
  ggtitle("Prediction of AUDPC using SNP markers")

