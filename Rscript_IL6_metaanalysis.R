##########################################################
# Script to perform meta-analysis on plasma IL-6 levels
# Gustav Nilsonne and Michael Ingre 2014
##########################################################

###########################
# READ AND PROCESS DATA
###########################

# Load packages
require(stats)
require(nlme)
require(plotrix)
require(foreign)
require(RColorBrewer)

# Define functions
fun_center <- function(x){ #Centers density plots with max at 0
  x$x <- x$x - x$x[x$y == max(x$y)]
  return(x)
}

# Read main data file
IL6_data <- read.csv2("IL6_data.csv", header = TRUE)
IL6_data$lnIL6 <- NA

# Read individual data and process into common format
Sothern1995 <- read.csv2("Sothern1995IndividualData.csv", header = FALSE)
Sothern1995ind <- data.frame("IL6_pgml.1" = as.vector(t(Sothern1995[3:10])))
Sothern1995ind[Sothern1995ind < 1] <- 1 # Values below detection limit set to detection limit
Sothern1995ind$Hour <- c(19, 22, 25, 28, 31, 34, 37, 40)
Sothern1995ind$Subject <- c(rep(1, 8), rep(2,8), rep(3, 8), rep(4, 8), rep(5, 8), rep(6, 8), rep(7, 8), rep(8, 8), rep(9, 8), rep(10, 8), rep(11, 8))
Sothern1995ind$Study = "Sothern1995"
Sothern1995ind$TimeFromCatheter_h <- c(0, 3, 6, 9, 12, 15, 18, 21)
Sothern1995ind$n <- 11
Sothern1995ind$TimeAsleep <- c(0, 0, 2, 5, 8, 0, 0, 0)
Sothern1995ind$TimeAwake <- c(12, 15, 0, 0, 0, 3, 6, 9)
Sothern1995ind$Day <- 1

Sothern1995Data <- data.frame(1:8)[, F]
Sothern1995Data$Study <- "Sothern1995"
#Sothern1995Data$Condition <- "Baseline"
#Sothern1995Data$PSD_Sleepduration <- NA
Sothern1995Data$Day <- c(1, 1, 2, 2, 2, 2, 2, 2)
Sothern1995Data$Hour <- c(19, 22, 25, 28, 31, 34, 37, 40)
Sothern1995Data$IL6_pgml.1 <- colMeans(Sothern1995)[3:10]
Sothern1995Data$lnIL6 <- log(colMeans(Sothern1995))[3:10]
Sothern1995Data$IL6_errorvalue <- NA
Sothern1995Data$IL6_sem <- (apply(Sothern1995, 2, sd)[3:10])/sqrt(11)
Sothern1995Data$n <- 11
Sothern1995Data$IL6_sd <- Sothern1995Data$IL6_sem * sqrt(Sothern1995Data$n)
Sothern1995Data$TimeFromCatheter_h <- c(0, 3, 6, 9, 12, 15, 18, 21)
Sothern1995Data$Sleep <- c(0, 0, 1, 1, 1, 0, 0, 0)
Sothern1995Data$TimeAsleep <- c(0, 0, 2, 5, 8, 0, 0, 0)
Sothern1995Data$PublicationYear <- 1995
Sothern1995Data$TimeAwake <- c(12, 15, 0, 0, 0, 3, 6, 9)
Sothern1995Data$Study <- "Sothern1995"

Knudsen2008ind <- read.csv2("Knudsen2008individual.csv", header = TRUE)
Knudsen2008ind$IL6_pgml.1[Knudsen2008ind$IL6_pgml.1 < 0.156] <- 0.156 # Values below detection limit set to detection limit - there were none
Knudsen2008ind$Study <- "Knudsen2008"
Knudsen2008ind$n <- 15
Knudsen2008ind$Day <- 1

Knudsen2008Data <- aggregate(Knudsen2008ind$IL6_pgml.1, list(Hour = Knudsen2008ind$Hour), mean)
names(Knudsen2008Data) <- c("Hour", "IL6_pgml.1")
Knudsen2008Data$lnIL6 <- aggregate(log(Knudsen2008ind$IL6_pgml.1), list(Hour = Knudsen2008ind$Hour), mean)$x
Knudsen2008Data$Study <- "Knudsen2008"
#Knudsen2008Data$Condition <- "Baseline"
#Knudsen2008Data$PSD_Sleepduration <- NA
Knudsen2008Data$Day <- c(1, 1, 1, 1, 1, 1, 2, 2)
Knudsen2008Data$IL6_errorvalue <- NA
Knudsen2008Data$n <- aggregate(Knudsen2008ind$IL6_pgml.1, list(Hour = Knudsen2008ind$Hour), length)$x
Knudsen2008Data$IL6_sd <- aggregate(Knudsen2008ind$IL6_pgml.1, list(Hour = Knudsen2008ind$Hour), sd)$x
Knudsen2008Data$IL6_sem <- Knudsen2008Data$IL6_sd/sqrt(Knudsen2008Data$n)
Knudsen2008Data$TimeFromCatheter_h <- c(0, 2, 3, 6, 9, 12, 21, 24)
Knudsen2008Data$Sleep <- c(0, 0, 0, 0, 0, 0, 0, 0)
Knudsen2008Data$TimeAsleep <- c(0, 0, 0, 0, 0, 0, 8, 0)
Knudsen2008Data$PublicationYear <- 2008
Knudsen2008Data$TimeAwake <- c(3, 5, 6, 9, 12, 15, 0, 3)

Lekander2013 <- read.dta("psd_srh_il6.dta")
Lekander2013ind <- data.frame("Subject" = Lekander2013$id)
Lekander2013ind$IL6_pgml.1 <- Lekander2013$il6m
Lekander2013ind$Hour <- Lekander2013$time
Lekander2013ind$Day <- Lekander2013$day
Lekander2013ind$Hour[Lekander2013ind$Day == 2] <- Lekander2013ind$Hour[Lekander2013ind$Day == 2] + 24
Lekander2013ind$Hour[Lekander2013ind$Day == 3] <- Lekander2013ind$Hour[Lekander2013ind$Day == 3] + 48
Lekander2013ind$Hour[Lekander2013ind$Day == 4] <- Lekander2013ind$Hour[Lekander2013ind$Day == 4] + 72
Lekander2013ind$Hour[Lekander2013ind$Day == 5] <- Lekander2013ind$Hour[Lekander2013ind$Day == 5] + 96
Lekander2013ind$Hour[Lekander2013ind$Day == 8] <- Lekander2013ind$Hour[Lekander2013ind$Day == 8] + 168
Lekander2013ind$Hour[Lekander2013ind$Day == 9] <- Lekander2013ind$Hour[Lekander2013ind$Day == 9] + 192
Lekander2013ind$Hour[Lekander2013ind$Day == 10] <- Lekander2013ind$Hour[Lekander2013ind$Day == 10] + 216
Lekander2013ind$Hour[Lekander2013ind$Day == 11] <- Lekander2013ind$Hour[Lekander2013ind$Day == 11] + 240
Lekander2013ind$TimeFromCatheter_h <- Lekander2013$time +1
Lekander2013ind$TimeAsleep <- NA ## TODO!
Lekander2013ind$TimeAwake <- NA ## TODO!
Lekander2013ind$n <- 9
Lekander2013ind$Study[Lekander2013ind$Day %in% c(2, 3, 9, 10, 11)] <- "Lekander2013full"
Lekander2013ind$Study[Lekander2013ind$Day %in% c(4, 5, 8)] <- "Lekander2013PSD"
Lekander2013ind <- Lekander2013ind[complete.cases(Lekander2013ind$IL6_pgml.1), ] # Remove NA values
Lekander2013ind$IL6_pgml.1[Lekander2013ind$IL6_pgml.1 < 0.039] <- 0.039 # Values below detection limit set to detection limit - there were none

Lekander2013ind_full <- Lekander2013ind[Lekander2013ind$Study == "Lekander2013full", ]
Lekander2013Data <- aggregate(Lekander2013ind_full$IL6_pgml.1, list(Hour = Lekander2013ind_full$Hour), mean)
names(Lekander2013Data) <- c("Hour", "IL6_pgml.1")
Lekander2013Data$lnIL6 <- aggregate(log(Lekander2013ind_full$IL6_pgml.1), list(Hour = Lekander2013ind_full$Hour), mean)$x
Lekander2013Data$Study <- "Lekander2013"
#Lekander2013Data$Condition <- "Baseline"
#Lekander2013Data$PSD_Sleepduration <- NA
Lekander2013Data$Day <- c(1, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 9, 9, 9, 9, 9, 9, 9, 9, 10, 10, 10, 10, 10, 10, 10, 10, 11, 11, 11, 11, 11, 11, 11, 11)
Lekander2013Data$IL6_errorvalue <- NA
Lekander2013Data$n <- aggregate(Lekander2013ind_full$IL6_pgml.1, list(Hour = Lekander2013ind_full$Hour), length)$x
Lekander2013Data$IL6_sd <- aggregate(Lekander2013ind_full$IL6_pgml.1, list(Hour = Lekander2013ind_full$Hour), sd)$x
Lekander2013Data$IL6_sem <- Lekander2013Data$IL6_sd/sqrt(Lekander2013Data$n)
Lekander2013Data$TimeFromCatheter_h <- seq(0, 21, 3)
Lekander2013Data$Sleep <- c(0, 1, 1, 1, 0, 0, 0, 0)
Lekander2013Data$TimeAsleep <- c(0, 0, 3, 6, 0, 0, 0, 0)
Lekander2013Data$PublicationYear <- 2013
Lekander2013Data$TimeAwake <- c(16, 0, 0, 1, 4, 7, 10, 13)

Axelsson2013ind <- read.csv2("Endotoxin_cytokine_data_140530.csv")
Axelsson2013ind$IL6_pgml.1[Axelsson2013ind$IL6_pgml.1 < 0.9] <- 0.9 # Values below detection limit set to detection limit

Axelsson2013Data <- aggregate(Axelsson2013ind$IL6_pgml.1, list(Hour = Axelsson2013ind$Hour), mean)
names(Axelsson2013Data) <- c("Hour", "IL6_pgml.1")
Axelsson2013Data$lnIL6 <- aggregate(log(Axelsson2013ind$IL6_pgml.1), list(Hour = Axelsson2013ind$Hour), mean)$x
Axelsson2013Data$Study <- "Karshikoff2015"
#Axelsson2013Data$Condition <- "Baseline"
#Axelsson2013Data$PSD_Sleepduration <- NA
Axelsson2013Data$Day <- c(1, 1, 1, 1)
Axelsson2013Data$IL6_errorvalue <- NA
Axelsson2013Data$n <- aggregate(Axelsson2013ind$IL6_pgml.1, list(Hour = Axelsson2013ind$Hour), length)$x
Axelsson2013Data$IL6_sd <- aggregate(Axelsson2013ind$IL6_pgml.1, list(Hour = Axelsson2013ind$Hour), sd)$x
Axelsson2013Data$IL6_sem <- Axelsson2013Data$IL6_sd/sqrt(Axelsson2013Data$n)
Axelsson2013Data$TimeFromCatheter_h <- c(1, 1.5, 3.5, 5)
Axelsson2013Data$Sleep <- c(0, 0, 0, 0)
Axelsson2013Data$TimeAsleep <- NA
Axelsson2013Data$PublicationYear <- 2013
Axelsson2013Data$TimeAwake <- c(3, 4.5, 6.5, 8.5)

IL6ind <- rbind(Sothern1995ind, Knudsen2008ind, Lekander2013ind, Axelsson2013ind)
IL6ind <- IL6ind[IL6ind$Study != "Lekander2013PSD", ]

IL6_data <- rbind(IL6_data, Sothern1995Data, Knudsen2008Data, Lekander2013Data, Axelsson2013Data)
IL6_data$Study[IL6_data$Study == "Axelsson2013"] <- "Karshikoff2015"

###########################
# ANALYSE GROUP DATA
###########################

# Exclude Born1997 on account of being an outlier
IL6_data <- IL6_data[IL6_data$Study != "Born1997", ]

# Construct weights
IL6_data$ivw <- NA # Inverse variance weight
IL6_data$nw <- NA # Weight by number of participants
IL6_data$sqw <- NA # Weight by number of participants and square root of number of observations
for (i in unique(IL6_data$Study)){
  mean_se <- mean(IL6_data[IL6_data$Study == i, ]$IL6_sem)
  study_ivw <- 1/(mean_se^2)
  IL6_data[IL6_data$Study == i, ]$ivw <- study_ivw/length(IL6_data[IL6_data$Study == i, ]$ivw)
  IL6_data[IL6_data$Study == i, ]$nw <- IL6_data[IL6_data$Study == i, ]$n/length(IL6_data[IL6_data$Study == i, ]$n)
  IL6_data[IL6_data$Study == i, ]$sqw <- (IL6_data[IL6_data$Study == i, ]$n * sqrt(length(IL6_data$n[IL6_data$Study == i])))/length(IL6_data[IL6_data$Study == i, ]$n)
}
IL6_data$sqw[IL6_data$Study == "Vgontzas2000"] <- IL6_data$sqw[IL6_data$Study == "Vgontzas2000"] * sqrt(3) # Increase weight of Vgontzas2000 because their data were anaverage of 3 nights 

# Make variates for modelling
IL6_data$lnIL6[IL6_data$Study == "Undar1999"] <- IL6_data$IL6_pgml.1[IL6_data$Study == "Undar1999"]
IL6_data$lnIL6[IL6_data$Study == "Agorastos2014"] <- IL6_data$IL6_pgml.1[IL6_data$Study == "Agorastos2014"]
IL6_data$lnIL6[IL6_data$Study == "Dugue1998"] <- IL6_data$IL6_pgml.1[IL6_data$Study == "Dugue1998"]
IL6_data$lnIL6[is.na(IL6_data$lnIL6)] <- log(IL6_data$IL6_pgml.1[is.na(IL6_data$lnIL6)]) - 0.5*log(((IL6_data$IL6_sd[is.na(IL6_data$lnIL6)])^2/(IL6_data$IL6_pgml.1[is.na(IL6_data$lnIL6)])^2)+1) # Following Higgins et al. 2008
IL6_data$SinHour <- sin(2*pi*IL6_data$Hour/24)
IL6_data$CosHour <- cos(2*pi*IL6_data$Hour/24)
IL6_data$SinHour12 <- sin(2*pi*IL6_data$Hour/12)
IL6_data$CosHour12 <- cos(2*pi*IL6_data$Hour/12)
IL6_data$SinHour6 <- sin(2*pi*IL6_data$Hour/6)
IL6_data$CosHour6 <- cos(2*pi*IL6_data$Hour/6)
IL6_data$SinHour3 <- sin(2*pi*IL6_data$Hour/3)
IL6_data$CosHour3 <- cos(2*pi*IL6_data$Hour/3)
IL6_data$SinHour1.5 <- sin(2*pi*IL6_data$Hour/1.5)
IL6_data$CosHour1.5 <- cos(2*pi*IL6_data$Hour/1.5)
IL6_data$IL6_variance <- IL6_data$IL6_sd^2
IL6_data$Hour_sq <- IL6_data$Hour^2
IL6_data$LogTimeFromCatheter_h <- log(1 + IL6_data$TimeFromCatheter_h)
IL6_data$Hour24 <- IL6_data$Hour
#IL6_data$Deprived <- IL6_data$Condition != "Baseline"
while (max(IL6_data$Hour24) > 24){
  IL6_data$Hour24[IL6_data$Hour24 > 24] <- IL6_data$Hour24[IL6_data$Hour24 > 24] - 24
}
IL6_data$TimeAsleep[is.na(IL6_data$TimeAsleep)] <- 0
IL6_data$TimeAwake[is.na(IL6_data$TimeAwake)] <- 0

# Build models
ModelNull <- lme(lnIL6 ~ TimeFromCatheter_h, data = IL6_data, weights = ~ 1/sqw, random = ~ TimeFromCatheter_h|Study, na.action = na.omit)
summary(ModelNull)
plot(ModelNull, resid(., type = "p") ~ Hour24, abline = 0)
residualsNull <- ModelNull$residuals

# Plot residuals from null model
pdf("Fig_residuals_null.pdf")
plot(residualsNull[, 2] ~ IL6_data$Hour24, frame.plot = F, xlab = "Time (h)", ylab = "ln IL-6, pg/ml, residual", type = "n", xaxt = "n", yaxt = "n", cex.lab = 1.5)
rect(-2, -5, 7, 40, col = "Misty Rose", border = NA)
rect(22, -5, 31, 40, col = "Misty Rose", border = NA)
axis(side = 1, at = c(0, 6, 12, 18, 24), labels = c("00:00", "06:00", "12:00", "18:00", "24:00"), cex.axis = 1.5)
axis(2, at = c(-3, -2, -1, 0, 1, 2, 3), cex.axis = 1.5)
abline(h = 0, col = "blue")
points(residualsNull[, 2] ~ IL6_data$Hour24, cex = 2*sqrt(IL6_data$sqw/pi))

dataforloess <- data.frame(residuals = residualsNull[, 2], Hour24 = IL6_data$Hour24, weights = IL6_data$sqw)
dataforloess2 <- rbind(dataforloess, dataforloess, dataforloess)
dataforloess2$Hour24[1:789] <- dataforloess2$Hour24[1:789] - 24
dataforloess2$Hour24[1779:2367] <- dataforloess2$Hour24[1779:2367] + 24
lo <- loess(residuals ~ Hour24, data = dataforloess2, weights = weights, span = 0.2)
pred <- predict(lo, seq(0, 24, 0.5), se=TRUE)
lines(pred$fit~ seq(0, 24, 0.5), col = "red", lwd = 2)

plot(residualsNull[, 2] ~ IL6_data$Hour24, frame.plot = F, xlab = "Time (h)", ylab = "ln IL-6, pg/ml, residual", type = "n", ylim = c(-0.5, 0.5), xaxt = "n", yaxt = "n", cex.lab = 1.5)
rect(-2, -5, 7, 40, col = "Misty Rose", border = NA)
rect(22, -5, 31, 40, col = "Misty Rose", border = NA)
axis(side = 1, at = c(0, 6, 12, 18, 24), labels = c("00:00", "06:00", "12:00", "18:00", "24:00"), cex.axis = 1.5)
axis(2, at = c(-0.5, 0, 0.5), cex.axis = 1.5)
abline(h = 0, col = "blue")
points(residualsNull[, 2] ~ IL6_data$Hour24, cex = 2*sqrt(IL6_data$sqw/pi))
lines(pred$fit~ seq(0, 24, 0.5), col = "red", lwd = 2)
dev.off()

Model24 <- lme(lnIL6 ~ SinHour + CosHour + TimeFromCatheter_h, data = IL6_data, weights = ~ 1/sqw, random = ~ TimeFromCatheter_h|Study)
summary(Model24)
plot(Model24, resid(., type = "p") ~ Hour, abline = 0)
plot(Model24, resid(., type = "p") ~ fitted(.), abline = 0)

Model12 <- lme(lnIL6 ~ SinHour + CosHour + SinHour12 + CosHour12 + TimeFromCatheter_h, data = IL6_data, weights = ~ 1/sqw, random = ~ TimeFromCatheter_h|Study, na.action = na.omit)
summary(Model12)

Model6 <- lme(lnIL6 ~ SinHour + CosHour + SinHour12 + CosHour12  + SinHour6 + CosHour6 + TimeFromCatheter_h, data = IL6_data, weights = ~ 1/sqw, random = ~ TimeFromCatheter_h|Study, na.action = na.omit)
summary(Model6)

# Test significance
#p_24vsNull <- 1 - pchisq(abs(Model24$logLik[1] - ModelNull$logLik[1])*2, 2)
#p_24vs12 <- 1 - pchisq(abs(Model24$logLik[1] - Model12$logLik[1])*2, 2)
#p_12vs6 <- 1 - pchisq(abs(Model12$logLik[1] - Model6$logLik[1])*2, 2)
# Use log likelihood instead of restricted log likelihood
logLik(ModelNull, REML = F)
logLik(Model24, REML = F)
logLik(Model12, REML = F)
logLik(Model6, REML = F)
p_24vsNull <- 1 - pchisq(abs(logLik(Model24, REML = F)[1] - logLik(ModelNull, REML = F)[1])*2, 2)
p_24vs12 <- 1 - pchisq(abs(logLik(Model24, REML = F)[1] - logLik(Model12, REML = F)[1])*2, 2)
p_12vs6 <- 1 - pchisq(abs(logLik(Model12, REML = F)[1] - logLik(Model6, REML = F)[1])*2, 2)

# Check model predictions
acroPhase <- function(sin, cos) {
  if (is.na(sin) | is.na(cos)) { return(NA) }
  if (cos >= 0 & sin  > 0) {K=0 }
  if (cos <  0 & sin >= 0) {K=12} 
  if (cos <= 0 & sin  < 0) {K=12} 
  if (cos >  0 & sin <= 0) {K=24} 
  return(atan(sin/cos)/(2*pi)*24 + K)
}

acroPhase(sin = Model24$coefficients$fixed["SinHour"], cos = Model24$coefficients$fixed["CosHour"])
amplitude <- sqrt(Model24$coefficients$fixed["SinHour"]^2 + Model24$coefficients$fixed["CosHour"]^2) # Amplitude

# Plot fixed-effects prediction
Intercept <- Model24$coefficients$fixed["(Intercept)"]
sinPred <- Model24$coefficients$fixed["SinHour"]*sin(2*pi*(1:24)/24)
cosPred <- Model24$coefficients$fixed["CosHour"]*cos(2*pi*(1:24)/24)

pdf("Fig_summary_pred.pdf")
plot(Intercept + sinPred + cosPred, 
     xaxt = "n", xaxs = "i", yaxt = "n",
     type = "n", bty = "n",
     ylim = c(0, 1),
     xlab = "Time of day, h",
     ylab = "ln IL-6, pg/ml",
     cex.lab = 1.5)
rect(-2, -5, 7, 40, col = "Misty Rose", border = NA)
rect(22, -5, 31, 40, col = "Misty Rose", border = NA)
axis(side = 1, at = c(1, 6, 12, 18, 24), labels = c("01:00", "06:00", "12:00", "18:00", "24:00"), cex.axis = 1.5)
axis(2, at = c(0, 0.5, 1), labels = c(0, 0.5, 1), cex.axis = 1.5)
lines(Intercept + sinPred + cosPred, lwd = 2, col = "red")
#legend("bottomleft", cex = 1.5, bty = "n", legend = c(paste("p", ifelse(p_24vsNull <0.001, "< 0.001", round(p_24vsNull, 3)), "for 24 vs null")))
#lines(Model12$coefficients$fixed["(Intercept)"] + Model12$coefficients$fixed["SinHour"]*sin(2*pi*(1:24)/24) + Model12$coefficients$fixed["CosHour"]*cos(2*pi*(1:24)/24) + Model12$coefficients$fixed["SinHour12"]*sin(2*pi*(1:12)/12) + Model12$coefficients$fixed["CosHour12"]*cos(2*pi*(1:12)/12), lwd = 2, col = "blue")
#lines(Model6$coefficients$fixed["(Intercept)"] + Model6$coefficients$fixed["SinHour"]*sin(2*pi*(1:24)/24) + Model6$coefficients$fixed["CosHour"]*cos(2*pi*(1:24)/24) + Model6$coefficients$fixed["SinHour12"]*sin(2*pi*(1:12)/12) + Model6$coefficients$fixed["CosHour12"]*cos(2*pi*(1:12)/12) + Model6$coefficients$fixed["SinHour6"]*sin(2*pi*(1:6)/6) + Model6$coefficients$fixed["CosHour12"]*cos(2*pi*(1:6)/6), lwd = 2, col = cols[3])
#legend("bottomleft", lty = c(2, 1), col = c("red", "blue"), legend = c("24 h cosinor", "24 + 12 h cosinor"), cex = 1.5, lwd = 2, bg = "white")

plot(Intercept + sinPred + cosPred, 
     xaxt = "n", xaxs = "i", yaxt = "n",
     type = "n", bty = "n",
     ylim = c(0, 1),
     xlab = "Time of day, h",
     ylab = "ln IL-6, pg/ml",
     cex.lab = 1.5)
rect(-2, -5, 7, 40, col = "Misty Rose", border = NA)
rect(22, -5, 31, 40, col = "Misty Rose", border = NA)
axis(side = 1, at = c(1, 6, 12, 18, 24), labels = c("01:00", "06:00", "12:00", "18:00", "24:00"), cex.axis = 1.5)
axis(2, at = c(0, 0.5, 1), labels = c(0, 0.5, 1), cex.axis = 1.5)
lines(Model12$coefficients$fixed["(Intercept)"] + Model12$coefficients$fixed["SinHour"]*sin(2*pi*(1:24)/24) + Model12$coefficients$fixed["CosHour"]*cos(2*pi*(1:24)/24) + Model12$coefficients$fixed["SinHour12"]*sin(2*pi*(1:12)/12) + Model12$coefficients$fixed["CosHour12"]*cos(2*pi*(1:12)/12), lwd = 2, col = "red")


plot(exp(Intercept + sinPred + cosPred), 
     xaxt = "n", xaxs = "i", yaxt = "n",
     type = "n", bty = "n",
     ylim = c(0, 3),
     xlab = "Time of day, h",
     ylab = "IL-6, pg/ml",
     cex.lab = 1.5)
rect(-2, -5, 7, 40, col = "Misty Rose", border = NA)
rect(22, -5, 31, 40, col = "Misty Rose", border = NA)
axis(side = 1, at = c(0, 6, 12, 18, 24), labels = c("00:00", "06:00", "12:00", "18:00", "24:00"), cex.axis = 1.5)
axis(2, at = c(0, 1, 2, 3), cex.axis = 1.5)
lines(exp(Intercept + sinPred + cosPred), lwd = 2, col = "red", lty = 2)
lines(exp(Model12$coefficients$fixed["(Intercept)"] + Model12$coefficients$fixed["SinHour"]*sin(2*pi*(1:24)/24) + Model12$coefficients$fixed["CosHour"]*cos(2*pi*(1:24)/24) + Model12$coefficients$fixed["SinHour12"]*sin(2*pi*(1:12)/12) + Model12$coefficients$fixed["CosHour12"]*cos(2*pi*(1:12)/12)), lwd = 2, col = "blue")
#lines(exp(Model6$coefficients$fixed["(Intercept)"] + Model6$coefficients$fixed["SinHour"]*sin(2*pi*(1:24)/24) + Model6$coefficients$fixed["CosHour"]*cos(2*pi*(1:24)/24) + Model6$coefficients$fixed["SinHour12"]*sin(2*pi*(1:12)/12) + Model6$coefficients$fixed["CosHour12"]*cos(2*pi*(1:12)/12) + Model6$coefficients$fixed["SinHour6"]*sin(2*pi*(1:6)/6) + Model6$coefficients$fixed["CosHour12"]*cos(2*pi*(1:6)/6)), lwd = 2, col = "green")
legend("bottomleft", lty = c(2, 1), col = c("red", "blue"), legend = c("24 h cosinor", "24 + 12 h cosinor"), cex = 1.5, lwd = 2)
dev.off()

# Plot prediction by dataset
pdf("Fig_eachstudy.pdf")
for (i in unique(IL6_data$Study)){
  with(IL6_data[IL6_data$Study == i, ],
       plotCI(Hour, lnIL6, IL6_sd, main = i, ylab = "", frame.plot = F, cex.axis=2, cex.lab=2, cex.main=2)
  )
  #lines(predict(Model24)[IL6_data$Study == i] ~ IL6_data[IL6_data$Study == i, ]$Hour, col = "red", lty = 2)
  lines(predict(Model12)[IL6_data$Study == i] ~ IL6_data[IL6_data$Study == i, ]$Hour, col = "red", lwd = 2)
}
dev.off()

# Get regression weights
weights <- aggregate(sqw ~ Study, data = IL6_data, FUN = "sum")
weights$sqw_percent <- 100*weights$sqw/sum(weights$sqw)
n_timepoints <- data.frame(table(IL6_data$Study))
weights <- merge(weights, n_timepoints, by.x = "Study", by.y = "Var1")
sum(aggregate(n ~ Study, data = IL6_data, FUN = "max")$n) # Total number of subjects
sum(weights$Freq) # Total number of time points

# Model variances
# I think the point is that higher variance may be confused for higher absolute values unless data are log-transformed
Model24sd <- lme(IL6_sd ~ SinHour + CosHour, data = IL6_data, weights = ~ 1/sqw, random = ~ 1|Study, na.action = "na.omit")
summary(Model24sd)
plot(predict(Model24sd) ~ IL6_data$Hour[!is.na(IL6_data$IL6_sd)])

# Investigate effect of sleep
modsleep24a <- update(Model24, ~ SinHour + CosHour + TimeFromCatheter_h + Sleep)
modsleep24b <- update(Model24, ~ SinHour + CosHour + TimeFromCatheter_h + TimeAsleep + TimeAwake)
modsleep24c <- update(Model24, ~ SinHour + CosHour + TimeFromCatheter_h + Sleep + TimeAsleep + TimeAwake)
summary(modsleep24a)
summary(modsleep24b)
summary(modsleep24c)
logLik(modsleep24a, REML = F)
logLik(modsleep24b, REML = F)
logLik(modsleep24c, REML = F)

modsleep12a <- update(Model24, ~ SinHour + CosHour + SinHour12 + CosHour12 + TimeFromCatheter_h + Sleep)
modsleep12b <- update(Model24, ~ SinHour + CosHour + SinHour12 + CosHour12 + TimeFromCatheter_h + TimeAsleep + TimeAwake)
modsleep12c <- update(Model24, ~ SinHour + CosHour + SinHour12 + CosHour12 + TimeFromCatheter_h + Sleep + TimeAsleep + TimeAwake)
summary(modsleep12a)
summary(modsleep12b)
summary(modsleep12c)
logLik(modsleep12a, REML = F)
logLik(modsleep12b, REML = F)
logLik(modsleep12c, REML = F)

p_12avs12 <- 1 - pchisq(abs(logLik(modsleep12a, REML = F)[1] - logLik(Model12, REML = F)[1])*2, 2)
p_12bvs12 <- 1 - pchisq(abs(logLik(modsleep12c, REML = F)[1] - logLik(Model12, REML = F)[1])*2, 2)
p_12cvs12 <- 1 - pchisq(abs(logLik(modsleep12c, REML = F)[1] - logLik(Model12, REML = F)[1])*2, 2)

p_24avs12 <- 1 - pchisq(abs(logLik(modsleep24a, REML = F)[1] - logLik(Model12, REML = F)[1])*2, 2)
p_24bvs12 <- 1 - pchisq(abs(logLik(modsleep24c, REML = F)[1] - logLik(Model12, REML = F)[1])*2, 2)
p_24cvs12 <- 1 - pchisq(abs(logLik(modsleep24c, REML = F)[1] - logLik(Model12, REML = F)[1])*2, 2)

p_24avs24 <- 1 - pchisq(abs(logLik(modsleep24a, REML = F)[1] - logLik(Model24, REML = F)[1])*2, 2)
p_24bvs24 <- 1 - pchisq(abs(logLik(modsleep24c, REML = F)[1] - logLik(Model24, REML = F)[1])*2, 2)
p_24cvs24 <- 1 - pchisq(abs(logLik(modsleep24c, REML = F)[1] - logLik(Model24, REML = F)[1])*2, 2)

#p_12avs12 <- 1 - pchisq(abs(modsleep12a$logLik[1] - Model12$logLik[1])*2, 2)
#p_12bvs12 <- 1 - pchisq(abs(modsleep12b$logLik[1] - Model12$logLik[1])*2, 2)
#p_12cvs12 <- 1 - pchisq(abs(modsleep12c$logLik[1] - Model12$logLik[1])*2, 2)

#p_24avs12 <- 1 - pchisq(abs(modsleep24a$logLik[1] - Model12$logLik[1])*2, 2)
#p_24bvs12 <- 1 - pchisq(abs(modsleep24b$logLik[1] - Model12$logLik[1])*2, 2)
#p_24cvs12 <- 1 - pchisq(abs(modsleep24c$logLik[1] - Model12$logLik[1])*2, 2)

#p_24avs24 <- 1 - pchisq(abs(modsleep24a$logLik[1] - Model24$logLik[1])*2, 2)
#p_24bvs24 <- 1 - pchisq(abs(modsleep24b$logLik[1] - Model24$logLik[1])*2, 2)
#p_24cvs24 <- 1 - pchisq(abs(modsleep24c$logLik[1] - Model24$logLik[1])*2, 2)


# Investigate time from catheterization effect
ranef <- ranef(Model12)
pdf("Fig_catheter.pdf")
boxplot(ranef$TimeFromCatheter_h[ranef$TimeFromCatheter_h < 0.1 & ranef$TimeFromCatheter_h > -0.1], frame.plot = F, ylim = c(-0.25, 0.2), ylab = "Estimated slope")
stripchart(ranef$TimeFromCatheter_h, add = T, method = "jitter", vertical = T, pch = 1, col = "red")
text(1.05, 0.2, labels = "Späth-Schwalbe 1998", pos = 4)
text(1.05, -0.24, labels = "Ündar 1999", pos = 4)
dev.off()
t.test(ranef$TimeFromCatheter_h)
# TODO: Add weights to t-test


###########################
# ANALYSE INDIVIDUAL PARTICIPANT DATA
###########################

# Transform data
IL6ind$lnIL6 <- log(IL6ind$IL6_pgml)
IL6ind$log2IL6 <- log2(IL6ind$IL6_pgml)
IL6ind$log10IL6 <- log10(IL6ind$IL6_pgml)
IL6ind$sqrtIL6 <- sqrt(IL6ind$IL6_pgml)

IL6ind$Study[IL6ind$Study == "Axelsson2013"] <- "Karshikoff2015"

# Inspect distributions
density_raw <- density(IL6ind$IL6_pgml)
density_raw <- fun_center(density_raw)
density_log2 <- density(IL6ind$log2IL6)
density_log2 <- fun_center(density_log2)
density_ln <- density(IL6ind$lnIL6)
density_ln <- fun_center(density_ln)
density_log10 <- density(IL6ind$log10IL6)
density_log10 <- fun_center(density_log10)
density_sqrt <- density(IL6ind$sqrtIL6)
density_sqrt <- fun_center(density_sqrt)

pdf("Fig_transformations.pdf")
plot(density_raw, xaxs="i", yaxs="i", xaxt="n", yaxt="n", xlim = c(-3, 4), ylim = c(0, 1.2), col = NULL, main = "", bty = "n", cex.axis = 1.5, cex.lab = 1.5)
axis(1, cex.axis = 1.5)
axis(2, at = c(0, 0.5, 1), labels = c(0, 0.5, 1), cex.axis = 1.5)
polygon(density_raw, col = "gray", border = NA)
lines(density_log2, col = "black", lwd = 2)
lines(density_ln, col = "green", lwd = 2)
lines(density_log10, col = "blue", lwd = 2)
lines(density_sqrt, col = "red", lwd = 2)
legend("topright", lty = 1, lwd = c(10, 2, 2, 2, 2), col = c("gray", "black", "green", "blue", "red"), legend = c("raw", expression('log'[2]), "ln", expression('log'[10]), expression(sqrt(phantom(0)))), bty = "n", cex = 1.5)
dev.off()

pdf("Fig_transformations2.pdf")
qqnorm(sort(IL6ind$IL6_pgml), main = "raw", bty = "n", cex.main = 1.5, cex.axis = 1.5, cex.lab = 1.5)
qqline(sort(IL6ind$IL6_pgml))

qqnorm(sort(IL6ind$sqrtIL6), main = expression(paste(sqrt(phantom(0)), "-transformed"), "-transformed"), bty = "n", cex.main = 1.5, cex.axis = 1.5, cex.lab = 1.5)
qqline(sort(IL6ind$sqrtIL6))

qqnorm(sort(IL6ind$lnIL6), main = "ln-transformed", bty = "n", cex.main = 1.5, cex.axis = 1.5, cex.lab = 1.5)
qqline(sort(IL6ind$lnIL6))
dev.off()

qqnorm(sort(IL6ind$log2IL6), main = "log2")
qqline(sort(IL6ind$log2IL6))
qqnorm(sort(IL6ind$log10IL6), main = "log10")
qqline(sort(IL6ind$log10IL6))

# Plot each set of individual participant data
pdf("Fig_ind.pdf")
with(IL6ind[IL6ind$Study == "Sothern1995", ],
     plot(lnIL6 ~ Hour, type = "n", main = "Sothern 1995", bty = "n", xaxt = "n", xlab = "Time, h", ylab = "ln(IL-6, pg/ml)", cex.main = 1.5, cex.axis = 1.5, cex.lab = 1.5, ylim = c(-3, 4)))
rect(23, -5, 31, 40, col = "Misty Rose", border = NA)
axis(1, at = c(19, 25, 31, 37, 43), labels = c("19:00", "01:00", "07:00", "13:00", "19:00"), cex.axis = 1.5, cex.lab = 1.5)
for (j in unique(IL6ind$Subject[IL6ind$Study == "Sothern1995"])){
  for (k in unique(IL6ind$Day[IL6ind$Study == "Sothern1995" & IL6ind$Subject == j])){
    with(IL6ind[IL6ind$Study == "Sothern1995" & IL6ind$Subject == j & IL6ind$Day == k, ],
         lines(lnIL6 ~ Hour, col = "gray", type = "o"))
  }
}
agg <- aggregate(IL6ind$lnIL6[IL6ind$Study == "Sothern1995"], list(h = IL6ind$Hour[IL6ind$Study == "Sothern1995"]), mean, na.action = na.omit)
n_vec <- NULL
for (j in unique(agg$h)){
  n <- sum(IL6ind$Hour[IL6ind$Study == "Sothern1995"] == j)
  n_vec <- c(n_vec, n)
}
agg$n <- n_vec
lo <- loess(x ~ h, weights = n, data = agg)
lines(predict(lo) ~ agg$h, lwd = 2)

with(IL6ind[IL6ind$Study == "Knudsen2008", ],
     plot(lnIL6 ~ Hour, type = "n", main = "Knudsen 2008", bty = "n", xaxt = "n", xlab = "Time, h", ylab = "ln(IL-6, pg/ml)", cex.main = 1.5, cex.axis = 1.5, cex.lab = 1.5, ylim = c(-3, 4)))
rect(23, -5, 31, 40, col = "Misty Rose", border = NA)
axis(1, at = c(10, 16, 22, 34), labels = c("10:00", "16:00", "22:00", "10:00"), cex.axis = 1.5, cex.lab = 1.5)
for (j in unique(IL6ind$Subject[IL6ind$Study == "Knudsen2008"])){
  for (k in unique(IL6ind$Day[IL6ind$Study == "Knudsen2008" & IL6ind$Subject == j])){
    with(IL6ind[IL6ind$Study == "Knudsen2008" & IL6ind$Subject == j & IL6ind$Day == k, ],
         lines(lnIL6 ~ Hour, col = "gray", type = "o"))
  }
}
agg <- aggregate(IL6ind$lnIL6[IL6ind$Study == "Knudsen2008"], list(h = IL6ind$Hour[IL6ind$Study == "Knudsen2008"]), mean, na.action = na.omit)
n_vec <- NULL
for (j in unique(agg$h)){
  n <- sum(IL6ind$Hour[IL6ind$Study == "Knudsen2008"] == j)
  n_vec <- c(n_vec, n)
}
agg$n <- n_vec
lo <- loess(x ~ h, weights = n, data = agg)
lines(predict(lo) ~ agg$h, lwd = 2)

with(IL6ind[IL6ind$Study == "Karshikoff2015", ],
     plot(lnIL6 ~ Hour, type = "n", main = "Karshikoff 2015", bty = "n", xaxt = "n", xlab = "Time, h", ylab = "ln(IL-6, pg/ml)", cex.main = 1.5, cex.axis = 1.5, cex.lab = 1.5, xlim = c(0, 24), ylim = c(-3, 4)))
rect(0, -5, 7, 40, col = "Misty Rose", border = NA)
rect(23, -5, 31, 40, col = "Misty Rose", border = NA)
axis(1, at = c(0, 6, 10, 15, 24), labels = c("00:00", "06:00", "10:00", "15:00", "24:00"), cex.axis = 1.5, cex.lab = 1.5)
for (j in unique(IL6ind$Subject[IL6ind$Study == "Karshikoff2015"])){
  for (k in unique(IL6ind$Day[IL6ind$Study == "Karshikoff2015" & IL6ind$Subject == j])){
    with(IL6ind[IL6ind$Study == "Karshikoff2015" & IL6ind$Subject == j & IL6ind$Day == k, ],
         lines(lnIL6 ~ Hour, col = "gray", type = "o"))
  }
}
agg <- aggregate(IL6ind$lnIL6[IL6ind$Study == "Karshikoff2015"], list(h = IL6ind$Hour[IL6ind$Study == "Karshikoff2015"]), mean, na.action = na.omit)
n_vec <- NULL
for (j in unique(agg$h)){
  n <- sum(IL6ind$Hour[IL6ind$Study == "Karshikoff2015"] == j)
  n_vec <- c(n_vec, n)
}
agg$n <- n_vec
lo <- loess(x ~ h, weights = n, data = agg)
lines(predict(lo) ~ agg$h, lwd = 2)

with(IL6ind[IL6ind$Study == "Lekander2013full", ],
     plot(lnIL6 ~ Hour, type = "n", main = "Lekander 2013", bty = "n", xaxt = "n", xlab = "Time, h", ylab = "ln(IL-6, pg/ml)", cex.main = 1.5, cex.axis = 1.5, cex.lab = 1.5, xlim = c(23, 96), ylim = c(-3, 4)))
rect(23, -5, 31, 40, col = "Misty Rose", border = NA)
rect(47, -5, 55, 40, col = "Misty Rose", border = NA)
rect(71, -5, 79, 40, col = "Misty Rose", border = NA)
rect(95, -5, 103, 40, col = "Misty Rose", border = NA)
axis(1, at = c(24, 36, 48, 60, 72, 84, 96), labels = c("00:00", "12:00", "00:00", "12:00", "00:00", "12:00", "00:00"), cex.axis = 1.5, cex.lab = 1.5)
for (j in unique(IL6ind$Subject[IL6ind$Study == "Lekander2013full"])){
  with(IL6ind[IL6ind$Study == "Lekander2013full" & IL6ind$Subject == j & IL6ind$Hour < 70, ],
       lines(lnIL6 ~ Hour, col = "tomato1", type = "o"))
}
for (j in unique(IL6ind$Subject[IL6ind$Study == "Lekander2013full"])){
  newdata <- IL6ind[IL6ind$Study == "Lekander2013full" & IL6ind$Subject == j & IL6ind$Hour > 70, ]
  newdata$Hour <- newdata$Hour - 168
  lines(lnIL6 ~ Hour, data = newdata, col = "lightblue", type = "o")
}
agg <- aggregate(IL6ind$lnIL6[IL6ind$Study == "Lekander2013full" & IL6ind$Hour < 70], list(h = IL6ind$Hour[IL6ind$Study == "Lekander2013full" & IL6ind$Hour < 70]), mean, na.action = na.omit)
n_vec <- NULL
for (j in unique(agg$h)){
  n <- sum(IL6ind$Hour[IL6ind$Study == "Lekander2013full" & IL6ind$Hour < 70] == j)
  n_vec <- c(n_vec, n)
}
agg$n <- n_vec
lo <- loess(x ~ h, weights = n, data = agg)
lines(predict(lo) ~ agg$h, lwd = 2, col = "red")
agg <- aggregate(IL6ind$lnIL6[IL6ind$Study == "Lekander2013full" & IL6ind$Hour > 70], list(h = IL6ind$Hour[IL6ind$Study == "Lekander2013full" & IL6ind$Hour > 70]), mean, na.action = na.omit)
n_vec <- NULL
for (j in unique(agg$h)){
  n <- sum(IL6ind$Hour[IL6ind$Study == "Lekander2013full" & IL6ind$Hour > 70] == j)
  n_vec <- c(n_vec, n)
}
agg$n <- n_vec
agg$h <- agg$h - 168
lo <- loess(x ~ h, weights = n, data = agg)
lines(predict(lo) ~ agg$h, lwd = 2, col = "blue")
dev.off()

# Define weights
IL6ind$nw <- NA
for (i in unique(IL6ind$Study)){
  IL6ind[IL6ind$Study == i, ]$nw <- IL6ind[IL6ind$Study == i, ]$n/length(IL6ind[IL6ind$Study == i, ]$n)
}

IL6ind$sqw <- NA
for (i in unique(IL6ind$Study)){
  IL6ind[IL6ind$Study == i, ]$sqw <- (IL6ind[IL6ind$Study == i, ]$n * sqrt(length(IL6ind$n[IL6ind$Study == i])))/length(IL6ind[IL6ind$Study == i, ]$n)
}

# Make variates for modelling
IL6ind$SinHour <- sin(2*pi*IL6ind$Hour/24)
IL6ind$CosHour <- cos(2*pi*IL6ind$Hour/24)
IL6ind$SinHour12 <- sin(2*pi*IL6ind$Hour/12)
IL6ind$CosHour12 <- cos(2*pi*IL6ind$Hour/12)
#IL6ind$IL6_variance <- IL6ind$IL6_sd^2
IL6ind$Hour_sq <- IL6ind$Hour^2
#IL6ind$LogTimeFromCatheter_h <- log(1 + IL6ind$TimeFromCatheter_h)
IL6ind$Hour24 <- IL6ind$Hour
#IL6ind$Deprived <- IL6ind$Condition != "Baseline"
while (max(IL6ind$Hour24) > 24){
  IL6ind$Hour24[IL6ind$Hour24 > 24] <- IL6ind$Hour24[IL6ind$Hour24 > 24] - 24
}

# Build models
Modelind_null <- lme(lnIL6 ~ TimeFromCatheter_h, data = IL6ind, weights = ~ 1/sqw, random = ~ TimeFromCatheter_h|Study/Subject)
Modelind_24 <- lme(lnIL6 ~ SinHour + CosHour + TimeFromCatheter_h, data = IL6ind, weights = ~ 1/sqw, random = ~ TimeFromCatheter_h|Study/Subject)
Modelind_12 <- lme(lnIL6 ~ SinHour + CosHour + SinHour12 + CosHour12 + TimeFromCatheter_h, data = IL6ind, weights = ~ 1/sqw, random = ~ TimeFromCatheter_h|Study/Subject)

summary(Modelind_null)
summary(Modelind_24)
summary(Modelind_12)

# Build loess function on residuals from null model
agg <- aggregate(rowSums(Modelind_null$residuals[, 1:3]), Modelind_null$data["Hour"], mean, na.action = na.omit)
w_vec <- NULL
for (i in unique(agg$Hour)){
  w <- sum((IL6ind$sqw[IL6ind$Study != "Lekander2013PSD"])[(IL6ind$Hour[IL6ind$Study != "Lekander2013PSD"] == i)])
  w_vec <- c(w_vec, w)
}
agg$w <- w_vec
while(any(agg$Hour < 0)){
  agg$Hour[agg$Hour < 0] <- agg$Hour[agg$Hour < 0] + 24
}
while(any(agg$Hour > 24)){
  agg$Hour[agg$Hour > 24] <- agg$Hour[agg$Hour > 24] - 24
}
agg <- agg[order(agg$Hour), ]
agg2 <- agg
agg3 <- agg
agg4 <- agg
agg5 <- agg
agg2$Hour <- agg2$Hour - 24
agg3$Hour <- agg3$Hour + 24
agg4$Hour <- agg4$Hour - 48
agg5$Hour <- agg5$Hour + 48
agg <- rbind(agg4, agg2, agg, agg3, agg5)
lo_ind_0.4 <- loess(x ~ Hour, weights = w, data = agg, span = 0.4)
lo_ind_0.2 <- loess(x ~ Hour, weights = w, data = agg, span = 0.2)
plot(agg$x ~ agg$Hour)
lines((Modelind_null$coefficients$fixed["(Intercept)"] + predict(lo_ind_0.2)) ~ agg$Hour, lwd = 2, col = "black")
lines((Modelind_null$coefficients$fixed["(Intercept)"] + predict(lo_ind_0.4)) ~ agg$Hour, lwd = 2, col = "gray")

# Test models for statistical significance
p_ind24vsNull <- 1 - pchisq(abs(Modelind_24$logLik[1] - Modelind_null$logLik[1])*2, 2)
p_ind24vs12 <- 1 - pchisq(abs(Modelind_24$logLik[1] - Modelind_12$logLik[1])*2, 2)
#p_ind12vs24 <- 1 - pchisq(abs(Modelind_12$logLik[1] - Modelind_24$logLik[1])*2, 2) #Same result as 24vs12
p_ind12vsNull <- 1 - pchisq(abs(Modelind_12$logLik[1] - Modelind_null$logLik[1])*2, 2)

# Plot fixed-effects prediction
Intercept24 <- Modelind_24$coefficients$fixed["(Intercept)"]
sinPred24 <- Modelind_24$coefficients$fixed["SinHour"]*sin(2*pi*(1:24)/24)
cosPred24 <- Modelind_24$coefficients$fixed["CosHour"]*cos(2*pi*(1:24)/24)

Intercept12 <- Modelind_12$coefficients$fixed["(Intercept)"]
sinPred12a <- Modelind_12$coefficients$fixed["SinHour"]*sin(2*pi*(1:24)/24)
cosPred12a <- Modelind_12$coefficients$fixed["CosHour"]*cos(2*pi*(1:24)/24)
sinPred12b <- Modelind_12$coefficients$fixed["SinHour12"]*sin(2*pi*(1:12)/12)
cosPred12b <- Modelind_12$coefficients$fixed["CosHour12"]*cos(2*pi*(1:12)/12)

pdf("Fig_ind_pred.pdf")
plot(Intercept12 + sinPred12a + cosPred12a + sinPred12b + cosPred12b, 
     bty = "n",
     xaxt= "n",
     yaxt= "n",
     xaxs= "i",
     type = "n", 
     #xlim = c(0, 2),
     ylim = c(0, 1),
     #main = "Predicted diurnal variation",
     xlab = "Time of day, h",
     ylab = "ln IL-6, pg/ml",
     cex.lab = 1.5)
rect(-2, -5, 7, 40, col = "Misty Rose", border = NA)
rect(22, -5, 31, 40, col = "Misty Rose", border = NA)
axis(side = 1, at = c(0, 6, 12, 18, 24), labels = c("00:00", "06:00", "12:00", "18:00", "24:00"), cex.axis = 1.5)
axis(2, at = c(0, 0.5, 1), cex.axis = 1.5)
#lines((Modelind_null$coefficients$fixed["(Intercept)"] + predict(lo_ind_0.3)) ~ agg$Hour, lwd = 2, col = "black")
#lines((Modelind_null$coefficients$fixed["(Intercept)"] + predict(lo_ind_0.4)) ~ agg$Hour, lwd = 2, col = "gray")
lines(Intercept12 + sinPred24 + cosPred24, lwd = 2, col = "red")
#lines(Intercept12 + sinPred12a + cosPred12a + sinPred12b + cosPred12b, lwd = 2, col = "blue")
#legend("bottomright", lwd = 2, col = c("black", "gray", "red", "blue"), legend = c("loess, span 0.2", "loess, span 0.4", "24 h period", "24 and 12 h period"), bty = "n", cex = 1.5)
legend("bottomleft", cex = 1.5, bty = "n", 
       legend = c(paste("p =", round(p_ind24vsNull, 3), "for 24 vs null"),
                  paste("p =", round(p_ind24vs12, 3), "for 24 vs 24 and 12")))

plot(exp(Intercept12 + sinPred12a + cosPred12a + sinPred12b + cosPred12b), 
     bty = "n",
     xaxt= "n",
     yaxt= "n",
     xaxs= "i",
     type = "n", 
     #xlim = c(0, 2),
     ylim = c(0, 3),
     #main = "Predicted diurnal variation",
     xlab = "Time of day, h",
     ylab = "IL-6, pg/ml",
     cex.lab = 1.5)
rect(-2, -5, 7, 40, col = "Misty Rose", border = NA)
rect(22, -5, 31, 40, col = "Misty Rose", border = NA)
axis(side = 1, at = c(0, 6, 12, 18, 24), labels = c("00:00", "06:00", "12:00", "18:00", "24:00"), cex.axis = 1.5)
axis(2, at = c(0, 1, 2, 3), cex.axis = 1.5)
#lines(exp(Modelind_null$coefficients$fixed["(Intercept)"] + predict(lo_ind_0.3)) ~ agg$Hour, lwd = 2, col = "black")
#lines(exp(Modelind_null$coefficients$fixed["(Intercept)"] + predict(lo_ind_0.4)) ~ agg$Hour, lwd = 2, col = "gray")
lines(exp(Intercept12 + sinPred24 + cosPred24), lwd = 2, col = "red")
#lines(exp(Intercept12 + sinPred12a + cosPred12a + sinPred12b + cosPred12b), lwd = 2, col = "blue")
dev.off()

# Plot predicted by observed for each dataset
getPalette = colorRampPalette(brewer.pal(9, "Set1"))
pdf("Fig_ind_pred_obs.pdf")
for (i in unique(IL6ind$Study[IL6ind$Study != "Lekander2013PSD"])){
  plot(lnIL6 ~ Hour, main = i, data = IL6ind[IL6ind$Study == i, ], ylab = "ln IL-6, pg/ml", bty = "n", cex.main = 1.5, cex.axis = 1.5, cex.lab = 1.5, xaxt = "n")
  axis(side = 1, at = c(-6, 0, 6, 12, 18, 24, 30, 36, 42), labels = c("18:00", "00:00", "06:00", "12:00", "18:00", "00:00", "06:00", "12:00", "18:00"), cex.axis = 1.5)
  cnt <- 1
  for (j in unique(IL6ind$Subject[IL6ind$Study == i])){
    lines(predict(Modelind_24)[IL6ind$Study == i & IL6ind$Subject == j] ~ IL6ind[IL6ind$Study == i & IL6ind$Subject == j, ]$Hour, col = getPalette(length(unique(IL6ind$Subject[IL6ind$Study == i])))[cnt])
    points(lnIL6 ~ Hour, main = i, data = IL6ind[IL6ind$Study == i & IL6ind$Subject == j, ], col = getPalette(length(unique(IL6ind$Subject[IL6ind$Study == i])))[cnt])
    cnt <- cnt + 1
  }
}
dev.off()

# Visualise weights for each study
Weights <- data.frame("Study" = unique(IL6ind$Study))
Weights <- merge(Weights, aggregate(nw ~ Study, data = IL6ind, sum), by = "Study")
Weights <- merge(Weights, aggregate(sqw ~ Study, data = IL6ind, sum), by = "Study")
nobs_vec <- NULL
for (i in sort(unique(IL6ind$Study))){
  nobs <- length(IL6ind$nw[IL6ind$Study == i])
  nobs_vec <- c(nobs_vec, nobs)
}
Weights$nobs <- nobs_vec
Weights_prop <- prop.table(as.matrix(Weights[, 2:4]), 2)

pdf("Weights_ind.pdf")
barplot(Weights_prop, col = brewer.pal(4, "Set1"), yaxt = "n", cex.names = 1.5, names = c("n_subj", "sqrt", "n_obs"))
text(0.2, 0.2, "Karshikoff2015", pos = 4, cex = 1.3)
text(0.2, 0.5, "Knudsen2008", pos = 4, cex = 1.3)
text(0.2, 0.7, "Lekander2013", pos = 4, cex = 1.3)
text(0.2, 0.9, "Sothern1995", pos = 4, cex = 1.3)
dev.off()