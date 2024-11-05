# Chargement des packages nécessaires
require(survival)
require(rms)
library(olsrr)
# Simulation de données
set.seed(123)  # pour la reproductibilité

n <- 500  # nombre d'observations

# Variables explicatives
CURSMOKE <- rbinom(n, 1, 0.3)  # 1 = fumeur, 0 = non-fumeur
SEX <- rbinom(n, 1, 0.5)       # 1 = homme, 0 = femme
BMI <- rnorm(n, mean = 25, sd = 4)  # indice de masse corporelle
AGE <- rpois(n, lambda = 60) + 20  # âge (en années)
ID <- 1:n
# Simulation du temps avec une distribution de Weibull
shape_param <- 2  # paramètre de forme pour la distribution Weibull
scale_param <- exp(-2 + 0.5 * CURSMOKE + 0.3 * SEX + 0.1 * BMI + 0.01 * AGE)  # échelle

TIMEDTH <- rweibull(n, shape = shape_param, scale = scale_param)
DEATH <- rbinom(n, 1, 0.7)  # indicateur d'événement observé (70% ont un événement observé)

# Création d'un data frame avec les données simulées
sim_data <- data.frame(ID,CURSMOKE, SEX, BMI, AGE, TIMEDTH, DEATH)
head(sim_data)
# Ajustement des modèles AFT avec différentes distributions

# Modèle exponentiel
fit.exp <- survreg(Surv(TIMEDTH, DEATH) ~ CURSMOKE + SEX + BMI + log(AGE), dist = "exponential", data = sim_data)
BIC.exp <- AIC(fit.exp, k = log(n))  # Calcul du BIC

# Modèle de Weibull
fit.weibull <- survreg(Surv(TIMEDTH, DEATH) ~ CURSMOKE + SEX + BMI + log(AGE), dist = "weibull", data = sim_data)
BIC.weibull <- AIC(fit.weibull, k = log(n))  # Calcul du BIC

# Modèle log-normal
fit.lognormal <- survreg(Surv(TIMEDTH, DEATH) ~ CURSMOKE + SEX + BMI + log(AGE), dist = "lognormal", data = sim_data)
BIC.lognormal <- AIC(fit.lognormal, k = log(n))  # Calcul du BIC

# Modèle log-logistique
fit.loglogistic <- survreg(Surv(TIMEDTH, DEATH) ~ CURSMOKE + SEX + BMI + log(AGE), dist = "loglogistic", data = sim_data)
BIC.loglogistic <- AIC(fit.loglogistic, k = log(n))  # Calcul du BIC

# Comparaison des BIC
BIC_values <- c(Exponential = BIC.exp, Weibull = BIC.weibull, Lognormal = BIC.lognormal, Loglogistic = BIC.loglogistic)
print(BIC_values)

# Choix du meilleur modèle (celui avec le BIC le plus faible)
best_model <- names(BIC_values)[which.min(BIC_values)]
cat("Le meilleur modèle selon le BIC est :", best_model, "\n")

  fit.weibull_psm <- psm(Surv(TIMEDTH, DEATH) ~ CURSMOKE + SEX + BMI + log(AGE), dist = "weibull", data = sim_data, y = TRUE)
  std.resid <- residuals(fit.weibull_psm, type = "censored.normalized")[, 1]
  cox.resid = -log(val.surv(fit.weibull_psm)$est.surv);

  par(mfrow = c(2, 3))
  # Absence de relation residuelle
  
  plot(x = sim_data$BMI, y = std.resid);
  lines(smooth.spline(x = sim_data$BMI, y = std.resid), col = "purple", lwd = 2);
  
  plot(x = sim_data$AGE, y = std.resid);
  lines(smooth.spline(x = sim_data$AGE, y = std.resid), col = "purple", lwd = 2);
  
  plot(x = sim_data$BMI, y = cox.resid);
  lines(smooth.spline(x = sim_data$BMI, y = cox.resid), col = "purple", lwd = 2);
  
  plot(x = sim_data$AGE, y = cox.resid);
  lines(smooth.spline(x = sim_data$AGE, y = cox.resid), col = "purple", lwd = 2);  
  

  # Absence d'observations influentes
  dfbetas = residuals(fit.weibull, type = "dfbetas");
  plot(y = dfbetas[,2], x = sim_data$ID, type = "l", xlab = "Obs",
       ylab = "DFBETA de CURSMOKE ", ylim = c(-0.25, 0.25));
  abline(h = c(-0.2, 0.2), col = "purple", lty = 2, lwd = 2);
  
  
  ## Hypothese distributionnelle
  plot(survfit(Surv(cox.resid, DEATH)~1, data = sim_data), fun = "cumhaz");
  abline(a = 0, b = 1, lwd = 2, lty = 3, col = "purple");

  
  fit.weibull;
  exp(fit.weibull$coef)[2];
  c(exp(fit.weibull$coef[2] - 1.96*sqrt(diag(fit.weibull$var)[2])),
    exp(fit.weibull$coef[2] + 1.96*sqrt(diag(fit.weibull$var)[2])));
  
