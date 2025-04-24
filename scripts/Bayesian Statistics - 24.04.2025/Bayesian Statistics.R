library(rethinking)

# Read data
my_data <- read.csv("db.csv", stringsAsFactors = FALSE)

# Process data
my_data$model <- as.factor(my_data$model)
i <- my_data$treatment %in% c("ED2", "ED1")
my_data$RxE <- i + 1L # 2 for ERI, 1 for else

# Create binary indicators for modeling
my_data$is_ERI <- as.integer(my_data$RxE == 2)
my_data$is_LM04 <- as.integer(my_data$model == "LM04") 

# Prepare data for modeling
d <- list(
  axis_length = my_data$axis_minor_length,
  is_ERI = my_data$is_ERI,
  is_LM04 = my_data$is_LM04,
  model_ERI = my_data$is_LM04 * my_data$is_ERI  # Interaction term
)

#########################  1  ##############################

# Model 1: Exponential distribution with full interaction model
# This allows for different effects of ERI for each model
m_full <- ulam(
  alist(
    axis_length ~ dexp(lambda),
    log(lambda) <- a + b_model*is_LM04 + b_ERI*is_ERI + b_int*model_ERI,
    a ~ dnorm(0, 1),
    b_model ~ dnorm(0, 0.5),
    b_ERI ~ dnorm(0, 0.5),
    b_int ~ dnorm(0, 0.5)
  ), data=d, chains=4, cores=4, iter=2000
)

# Print model results
precis(m_full)

# Plot posterior distributions
post <- extract.samples(m_full)

# Calculate lambda for each group
lambda_LM02_else <- exp(post$a)
lambda_LM02_ERI <- exp(post$a + post$b_ERI)
lambda_LM04_else <- exp(post$a + post$b_model)
lambda_LM04_ERI <- exp(post$a + post$b_model + post$b_ERI + post$b_int)

# Calculate expected values (means) for each group (1/lambda for exponential)
mean_LM02_else <- 1/lambda_LM02_else
mean_LM02_ERI <- 1/lambda_LM02_ERI
mean_LM04_else <- 1/lambda_LM04_else
mean_LM04_ERI <- 1/lambda_LM04_ERI

# Calculate mean differences due to ERI for each model
diff_LM02 <- mean_LM02_ERI - mean_LM02_else
diff_LM04 <- mean_LM04_ERI - mean_LM04_else

# Plot the differences
png("ERI_effect_sizes.png", width=800, height=600)
par(mfrow=c(2,2))

# Density plots of posterior differences
plot(density(diff_LM02), main="ERI Effect on LM02", 
     xlab="Difference in Mean Minor Axis Length", lwd=2)
abline(v=0, lty=2)
abline(v=mean(diff_LM02), col="red", lwd=2)

plot(density(diff_LM04), main="ERI Effect on LM04", 
     xlab="Difference in Mean Minor Axis Length", lwd=2)
abline(v=0, lty=2)
abline(v=mean(diff_LM04), col="red", lwd=2)

# Plot expected means for each group
means <- data.frame(
  group = c("LM02_else", "LM02_ERI", "LM04_else", "LM04_ERI"),
  mean = c(mean(mean_LM02_else), mean(mean_LM02_ERI), 
           mean(mean_LM04_else), mean(mean_LM04_ERI)),
  lower = c(quantile(mean_LM02_else, 0.025), quantile(mean_LM02_ERI, 0.025),
            quantile(mean_LM04_else, 0.025), quantile(mean_LM04_ERI, 0.025)),
  upper = c(quantile(mean_LM02_else, 0.975), quantile(mean_LM02_ERI, 0.975),
            quantile(mean_LM04_else, 0.975), quantile(mean_LM04_ERI, 0.975))
)

# Calculate probability of directional effect
prob_positive_LM02 <- mean(diff_LM02 > 0)
prob_positive_LM04 <- mean(diff_LM04 > 0)

# Print summary of effects
cat("\nEffect of ERI on LM02:\n")
cat("Mean difference:", mean(diff_LM02), "\n")
cat("95% CI:", quantile(diff_LM02, c(0.025, 0.975)), "\n")
cat("Probability of positive effect:", prob_positive_LM02, "\n\n")

cat("Effect of ERI on LM04:\n")
cat("Mean difference:", mean(diff_LM04), "\n")
cat("95% CI:", quantile(diff_LM04, c(0.025, 0.975)), "\n")
cat("Probability of positive effect:", prob_positive_LM04, "\n")

# Visualize posterior distributions of means
plotIdx <- 1:4
plot(plotIdx, means$mean, ylim=c(min(means$lower), max(means$upper)),
     pch=19, xlab="Group", ylab="Expected Minor Axis Length",
     main="Posterior Means with 95% CI", xaxt="n")
axis(1, at=plotIdx, labels=means$group)
segments(plotIdx, means$lower, plotIdx, means$upper)
dev.off()


### Saving of data su file

# Crea un file per salvare i risultati numerici dettagliati
sink("Results_ERI.txt")

# Stampa intestazione
cat("==========================================\n")
cat("ANALISI BAYESIANA DELL'EFFETTO ERI\n")
cat("==========================================\n\n")

# Effetto di ERI su LM02
cat("EFFETTO ERI SU LM02:\n")
cat("-------------------\n")
cat("Differenza media:", round(mean(diff_LM02), 2), "\n")
cat("Mediana della differenza:", round(median(diff_LM02), 2), "\n")
cat("Deviazione standard:", round(sd(diff_LM02), 2), "\n")
cat("Intervallo di credibilità al 95%:", round(quantile(diff_LM02, c(0.025, 0.975)), 2), "\n")
cat("Intervallo di credibilità al 89%:", round(quantile(diff_LM02, c(0.055, 0.945)), 2), "\n")
cat("Intervallo di credibilità al 50%:", round(quantile(diff_LM02, c(0.25, 0.75)), 2), "\n")
cat("Probabilità che l'effetto sia positivo:", round(mean(diff_LM02 > 0) * 100, 2), "%\n")
cat("Probabilità che l'effetto sia > 2:", round(mean(diff_LM02 > 2) * 100, 2), "%\n")
cat("Probabilità che l'effetto sia > 4:", round(mean(diff_LM02 > 4) * 100, 2), "%\n\n")

# Effetto di ERI su LM04
cat("EFFETTO ERI SU LM04:\n")
cat("-------------------\n")
cat("Differenza media:", round(mean(diff_LM04), 2), "\n")
cat("Mediana della differenza:", round(median(diff_LM04), 2), "\n")
cat("Deviazione standard:", round(sd(diff_LM04), 2), "\n")
cat("Intervallo di credibilità al 95%:", round(quantile(diff_LM04, c(0.025, 0.975)), 2), "\n")
cat("Intervallo di credibilità al 89%:", round(quantile(diff_LM04, c(0.055, 0.945)), 2), "\n")
cat("Intervallo di credibilità al 50%:", round(quantile(diff_LM04, c(0.25, 0.75)), 2), "\n")
cat("Probabilità che l'effetto sia positivo:", round(mean(diff_LM04 > 0) * 100, 2), "%\n")
cat("Probabilità che l'effetto sia > 2:", round(mean(diff_LM04 > 2) * 100, 2), "%\n")
cat("Probabilità che l'effetto sia > 4:", round(mean(diff_LM04 > 4) * 100, 2), "%\n\n")

# Confronto tra modelli
cat("CONFRONTO TRA EFFETTI IN LM02 E LM04:\n")
cat("-----------------------------------\n")
cat("Differenza tra effetti (LM02 - LM04):", round(mean(diff_LM02 - diff_LM04), 2), "\n")
cat("Intervallo di credibilità al 95%:", round(quantile(diff_LM02 - diff_LM04, c(0.025, 0.975)), 2), "\n")
cat("Probabilità che l'effetto sia maggiore in LM02:", round(mean(diff_LM02 > diff_LM04) * 100, 2), "%\n\n")

# Statistiche descrittive per ogni gruppo
cat("VALORI MEDI PER GRUPPO:\n")
cat("----------------------\n")
cat("Media LM02_else:", round(mean(mean_LM02_else), 2), 
    "(95% CI:", round(quantile(mean_LM02_else, 0.025), 2), "-", 
    round(quantile(mean_LM02_else, 0.975), 2), ")\n")
cat("Media LM02_ERI:", round(mean(mean_LM02_ERI), 2), 
    "(95% CI:", round(quantile(mean_LM02_ERI, 0.025), 2), "-", 
    round(quantile(mean_LM02_ERI, 0.975), 2), ")\n")
cat("Media LM04_else:", round(mean(mean_LM04_else), 2), 
    "(95% CI:", round(quantile(mean_LM04_else, 0.025), 2), "-", 
    round(quantile(mean_LM04_else, 0.975), 2), ")\n")
cat("Media LM04_ERI:", round(mean(mean_LM04_ERI), 2), 
    "(95% CI:", round(quantile(mean_LM04_ERI, 0.025), 2), "-", 
    round(quantile(mean_LM04_ERI, 0.975), 2), ")\n\n")

# Parametri del modello
cat("PARAMETRI DEL MODELLO:\n")
cat("--------------------\n")
cat("a (intercetta):", round(mean(post$a), 3), 
    "(95% CI:", round(quantile(post$a, 0.025), 3), "-", 
    round(quantile(post$a, 0.975), 3), ")\n")
cat("b_model (effetto LM04):", round(mean(post$b_model), 3), 
    "(95% CI:", round(quantile(post$b_model, 0.025), 3), "-", 
    round(quantile(post$b_model, 0.975), 3), ")\n")
cat("b_ERI (effetto ERI):", round(mean(post$b_ERI), 3), 
    "(95% CI:", round(quantile(post$b_ERI, 0.025), 3), "-", 
    round(quantile(post$b_ERI, 0.975), 3), ")\n")
cat("b_int (interazione):", round(mean(post$b_int), 3), 
    "(95% CI:", round(quantile(post$b_int, 0.025), 3), "-", 
    round(quantile(post$b_int, 0.975), 3), ")\n\n")

# Chiudi il file di testo
sink()

# Stampa anche a schermo
cat("Risultati dettagliati salvati nel file 'risultati_dettagliati_ERI.txt'\n\n")

# Creazione di un grafico aggiuntivo che mostra la probabilità cumulativa
png("probabilita_cumulativa_effetti.png", width=800, height=400)
par(mfrow=c(1,2))

# Funzione per calcolare la CDF empirica
plot_ecdf <- function(x, main, col="blue") {
  plot(ecdf(x), main=main, xlab="Dimensione dell'Effetto", 
       ylab="Probabilità Cumulativa", col=col, lwd=2)
  abline(v=0, lty=2)
  abline(h=0.5, lty=3)
  text(max(x)*0.7, 0.2, 
       paste("P(effetto > 0) =", round(mean(x > 0)*100, 1), "%"), 
       cex=1.2)
}

plot_ecdf(diff_LM02, "CDF Effetto ERI su LM02", "darkblue")
plot_ecdf(diff_LM04, "CDF Effetto ERI su LM04", "darkred")

dev.off()
