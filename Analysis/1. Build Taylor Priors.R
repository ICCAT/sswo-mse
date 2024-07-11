################################################################################




################################################################################


library(SSWO)

# ----- Sample the LH parameters from a truncated log-normal distribution -----

# Changes:
# - samples from log-normal distribution using CV
# - samples L50, L50_L95, and L50/Linf (assumed CV of 0.2) (don't use L50/Linf atm)
# - changed vtto CV to 0.2 to avoid sampling extreme negative values
# - changed vonK CV to an assumed 0.2 to avoid sampling extreme values (e.g < 0.02)
# - truncates at 1.96 SDs

simulated.data.trunc <- Sample_Trunc_LogNormal_Taylor()


# ----- Generate the Correlation Matrix from FishLife -----

# Changes:
# - returns 2nd rather than 4th list level from `FishLife::Plot_taxa`

LF_preds <- get_FishLife_Predictions()

cov_pred <- LF_preds$Cov_pred
cor_matrix <- calc_Correlation_Matrix(cov_pred)

# ---- Add correlation structure to M, Linf, K, and Lm ----

# Changes:
# - modified lower and upper bounds of truncated values to 1.96 SDs as above
# - sampled in log space

Var_Names <- data.frame(FishLife=c('M','Loo','K', 'Lm'),
                     TaylorCSV=c('mort', 'vLinf', 'vonK', 'L50')
)

simulated.data.trunc_updated <- Add_Correlation_Taylor(simulated.data.trunc,
                                            cor_matrix,
                                            Var_Names)



# ---- Calculate steepness -----

# Changes:
# - used L50 and L50_95 to calculate maturity-at-age

steepness <- sapply(1:nrow(simulated.data.trunc_updated),
                    function(x)
                      Calculate_Steepness(x, simulated.data.trunc_updated))

simulated.data.trunc_updated$steepness <- steepness
ind <- which(is.finite(simulated.data.trunc_updated$steepness))
simulated.data.trunc_updated <- simulated.data.trunc_updated[ind,]

head(simulated.data.trunc_updated)

###### CHECKS #####
dev.off()
plot(simulated.data.trunc_updated$mort/simulated.data.trunc_updated$vonK,
     simulated.data.trunc_updated$L50/simulated.data.trunc_updated$vLinf,
     xlim=c(0,8), ylim=c(0,1))

# Predicted M/K ratios appear too high
# M/K values > 3 are extremely rare

mean(simulated.data.trunc$mort)
mean(simulated.data.trunc$vonK)

# von K too low??

# predicted vonK from FishLife
exp(LF_preds$Mean_pred[2])

# Predicted M/K from FishLife
exp(LF_preds$Mean_pred[6])/exp(LF_preds$Mean_pred[2])

# Predicted L50/Linf from FishLife
exp(LF_preds$Mean_pred[7])/exp(LF_preds$Mean_pred[1])
