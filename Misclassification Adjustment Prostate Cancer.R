#Set the distributions for the bias parameters based on Validation Study Data

SE_SD <- (0.99-0.90)/(2*1.96)
SP_SD <- (1.00-0.95)/(2*1.96)

SEmean <- 0.945
SPmean <- 0.995

RRvector <- c()

alphaSE <- SEmean * ( ( ( SEmean*(1-SEmean) )/(SE_SD^2) ) -1 )
betaSE <- (1-SEmean) * ( ( ( SEmean*(1-SEmean) )/(SE_SD^2) ) -1 )

alphaSE <- alphaSE*1.5

alphaSP <- SPmean * ( ( ( SPmean*(1-SPmean) )/(SP_SD^2) ) -1 )
betaSP <- (1-SPmean) * ( ( ( SPmean*(1-SPmean) )/(SP_SD^2) ) -1 )

#Get the 2.5th, 50th, and 97.5th percentile values and plot distributions to check that our distribution matches the validation data

SEverif <- rbeta(1000,alphaSE, betaSE)
SEquant <- quantile(SEverif, c(0.025, 0.50, 0.975))

SPverif <- rbeta(1000,alphaSP, betaSP)
SPquant <- quantile(SPverif, c(0.025, 0.50, 0.975))

range <- seq(0, 1, length=100)

#plot(range, dbeta(range, alphaSE, betaSE), type='l')
#plot(range, dbeta(range, alphaSP, betaSP), type='l')

for (i in 1:250) {

#Non-differentially sample the bias parameters from our calculated distributions

SEsampled <- rbeta(1,alphaSE, betaSE)
SPsampled <- rbeta(1,alphaSP, betaSP)

#Generate an estimate of the true prevalence of the exposure while also propagating the conventional random error

D1_total <- 5448373
D0_total <- 296495347

A_0 <- ((0.813*D1_total) - (D1_total*(1 - SPsampled)))/(SEsampled - 1 + SPsampled)
B_0 <- D1_total - A_0

C_0 <- ((0.854*D0_total) - (D0_total*(1 - SPsampled)))/(SEsampled - 1 + SPsampled)
D_0 <- D0_total-C_0

P_ED1 <- rbeta(1, A_0, B_0)
P_ED0 <- rbeta(1, C_0, D_0)

#Calculate PPV and NPV and use it to perform misclassification adjustment

PPV1 <- (SEsampled*P_ED1)/((SEsampled*P_ED1) + (1-SPsampled)*(1-P_ED1))
NPV1 <- (SPsampled*(1-P_ED1))/(((1-SEsampled)*P_ED1) + SPsampled*(1-P_ED1))

W1 <- sum(rbinom(0.813*D1_total,1,PPV1))
Y1 <- sum(rbinom(0.187*D1_total,1,(1-NPV1)))
A <- W1 + Y1
B <- D1_total - A

PPV0 <- (SEsampled*P_ED0)/((SEsampled*P_ED0) + (1-SPsampled)*(1-P_ED0))
NPV0 <- (SPsampled*(1-P_ED0))/(((1-SEsampled)*P_ED0) + SPsampled*(1-P_ED0))

W0 <- sum(rbinom(0.854*D0_total,1,PPV0))
Y0 <- sum(rbinom(0.146*D0_total,1,(1-NPV0)))
C <- W0 + Y0
D <- D0_total - C

#Incorporate conventional random error

z <- rnorm(1,0,1)
logRR <- log((C/D)/(A/B))
stderror <- sqrt((1/A) + (1/B) + (1/C) + (1/D))
totalRR <- exp(logRR - (z*stderror))

RRvector <- append(RRvector,totalRR)
}
#plot histogram of ORs

RRquant <- quantile(RRvector, c(0.025, 0.50, 0.975), na.rm = TRUE)

hist(RRvector, breaks = 20, xlim = c(1.2,1.8), xlab = 'Odds Ratio', 
     ylab = 'Frequency', main = 'Odds Ratio Adjusted For Misclassification Bias')
abline(v= RRquant[2],col="blue",lwd=2)
abline(v= c(RRquant[1],RRquant[3]),col="blue",lwd=1, lty = 2)
abline(v=1.33,col="red",lwd=2)
abline(v=1.324,col="red",lwd=1, lty = 2)
abline(v=1.336,col="red",lwd=1, lty = 2)