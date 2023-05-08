#set distributions for bias parameters
alpha_e_d <- 12*8
beta_e_d <- 88*8

alpha_e_nd <- 11*8
beta_e_nd <- 89*8
    
alpha_ne_d <- 13*8
beta_ne_d <- 87*8
    
alpha_ne_nd <- 10*8
beta_ne_nd <- 90*8

RRvector <- c()

for (i in 1:10000) {
#sample from bias parameter distributions
s_e_d_sampled <- rbeta(1, alpha_e_d, beta_e_d)
s_e_nd_sampled <- rbeta(1, alpha_e_nd, beta_e_nd)
s_ne_d_sampled <- rbeta(1, alpha_ne_d, beta_ne_d)
s_ne_nd_sampled <- rbeta(1, alpha_ne_nd, beta_ne_nd)

# s_e_d_verif <- rbeta(1000, alpha_e_d, beta_e_d)
# s_e_d_quants <- quantile(s_e_d_verif, c(0.025, 0.50, 0.975))
# 
# s_e_nd_verif <- rbeta(1000, alpha_e_nd, beta_e_nd)
# s_e_nd_quants <- quantile(s_e_nd_verif, c(0.025, 0.50, 0.975))
# 
# s_ne_d_verif <- rbeta(1000, alpha_ne_d, beta_ne_d)
# s_ne_d_quants <- quantile(s_ne_d_verif, c(0.025, 0.50, 0.975))
# 
# s_ne_nd_verif <- rbeta(1000, alpha_ne_nd, beta_ne_nd)
# s_ne_nd_quants <- quantile(s_ne_nd_verif, c(0.025, 0.50, 0.975))

n <- 301943720

#perform adjustment for selection bias
A <- (0.813*n)/s_e_d_sampled
B <- (0.854*n)/s_e_nd_sampled
C <- (0.187*n)/s_ne_d_sampled
D <- (0.146*n)/s_ne_nd_sampled

#incorporate conventional random error
z <- rnorm(1,0,1)
logRR <- log((C/D)/(A/B))
stderror <- sqrt((1/A) + (1/B) + (1/C) + (1/D))
totalRR <- exp(logRR - (z*stderror))
RRvector <- append(RRvector,totalRR)
}

#plot histogram of adjusted ORs
RRquant <- quantile(RRvector, c(0.025, 0.50, 0.975), na.rm = TRUE)
hist(RRvector, breaks = 20, xlab = 'Odds Ratio', 
     ylab = 'Frequency', main = 'Odds Ratio Adjusted For Selection Bias')
abline(v= RRquant[2],col="blue",lwd=2)
abline(v= c(RRquant[1],RRquant[3]),col="blue",lwd=1, lty = 2)
abline(v=1.33,col="red",lwd=2)
abline(v=1.324,col="red",lwd=1, lty = 2)
abline(v=1.336,col="red",lwd=1, lty = 2)