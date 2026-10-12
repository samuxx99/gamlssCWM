#######################
## TUTORIALS SNIPPET ##
#######################


#----------------#
# ENRON tutorial #
#----------------#

# library("gamlssCWM")
devtools::load_all()
data("ENRON", package = "gamlssCWM")
attach(ENRON)
str(ENRON)

plot(ENRON, pch = 16, cex = .5)

fit <- gamlssCWM(
    formula       = UTILITY ~ poly(ENRON, 2, raw = TRUE),
    sigma_formula = ~ 1,
    data          = ENRON,
    famX          = gamlss.dist::PE(),
    famY          = gamlss.dist::NO(),
    init          = "kmeans",
    nstart        = 10,
    max_iter      = 100,
    tol           = 1e-4,
    k             = 3,
    verbose       = FALSE,
    seed          = 123
)
print(fit)

summary(fit)

plot(fit, dataX = ENRON, dataY = UTILITY, pch = 16, contours = TRUE, legend = "topleft")

persp(fit, dataX = ENRON, dataY = UTILITY, theta = 50, shade = 0.2)


#-------------------#
# ELECTION tutorial #
#-------------------#

data("ELECTION", package = "gamlssCWM")
str(ELECTION)
attach(ELECTION)
plot(age.tweetgenie ~ age, pch = 16, cex = .5)

k <- c(1, 2, 3, 4)
ICs <- vector(mode = "list", length = length(k))
fits <- vector(mode = "list", length = length(k))

for (j in 1:length(k)) {
    fits[[j]] <- gamlssCWM(
      formula       = age.tweetgenie ~ age,
      sigma_formula = ~ 1,
      data             = ELECTION,
      famX             = gamlss.dist::PO(),
      famY             = gamlss.dist::GA(),
      init             = "rand.soft",
      nstart           = 10,
      short_iterations = 5, 
      max_iter         = 550, 
      tol              = 1e-3,
      k                = k[j], 
      verbose          = FALSE,
      seed             = 123
    )
    ICs[[j]] <- getIC(fits[[j]])
}

par(mfrow=c(1,4), mar=c(3.5,2,2,1), oma=c(3,3,1,0), mgp=c(1.5,.4,0))

for(j in seq_along(fits)){
  plot(fits[[j]], dataX=age, dataY=age.tweetgenie,
       pch=16, cex=.4, xlab="", ylab="", main="")
  mtext(paste0("k = ", k[j]), 3, .5, font=2, cex=.9)
  mtext(paste0("(", letters[j], ")"), 1, 2.2, adj=.5, cex=.8)
}

mtext("age", 1, outer=TRUE, line=1.5)
mtext("age tweetgenie", 2, outer=TRUE, line=1.5)
par(mfrow=c(1,1))

ICs_mat <- t(as.data.frame(ICs)); n <- ncol(ICs_mat)
par(mar = c(4, 5, 2, 1))
matplot(k, ICs_mat, type = "b", lty = 1, pch = 1:n, col = 1:n, 
xlab = "k", ylab = "", las = 1, cex = .8, xaxt = "n")
best_k <- k[which.max(tabulate(apply(ICs_mat, 2, which.min), length(k)))]
axis(1, at = k); abline(v = best_k, lwd = 1.5, col = "red", lty = 2)
k_icl <- k[which.min(ICs_mat[, "ICL"])]
abline(v = k_icl, lwd = 1.5, col = "blue", lty = 2)
mtext("Information Criterion", 2, 4)
legend(1.6, 35900, colnames(ICs_mat),
 col = 1:n, lty = 1, pch = 1:n,
 bty = "n", cex = .8, ncol = 2, x.intersp = .3, y.intersp = .52,
 seg.len = .9)

fit <- gamlssCWM(
    formula       = age.tweetgenie ~ age,
    sigma_formula = ~ age,
    data          = ELECTION,
    famX          = gamlss.dist::PO(),
    famY          = gamlss.dist::GA(),
    init          = "manual",
    start.z       = fits[[3]]$posteriors,
    max_iter      = 120,
    tol           = 1e-3,
    k             = 3,
    verbose       = FALSE,
    seed          = 123
)

print(fit)

summary(fit)

plot(fit, dataX = age, dataY = age.tweetgenie, pch = 16, cex = .4, legend = "topleft")


sessionInfo()
