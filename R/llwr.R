llwr <- function(theta, psi, Bd, Pr,
	particle.density, air.porosity, critical.PR,
	psi.FC, psi.WP,
	pars.water = NULL, pars.busscher = NULL,
	graph = TRUE, xlab = "Bulk density", ylab = "Water content", 
	main = "Least Limiting Water Range", ...)
{
   n <- length(theta)
   if (length(psi) != n || length(Bd) != n ||
      length(Pr) != n)
      stop("incompatible dimensions!")
   dat <- cbind(theta, psi, Bd, Pr)
   if (!is.numeric(dat))
      stop("non-numeric data!")
   limits <- c(particle.density, air.porosity,
      critical.PR, psi.FC, psi.WP)
   if (length(limits) > 5)
      stop("each limiting value must be a single value!")

   if (is.null(pars.water)) {
      fit1. <- lm(log(theta) ~ Bd + log(psi))
      a. <- coef(fit1.)
      fit1 <- nls(theta ~ exp(a + b*Bd) * psi^c,
         start = list(a = a.[1], b = a.[2], c = a.[3]))
      a <- coef(fit1)
   } else if (length(pars.water) == 3) {
      a <- pars.water
   }
   if (is.null(pars.busscher)) {
      fit2 <- fitbusscher(Pr, theta, Bd)
      b <- coef(fit2)
   } else if (length(pars.busscher) == 3) {
      b <- pars.busscher
   }
   names(b) <- letters[4:6]

   Dp <- particle.density
   thetaAFP <- 1 - Bd/Dp - air.porosity
   PRc <- critical.PR
   thetaPR <- (PRc / (b[1] * Bd ^ b[3])) ^ (1/b[2])
   thetaFC <- exp(a[1] + a[2] * Bd) * psi.FC ^ a[3]
   thetaWP <- exp(a[1] + a[2] * Bd) * psi.WP ^ a[3]
   theta. <- cbind(thetaAFP, thetaPR, thetaFC, thetaWP)

   if (graph) {
      plot(range(Bd), range(theta.)*c(1, 1.1),
         pch = "",
         xlab = xlab,
         ylab = ylab,
         main = main, ...)
      x. <- seq(range(Bd)[1], range(Bd)[2], length.out = 100)
      thetaAFP. <- 1 - x./Dp - air.porosity
      thetaPR. <- (PRc / (b[1] * x. ^ b[3])) ^ (1/b[2])
      mi <- which.min((thetaAFP. - thetaPR.)^2)
      x <- seq(range(Bd)[1], x.[mi], length.out = 100)
      yUp. <- cbind(1 - x/Dp - air.porosity,
         exp(a[1] + a[2] * x) * psi.FC ^ a[3])
      yUp <- apply(yUp., 1, min)
      yLow. <- cbind((PRc / (b[1] * x ^ b[3])) ^ (1/b[2]),
         exp(a[1] + a[2] * x) * psi.WP ^ a[3])
      yLow <- apply(yLow., 1, max)
      polygon(c(x[1], x, x[100]), c(yLow[1], yUp, yLow[100]),
         col = gray(0.8, alpha = 1), border = FALSE)
      polygon(c(x[1], x, x[100]), c(yUp[1], yLow, yUp[100]),
         col = gray(0.8, alpha = 1), border = FALSE)
      curve(1 - x/Dp - air.porosity, add = TRUE) # AFP
      curve(exp(a[1] + a[2] * x) * psi.FC ^ a[3], add = TRUE) # FC
      curve((PRc / (b[1] * x ^ b[3])) ^ (1/b[2]),
         add = TRUE, col = "blue", lty = 2) # PR
      curve(exp(a[1] + a[2] * x) * psi.WP ^ a[3],
         add = TRUE, col = "blue", lty = 2) # WP
      points(Bd, thetaAFP)
      points(Bd, thetaFC, pch = 16)
      points(Bd, thetaPR, col = "blue")
      points(Bd, thetaWP, col = "blue", pch = 16)
      tex <- c(expression(theta[AFP]), expression(theta[FC]),
         expression(theta[PR]), expression(theta[WP]))
      legend('topright', tex, lty = c(1, 1, 2, 2),
         col = c(1, 1, 4, 4), pch = c(1, 16, 1, 16),
         cex = 0.8, bg = "white")
   }

   out <- list(limiting.theta = theta., pars.water = a,
      pars.busscher = b)
   class(out) <- "llwr"
   return(out)
}