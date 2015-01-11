llwr <- 
function (theta, psi, Bd, Pr, particle.density, air.porosity, 
    critical.PR, psi.FC, psi.WP, water.model = c("Silva", "Ross"), 
    pars.water = NULL, pars.busscher = NULL, graph = TRUE, graph2 = TRUE, 
    xlab = expression(Bulk~Density~(Mg~m^{-3})), 
    ylab = expression(theta~(m^{3}~m^{-3})), 
    main = "Least Limiting Water Range", ...) 
{
    n <- length(theta)
    if (length(psi) != n || length(Bd) != n || length(Pr) != 
        n) 
        stop("incompatible dimensions!")
    dat <- cbind(theta, psi, Bd, Pr)
    if (!is.numeric(dat)) 
        stop("non-numeric data!")
    limits <- c(particle.density, air.porosity, critical.PR, 
        psi.FC, psi.WP)
    if (length(limits) > 5) 
        stop("each limiting value must be a single value!")
    water.model <- match.arg(water.model)
    if (is.null(pars.water)) {
        if (water.model == "Silva") {
            fit1. <- lm(log(theta) ~ Bd + log(psi))
            a. <- coef(fit1.)
            fit1 <- nls(theta ~ exp(a + b * Bd) * psi^c, start = list(a = a.[1], 
                b = a.[2], c = a.[3]))
            a <- coef(fit1)
        }
        else {
            fit1. <- lm(log(theta) ~ log(psi))
            a. <- coef(fit1.)
            fit1 <- nls(theta ~ a * psi^b, start = list(a = exp(a.[1]), 
                b = a.[2]))
            a <- c(log(coef(fit1)[1]), 0, coef(fit1)[2])
        }
        rsq1 <- Rsq(fit1)
    }
    else if (length(pars.water) == 3) {
        a <- pars.water
    }
    if (is.null(pars.busscher)) {
        fit2 <- fitbusscher(Pr, theta, Bd)
        b <- coef(fit2)
        rsq2 <- Rsq(fit2)
    }
    else if (length(pars.busscher) == 3) {
        b <- pars.busscher
    }
    Dp <- particle.density
    thetaAFP <- 1 - Bd/Dp - air.porosity
    PRc <- critical.PR
    thetaPR <- (PRc/(b[1] * Bd^b[3]))^(1/b[2])
    thetaFC <- exp(a[1] + a[2] * Bd) * psi.FC^a[3]
    thetaWP <- exp(a[1] + a[2] * Bd) * psi.WP^a[3]
    theta. <- cbind(thetaAFP, thetaPR, thetaFC, thetaWP)
    x. <- seq(range(Bd)[1], range(Bd)[2], length.out = 100)
    thetaAFP. <- 1 - x./Dp - air.porosity
    thetaPR. <- (PRc/(b[1] * x.^b[3]))^(1/b[2])
    mi <- which.min((thetaAFP. - thetaPR.)^2)
    x <- seq(range(Bd)[1], x.[mi], length.out = 100)
    yUp. <- cbind(1 - x/Dp - air.porosity, exp(a[1] + a[2] * 
        x) * psi.FC^a[3])
    yUp <- apply(yUp., 1, min)
    yLow. <- cbind((PRc/(b[1] * x^b[3]))^(1/b[2]), exp(a[1] + 
        a[2] * x) * psi.WP^a[3])
    yLow <- apply(yLow., 1, max)
    if (graph) {
        plot(range(Bd), range(theta.) * c(1, 1.1), pch = "", 
            xlab = xlab, ylab = ylab, main = main, ...)
        polygon(c(x[1], x, x[100]), c(yLow[1], yUp, yLow[100]), 
            col = gray(0.8, alpha = 1), border = FALSE)
        polygon(c(x[1], x, x[100]), c(yUp[1], yLow, yUp[100]), 
            col = gray(0.8, alpha = 1), border = FALSE)
        curve(1 - x/Dp - air.porosity, add = TRUE)
        curve(exp(a[1] + a[2] * x) * psi.FC^a[3], add = TRUE)
        curve((PRc/(b[1] * x^b[3]))^(1/b[2]), add = TRUE, col = "blue", 
            lty = 2)
        curve(exp(a[1] + a[2] * x) * psi.WP^a[3], add = TRUE, 
            col = "blue", lty = 2)
        points(Bd, thetaAFP)
        points(Bd, thetaFC, pch = 16)
        points(Bd, thetaPR, col = "blue")
        points(Bd, thetaWP, col = "blue", pch = 16)
        tex <- c(expression(theta[AFP]), expression(theta[FC]), 
            expression(theta[PR]), expression(theta[WP]))
        legend("topright", tex, lty = c(1, 1, 2, 2), col = c(1, 
            1, 4, 4), pch = c(1, 16, 1, 16), cex = 0.8, bg = "white")
        if (graph2) {
            dev.new(width = 3, height = 3)
            plot(x, yUp - yLow, type = "l", xlab = xlab, ylab = "LLWR", 
                ...)
        }
    }
    area <- trapez(x, yUp) - trapez(x, yLow)
    out <- list(limiting.theta = theta., 
        pars.water = if (is.null(pars.water)) fit1 else a, 
        r.squared.water = if (is.null(pars.water)) rsq1,
        pars.busscher = if (is.null(pars.busscher)) fit2 else b,
        r.squared.busscher = if (is.null(pars.busscher)) rsq2, 
        area = area)
    class(out) <- "llwr"
    return(out)
}

# -----------------------------------------------
# print method
print.llwr <- 
function (x, ...) 
{
    cat("\n          Least Limiting Water Range \n")
    cat("\n---------- \nLimiting theta (6 first rows):\n")
    print(head(x$limiting.theta))

    cat("\n\n---------- \nEstimates of the soil water content model: \n")
    if (inherits(x$pars.water, "nls")) {
        print(summary(x$pars.water))
        rsq1 <- x$r.squared.water$pseudo.R.squared
        cat("pseudo R-squared:", rsq1, "\n")
        adj.rsq1 <- x$r.squared.water$adj.R.squared
        cat("adjusted R-squared:", adj.rsq1, "\n")
    }
    else if (inherits(x$pars.water, "numeric")) {
        print(x$pars.water)
    }

    cat("\n\n---------- \nEstimates of the soil penetration resistance model (Busscher, 1990): \n")
    if (inherits(x$pars.busscher, "nls")) {
        print(summary(x$pars.busscher))
        rsq2 <- x$r.squared.busscher$pseudo.R.squared
        cat("pseudo R-squared:", rsq2, "\n")
        adj.rsq2 <- x$r.squared.busscher$adj.R.squared
        cat("adjusted R-squared:", adj.rsq2, "\n")
    }
    else if (inherits(x$pars.busscher, "numeric")) {
        print(x$pars.busscher)
    }

    cat("\n\n---------- \nShaded area:", x$area, "\n")
    invisible(x)
}

# ---------------------------------------------------
# numerical integration: trapezoidal rule
trapez <- 
function(x, y)
{
   if (length(x) != length(y))
      stop("incompatible dimensions!")
   xy <- sortedXyData(x, y)
   x <- xy[["x"]]
   y <- xy[["y"]]
   n <- nrow(xy)
   out <- 0.5 * sum((x[2:n] - x[2:n - 1]) * 
      (y[2:n] + y[2:n - 1]))
   return(out)
}
