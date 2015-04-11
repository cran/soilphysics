iwc <- 
function(theta_R, theta_S, alpha, n, a, b, hos = 0, 
      graph = TRUE, xlab = "Matric potential (hPa)", 
      ylab = "Water content", ...)
{
    stopifnot(hos >= 0)

    # Ch (van Genuchten's first derivative)
    x <- NULL
    der. <- function(x, theta_R, theta_S, alpha, n) 
    {
        out <- abs((theta_S - theta_R) * 
           (1 + (alpha * x)^n)^(1/n - 1) * 
           (1/n - 1) * (alpha * x)^n * 
           (n/(x * (1 + (alpha * x)^n))))
        return(out)
    }

    # weight limiting function for salinity
    ws <- function(x, theta_R, theta_S, alpha, n)
    {
        theta <- soilwater(x, theta_R, theta_S, alpha, n)
        out <- ( 1 + hos * theta_S * theta^-2 *
           der.(x, theta_R, theta_S, alpha, n) )^-1
        return(out)
    }

    # Es
    der <- function(x, theta_R, theta_S, alpha, n)
    {
        out <- ws(x, theta_R, theta_S, alpha, n) *
           der.(x, theta_R, theta_S, alpha, n)
        return(out)
    }

    # Hydraulic conductivity limiting function
    Kr <- function(x, theta_R, theta_S, alpha, n) 
    {
       out <- (1 - (alpha*x)^(n-1) * 
          (1 + (alpha*x)^n)^(1/n - 1) )^2 *
          (1 + (alpha*x)^n)^( (1 - n) / 2*n)
       return(out)
    }

    # Kr weight function
    wK <- function(x, theta_R, theta_S, alpha, n) 
    {
        out <- ( Kr(x = 330, theta_R, theta_S, alpha, n) / 
           Kr(x, theta_R, theta_S, alpha, n) )^0.08
        return(out)
    }

    EK <- function(x, theta_R, theta_S, alpha, n) 
    {
        out <- wK(x, theta_R, theta_S, alpha, n) * 
           der(x, theta_R, theta_S, alpha, n)
        return(out)
    }

    # inverted van Genuchten's function for potential determination
    potential <- function(theta, theta_R, theta_S, alpha, n)
    {
        m <- 1 - 1/n
        out <- (( (theta_S - theta_R) / 
           (theta - theta_R) )^(1/m) - 1 )^(1/n) * 1/alpha
        return(out)
    }
    h10 <- potential(theta_S - 0.10, theta_R, theta_S, alpha, n)
    h15 <- potential(theta_S - 0.15, theta_R, theta_S, alpha, n)

    # weight function for air porosity limiting function
    wa <- function(x) 
    {
        A <- 1 / ( log10( h15 / h10 ) )
        out <- A * log10(x / h10 )
        for(i in 1:length(x)) {
        if (out[i] < 0) out[i] <- 0
        if (out[i] > 1) out[i] <- 1
        }
        return(out)
    }

    EKa <- function(x, theta_R, theta_S, alpha, n) 
    {
        out <- wK(x, theta_R, theta_S, alpha, n) * 
           wa(x) * 
           der(x, theta_R, theta_S, alpha, n)
        return(out)
    }

    # weight function for PR limiting function
    wR <- function(x, a, b)
    {
        out <- 2.5 - a * x ^ b
        return(out)
    }

    ER <- function(x, a, b, theta_R, theta_S, alpha, n) 
    {
        out <- wR(x, a, b) *
           der(x, theta_R, theta_S, alpha, n)
        return(out)
    }

    # integration interval for PR = 1.5 and PR = 2.5 MPa
    h1.5 <- (1.5 / a) ^ (1/b)
    h2.5 <- (2.5 / a) ^ (1/b)

    # dry hydraulic conductivity
    Kdry <- Kr(x = seq(12000, 15000, by = 50), 
        theta_R, theta_S, alpha, n)
    fit2 <- fitsoilwater4(theta = Kdry, 
        psi = seq(12000, 15000, by = 50), 
        model = "Ross")

    # weight limiting function for Kdry
    wKdry <- function(x, d)
    {
        out <- (12000 / x) ^ (-d)
        return(out)
    }

    ERKdry <- function(x, a, b, d, theta_R, theta_S, alpha, n) 
    {
        out <- wR(x, a, b) * wKdry(x, d) *
           der(x, theta_R, theta_S, alpha, n)
        return(out)
    }

    # ----------------------------------------------
    # IWC integrations
    # EKa
    ave1 <- trapez(x = seq(h10, h15, length = 1000), 
        y = EKa(x = seq(h10, h15, length = 1000), 
           theta_R, theta_S, alpha, n))

    # EK
    ave2 <- trapez(x = seq(h15, 330, length = 1000), 
        y = EK(x = seq(h15, 330, length = 1000), 
           theta_R, theta_S, alpha, n))

    # Ch
    ave3 <- trapez(x = seq(330, h1.5, length = 1000), 
        y = der(x = seq(330, h1.5, length = 1000), 
           theta_R, theta_S, alpha, n))

    # ER
    if (h2.5 > 15000) h2.5. <- 15000 else h2.5. <- h2.5
    ave4 <- trapez(x = seq(h1.5, h2.5., length = 1000), 
        y = ER(x = seq(h1.5, h2.5., length = 1000), 
            a, b, theta_R, theta_S, alpha, n))

    # ERKdry
    ave5 <- trapez(x = seq(12000, 15000, length = 1000), 
        y = ERKdry(x = seq(12000, 15000, length = 1000), 
            a, b,	d = coef(fit2)[2], theta_R, theta_S, alpha, n))

    # ------------------------------------------------
    # integral energy calculation
    # hEKa
    hEKa <- function(x, theta_R, theta_S, alpha, n) 
    {
        out <- x * EKa(x, theta_R, theta_S, alpha, n)
        return(out)
    }
    ave1h <- trapez(x = seq(h10, h15, length = 1000), 
        y = hEKa(x = seq(h10, h15, length = 1000), 
           theta_R, theta_S, alpha, n))

    # hEK
    hEK <- function(x, theta_R, theta_S, alpha, n) 
    {
        out <- x * EK(x, theta_R, theta_S, alpha, n)
        return(out)
    }
    ave2h <- trapez(x = seq(h15, 330, length = 1000), 
        y = hEK(x = seq(h15, 330, length = 1000), 
           theta_R, theta_S, alpha, n))

    # hCh
    hder <- function(x, theta_R, theta_S, alpha, n)
    {
        out <- x * der(x, theta_R, theta_S, alpha, n)
        return(out)
    }
    ave3h <- trapez(x = seq(330, h1.5, length = 1000), 
        y = hder(x = seq(330, h1.5, length = 1000), 
           theta_R, theta_S, alpha, n))

    # hER
    hER <- function(x, a, b, theta_R, theta_S, alpha, n) 
    {
        out <- x * ER(x, a, b, theta_R, theta_S, alpha, n)
        return(out)
    }
    ave4h <- trapez(x = seq(h1.5, h2.5., length = 1000), 
        y = hER(x = seq(h1.5, h2.5., length = 1000), 
            a, b, theta_R, theta_S, alpha, n))

    # hERKdry
    hERKdry <- function(x, a, b, d, theta_R, theta_S, alpha, n) 
    {
        out <- x * ERKdry(x, a, b,	
            d, theta_R, theta_S, alpha, n)
        return(out)
    }
    ave5h <- trapez(x = seq(12000, 15000, length = 1000), 
        y = hERKdry(x = seq(12000, 15000, length = 1000), 
            a, b,	d = coef(fit2)[2], theta_R, theta_S, alpha, n))

    # -------------------------------------------------
    # graphic
    if (graph) {
        curve(der(x, theta_R, theta_S, alpha, n), 
            from = 0, to = 350, 
            xlab = xlab, ylab = ylab, 
            main = "Wet range", ...)
        curve(EK(x, theta_R, theta_S, alpha, n), 
           add = TRUE, lty = 3, col = 4)
        curve(EKa(x, theta_R, theta_S, alpha, n), 
           add = TRUE, lty = 2, col = 2)
        legend("topright",  c(expression(C(h)), 
             expression(E[R](h)), expression(E[RKdry](h))), 
             lty = c(1, 3, 2), col = c(1, 4, 2), 
             bg = "lightyellow")
    
        dev.new(width = 5, height = 5)
        curve(der(x, theta_R, theta_S, alpha, n), 
            from = h1.5, to = 15000, 
            xlab = xlab, ylab = ylab, 
            main = "Dry range", ...)
        curve(ER(x, a, b, theta_R, theta_S, alpha, n), 
            add = TRUE, from = h1.5, to = h2.5., 
            lty = 3, col = 4)
        curve(ERKdry(x, a, b, d = coef(fit2)[2], 
            theta_R, theta_S, alpha, n), 
            add = TRUE, from = 12000, to = 15000, 
            lty = 2, col = 2)
        legend("topright",  c(expression(C(h)), 
             expression(E[R](h)), expression(E[RKdry](h))), 
             lty = c(1, 3, 2), col = c(1, 4, 2), 
             bg = "lightyellow")

        if (hos != 0) {
           dev.new(width = 5, height = 5)
           curve(der.(x, theta_R, theta_S, alpha, n), 
              from = 1, to = 15000, log = "x",
              xlab = xlab, ylab = ylab, 
              main = "Salinity effect", ...)
           curve(der(x, theta_R, theta_S, alpha, n), 
               add = TRUE, from = 1, to = 15000, 
               log = "x", col = 2, lty = 2)
           legend("topright", c(expression(C(h)), 
               expression(E[S](h))), lty = c(1, 2), 
               col = c(1, 2), bg = "lightyellow")
        }        
    }

    # ----------------------------------------
    # output
    IWC <- ave1 + ave2 + ave3 + ave4 + ave5
    EI <- (ave1h + ave2h + ave3h + ave4h + ave5h) / (10 * IWC)
    return(list(IWC = IWC, EI = EI))
}


# --------------------------------------------
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
