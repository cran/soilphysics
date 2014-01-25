sigmaP <- 
function (voidratio, stress, n4VCL = 2, method = c("casagrande", 
    "VCLzero", "reg1", "reg2", "reg3", "reg4", "pacheco"), mcp = NULL, 
    graph = TRUE, ...) 
{
    if (length(voidratio) != length(stress)) 
        stop("incompatible dimensions!")
    stopifnot(is.numeric(voidratio))
    stopifnot(is.numeric(stress))
    method <- match.arg(method)
    x <- NULL
    xy <- sortedXyData(log10(stress), voidratio)
    if (n4VCL < 2) 
        stop("'n4VCL' must be a positive integer >= 2!")
    b. <- coef(lm(y ~ x, data = tail(xy, n4VCL)))
    f.plot <- function(...) {
        plot(y ~ x, data = xy, xaxt = "n", type = "b", las = 1, 
            ylab = "Void ratio", xlab = expression(Log[10] ~ 
                Applied ~ stress), main = "Compression curve", 
            ...)
        xval <- pretty(par("usr")[1:2])
        axis(side = 1, at = xval, labels = 10^xval)
    }
    if (method == "casagrande") {
        if (is.null(mcp)) 
            stop("'mcp' must be defined for using 'casagrande' method!")
        fit <- lm(y ~ x + I(x^2) + I(x^3) + I(x^4), data = xy)
        est <- coef(fit)
        Xmax <- mcp
        Ymax <- predict(fit, newdata = data.frame(x = Xmax))
        b1.tan <- est[2] + 2 * est[3] * Xmax + 3 * est[4] * Xmax^2 + 
            4 * est[5] * Xmax^3
        b0.tan <- Ymax - b1.tan * Xmax
        b0.bis <- Ymax - b1.tan/2 * Xmax
        b1.bis <- b1.tan/2
        x0 <- as.vector((b.[1] - b0.bis)/(b1.bis - b.[2]))
        y0 <- b.[1] + b.[2] * x0
        if (graph) {
            f.plot(...)
            abline(b., lty = 3)
            lines(x = c(Xmax, Xmax), y = c(-1e+09, Ymax), lty = 3)
            abline(b0.tan, b1.tan, lty = 3)
            lines(x = c(Xmax, 1e+09), y = c(Ymax, Ymax), lty = 3)
            curve(b0.bis + b1.bis * x, from = Xmax, to = 1e+09, 
                add = TRUE, lty = 3)
            lines(x = c(x0, x0), y = c(-1e+09, y0), col = "red")
        }
    }
    else if (method == "VCLzero") {
        b <- c(xy$y[1], 0)
        x0 <- as.vector((b.[1] - b[1])/(b[2] - b.[2]))
        if (graph) {
            f.plot(...)
            abline(b., lty = 3)
            abline(b, lty = 3)
            lines(x = c(x0, x0), y = c(-1e+09, b[1] + b[2] * 
                x0), lty = 1, col = "red")
        }
    }
    else if (method == "reg1") {
        b <- coef(lm(y ~ x, data = head(xy, 2)))
        x0 <- as.vector((b.[1] - b[1])/(b[2] - b.[2]))
        if (graph) {
            f.plot(...)
            abline(b., lty = 3)
            abline(b, lty = 3)
            lines(x = c(x0, x0), y = c(-1e+09, b[1] + b[2] * 
                x0), lty = 1, col = "red")
        }
    }
    else if (method == "reg2") {
        b <- coef(lm(y ~ x, data = head(xy, 3)))
        x0 <- as.vector((b.[1] - b[1])/(b[2] - b.[2]))
        if (graph) {
            f.plot(...)
            abline(b., lty = 3)
            abline(b, lty = 3)
            lines(x = c(x0, x0), y = c(-1e+09, b[1] + b[2] * 
                x0), lty = 1, col = "red")
        }
    }
    else if (method == "reg3") {
        b <- coef(lm(y ~ x, data = head(xy, 4)))
        x0 <- as.vector((b.[1] - b[1])/(b[2] - b.[2]))
        if (graph) {
            f.plot(...)
            abline(b., lty = 3)
            abline(b, lty = 3)
            lines(x = c(x0, x0), y = c(-1e+09, b[1] + b[2] * 
                x0), lty = 1, col = "red")
        }
    }
    else if (method == "reg4") {
        b <- coef(lm(y ~ x, data = head(xy, 5)))
        x0 <- as.vector((b.[1] - b[1])/(b[2] - b.[2]))
        if (graph) {
            f.plot(...)
            abline(b., lty = 3)
            abline(b, lty = 3)
            lines(x = c(x0, x0), y = c(-1e+09, b[1] + b[2] * 
                x0), lty = 1, col = "red")
        }
    }
    else if (method == "pacheco") {
        fit <- lm(y ~ x + I(x^2) + I(x^3) + I(x^4), data = xy)
        est <- coef(fit)
        x0. <- (xy$y[1] - b.[1])/b.[2]
        y0. <- predict(fit, newdata = data.frame(x = x0.))
        x0 <- (y0. - b.[1])/b.[2]
        if (graph) {
            f.plot(...)
            abline(b., lty = 3)
            abline(h = xy$y[1], lty = 3)
            lines(x = c(x0., x0.), y = c(xy$y[1], y0.), lty = 3, 
                col = "red")
            lines(x = c(x0., x0), y = c(y0., y0.), lty = 3, col = "red")
            lines(x = c(x0, x0), y = c(-1e+09, y0.), lty = 1, 
                col = "red")
        }
    }
    sigmaP <- as.vector(10^x0)
    cat("\nPreconsolidation stress:", sigmaP, "\n")
    cat("Method:", method, "\n")
    invisible(sigmaP)
}
