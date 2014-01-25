maxcurv <-
function (x.range, fun, graph = TRUE, ...) 
{
    stopifnot(is.atomic(x.range))
    if (!is.numeric(x.range)) 
        stop("'x.range' must be a numeric vector!")
    if (length(x.range) != 2) 
        stop("'x.range' must be a vector of length two!")
    if (diff(x.range) < 0) 
        stop("please, reorder 'x.range'.")
    if (!inherits(fun, "function")) 
        stop("'fun' must be a 'function' of x!")
    dfun <- deriv3(fun2form(fun), "x", func = TRUE)
    if (attr(dfun(x.range[1]), "gradient") == attr(dfun(x.range[2]), 
        "gradient")) 
        stop("'fun' should not be a linar function of x!")
    b <- lm(range(fun(x.range)) ~ x.range)$coef
    newx <- seq(x.range[1], x.range[2], length.out = 2000)
    newy <- fun(newx)
    if (fun(x.range[1]) > fun(x.range[2])) {
        si <- -1
    } else {
        si <- 1
    }
    pred.lm <- mean(fun(x.range)) - si * b[2] * mean(x.range) + 
        si * b[2] * newx
    delta <- NULL
    for (i in 1:length(newy)) {
        delta[i] <- abs(newy[i] - pred.lm[i])
    }
    ind <- which.max(delta)
    if (graph) {
        x <- NULL
        curve(fun(x), from = x.range[1], to = x.range[2],
            main = "Maximum curvature point", las = 1, ...)
        lines(x = c(newx[ind], newx[ind]),
            y = c(-1e+20, newy[ind]), lty = 3)
        lines(x = c(-1e+20, newx[ind]),
            y = c(newy[ind], newy[ind]), lty = 3)
    }
    cat("\n          Maximum curvature point \n",
        "\nfunction:", deparse(fun)[2],
        "\ncritical x: ", newx[ind],
        "\ncritical y: ", newy[ind], "\n")
    out <- list(x0 = newx[ind], y0 = newy[ind])
    class(out) <- "maxcurv"
    invisible(out)
}
