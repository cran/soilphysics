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

    # first derivative
    dfun <- deriv3(fun2form(fun), "x", func = TRUE)
    if (attr(dfun(x.range[1]), "gradient") == attr(dfun(x.range[2]), 
        "gradient")) 
        stop("'fun' should not be a linar function of x!")

    # gradient and hessian
    x <- seq(x.range[1], x.range[2], length.out = 2000)
    gr <- attr(dfun(x), "gradient")
    he <- attr(dfun(x), "hessian")[,, "x"]
    k <- abs(he) / (1 + gr^2)^(3/2)
    mcp <- x[which.max(k)]

    # graph
    if (graph) {
        curve(fun, from = x.range[1], to = x.range[2], ...)
        lines(x = c(mcp, mcp, -9e9), 
           y = c(-9e9, fun(mcp), fun(mcp)), lty = 3)
        devAskNewPage(ask = TRUE)
        plot(x, k, type = "l", ...)
        lines(x = c(mcp, mcp, -9e9), 
           y = c(-9e9, max(k), max(k)), lty = 3)
    }

    # output
    out <- list(fun = deparse(fun)[2], 
       x0 = mcp, y0 = fun(mcp))
    class(out) <- "maxcurv"
    return(out)
}

# -------------------------------------------
# print method
print.maxcurv <- 
function (x, digits = 4L, quote = TRUE, ...) 
{
    cat("\n          Maximum curvature point \n",
        "\nfunction:", x$fun,
        "\ncritical x: ", x$x0,
        "\ncritical y: ", x$y0, "\n")
    invisible(x)
}