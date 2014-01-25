Sindex <-
function (theta_R, theta_S, alpha, n, m = 1 - 1/n, graph = TRUE, 
    ...) 
{
    x <- NULL
    h_i <- 1/alpha * (1/m)^(1/n)
    theta_i <- theta_R + (theta_S - theta_R) * (1 + 1/m)^(-m)
    if (graph) {
        curve(soilwater(x, theta_R, theta_S, alpha, n, m),
            xlab = "Matric potential", 
            ylab = "Soil water content",
            main = "Soil Water Retention Curve", ...)
        lines(x = c(h_i, h_i), y = c(-1e+09, theta_i), lty = 3)
        lines(x = c(-1e+09, h_i), y = c(theta_i, theta_i), lty = 3)
    }
    S <- abs(-n * (theta_S - theta_R) * (1 + 1/m)^(-(1 + m)))
    if (S >= 0.05) {
        clas <- "Very good"
    }
    else if (S < 0.05 & S >= 0.035) {
        clas <- "Good"
    }
    else if (S < 0.035 & S > 0.02) {
        clas <- "Poor"
    }
    else {
        clas <- "Very poor"
    }
    out <- list(h_i = h_i, theta_i = theta_i,
        S.index = S, PhysicalQuality = clas)
    cat("\n          The S Index \n")
    cat("\nh_i :", round(h_i, 4),
        "\ntheta_i :", round(theta_i, 4),
        "\n|S| :", round(S, 4),
        "\nSoil physical quality :", clas, "\n")
    class(out) <- "S.index"
    invisible(out)
}