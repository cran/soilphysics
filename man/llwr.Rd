\name{llwr}
\alias{llwr}

\title{Least Limiting Water Range}

\description{A graphical solution for the Least Limiting Water Range.}

\usage{
llwr(theta, psi, Bd, Pr, particle.density, air.porosity, 
    critical.PR, psi.FC, psi.WP, water.model = c("Silva", "Ross"), 
    pars.water = NULL, pars.busscher = NULL, graph = TRUE, graph2 = TRUE, 
    xlab = expression(Bulk~Density~(Mg~m^{-3})), 
    ylab = expression(theta~(m^{3}~m^{-3})), 
    main = "Least Limiting Water Range", ...) 
}

\arguments{
  \item{theta}{a numeric vector containing values of volumetric water content (\eqn{m^3~m^{-3}}).}
  \item{psi}{a numeric vector containing values of water potential (Psi, MPa, kPa, etc.).}
  \item{Bd}{a numeric vector containing values of dry bulk density (\eqn{Mg~m^{-3}}).}
  \item{Pr}{a numeric vector containing values of penetration resistance (MPa, kPa, etc.).}
  \item{particle.density}{the value of the soil particle density.}
  \item{air.porosity}{the value of the limiting soil air porosity.}
  \item{critical.PR}{the value of the critical soil penetration resistance.}
  \item{psi.FC}{the value of the water potential at the field capacity.}
  \item{psi.WP}{the value of the water potential at the wilting point.}
  \item{water.model}{a character; the model to be used for calculating the soil water content. It must be one of the
	two: \code{"Silva"} (default) or \code{"Ross"}. See details.}
  \item{pars.water}{optional; a numeric vector containing the estimates of the parameters of the soil water retention 
	model employed. If \code{NULL} (default), \code{llwr()} estimates them using a Newton-Raphson algorithm. See details.}
  \item{pars.busscher}{optional; a numeric vector containing estimates of the three parameters of the model proposed by
	 Busscher (1990) for the functional relationship among \code{Pr}, \code{theta} and \code{Bd}.
	If \code{NULL} (default), \code{llwr()} estimates them using a Newton-Raphson algorithm.}
  \item{graph}{logical; if TRUE (default) a graphical solution for the Least Limiting Water Range is plotted.}
  \item{graph2}{logical; if TRUE (default) a line of the Least Limiting Water Range as a function of bulk density is plotted. 
	If \code{graph = FALSE}, then \code{llwr()} will automatically consider \code{graph2 = FALSE} too.}
  \item{xlab}{a title for the x axis; the default is \eqn{Bulk~Density~(Mg~m^{-3})}.}
  \item{ylab}{a title for the y axis; the default is \eqn{\theta~(m^{3}~m^{-3})}.}
  \item{main}{a main title for the graphic; the default is "Least Limiting Water Range"}
  \item{\dots}{further graphical arguments.}
}

\details{
The numeric vectors \code{theta}, \code{psi}, \code{Bd} and \code{Pr} are supposed to have the same length,
and their values should have appropriate unit of measurement. For fitting purposes, it is not advisable to use
vectors with less than five values.

The model employed by Silva et al. (1994) for the soil water content (\eqn{\theta}) as a function of the soil bulk density (\eqn{\rho}) 
and the water potential (\eqn{\psi}) is:

\deqn{\theta = exp(a + b\rho)\psi^c}

The model proposed by Ross et al. (1991) for the soil water content (\eqn{\theta}) as a function of the water potential (\eqn{\psi}) is:

\deqn{\theta = a \psi^c}

}

\value{
A list of
  \item{limiting.theta }{a \eqn{n x 4} matrix containing the limiting values of water content for each input 
	value of bulk density at the air field porosity (\code{thetaAFP}), penetration resistance (\code{thetaPR}),
	field capacity (\code{thetaFC}) and wilting point (\code{thetaWP}).}
  \item{pars.water }{a "nls" object or a numeric vector containing estimates of the three parameters of the model employed by 
	Silva et al. (1994) for	the functional relationship among \code{theta}, \code{Bd} and \code{psi}.}
  \item{r.squared.water }{a "Rsq" object containing the pseudo and the adjusted R-squared for the water model.}
  \item{pars.busscher }{a "nls" object or a numeric vector containing estimates of the three parameters of the model proposed by
	 Busscher (1990) for the functional relationship among \code{Pr}, \code{theta} and \code{Bd}.}
  \item{r.squared.busscher }{a "Rsq" object containing the pseudo and the adjusted R-squared for the penetration resistance model.}
  \item{area }{numeric; the shaded (LLWR) area.}
}

\references{
Busscher, W. J. (1990). Adjustment of flat-tipped penetrometer resistance data to common water content. 
	\emph{Transactions of the ASAE}, 3:519-524.

Leao et al. (2005). An Algorithm for Calculating the Least Limiting Water Range of Soils. 
	\emph{Agronomy Journal}, 97:1210-1215.

Leao et al. (2006). Least limiting water range: A potential indicator of changes in near-surface 
	soil physical quality after the conversion of Brazilian Savanna into pasture. 
	\emph{Soil & Tillage Research}, 88:279-285.

Ross et al. (1991). Equation for extending water-retention curves to dryness.
	\emph{Soil Science Society of America Journal}, 55:923-927. 

Silva et al. (1994). Characterization of the least limiting water range of soils. 
	\emph{Soil Science Society of America Journal}, 58:1775-1781.
}

\author{
Anderson Rodrigo da Silva <anderson.agro@hotmail.com>
}

\seealso{
\code{\link{fitbusscher}}
}

\examples{
# Example 1 - part of the data set used by Leao et al. (2005)
data(skp1994)
ex1 <- with(skp1994,
	llwr(theta = Water, psi = Potential, Bd = BD, Pr = PR,
		particle.density = 2.65, air.porosity = 0.1,
		critical.PR = 2, psi.FC = 0.01, psi.WP = 1.5))
ex1

# Example 2 - specifying the parameters (Leao et al., 2005)
a <- c(-0.9175, -0.3027, -0.0835)  # Silva et al. model of water content
b <- c(0.0827, -1.6087, 3.0570)   # Busscher's model
ex2 <- with(skp1994,
	llwr(theta = Water, psi = Potential, Bd = BD, Pr = PR,
		particle.density = 2.65, air.porosity = 0.1,
		critical.PR = 2, psi.FC = 0.01, psi.WP = 1.5,
		pars.water = a, pars.busscher = b))
ex2

# End (not run)
}
