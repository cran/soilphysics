\name{llwr}
\alias{llwr}

\title{Least Limiting Water Range}

\description{A graphical solution for the Least Limiting Water Range.}

\usage{
llwr(theta, psi, Bd, Pr, particle.density, air.porosity, 
	critical.PR, psi.FC, psi.WP, pars.water = NULL, 
	pars.busscher = NULL, graph = TRUE, 
	xlab = "Bulk density", ylab = "Water content", 
	main = "Least Limiting Water Range", ...)
}

\arguments{
  \item{theta}{a numeric vector containing values of volumetric water content.}
  \item{psi}{a numeric vector containing values of water potential (Psi).}
  \item{Bd}{a numeric vector containing values of dry bulk density.}
  \item{Pr}{a numeric vector containing values of penetration resistance.}
  \item{particle.density}{the value of the soil particle density.}
  \item{air.porosity}{the value of the limiting soil air porosity.}
  \item{critical.PR}{the value of the critical soil penetration resistance.}
  \item{psi.FC}{the value of the water potential at the field capacity.}
  \item{psi.WP}{the value of the water potential at the wilting point.}
  \item{pars.water}{optional; a numeric vector containing estimates of the three parameters of the model employed by 
	Silva et al. (1994) for	the functional relationship among \code{theta}, \code{Bd} and \code{psi}. 
	If \code{NULL} (default), \code{llwr()} estimates them using a Newton-Raphson algorithm.}
  \item{pars.busscher}{optional; a numeric vector containing estimates of the three parameters of the model proposed by
	 Busscher (1990) for the functional relationship among \code{Pr}, \code{theta} and \code{Bd}.
	If \code{NULL} (default), \code{llwr()} estimates them using a Newton-Raphson algorithm.}
  \item{graph}{logical; if TRUE (default) a graphical solution for the Least Limiting Water Range is plotted.}
  \item{xlab}{a title for the x axis; the default is "Bulk density".}
  \item{ylab}{a title for the x axis; the default is "Water content".}
  \item{main}{a main title for the graphic; the default is "Least Limiting Water Range"}
  \item{\dots}{further graphical arguments.}
}

\details{
The numeric vectors \code{theta}, \code{psi}, \code{Bd} and \code{Pr} are supposed to have the same length,
and their values should have appropriate unit of measurement. For fitting purposes, it is not advisable to use
vectors with less than five values.

The model for the water content as a function of the soil bulk density (\eqn{\rho}) and the water potential (\eqn{\psi}),
	 employed by Silva et al. (1994), is:

\deqn{\theta = exp(a + b\rho)\psi^c} 

}

\value{
A list of
  \item{limiting.theta }{a \eqn{n x 4} matrix containing the limiting values of water content for each input 
	value of bulk density at the air field porosity (\code{thetaAFP}), penetration resistance (\code{thetaPR}),
	field capacity (\code{thetaFC}) and wilting point (\code{thetaWP}).}
  \item{pars.water }{a numeric vector containing estimates of the three parameters of the model employed by 
	Silva et al. (1994) for	the functional relationship among \code{theta}, \code{Bd} and \code{psi}.}
  \item{pars.busscher }{a numeric vector containing estimates of the three parameters of the model proposed by
	 Busscher (1990) for the functional relationship among \code{Pr}, \code{theta} and \code{Bd}.}
}

\references{
Busscher, W. J. (1990). Adjustment of flat-tipped penetrometer resistance data to common water content. 
	\emph{Transactions of the ASAE}, 3:519-524.

Leao et al. (2005). An Algorithm for Calculating the Least Limiting Water Range of Soils. 
	\emph{Agronomy Journal}, 97:1210-1215.

Leao et al. (2006). Least limiting water range: A potential indicator of changes in near-surface 
	soil physical quality after the conversion of Brazilian Savanna into pasture. 
	\emph{Soil & Tillage Research}, 88:279-285.

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
with(skp1994,
	llwr(theta = Water, psi = Potential, Bd = BD, Pr = PR,
		particle.density = 2.65, air.porosity = 0.1,
		critical.PR = 2, psi.FC = 0.01, psi.WP = 1.5))

# Example 2 - specifying the parameters (Leao et al., 2005)
a <- c(-0.9175, -0.3027, -0.0835)  # model of water content
b <- c(0.0827, -1.6087, 3.0570)   # Busscher's model
with(skp1994,
	llwr(theta = Water, psi = Potential, Bd = BD, Pr = PR,
		particle.density = 2.65, air.porosity = 0.1,
		critical.PR = 2, psi.FC = 0.01, psi.WP = 1.5,
		pars.water = a, pars.busscher = b))

}
