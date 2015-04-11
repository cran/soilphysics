\name{iwc}
\alias{iwc}

\title{ Integral Water Capacity }

\description{
Quantifying the soil water availability for plants through the IWC approach.
The theory was based on the work of Asgarzadeh et al. (2014), using the van Genuchten-Mualem Model 
for estimation of the water retention curve and a simple power model for penetration resistance.
}

\usage{
iwc(theta_R, theta_S, alpha, n, a, b, hos = 0, 
	graph = TRUE, xlab = "Matric potential (hPa)", 
	ylab = "Water content", ...)
}

\arguments{
  \item{theta_R}{ the residual water content; a numeric parameter from van Genuchten's model; see details.}
  \item{theta_S}{ the water content at saturation; a numeric parameter from van Genuchten's model; see details.}
  \item{alpha}{ a scale parameter from van Genuchten's model; see details.}
  \item{n}{ a shape parameter from van Genuchten's model; see details.}
  \item{a}{ a parameter of the soil penetration resistance model; see details.}
  \item{b}{ a parameter of the soil penetration resistance model; see details.}
  \item{hos}{ optional; the value of osmotic head at which the soil pores become saturated. Used only if one is 
	concerned about the salinity effects on the water available for plants. Default is zero. 
	See Groenevelt et al. (2004) for more details.}
  \item{graph}{logical; if TRUE (default), graphics for both dry and wet range are built.}
  \item{xlab}{ a label for x-axis.}
  \item{ylab}{ a label for y-axis.}
  \item{\dots}{ further graphical parameters.}
}

\details{
The parameters of the van Genuchten-Mualem Model can be estimated through the function \code{fitsoilwater()}. 
The soil penetration resistance model: \eqn{PR = a*h^b}.
}

\value{
A list of
   \item{IWC}{ integral water capacity}
   \item{EI}{ integral energy calculation}
}

\references{
Asgarzadeh, H.; Mosaddeghi, M.R.; Nikbakht, A.M. (2014) SAWCal: A user-friendly program for 
	calculating soil available water quantities and physical quality indices. 
	\emph{Computers and Electronics in Agriculture}, 109:86--93.

Groenevelt, P.H.; Grant, C.D.; Semetsa, S. (2001) A new procedure to determine soil water availability.
	\emph{Australian Journal Soil Research}, 39:577--598.

Groenevelt, P.H., Grant, C.D., Murray, R.S. (2004) On water availability in saline soils. 
	\emph{Australian Journal Soil Research}, 42:833--840.
}

\author{
Anderson Rodrigo da Silva <anderson.agro@hotmail.com>
}

\seealso{
\code{\link{soilwater}}, \code{\link{fitsoilwater}}, \code{\link{llwr}}
}

\examples{
# example 1 (Fig 1b, Asgarzadeh et al., 2014)
iwc(theta_R = 0.0160, theta_S = 0.4828, alpha = 0.0471, n = 1.2982, 
	a = 0.2038, b = 0.2558, graph = TRUE) 

# example 2 (Table 1, Asgarzadeh et al., 2014)
iwc(theta_R = 0.166, theta_S = 0.569, alpha = 0.029, n = 1.308, 
	a = 0.203, b = 0.256, graph = TRUE) 

# example 3: evaluating the salinity effect
iwc(theta_R = 0.166, theta_S = 0.569, alpha = 0.029, n = 1.308, 
	a = 0.203, b = 0.256, hos = 200, graph = TRUE)

# End (Not run)
}
