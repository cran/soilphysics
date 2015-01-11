\name{voidratio}
\alias{voidratio}

\title{ Void Ratio }

\description{ A function to calculate the soil void ratio. }

\usage{
voidratio(wetsoil, drysoil, diam.cylinder, height.cylinder,
    dens.particle, deformation)
}

\arguments{
  \item{wetsoil}{ the weight of wet soil. }
  \item{drysoil}{ the weight of dry soil. }
  \item{diam.cylinder}{ the diameter of the cylinder. }
  \item{height.cylinder}{ the heigth of the cylinder. }
  \item{dens.particle}{ the particle density. }
  \item{deformation}{ a numeric vector containing soil deformation values. }
}

\value{ A numeric vector with same length of \code{deformation} containig void ratio values. }

\author{ Anderson Rodrigo da Silva <anderson.agro@hotmail.com> }

\seealso{ \code{\link{sigmaP}} }

\examples{
def <- c(0.10, 0.12, 0.13, 0.15, 0.20, 0.26, 0.34, 0.43, 0.52)
pres <- c(1, 12.5, 25, 50, 100, 200, 400, 800, 1600)
VR <- voidratio(147.66, 118.13, 6.4, 3.4, 2.7, def)
VR
plot(VR ~ log10(pres), xaxt = "n", type = "b",
	ylab = "Void ratio",
	xlab = expression(Log[10]~Applied~stress~(kPa)),
	main = "Compression curve")
axis(side = 1, at = pretty(log10(pres)),
	labels = 10^pretty(log10(pres)))

# End (not run)
}
