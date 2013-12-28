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
def <- c(0.1, 17, 30, 51, 97, 145, 201, 260, 325) / 100
pres <- c(1, 12.2, 25, 50, 100, 200, 400, 800, 1600)
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
