## ----setup, echo = FALSE-------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----eval=FALSE----------------------------------------------------------
#  install.packages("soilphysics")

## ---- eval=FALSE---------------------------------------------------------
#  install.packages("devtools")
#  devtools::install_github("arsilva87/soilphysics")

## ----warning=FALSE-------------------------------------------------------
library(soilphysics) 

## ------------------------------------------------------------------------
stress <- stressTraffic(inflation.pressure=200, 
	   recommended.pressure=200, 
	   tyre.diameter=1.8, 
	   tyre.width=0.4, 
	   wheel.load=4000, 
	   conc.factor=c(4,5,5,5,5,5),
           layers=c(0.05,0.1,0.3,0.5,0.7,1), 
	   plot.contact.area = TRUE)


## ------------------------------------------------------------------------
soilDeformation(stress = 300,
                p.density = 2.67,
	        iBD = 1.55, 
                N = 1.9392, 
                CI = 0.06037, 
                k = 0.00608, 
                k2 = 0.01916,
                m = 1.3,graph=TRUE,ylim=c(1.4,2.0))


## ------------------------------------------------------------------------
mois <- c(0.083, 0.092, 0.108, 0.126, 0.135)
bulk <- c(1.86, 1.92, 1.95, 1.90, 1.87)
criticalmoisture(theta = mois, Bd = bulk)

## ------------------------------------------------------------------------
iwc(theta_R = 0.166, theta_S = 0.569, alpha = 0.029, n = 1.308, 
    a = 0.203, b = 0.256, hos = 200, graph = TRUE)

## ------------------------------------------------------------------------
# Usage
data(skp1994)
with(skp1994,
	llwr(theta = W, h = h, Bd = BD, Pr = PR,
		particle.density = 2.65, air = 0.1,
		critical.PR = 2, h.FC = 100, h.WP = 15000))

## ------------------------------------------------------------------------
par(mfrow=c(1,2))
llwr_llmpr(thetaR=0.1180, thetaS=0.36, alpha=0.133, n=1.30, 
        d=0.005, e=-2.93, f=3.54, PD=2.65,
        critical.PR=4, h.FC=100, h.PWP=15000, air.porosity=0.1,
        labels=c("AFP", "FC","PWP", "PR"),
        graph1=TRUE,graph2=FALSE, ylab=expression(LLMPR~(hPa)), ylim=c(15000,1))
mtext(expression("Bulk density"~(Mg~m^-3)),1,line=2.2, cex=0.8)

llwr_llmpr(thetaR=0.1180, thetaS=0.36, alpha=0.133, n=1.30, 
        d=0.005, e=-2.93, f=3.54, PD=2.65,
        critical.PR=4, h.FC=100, h.PWP=15000, air.porosity=0.1,
        labels=c("AFP", "FC","PWP", "PR"),
        graph1=FALSE,graph2=TRUE, ylab=expression(LLMPR~(hPa)), ylim=c(0.1,0.5))
mtext(expression("Bulk density"~(Mg~m^-3)),1,line=2.2, cex=0.8)

## ------------------------------------------------------------------------
pres <- c(1, 12.5, 25, 50, 100, 200, 400, 800, 1600)
VR <- c(0.846, 0.829, 0.820, 0.802, 0.767, 0.717, 0.660, 0.595, 0.532)
sigmaP(VR, pres, method = "casagrande", n4VCL = 2)

## ---- echo=FALSE---------------------------------------------------------
h <- c(0.001, 50.65, 293.77, 790.14, 992.74, 5065, 10130, 15195)
w <- c(0.5650, 0.4013, 0.2502, 0.2324, 0.2307, 0.1926, 0.1812, 0.1730)
fitsoilwater(theta=w, x=h, ylim=c(0.1,0.6))

## ------------------------------------------------------------------------
Sindex(theta_R=0, theta_S=0.395, alpha=0.0217, n=1.103, xlim = c(0, 1000))

## ------------------------------------------------------------------------
data(SoilAggregate)
head(SoilAggregate)
classes <- c(3, 1.5, 0.75, 0.375, 0.178, 0.053)
out <- aggreg.stability(sample.id = SoilAggregate[ ,1], 
                        dm.classes = classes, 
                        aggre.mass = SoilAggregate[ ,-1])
head(out)

