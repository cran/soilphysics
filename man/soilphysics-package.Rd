\name{soilphysics-package}
\alias{soilphysics-package}
\alias{soilphysics}
\docType{package}

\title{ Soil Physical Analysis }

\description{
Basic and model-based soil physical analyses. 
}

\details{
\tabular{ll}{
Package: \tab soilphysics\cr
Type: \tab Package\cr
Version: \tab 2.3\cr
Date: \tab 2015-05-09\cr
License: \tab GPL (>= 2)\cr
}
Functions for modelling the load bearing capacity and the penetration resistance. The package
allows one to model the soil water retention through six different models. There are some useful
and easy-to-use functions to perform parameter estimation of these models. Methods to obtain the 
preconsolidation stress are available, such as the standard of Casagrande (1936) and so on. It is 
possible to quantify soil water availability for plants through the Least Limiting Water Range 
approach as well as the Integral Water Capacity. Also, it is possible to determine the water suction 
at the point of hydraulic cut-off. Moreover, users can deal with the high-energy-moisture-characteristics 
(HEMC) methodology proposed by Pierson & Mulla (1989), which is used to analyze the aggregate stability. 
There is a function to determine the soil critical moisture and the maximum bulk density for one or more 
samples, based on the Proctor (1933) compaction test. Other utilities like a function to calculate the 
soil liquid limit, the void ratio and to determine the maximum curvature point are available.
}

\note{
\emph{soilphysics} is an ongoing project.
We welcome any and all criticism, comments and suggestions.
}

\author{
Anderson Rodrigo da Silva, Renato Paiva de Lima

Maintainer: Anderson Rodrigo da Silva <anderson.agro@hotmail.com>
}

\references{
ABNT - Associacao Brasileira de Normas Tecnicas. (1990). \emph{Ensaio de
    adensamento unidimensional}: NBR 12007. Rio de Janeiro. 13p.

Arvidsson, J.; Keller, T. (2004). Soil precompression stress
    I. A survey of Swedish arable soils. \emph{Soil & Tillage Research}, 77:85-95.

Bowles, J. A. (1986). \emph{Engineering Properties of Soils and their
    Measurements}. 3rd ed. McGraw-Hill Book Company, Inc. NY, 218pp.

Busscher, W.J. (1990) Adjustment of flat-tipped penetrometer resistance
   data to common water content. \emph{Transactions of the ASAE}, 3:519-524.

Casagrande, A. (1936). \emph{The determination of the pre-consolidation
   load and its practical significance}. In: Proceedings of the
   International Conference on Soil Mech. and Found. Eng.
   (ICSMFE), Cambridge, MA, 22-26 June 1936, vol. 3. Harvard
   University, Cambridge, MA, USA, pp. 60-64.

Dexter, A. R. (2004). Soil physical quality Part I. Theory,
   effects of soil texture, density, and organic matter,
   and effects on root growth. \emph{Geoderma}, 120:201-214.

Dexter et al. (2008). A user-friendly water retention function that takes account of the textural
	and structural pore spaces in soil. \emph{Geoderma}, 143:243-253.

Dexter et al. (2012). Equilibrium, non-equilibrium and residual water: Consequences for soil
	water retention. Geoderma, 177-178: 63-71.

Dias Junior, M. S. (1994). \emph{Compression of three soils under longterm
   tillage and wheel traffic}. 1994. 114p. Ph.D. Thesis - Michigan
   State University, East Lansing.

Dias Junior, M. S.; Pierce, F. J. (1995). A simple procedure for
   estimating preconsolidation pressure from soil compression curves.
   \emph{Soil Technology}, 8:139-151.

Genuchten, M.T. van. (1980) A closed form equation for predicting the
   hydraulic conductivity of unsaturated soils. \emph{Soil Science Society of
   America Journal}, 44:892-898.

Groenevelt & Grant (2004). A newmodel for the soil-water retention curve that solves the problem 
	of residualwater contents. \emph{European Journal of Soil Science}, 55:479-485.

Leao et al. (2005). An Algorithm for Calculating the Least Limiting Water Range of Soils. 
	\emph{Agronomy Journal}, 97:1210-1215.

Leao et al. (2006). Least limiting water range: A potential indicator of changes in near-surface 
	soil physical quality after the conversion of Brazilian Savanna into pasture. 
	\emph{Soil & Tillage Research}, 88:279-285.

Mualem, Y. (1976) A new model for predicting the hydraulic conductivity
   of unsaturated porous media. \emph{Water Resources Research}, 12:513-522.

Pierson, F.B.; Mulla, D.J. (1989) An Improved Method for Measuring Aggregate Stability
of a Weakly Aggregated Loessial Soil. \emph{Soil Sci. Soc. Am. J.}, 53:1825--1831.

Proctor, E. R. (1933). Design and construction of rolled earth dams.
   \emph{Eng. News Record}, 3:245-284, 286-289, 348-351, 372-376.

Ross et al. (1991). Equation for extending water-retention curves to dryness.
	\emph{Soil Science Society of America Journal}, 55:923-927. 

Silva et al. (2010). Indicadores da qualidade fisica do solo. In: Jong Van Lier, Q. (Ed).
   \emph{Fisica do solo}. Vicosa (MG): Sociedade Brasileira de Ciencia do Solo. p.541-281.

Silva et al. (1994). Characterization of the least limiting water range of soils. 
	\emph{Soil Science Society of America Journal}, 58:1775-1781.

Sowers, G. F. (1965). Consistency. In: BLACK, C.A. (Ed.). \emph{Methods of soil analysis}.
   Madison: American Society of Agronomy. Part 1, p.391-399. (Agronomy, 9).

Sowers, G. F. (1965). Consistency. In: KLUTE, A. (Ed.). 2nd ed. \emph{Methods of soil analysis}.
   Madison: American Society of Agronomy. Part 1, p.545-566.
}
