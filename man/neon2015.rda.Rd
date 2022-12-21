\name{neon2015}

\docType{data}

\alias{neon2015}

\title{Detection-nondetection data of 12 foliage gleaning bird species
  in 2015 in Bartlett Experimental Forest in New Hampshire, USA}

\description{
Detection-nondetection data of 12 foliage gleaning bird species in 2015 in the 
Bartlett Experimental Forest in New Hampshire, USA. These data were collected as 
part of the National Ecological Observatory Network (NEON). Data were collected
at 80 sites where observers recorded the number of all bird species observed during a
six minute, 125m radius point count survey once during the breeding season. The six
minute survey was split into three two-minute intervals following a removal design
where the observer recorded the interval during which a species was first observed
(if any) with a 1, intervals prior to observation with a 0, and then mentally removed
the species from subsequent intervals (marked with NA), which enables modeling of data 
in an occupancy modeling framework. The 12 species included in the data set are as follows: 
(1) AMRE: American Redstart; (2) BAWW: Black-and-white Warbler; 
(3) BHVI: Blue-headed Vireo; (4) BLBW: Blackburnian Warbler; 
(5) BLPW: Blackpoll Warbler; (6) BTBW: Black-throated Blue Warbler; (7) BTNW: 
BLack-throated Green Warbler; (8) CAWA: Canada Warbler; (9) MAWA: Magnolia Warbler; 
(10) NAWA: Nashville Warbler; (11) OVEN: Ovenbird; (12) REVI: Red-eyed Vireo. 
}

\usage{
data(neon2015)
}

\source{NEON (National Ecological Observatory Network). Breeding landbird point counts, 
	RELEASE-2021 (DP1.10003.001). https://doi.org/10.48443/s730-dy13. 
	Dataset accessed from https://data.neonscience.org on October 10, 2021}
\references{
  Doser, J. W., Leuenberger, W., Sillett, T. S., Hallworth, M. T. & Zipkin, E. F. (2022). 
  Integrated community occupancy models: A framework to assess occurrence and biodiversity 
  dynamics using multiple data sources. Methods in Ecology and Evolution, 00, 1-14. \doi{10.1111/2041-210X.13811}

  Barnett, D. T., Duffy, P. A., Schimel, D. S., Krauss, R. E., Irvine, K. M., 
  Davis, F. W.,Gross, J. E., Azuaje, E. I., Thorpe, A. S., Gudex-Cross, D., et al. (2019). 
  The terrestrial organism and biogeochemistry spatial sampling design for 
  the national ecological observatory network. Ecosphere, 10(2):e02540.
}

\format{
  \code{neon2015} is a list with four elements: 

   \code{y}: a three-dimensional array of detection-nondetection data with 
     dimensions of species (12), sites (80) and replicates (3). 

   \code{occ.covs}: a numeric matrix with 80 rows and one column consisting of the 
     elevation at each site.

   \code{det.covs}: a list of two numeric vectors with 80 elements. The 
     first element is the day of year when the survey was conducted for a given 
     site. The second element is the time of day when the survey began.

   \code{coords}: a numeric matrix with 80 rows and two columns containing the 
     site coordinates (Easting and Northing) in UTM Zone 19. The proj4string is
     "+proj=utm +zone=19 +units=m +datum=NAD83".
}

\keyword{datasets}


