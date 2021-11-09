\name{hbef2015}

\docType{data}

\alias{hbef2015}

\title{Detection-nondetection data of 12 foliage gleaning bird species
  in 2015 in the Hubbard Brook Experimental Forest}

\description{
Detection-nondetection data of 12 foliage gleaning bird species in 2015 in the 
Hubbard Brook Experimental Forest (HBEF) in New Hampshire, USA. Data were collected at 373
sites over three replicate point counts each of 10 minutes in length, with a detection 
radius of 100m. Some sites were not visited for all three replicates. The 12 species
included in the data set are as follows: (1) AMRE: American Redstart; (2) BAWW: 
Black-and-white Warbler; (3) BHVI: Blue-headed Vireo; (4) BLBW: Blackburnian Warbler; 
(5) BLPW: Blackpoll Warbler; (6) BTBW: Black-throated Blue Warbler; (7) BTNW: 
BLack-throated Green Warbler; (8) CAWA: Canada Warbler; (9) MAWA: Magnolia Warbler; 
(10) NAWA: Nashville Warbler; (11) OVEN: Ovenbird; (12) REVI: Red-eyed Vireo. 
}

\usage{
data(hbef2015)
}

\source{Rodenhouse, N. and S. Sillett. 2019. Valleywide Bird Survey, Hubbard Brook 
	Experimental Forest, 1999-2016 (ongoing) ver 3. Environmental Data Initiative. 
	\doi{10.6073/pasta/faca2b2cf2db9d415c39b695cc7fc217} 
	(Accessed 2021-09-07)}

\references{
  Doser, J.W., Leuenberger, W., Sillett, T.S., Hallworth, M.T., Zipkin, E.F. 
  Integrated community occupancy models:  A framework to assess occurrence 
  and biodiversity dynamics using multiple data sources.
  https://arxiv.org/abs/2109.01894
}

\format{
  \code{hbef2015} is a list with four elements: 

   \code{y}: a three-dimensional array of detection-nondetection data with 
     dimensions of species (12), sites (373) and replicates (3). 

   \code{occ.covs}: a numeric matrix with 373 rows and one column consisting of the 
     elevation at each site.

   \code{det.covs}: a list of two numeric matrices with 373 rows and 3 columns. 
     The first element is the day of year when the survey was 
     conducted for a given site and replicate. The second element is the
     time of day when the survey was conducted. 

   \code{coords}: a numeric matrix with 373 rows and two columns containing the 
     site coordinates (Easting and Northing) in UTM Zone 19. The proj4string is
     "+proj=utm +zone=19 +units=m +datum=NAD83".
}

\keyword{datasets}

