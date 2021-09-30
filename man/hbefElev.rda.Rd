\name{hbefElev}

\docType{data}

\alias{hbefElev}

\title{Elevation in meters extracted at a 30m resolution of the Hubbard Brook Experimental Forest}

\description{
Elevation in meters extracted at a 30m resolution of the Hubbard Brook Experimental Forest. 
Data come from the National Elevation Dataset. 
}

\usage{
data(hbefElev)
}

\format{
  \code{hbefElev} is a data frame with three columns:

   \code{val}: the elevation value. The value is standardized using the mean 
   and standard deviation of the elevation values in the \code{hbef2015} that 
   are used in the model fitting process.

   \code{Easting}: the x coordinate of the point. The proj4string is 
   "+proj=utm +zone=19 +units=m +datum=NAD83".

   \code{Northing}: the y coordinate of the point. The proj4string is 
   "+proj=utm +zone=19 +units=m +datum=NAD83".
}

\keyword{datasets}
