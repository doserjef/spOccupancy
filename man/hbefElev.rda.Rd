\name{hbefElev}

\docType{data}

\alias{hbefElev}

\title{Elevation in meters extracted at a 30m resolution across the Hubbard Brook Experimental Forest}

\description{
Elevation in meters extracted at a 30m resolution of the Hubbard Brook Experimental Forest. 
Data come from the National Elevation Dataset. 
}

\source{Gesch, D., Oimoen, M., Greenlee, S., Nelson, C., Steuck, M., & Tyler, D. (2002). 
	The national elevation dataset. Photogrammetric engineering and remote sensing, 
	68(1), 5-32.
}

\references{Gesch, D., Oimoen, M., Greenlee, S., Nelson, C., Steuck, M., & Tyler, D. (2002). 
	The national elevation dataset. Photogrammetric engineering and remote sensing, 
	68(1), 5-32.
}

\usage{
data(hbefElev)
}

\format{
  \code{hbefElev} is a data frame with three columns:

   \code{val}: the elevation value in meters.

   \code{Easting}: the x coordinate of the point. The proj4string is 
   "+proj=utm +zone=19 +units=m +datum=NAD83".

   \code{Northing}: the y coordinate of the point. The proj4string is 
   "+proj=utm +zone=19 +units=m +datum=NAD83".
}

\keyword{datasets}

