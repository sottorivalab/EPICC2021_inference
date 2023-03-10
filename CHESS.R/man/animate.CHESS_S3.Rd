% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/chess_s3_object-plotting_methods.R
\name{animate.CHESS_S3}
\alias{animate.CHESS_S3}
\title{Function to plot a animation of a CHESS S3 simulation}
\usage{
\method{animate}{CHESS_S3}(
  x,
  mm_per_dot = 0.1,
  background_color = NULL,
  color_by = "type",
  alpha_burden = TRUE
)
}
\arguments{
\item{x}{CHESS S3 object. Must include snapshots (see 'snapshot_times' argument of new_simulation function (?new_simulation))}

\item{mm_per_dot}{(optional) value that scales the grid points to a pseudo mm value (default: 0.1).}

\item{background_color}{(optional) color by which to color the background of the simulated space (default: NULL, no color).}

\item{color_by}{(optional) column in snapshot data by which to color cells (default: type, the cell type).}
}
\value{
A 'gganim' object containg the animation.
}
\description{
Function to plot a animation of a CHESS S3 simulation
}
\examples{
anim = animate.CHESS_S3(new_simulation(20, snapshot_times = 1:100))
}
