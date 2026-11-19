\name{harvest_read}
\alias{harvest_read}
\title{Import LANDIS-II Biomass Harvest inputs}
\usage{
harvest_read(path, project, scenario_directory = dirname(path),
  raster_template = NULL)
}
\arguments{
  \item{path}{Path to the LANDIS-II Biomass Harvest parameter file.}
  \item{project}{A landisR project containing the landscape grid and duration.}
  \item{scenario_directory}{LANDIS-II scenario directory used to resolve
  relative \code{ManagementAreas} and \code{Stands} map paths.}
  \item{raster_template}{Optional single-layer \code{terra} raster or path
  identifying the reference grid, usually the ecoregion map. Defaults to the
  project ecoregion map when its path is known, otherwise the management-area
  map.}
}
\description{
Read a Biomass Harvest parameter file and the two spatial inputs it
references. Join the maps to the active cell IDs of a landisR project and
prepare the tables accepted by \code{\link{harvest_extension}}.
}
\details{
The reader supports the LANDIS-II syntax used by the 6est scenario:
\code{Economic} stand ranking, \code{MinimumAge}, optional
\code{MinimumTimeSinceLastHarvest}, \code{CompleteStandSpread},
\code{SpeciesList} removal rules and percentage-of-area assignments in
\code{HarvestImplementations}. An age range without a percentage specifies
100 percent removal. The harvest engine also accepts other prescription
types through the table interface. Output-file settings in the input are for
LANDIS-II. Set the output paths and formats in the R simulation script.

The two input rasters must align with the reference grid. Zero or missing
management-area values on active cells are assigned an unmanaged sector
code. Each active cell requires a positive stand code.
}
\value{
A list with \code{management} (\code{cell_id}, \code{management_area},
\code{stand}), \code{prescriptions} (including a list column of economic
ranking tables), \code{rules_by_prescription} (named cohort-removal
tables), and \code{timestep} (years).
}
\examples{
\dontrun{
harvest <- harvest_read("inputs/biomass-harvest.txt", project,
  scenario_directory = "scenario", raster_template = "inputs/ecoregions.tif")
extension_harvest <- harvest_extension(
  management = harvest$management,
  prescriptions = harvest$prescriptions,
  rules_by_prescription = harvest$rules_by_prescription,
  timestep = harvest$timestep,
  allocation_mode = "landis_compatible"
)
}
}
