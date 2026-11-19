\name{landisR.harvest-package}
\alias{landisR.harvest-package}
\alias{landisR.harvest}
\title{landisR.harvest: Biomass Harvest for LANDIS-R}
\description{
Plan and apply harvest prescriptions to stands on the shared living-cohort
state, and track harvest events and history.
}
\details{
Use \code{\link{harvest_read}} to import supported LANDIS-II inputs.
For coupled simulations, register \code{\link{harvest_extension}}.
For a manual loop, prepare the management topology with
\code{\link{harvest_prepare}} and carry the updated simulation and history
between calls to \code{\link{harvest_step}}. Use
\code{\link{harvest_plan_prescriptions}} to inspect an allocation before
applying it.
}
\seealso{\code{\link{harvest_apply_cells}}, \code{\link{harvest_legacy}}}
\docType{package}
