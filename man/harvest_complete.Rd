\name{harvest_complete}
\alias{harvest_complete}
\alias{harvest_prepare_management}
\alias{harvest_plan_complete}
\title{Experimental whole-stand harvest allocation}
\usage{
harvest_prepare_management(simulation, management)
harvest_plan_complete(simulation, topology, targets)
harvest_complete(simulation, topology, targets)
}
\arguments{
\item{simulation}{An initialized landisR simulation.}
\item{management}{Data frame with cell_id, management_area, stand; one row per active cell.}
\item{topology}{Prepared management topology for the same active grid.}
\item{targets}{Data frame with management_area, fraction (0 to 1), minimum_age; one row per area.}
}
\description{Pilot MaxCohortAge, inclusive MinimumAge, Complete and ClearCut. Computes mean of cell maximum ages, truncated, using native CSR access. Allocates full stands up to or beyond the target; stops if candidates exhausted. Eligible ties and selected stands with empty cells are rejected pending qualification. No time advance, RNG, harvest history, or regeneration specific to harvest.}
\value{Prepared topology, allocation plan, or harvest_apply_cells result with plan.}

\examples{
grid <- data.frame(cell_id = 1:2, row = 1L, col = 1:2, active = TRUE)
project <- landisR::landis_build_landscape(grid, duration = 20L,
  cell_length = 100, species = data.frame(species = "PINE"))
sim <- landisR::landis_initialize(project)
cohorts <- data.frame(cell_id = 1:2, species = "PINE",
  age = c(40L, 20L), biomass = 100)
sim <- landisR::landis_apply(sim, landisR::landis_update(cohorts = cohorts))
management <- data.frame(cell_id = 1:2, management_area = 1L, stand = 1:2)
topology <- harvest_prepare_management(sim, management)
target <- data.frame(management_area = 1L, fraction = 0.5,
  minimum_age = 30L)
plan <- harvest_plan_complete(sim, topology, target)
plan$cells
cut <- harvest_complete(sim, topology, target)
cut$sites
}
