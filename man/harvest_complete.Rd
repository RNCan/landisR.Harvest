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
\item{simulation}{landisR 0.1.0.9007 simulation.}
\item{management}{Data frame with cell_id, management_area, stand; one row per active cell.}
\item{topology}{Prepared management topology for the same active grid.}
\item{targets}{Data frame with management_area, fraction (0 to 1), minimum_age; one row per area.}
}
\description{Pilot MaxCohortAge, inclusive MinimumAge, Complete and ClearCut. Computes mean of cell maximum ages, truncated, using native CSR access. Allocates full stands up to or beyond the target; stops if candidates exhausted. Eligible ties and selected stands with empty cells are rejected pending qualification. No time advance, RNG, harvest history, or regeneration specific to harvest.}
\value{Prepared topology, allocation plan, or harvest_apply_cells result with plan.}
