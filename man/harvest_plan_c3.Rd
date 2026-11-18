\name{harvest_plan_c3}
\alias{harvest_plan_c3}
\alias{harvest_complete_c3}
\title{Pilot LANDIS-II v7 Random stand ranking with an explicit MT19937 stream}
\usage{
harvest_plan_c3(simulation, topology, targets, economic = NULL,
  forest_type = NULL, allocation_mode = c("corrected", "landis_compatible"),
  landis_seed, landis_rng_offset = 0)

harvest_complete_c3(simulation, topology, targets, economic = NULL,
  forest_type = NULL, allocation_mode = c("corrected", "landis_compatible"),
  landis_seed, landis_rng_offset = 0)
}
\arguments{
  \item{simulation}{An existing landisR simulation with a shared live CSR state.}
  \item{topology}{Topology created by \code{harvest_prepare_management}.}
  \item{targets}{A data frame with one row per management area; columns \code{management_area}, \code{fraction}, \code{ranking}, \code{target_type} and optional eligibility criteria. Ranking is one of \code{Random}, \code{MaxCohortAge}, \code{Economic}, \code{RegulateAges}.}
  \item{economic}{Optional table for the Economic ranking.}
  \item{forest_type}{Optional table for ForestTypeTable eligibility.}
  \item{allocation_mode}{Allocation and ForestTypeTable interpretation.}
  \item{landis_seed}{Seed for the LANDIS-II v7 32-bit MT19937 generator.}
  \item{landis_rng_offset}{Number of draws already consumed by the \emph{same} external LANDIS-II generator. It must be accurately supplied for coupled simulations.}
}
\description{
  Reproduces the Random stand-ranking and selection sequence in H04-C3a benchmark scenarios using LANDIS-II v7 MT19937 and the core's shuffle-with-removal algorithm.
}
\details{
  This is an isolated pilot for SiteSelection Complete and CohortsRemoved ClearCut. The external RNG is NOT automatically connected to Biomass Succession, Base Fire or other perturbations, and the returned \code{landis_rng$next_offset} does not mutate the landisR simulation RNG. Exact parity of full coupled simulations requires a single shared RNG state and the exact order of draws across components.
  One random draw per harvested stand is modeled for the H04-C3a one-cell Complete fixtures, as observed in the LANDIS-II binaries. More complicated spatial selection methods remain to be qualified.
}
\value{\code{harvest_plan_c3} returns a plan with stand ranks, selected cells, management-area statistics and \code{landis_rng} offsets. \code{harvest_complete_c3} also applies the harvest to a new live CSR state and returns the plan and effects.}

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
  minimum_age = 0L, ranking = "MaxCohortAge", target_type = "area")
target$ranking <- "Random"
plan <- harvest_plan_c3(sim, topology, target, landis_seed = 42)
plan$cells
cut <- harvest_complete_c3(sim, topology, target, landis_seed = 42)
cut$sites
}
