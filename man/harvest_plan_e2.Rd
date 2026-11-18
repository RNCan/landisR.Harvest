\name{harvest_plan_e2}
\alias{harvest_plan_e2}
\alias{harvest_step_e2}
\title{Multiple-prescription experimental LANDIS-II Harvest selection}
\description{Allocates site-selection Complete harvest stands across competing prescriptions in each management area, using remaining-target probability weights and an explicit MT19937-compatible draw stream.}\usage{
harvest_plan_e2(simulation, topology, prescriptions, year=simulation$year,
  landis_seed, landis_rng_offset=0,
  allocation_mode=c("corrected","landis_compatible"))
harvest_step_e2(simulation, topology, prescriptions, rules_by_prescription,
  year=simulation$year, landis_seed, landis_rng_offset=0,
  allocation_mode=c("corrected","landis_compatible"))
}
\arguments{
\item{simulation}{A landisR simulation sharing the current CSR forest state.}
\item{topology}{The precomputed stand and management-area topology.}
\item{prescriptions}{A data frame with management_area, prescription, fraction (0..1), ranking (MaxCohortAge), target_type (area), minimum_age, maximum_age, begin_time, end_time, site_selection (Complete). The row order is LANDIS-II implementation order within management areas.}
\item{year}{Current LANDIS-II harvest time.}
\item{landis_seed}{UInt32 LANDIS-II MT19937 seed.}
\item{landis_rng_offset}{Number of draws already consumed by the LANDIS-II stream.}
\item{allocation_mode}{Corrected semantics, or legacy LANDIS-compatible input validation.}
\item{rules_by_prescription}{Named list of cohort-removal rule data frames (species, age_min, age_max, fraction), with one entry for each active prescription.}
}
\details{H04-E2 is an experimental Complete and MaxCohortAge area-target pilot. It is not a full simulation-level RNG coupling. The C++ kernel executes one prescription-choice draw per harvested stand. The site effects are applied per prescription to disjoint selected cells using the H03 cohort engine. No dead pools are created. Spatial adjacency constraints, Random rankings, spreading, and repeat prescriptions are deliberately outside this pilot.}
\value{The planner returns exact-ordered events, selected cells, per-prescription allocations, and RNG draw offsets. The step function also returns the updated shared living-cohort simulation and sparse effects.}

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
prescriptions <- transform(target, prescription = "CUT", maximum_age = 120L,
  begin_time = 10L, end_time = 20L, site_selection = "Complete")
rules <- list(CUT = data.frame(species = "PINE",
  age_min = 1L, age_max = 120L, fraction = 1))
plan <- harvest_plan_e2(sim, topology, prescriptions, year = 10L,
  landis_seed = 42)
cut <- harvest_step_e2(sim, topology, prescriptions, rules, year = 10L,
  landis_seed = 42)
cut$sites
}
