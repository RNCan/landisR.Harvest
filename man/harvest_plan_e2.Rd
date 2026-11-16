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
