\name{harvest_apply_cells}
\alias{harvest_apply_cells}
\title{Experimental deterministic cohort harvest on specified cells}
\usage{harvest_apply_cells(simulation, cells, rules = NULL, clearcut = FALSE,
  preserve_cohort_order = FALSE)}
\arguments{
\item{preserve_cohort_order}{Preserve the order of surviving cohorts.}
\item{simulation}{An initialized landisR simulation.}
\item{cells}{Unique active canonical cell IDs.}
\item{rules}{Data frame with species, age_min, age_max, fraction. Inclusive age limits; disjoint ranges per species. Fraction 1 means whole removal.}
\item{clearcut}{Logical; remove all cohorts. Mutually exclusive with rules.}
}
\description{H03-A prototype. Reads native compressed living biomass states, calculates deterministic removals and applies a single batch through landis_apply. No spatial allocation, regeneration, dead pools, or clock advancement.}
\value{List with simulation, replacement, effects, sites, changed_cells, clear_cells, states_evaluated. Biomass units are g/m2.}

\examples{
grid <- data.frame(cell_id = 1:2, row = 1L, col = 1:2, active = TRUE)
project <- landisR::landis_build_landscape(grid, duration = 20L,
  cell_length = 100, species = data.frame(species = "PINE"))
sim <- landisR::landis_initialize(project)
cohorts <- data.frame(cell_id = 1:2, species = "PINE",
  age = c(40L, 20L), biomass = 100)
sim <- landisR::landis_apply(sim, landisR::landis_update(cohorts = cohorts))
rule <- data.frame(species = "PINE",
  age_min = 30L, age_max = 120L, fraction = 0.25)
cut <- harvest_apply_cells(sim, cells = 1L, rules = rule)
cut$effects
landisR::landis_state(cut$simulation, representation = "cohorts")
}
