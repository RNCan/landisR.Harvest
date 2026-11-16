\name{harvest_apply_cells}
\alias{harvest_apply_cells}
\title{Experimental deterministic cohort harvest on specified cells}
\usage{harvest_apply_cells(simulation, cells, rules = NULL, clearcut = FALSE,
  preserve_cohort_order = FALSE)}
\arguments{
\item{preserve_cohort_order}{Preserve the order of surviving cohorts.}
\item{simulation}{A landisR 0.1.0.9007 interactive simulation.}
\item{cells}{Unique active canonical cell IDs.}
\item{rules}{Data frame with species, age_min, age_max, fraction. Inclusive age limits; disjoint ranges per species. Fraction 1 means whole removal.}
\item{clearcut}{Logical; remove all cohorts. Mutually exclusive with rules.}
}
\description{H03-A prototype. Reads native compressed living biomass states, calculates deterministic removals and applies a single batch through landis_apply. No spatial allocation, regeneration, dead pools, or clock advancement.}
\value{List with simulation, replacement, effects, sites, changed_cells, clear_cells, states_evaluated. Biomass units are g/m2.}
