\name{harvest_plan}
\alias{harvest_plan}
\alias{harvest_complete_b2}
\title{H04-B2 native whole-stand allocation with Economic ranking and stands targets}
\usage{
harvest_plan(simulation, topology, targets, economic = NULL,
             allocation_mode = c("corrected", "landis_compatible"))
harvest_complete_b2(simulation, topology, targets, economic = NULL,
                    allocation_mode = c("corrected", "landis_compatible"))
}
\arguments{
\item{simulation}{Version-pinned An initialized landisR simulation.}
\item{topology}{Object returned by harvest_prepare_management.}
\item{targets}{One row per management area, with columns management_area, fraction (0..1), minimum_age (inclusive), ranking (MaxCohortAge or Economic), target_type (area or stands). Optional site_selection=Complete, cohorts_removed=ClearCut.}
\item{economic}{For Economic ranking, table with species, rank (0..100), merchantable_age (positive integer); optional management_area supports distinct economic tables by unit. Omitted species contribute zero rank.}
\item{allocation_mode}{corrected (default), or landis_compatible for observed LANDIS-II Biomass Harvest 4.4 behavior of percentage-of-stands objectives with Complete.}
}
\description{Prepare or apply a whole-stand harvest allocation using a deterministic ranking.}
\details{The forest state is not mutated by harvest_plan. harvest_complete_b2 applies a single group CSR transaction via the validated harvest_apply_cells kernel. Stand ages are PRE-HARVEST ages; LANDIS-II Economic event logs may instead report 0 after cutting. Cohort ages must already reflect the succession engine's age conventions; the Harvest module never rounds them. Negative ranks, eligible rank ties, selected stands with empty sites, multiple prescriptions, other selectors and zero merchantable age remain out of scope in H04-B2. LANDIS-compatible stands allocation was experimentally qualified for positive fractions .10/.30/.50/1.00 with Complete; the mode deliberately reproduces the legacy anomalous behavior in these cases.}
\value{harvest_plan returns a plan. harvest_complete_b2 returns the same effect list as harvest_apply_cells, with plan attached.}

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
plan <- harvest_plan(sim, topology, target)
plan$cells
cut <- harvest_complete_b2(sim, topology, target)
cut$sites
}
