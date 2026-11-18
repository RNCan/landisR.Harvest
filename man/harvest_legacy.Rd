\name{harvest_legacy}
\alias{harvest_legacy}
\alias{harvest_plan_d2}
\alias{harvest_prepare_adjacency_d2}
\alias{harvest_step_d2}
\alias{harvest_plan_e3}
\alias{harvest_step_e3}
\alias{harvest_plan_e4d}
\alias{harvest_step_e4d}
\title{Compatibility interfaces for earlier Harvest configurations}
\description{
Prepare stand adjacency, plan harvest on current living cohorts, or apply
the plan. These interfaces preserve earlier calling conventions for
deterministic single-target configurations and multi-prescription planning.
}
\usage{
harvest_prepare_adjacency_d2(simulation, topology)

harvest_plan_d2(simulation, topology, targets, history = NULL,
  year = simulation$year, economic = NULL, forest_type = NULL,
  allocation_mode = c("corrected", "landis_compatible"))
harvest_step_d2(simulation, topology, targets, history = NULL,
  year = simulation$year, rules = NULL, clearcut = FALSE, economic = NULL,
  forest_type = NULL, allocation_mode = c("corrected", "landis_compatible"))

harvest_plan_e3(simulation, topology, prescriptions, history = NULL,
  year = simulation$year, landis_seed, landis_rng_offset = 0,
  allocation_mode = c("corrected", "landis_compatible"))
harvest_step_e3(simulation, topology, prescriptions, rules_by_prescription,
  history = NULL, year = simulation$year, landis_seed, landis_rng_offset = 0,
  allocation_mode = c("corrected", "landis_compatible"))

harvest_plan_e4d(simulation, topology, prescriptions, history = NULL,
  year = simulation$year, landis_seed, landis_rng_offset = 0,
  allocation_mode = c("corrected", "landis_compatible"),
  rules_by_prescription = NULL,
  harvest_timestep = simulation$project$canonical$timestep,
  damage_history = NULL, report = getOption("landisR.progress.report"))
harvest_step_e4d(simulation, topology, prescriptions, rules_by_prescription,
  history = NULL, year = simulation$year, landis_seed, landis_rng_offset = 0,
  allocation_mode = c("corrected", "landis_compatible"),
  harvest_timestep = simulation$project$canonical$timestep,
  damage_history = NULL, progress = getOption("landisR.progress", interactive()))
}
\arguments{
  \item{simulation}{An initialized landisR simulation.}
  \item{topology}{Stand and management topology for the same grid. Start
  with \code{harvest_prepare_management()} and add rook adjacency with
  \code{harvest_prepare_adjacency_d2()}. The E4d functions can also use the
  result of \code{harvest_prepare()}.}
  \item{targets}{One row per management area, with \code{management_area},
  \code{ranking}, \code{fraction}, and \code{target_type}. D2 uses Complete
  selection and deterministic ranking. Optional constraints include
  \code{minimum_age}, \code{maximum_age}, \code{stand_adjacency},
  \code{adjacency_type} and \code{minimum_time_since_last_harvest}.}
  \item{history}{Previous stand history with \code{stand} and
  \code{last_harvest_year}. When continuing E4d, keep all history attributes.}
  \item{year}{Nonnegative integer simulation year for scheduling and history.}
  \item{economic, forest_type}{Optional species ranking and forest-type
  rule tables, as in \code{harvest_plan_c2()}.}
  \item{allocation_mode}{Use the default corrected mode, or specify
  \code{"landis_compatible"} for the compatible allocation rules available
  to the selected planner.}
  \item{rules}{Cohort-removal rule data frame for D2, as in
  \code{harvest_apply_cells()}.}
  \item{clearcut}{For D2, remove all cohorts from the selected cells.
  Supply \code{rules} or set \code{clearcut = TRUE}.}
  \item{prescriptions}{Prescription table. E3 accepts \code{Complete}
  selection, \code{MaxCohortAge} ranking and area quotas. E4d accepts the
  prescriptions documented in \code{\link{harvest_step}}.}
  \item{rules_by_prescription}{Named list keyed by prescription. E3
  takes cohort-removal tables. For E4d, see the treatment formats in
  \code{\link{harvest_step}}.}
  \item{landis_seed, landis_rng_offset}{Unsigned 32-bit MT19937 seed and
  number of prior draws. E3 uses an explicit standalone stream. With E4d
  and a shared LANDIS-II stream, use the simulation seed and
  \code{allocation_mode = "landis_compatible"}. The offset must be zero
  or match the current offset in the simulation.}
  \item{harvest_timestep}{Harvest interval in years for E4d.}
  \item{damage_history}{Optional external cell damage records for E4d,
  with \code{cell_id} and \code{year} columns.}
  \item{report}{Optional callback taking a progress fraction and
  phase label.}
  \item{progress}{Show a console progress bar for the E4d step.}
}
\details{
D2 can apply adjacency and elapsed-time restrictions with one target per
management area. Its Complete-only planner does not activate neighboring
stand reservations. For new multi-prescription runs, use
\code{\link{harvest_step}} or \code{\link{harvest_extension}}.
}
\value{
\code{harvest_prepare_adjacency_d2()} returns the topology with its
adjacency cache. \code{harvest_plan_d2()} returns stand and area
summaries, selected cells, and planned history. E3 and E4d planners
return prescription allocations, sparse event records, selected cells,
history and random-stream offsets. Step functions also return the
simulation after removal, the plan and sparse treatment effects.
}
\seealso{\code{\link{harvest_prepare}}, \code{\link{harvest_step}}}
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
topology <- harvest_prepare_adjacency_d2(sim, topology)
target <- data.frame(management_area = 1L, fraction = 0.5,
  minimum_age = 0L, ranking = "MaxCohortAge", target_type = "area")
d2 <- harvest_plan_d2(sim, topology, target, year = 10L)
cut_d2 <- harvest_step_d2(sim, topology, target, year = 10L, clearcut = TRUE)

prescriptions <- transform(target, prescription = "CUT", maximum_age = 120L,
  begin_time = 10L, end_time = 20L, site_selection = "Complete")
rules <- list(CUT = data.frame(species = "PINE",
  age_min = 1L, age_max = 120L, fraction = 1))
e3 <- harvest_plan_e3(sim, topology, prescriptions, year = 10L, landis_seed = 42)
cut_e3 <- harvest_step_e3(sim, topology, prescriptions, rules,
  year = 10L, landis_seed = 42)
e4d <- harvest_plan_e4d(sim, topology, prescriptions, year = 10L,
  landis_seed = 42, rules_by_prescription = rules)
cut_e4d <- harvest_step_e4d(sim, topology, prescriptions, rules,
  year = 10L, landis_seed = 42, progress = FALSE)
}
