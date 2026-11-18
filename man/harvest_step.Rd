\name{harvest_step}
\alias{harvest_prepare}
\alias{harvest_plan_prescriptions}
\alias{harvest_step}
\alias{harvest_extension}
\title{Prepare, plan and apply biomass harvest prescriptions}
\usage{
harvest_prepare(simulation, management, prescriptions = NULL)

harvest_plan_prescriptions(simulation, topology, prescriptions,
  rules_by_prescription, history = NULL, year = simulation$year,
  seed = simulation$seed, rng_offset = 0,
  allocation_mode = c("corrected", "landis_compatible"),
  harvest_timestep = simulation$project$canonical$timestep, damage_history = NULL)

harvest_step(simulation, topology, prescriptions, rules_by_prescription,
  history = NULL, year = simulation$year, seed = simulation$seed,
  rng_offset = 0, allocation_mode = c("corrected", "landis_compatible"),
  harvest_timestep = simulation$project$canonical$timestep, damage_history = NULL)

harvest_extension(management, prescriptions, rules_by_prescription,
  timestep = NULL, seed = NULL, rng_offset = 0,
  allocation_mode = c("corrected", "landis_compatible"), damage_history = NULL)
}
\arguments{
  \item{simulation}{Initialized landisR simulation.}
  \item{management}{Data frame with \code{cell_id}, \code{management_area}
  and \code{stand}, with one row per active cell.}
  \item{topology}{Result of \code{harvest_prepare()}.}
  \item{prescriptions}{Data frame with \code{management_area},
  \code{prescription}, \code{fraction}, \code{ranking}, \code{target_type},
  \code{minimum_age}, \code{maximum_age}, \code{begin_time}, \code{end_time}
  and \code{site_selection}. Areas without prescriptions are left unharvested.
  With \code{CompleteStandSpread} or \code{PartialStandSpread},
  include numeric \code{min_block_ha} and \code{max_block_ha}. The aliases
  \code{CompleteStandSpreading} and \code{PartialStandSpreading} are accepted.
  The minimum may be zero and the maximum must be positive and at least
  the minimum. With \code{PatchCutting}, include numeric \code{patch_fraction}
  (greater than zero and at most one, as a proportion of the stand) and
  \code{patch_size_ha} (positive). The optional \code{patch_option} is
  \code{AllowOverlap}, \code{PatchSize}, \code{PercentCut}, or empty.
  \code{PatchSize} and \code{PercentCut} are accepted as in the tested
  LANDIS-II v7 parser and have the same effect as the empty option.
  Supplying prescriptions to \code{harvest_prepare()} caches the cell
  topology for partial spreading and patch cutting.
  Optional columns specify adjacency, minimum
  time since harvest or a required species and age range. For adjacency,
  use \code{stand_adjacency} and \code{adjacency_type} (\code{StandAge} or
  \code{TimeSinceLastHarvested}). The optional
  \code{adjacency_neighbor_set_aside} sets the number of years to reserve
  neighboring stands after a harvest. Supply a nonnegative whole number
  for each prescription. The default is zero. For stand history, use
  \code{minimum_time_since_last_harvest}. For a required forest type, use
  \code{forest_type_species}, \code{forest_type_age_min},
  \code{forest_type_age_max} and \code{forest_type_percent}. Set
  \code{repeat_type} to \code{"SingleRepeat"} or \code{"MultipleRepeat"}
  and supply a positive integer \code{repeat_interval} in years. Leave
  \code{repeat_type} empty for a prescription with no scheduled return.
  For \code{MultipleRepeat}, \code{times_to_repeat} limits the number of
  additional entries after the initial harvest. Use an integer of at least
  two, or leave it unset or zero for repeats through the simulation.
  Set \code{min_time_since_damage} to a nonnegative integer in years
  to set the cell-level minimum. The default is zero.
  The optional numeric column \code{minimum_biomass_mg_ha} sets an inclusive
  minimum living biomass density for each prescription. Omit the column or
  supply \code{NA} or zero to disable the criterion for that prescription.}
  \item{rules_by_prescription}{Named list of treatments. Each treatment can
  be a removal data frame with columns \code{species}, \code{age_min},
  \code{age_max} and \code{fraction}, the string \code{"PlantOnly"}, or a list
  with optional \code{rules}, \code{plant} and \code{prevent_establishment}.
  Supply a removal table or \code{"PlantOnly"} as \code{rules}, a character
  vector of species as \code{plant}, and a logical flag as
  \code{prevent_establishment}. Omitted rules give zero removal.
  A SingleRepeat for \code{CUT} uses treatment \code{CUT} initially and
  \code{CUT(1)} on return. MultipleRepeat reuses the initial treatment.}
  \item{history}{Previous stand history returned by \code{harvest_step()}.
  Preserve the object with its attributes to carry pending repeat entries
  into the next call.}
  \item{year}{Year used to evaluate prescription dates and harvest history.}
  \item{seed}{Harvest RNG seed. The extension uses the simulation seed when
  \code{NULL}.}
  \item{rng_offset}{Number of draws already consumed in the standalone Harvest
  stream. With \code{recruitment_rng = "landis_v7"} the shared stream is read
  from the simulation. Use the simulation seed and \code{landis_compatible}.}
  \item{allocation_mode}{Compatibility mode for the planner. The compatible
  mode rejects explicitly specified zero quotas. For PartialStandSpread,
  and PatchCutting, \code{landis_compatible} uses the LANDIS-II v7 start-cell
  index \code{floor(u * (n - 1))} (zero-based). The default \code{corrected}
  mode allows every candidate cell to be drawn. With PatchCutting, it also
  completes eligible queued cells when the candidate list is exhausted and
  stops drawing starts once the stand target is met. The compatible mode reproduces the observed
  absence of a reservation effect in LANDIS-II v7. With \code{Economic} ranking and
  \code{CompleteStandSpread}, it also reproduces the LANDIS-II v7 ordering for
  equal ranks and adds equally ranked adjacent stands to the block in discovery order.}
  \item{harvest_timestep}{Harvest interval for a manual planning or harvest
  loop. By default it is the succession interval. Use the same interval
  in all calls.}
  \item{damage_history}{Optional data frame with active \code{cell_id} and
  \code{year} columns for recorded external damage, such as fire or wind.
  When a cell has several records, the most recent applies. Supply the latest
  report at each manual step. The extension also accepts a function of
  \code{simulation} and \code{year} that returns such a table at each Harvest
  step. Years greater than the proposed harvest year are invalid.}
  \item{timestep}{Harvest interval in years. \code{NULL} uses the succession
  interval. A supplied interval must be a positive multiple of it.}
}
\description{
The planner supports area quotas with Complete, CompleteStandSpread and
PartialStandSpread or PatchCutting, with Random, MaxCohortAge or Economic ranking. CompleteStandSpread
adds whole eligible stands in the same management area, in rank order
from the neighbors of any stand in the block. It stops adding stands
when the maximum area is reached or no eligible neighbor remains. An
added stand can take the block above the maximum. PartialStandSpread
grows cell by cell in eight directions from a random start in the initial
stand. Subsequent stands enter the block from their first cell in raster
order that touches the block. It stops when the maximum is reached or
the available neighbors are exhausted. A final cell can take the area
above a fractional maximum. If the accumulated area falls below the
minimum, the block is discarded. PatchCutting starts successive patches
within one stand until it reaches the stand fraction or exhausts candidates.
Each patch grows in eight directions. Adjacent cells are normally excluded
from later patches. \code{AllowOverlap} permits adjacent patches to meet.
A cell can be selected only once in the harvest entry. A site contributes to the quota when it
has positive biomass removal. Attempts without removal appear in the
event journal. The extension runs before Succession on the same forest
state and carries the history and Harvest RNG offset between steps.

In \code{corrected} mode, each stand with positive biomass removal in
an accepted harvest reserves the available stands that share a cell edge
with it within the same management area. The reserve takes effect as
soon as the harvest is accepted, for all subsequent prescriptions and
blocks, including within the same time step. The stand becomes eligible
again at the harvest year plus the specified duration. If multiple
harvests reserve a stand, the latest expiry applies. Reserved stands
cannot start or join a harvest block while the reserve remains in effect.

The biomass criterion applies to candidate stands for new entries,
including neighboring stands in a spread block. The planner averages
living biomass over all active managed cells in the stand, including
cells with zero biomass. It sums all living cohorts, regardless of age,
species or treatment rules. The density is in Mg/ha (100 g/m2 equals
1 Mg/ha). It is evaluated on the shared forest state at the start of
the Harvest step, after any extensions that already ran in that step.
Age and other eligibility criteria still apply. Scheduled returns keep
their existing eligibility rules. The biomass criterion is an optional
LANDIS-R extension available in either allocation mode.

A successful initial entry can schedule a return. The interval is rounded
up to a multiple of the Harvest timestep. Reserved stands are unavailable
to new entries until the return sequence ends. Due entries are processed
before new entries in their management area and are independent of the
initial prescription's dates and quota. A return can have zero biomass
removal and still counts toward the configured number of entries.
In \code{corrected} mode, positive removal on the return updates the
last-harvest year. Single returns use cells with positive removal in each
member stand. Periodic spread returns use the original block footprint.
A periodic complete or patch-cut return uses the same selector on its
reserved stand.

The damage filter admits cells when the elapsed time is at least the
minimum. The date is the latest of the harvest and supplied external
damage records. In \code{corrected} mode, harvest updates the date
after positive biomass removal. All selectors apply the filter and
spread does not cross stands with no eligible cells. In compatible
mode, all visited cells update the date, including visits with zero
removal. The Complete selector bypasses the filter, a spread can pass
through an entirely filtered stand and an unrecorded cell uses calendar
year -100 as its initial damage date.
Pass the history object unchanged between steps to retain sparse
Harvest dates. Fire dates recorded on the simulation are read automatically
and combined with any caller-supplied external dates.

In \code{landis_compatible} mode, a single return leaves the original
harvest-history date in place. The compatible single return to a
spread block uses the initial stand's recorded cells and retains
journal entries for the neighboring stands, as in the tested LANDIS-II
v7 runs. Compatible periodic spread uses the historical reserved-stand
selector behavior inferred from the HarvestManagement implementation.
}
\value{
\code{harvest_prepare()} returns the topology. \code{harvest_plan_prescriptions()}
returns a plan. \code{harvest_step()} returns the updated simulation, plan,
events, history, effects, sites and next Harvest RNG offset. When the
reservation feature is enabled, stand history includes a numeric
\code{set_aside_until} column with the exclusive expiry year or \code{NA}
if no reserve is active. Each row in
the event journal describes one stand. Stands in the same harvest block
share an \code{EventID} within the time step. In compatible mode, structured
treatments use continuous event IDs saved in the returned history. Treatments
specified only as removal tables retain per-step IDs for schedules without returns. In simulations with
scheduled returns, \code{RepeatNumber} is zero at the initial entry and
positive at a return. The journal adds the number in parentheses to the
prescription name at a return. A single return reuses the initial event
ID. \code{treatment} gives the name of the rules table applied. \code{NumberOfSites} reports
the full stand size, \code{cells} lists the selected cells, and
\code{HarvestedSites} counts cells with positive removal. Planting is
scheduled on selected cells, including visits with zero removal. With
\code{prevent_establishment = TRUE}, cells with positive removal are
delisted from future harvest and excluded from establishment, including
planting. Excluded cells retain their living cohorts and growth. With establishment directives, events also include
\code{plant_species} and \code{establishment_excluded}. The returned
simulation carries these directives into the next \code{landisR::landis_step()}.

\code{harvest_extension()} returns a descriptor for the \code{extensions}
argument of \code{landisR::landis_sim()}.
}
\seealso{\code{harvest_plan()}, \code{harvest_plan_c2()}, \code{harvest_plan_c3()},
\code{harvest_plan_d2()} for earlier qualified configurations.}
