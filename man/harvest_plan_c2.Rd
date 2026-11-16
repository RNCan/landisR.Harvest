\name{harvest_plan_c2}
\alias{harvest_plan_c2}
\alias{harvest_complete_c2}
\title{H04-C2 deterministic qualified harvesting planner}
\description{Pilot for one Complete/ClearCut prescription per management area. Supports MaxCohortAge, Economic, RegulateAges, MinimumAge, MaximumAge and ForestTypeTable. Random ranking is not yet qualified.}
\usage{
harvest_plan_c2(simulation, topology, targets, economic = NULL,
                forest_type = NULL, allocation_mode = c("corrected", "landis_compatible"))
harvest_complete_c2(simulation, topology, targets, economic = NULL,
                    forest_type = NULL, allocation_mode = c("corrected", "landis_compatible"))
}
\arguments{
\item{simulation}{An initialized landisR forest simulation.}
\item{topology}{Topology returned by harvest_prepare_management.}
\item{targets}{One row per management area with fraction, ranking, target_type, and optional minimum_age and maximum_age.}
\item{economic}{Optional species economic ranking data frame.}
\item{forest_type}{Optional data frame with management_area, rule_id, inclusion, species, age_min, age_max, percent, and highest. One row per rule/species; repeated rule_id rows form a species union.}
\item{allocation_mode}{corrected enforces true ForestType percentages and correct eligible-stands targets; landis_compatible models the observed presence-based ForestType behaviour of the tested small stands and the legacy Complete-standoverharvest issue. Exact compatibility outside tested fixtures is not claimed.}
}
\value{harvest_plan_c2 returns stand and area tables plus selected cell IDs. harvest_complete_c2 returns a harvest removal result with the plan attached.}
\details{No forest-state duplication per event. Cohorts are read directly from the shared CSR representation. No age normalization is performed by Harvest.}
