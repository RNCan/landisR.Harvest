\name{harvest_legacy}
\alias{harvest_plan_d2}
\alias{harvest_prepare_adjacency_d2}
\alias{harvest_step_d2}
\alias{harvest_plan_e3}
\alias{harvest_step_e3}
\alias{harvest_plan_e4d}
\alias{harvest_step_e4d}
\title{Compatibility interfaces for earlier Harvest configurations}
\description{
\code{harvest_prepare_adjacency_d2()} adds rook adjacency to prepared
management topology. \code{harvest_plan_d2()} and \code{harvest_step_d2()}
support deterministic single-prescription configurations with adjacency
or harvest-history rules. The E3 functions retain the earlier
MaxCohortAge multi-prescription planner, and the E4d functions retain the
Random and MaxCohortAge planner and its seed-plus-offset arguments.
For the public multi-prescription interface see \code{harvest_step()}.
}
