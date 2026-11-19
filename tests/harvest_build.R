library(landisR.harvest)

grid <- data.frame(cell_id = 1:2, row = 1L, col = 1:2, active = TRUE)
project <- landisR::landis_build_landscape(grid, duration = 20L,
  cell_length = 100, species = data.frame(species = "PINE"))
sim <- landisR::landis_initialize(project)
cohorts <- data.frame(cell_id = 1:2, species = "PINE", age = c(40L, 20L),
                      biomass = c(1000, 1500))
sim <- landisR::landis_apply(sim, landisR::landis_update(cohorts = cohorts))

management <- data.frame(cell_id = 1:2, management_area = 1L, stand = 1:2)
prescriptions <- data.frame(management_area = 1L, prescription = "CUT",
  fraction = 0.5, ranking = "Random", target_type = "area", minimum_age = 0L,
  maximum_age = 120L, begin_time = 10L, end_time = 20L, site_selection = "Complete")
rules <- list(CUT = data.frame(species = "PINE", age_min = 1L, age_max = 120L,
                               fraction = 0.5))
parameters <- harvest_build(10L, management, prescriptions, rules)
stopifnot(inherits(parameters, "landisR_harvest_parameters"))

run <- function(extension) {
  initial <- landisR::landis_sim(sim, steps = 0L, extensions = list(harvest = extension),
                                 progress = FALSE)
  landisR::landis_sim(initial, steps = 1L, progress = FALSE)
}
from_tables <- run(harvest_extension(management, prescriptions, rules, 10L))
from_parameters <- run(harvest_extension(parameters = parameters))
stopifnot(
  nrow(from_tables$events$harvest[["10"]]) == 1L,
  identical(from_tables$events, from_parameters$events),
  identical(from_tables$extension_states, from_parameters$extension_states),
  identical(landisR::landis_csr(from_tables$simulation),
            landisR::landis_csr(from_parameters$simulation))
)

expect_error <- function(expr, message) {
  result <- tryCatch(force(expr), error = conditionMessage)
  stopifnot(is.character(result), grepl(message, result, fixed = TRUE))
}
expect_error(harvest_build(0L, management, prescriptions, rules), "Invalid timestep")
expect_error(harvest_build(10L, management, transform(prescriptions, maximum_age = -1L),
                           rules), "Invalid maximum_age")
expect_error(harvest_build(10L, management, prescriptions, list(OTHER = rules$CUT)),
             "Traitement manquant")
expect_error(harvest_extension(management, parameters = parameters), "Supply either")
wrong_species <- harvest_build(10L, management, prescriptions,
  list(CUT = transform(rules$CUT, species = "UNKNOWN")))
expect_error(run(harvest_extension(parameters = wrong_species)), "Unknown rule species")
