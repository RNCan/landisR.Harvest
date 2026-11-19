if (requireNamespace("terra", quietly = TRUE)) local({
  root <- tempfile("harvest-read-")
  scenario <- file.path(root, "scenario")
  dir.create(scenario, recursive = TRUE)
  on.exit(unlink(root, recursive = TRUE))

  grid <- data.frame(cell_id = 1:6, row = rep(1:2, 3), col = rep(1:3, each = 2),
                     active = c(rep(TRUE, 5), FALSE))
  project <- landisR::landis_build_landscape(grid, duration = 20L, cell_length = 100,
    species = data.frame(species = "PINE"))
  template <- terra::rast(nrow = 2, ncol = 3, xmin = 0, xmax = 300,
                          ymin = 0, ymax = 200, crs = "EPSG:32618")
  areas <- terra::setValues(template, c(1, 2, 2, NA, 0, NA))
  stands <- terra::setValues(template, c(11, 21, 22, 12, 23, NA))
  terra::writeRaster(areas, file.path(root, "management areas.tif"))
  terra::writeRaster(stands, file.path(root, "stands.tif"))

  path <- file.path(root, "harvest.txt")
  lines <- c(
    'LandisData "Biomass Harvest"', "Timestep 10",
    'ManagementAreas "..\\management areas.tif"', 'Stands "..\\stands.tif"',
    "Prescription CUT", "StandRanking Economic", "PINE 10 40",
    "MinimumAge 25", "MinimumTimeSinceLastHarvest 30", "SiteSelection CompleteStandSpread 0 50",
    "CohortsRemoved SpeciesList", "PINE 40-150(50%)",
    "HarvestImplementations", "1 CUT 100%", "2 CUT 2.5%",
    'EventLog "harvest.log"'
  )
  writeLines(lines, path)
  harvest <- landisR.harvest::harvest_read(path, project, scenario_directory = scenario)
  stopifnot(inherits(harvest, "landisR_harvest_parameters"),
            identical(names(harvest), c("management", "prescriptions",
                                      "rules_by_prescription", "timestep")),
            identical(harvest$timestep, 10L),
            identical(harvest$management$management_area, c(1L, 3L, 2L, 3L, 2L)),
            identical(harvest$management$stand, c(11L, 12L, 21L, 23L, 22L)),
            identical(harvest$prescriptions$fraction, c(1, 0.025)),
            harvest$rules_by_prescription$CUT$fraction == 0.5)

  sim <- landisR::landis_initialize(project)
  cohorts <- data.frame(cell_id = 1:5, species = "PINE", age = 70L, biomass = 1000)
  sim <- landisR::landis_apply(sim, landisR::landis_update(cohorts = cohorts))
  run <- function(extension) {
    initial <- landisR::landis_sim(sim, steps = 0L,
      extensions = list(harvest = extension), progress = FALSE)
    landisR::landis_sim(initial, steps = 1L, progress = FALSE)
  }
  imported <- run(landisR.harvest::harvest_extension(
    parameters = harvest, allocation_mode = "landis_compatible"))
  direct <- run(landisR.harvest::harvest_extension(
    management = harvest$management, prescriptions = harvest$prescriptions,
    rules_by_prescription = harvest$rules_by_prescription, timestep = harvest$timestep,
    allocation_mode = "landis_compatible"))
  stopifnot(any(imported$events$harvest[["10"]]$HarvestedSites > 0L),
            identical(imported$events, direct$events),
            identical(landisR::landis_csr(imported$simulation),
                      landisR::landis_csr(direct$simulation)))

  terra::writeRaster(terra::shift(stands, dx = 100), file.path(root, "stands.tif"),
                     overwrite = TRUE)
  fail <- tryCatch(landisR.harvest::harvest_read(path, project, scenario),
                   error = conditionMessage)
  stopifnot(is.character(fail), grepl("align with the reference grid", fail))

  writeLines(sub("StandRanking Economic", "StandRanking Random", lines, fixed = TRUE), path)
  fail <- tryCatch(landisR.harvest::harvest_read(path, project, scenario),
                   error = conditionMessage)
  stopifnot(is.character(fail), grepl("Unsupported Harvest prescription format", fail))
})
