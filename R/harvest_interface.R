# Public interfaces for area-based harvest prescriptions.
harvest_prepare <- function(simulation, management, prescriptions = NULL) {
  topology <- harvest_prepare_adjacency_d2(simulation, harvest_prepare_management(simulation, management))
  if (any(.harvest_is_partial(prescriptions$site_selection) |
          prescriptions$site_selection %in% "PatchCutting"))
    topology <- .harvest_prepare_cell_topology(simulation, topology)
  topology
}

harvest_plan_prescriptions <- function(simulation, topology, prescriptions, rules_by_prescription,
                                       history = NULL, year = simulation$year,
                                       seed = simulation$seed, rng_offset = 0,
                                       allocation_mode = c("corrected", "landis_compatible"),
                                       harvest_timestep = simulation$project$canonical$timestep,
                                       damage_history = NULL) {
  if (is.null(rules_by_prescription))
    stop("Fournir les r\xe8gles de pr\xe9l\xe8vement des prescriptions.", call. = FALSE)
  harvest_plan_e4d(simulation, topology, prescriptions, history = history, year = year,
                   landis_seed = seed, landis_rng_offset = rng_offset,
                   allocation_mode = match.arg(allocation_mode),
                   rules_by_prescription = rules_by_prescription, harvest_timestep = harvest_timestep,
                   damage_history = damage_history)
}

harvest_step <- function(simulation, topology, prescriptions, rules_by_prescription,
                         history = NULL, year = simulation$year,
                         seed = simulation$seed, rng_offset = 0,
                         allocation_mode = c("corrected", "landis_compatible"),
                         harvest_timestep = simulation$project$canonical$timestep,
                                       damage_history = NULL) {
  harvest_step_e4d(simulation, topology, prescriptions, rules_by_prescription,
                   history = history, year = year, landis_seed = seed,
                   landis_rng_offset = rng_offset, allocation_mode = match.arg(allocation_mode),
                   harvest_timestep = harvest_timestep,
                   damage_history = damage_history)
}

harvest_extension <- function(management, prescriptions, rules_by_prescription,
                              timestep = NULL, seed = NULL, rng_offset = 0,
                              allocation_mode = c("corrected", "landis_compatible"),
                              damage_history = NULL) {
  allocation_mode <- match.arg(allocation_mode)
  force(management)
  force(prescriptions)
  force(rules_by_prescription)
  force(timestep)
  force(seed)
  force(rng_offset)
  force(damage_history)
  if (!is.list(rules_by_prescription) || is.null(names(rules_by_prescription)) ||
      anyNA(names(rules_by_prescription)) || any(!nzchar(names(rules_by_prescription))) ||
      anyDuplicated(names(rules_by_prescription)))
    stop("rules_by_prescription doit \xeatre une liste de traitements nomm\xe9s.", call. = FALSE)

  initialize <- function(simulation) {
    core_step <- simulation$project$canonical$timestep
    interval <- if (is.null(timestep)) core_step else .h04_integer(timestep, "timestep")
    if (length(interval) != 1L || interval %% core_step != 0L)
      stop("timestep Harvest doit \xeatre un multiple positif du pas Succession.", call. = FALSE)
    topology <- harvest_prepare(simulation, management, prescriptions)
    harvest_seed <- if (is.null(seed)) simulation$seed else seed
    # The preview validates planner parameters without updating either RNG state.
    harvest_plan_prescriptions(simulation, topology, prescriptions, rules_by_prescription,
                               year = simulation$year, seed = harvest_seed, rng_offset = rng_offset,
                               allocation_mode = allocation_mode, harvest_timestep = interval)
    active <- prescriptions$fraction > 0
    single <- prescriptions$repeat_type %in% "SingleRepeat"
    names_to_validate <- unique(c(as.character(prescriptions$prescription[active]),
      if (any(single)) paste0(prescriptions$prescription[single], "(1)")))
    species <- as.character(landisR::landis_csr(simulation)$species)
    for (name in names_to_validate) {
      if (!name %in% names(rules_by_prescription))
        stop("Traitement manquant pour la prescription ", name, ".", call. = FALSE)
      .harvest_treatment(rules_by_prescription[[name]], species)
    }
    list(topology = topology, history = NULL, seed = harvest_seed,
         rng_offset = rng_offset, timestep = interval)
  }

  step <- function(simulation, state, year) {
    if (year %% state$timestep != 0L)
      return(list(simulation = simulation, state = state))
    result <- harvest_step(simulation, state$topology, prescriptions, rules_by_prescription,
                           history = state$history, year = year, seed = state$seed,
                           rng_offset = if (is.null(simulation$landis_rng)) state$rng_offset else 0,
                           allocation_mode = allocation_mode, harvest_timestep = state$timestep,
                           damage_history = if (is.function(damage_history))
                             damage_history(simulation, year) else damage_history)
    state$history <- result$history
    state$rng_offset <- result$landis_rng_offset_end
    list(simulation = result$simulation, state = state, events = result$events)
  }

  landisR::landis_extension(step, initialize, phase = "before_succession")
}
