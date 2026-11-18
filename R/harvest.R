# H03-A: one deterministic treatment, on a supplied set of unique cells.
# Read the shared CSR through the core and apply one sparse cohort edit.
harvest_apply_cells <- function(simulation, cells, rules=NULL, clearcut=FALSE,
                                preserve_cohort_order = FALSE) {
  .harvest_apply_cells(simulation, cells, rules, clearcut, preserve_cohort_order)
}

.harvest_apply_cells <- function(simulation, cells, rules, clearcut,
                                 preserve_cohort_order, timer = NULL) {
  .h04_sim(simulation)
  if (!is.logical(clearcut) || length(clearcut)!=1L || is.na(clearcut)) stop("Invalid clearcut")
  if (!is.numeric(cells) || anyNA(cells) || any(!is.finite(cells)) || any(cells!=trunc(cells)) || any(cells<1 | cells>.Machine$integer.max) || anyDuplicated(cells)) stop("cells must be unique integer IDs")
  forest <- landisR::landis_csr(simulation)
  ids <- as.integer(forest$cell_id)
  positions <- match(as.integer(cells),ids)
  if (anyNA(positions)) stop("Unknown or inactive cell")
  sp <- as.character(forest$species)
  if (clearcut && !is.null(rules)) stop("clearcut and rules are mutually exclusive")
  rules <- .harvest_validate_rules(rules, sp)
  if (!is.null(timer)) timer("validation_prelevements")
  ans <- .harvest_remove_at_positions(forest, ids, positions, sp, rules, clearcut)
  if (!is.null(timer)) timer("prelevements_cpp")
  # Single commit, no expansion of unaffected cohorts; no RNG or time advance.
  if (length(ans$changed_cells)) {
    update <- landisR::landis_update(
      cohorts = ans$replacement,
      clear_cells = ans$clear_cells,
      preserve_cohort_order = preserve_cohort_order,
      metadata = .harvest_apply_metadata(ans)
    )
    if (!is.null(timer)) timer("preparation_ecriture_CSR")
    simulation <- landisR::landis_apply(simulation, update)
    if (!is.null(timer)) timer("ecriture_CSR")
    deaths <- ans$effects
    simulation <- landisR::landis_schedule_regeneration(simulation, deaths[deaths$after == 0, ])
    if (!is.null(timer)) timer("regeneration")
  }
  ans$simulation <- simulation
  ans
}

.harvest_remove_at_positions <- function(forest, ids, positions, species, rules,
                                         clearcut = FALSE) {
  .Call("h03_apply", forest$states, as.integer(forest$cell_state),
        ids, as.integer(positions), species,
        as.integer(match(rules$species, species) - 1L),
        as.numeric(rules$age_min), as.numeric(rules$age_max), as.numeric(rules$fraction),
        clearcut, PACKAGE = "landisR.harvest")
}

.harvest_apply_metadata <- function(result) {
  list(source = "harvest_h03a", changed_cells = length(result$changed_cells),
       removed_g_m2_sum = sum(result$sites$removed_g_m2))
}

.harvest_validate_rules <- function(rules, sp) {
  if (is.null(rules)) rules <- data.frame(species=character(),age_min=numeric(),age_max=numeric(),fraction=numeric())
  if (!is.data.frame(rules) || !all(c("species","age_min","age_max","fraction") %in% names(rules))) stop("Invalid rules table")
  if (!is.character(rules$species) || anyNA(rules$species) || any(!rules$species %in% sp)) stop("Unknown rule species")
  for (k in c("age_min","age_max","fraction")) if (!is.numeric(rules[[k]]) || any(!is.finite(rules[[k]]))) stop("Invalid rule values: ",k)
  if (any(rules$age_min<1 | rules$age_min!=trunc(rules$age_min) | rules$age_max<rules$age_min | rules$age_max!=trunc(rules$age_max) | rules$fraction<0 | rules$fraction>1)) stop("Invalid age range or fraction")
  for (s in unique(rules$species)) {
    r <- rules[rules$species==s,,drop=FALSE]; r <- r[order(r$age_min),,drop=FALSE]
    if (nrow(r)>1L && any(r$age_min[-1L] <= r$age_max[-nrow(r)])) stop("Overlapping rules for ",s)
  }
  rules
}
