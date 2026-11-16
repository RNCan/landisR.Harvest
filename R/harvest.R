# H03-A: one deterministic treatment, on a supplied set of unique cells.
# Internal layout is version-pinned; public application is through landis_apply.
harvest_apply_cells <- function(simulation, cells, rules=NULL, clearcut=FALSE,
                                preserve_cohort_order = FALSE) {
  .h04_sim(simulation)
  if (!is.logical(clearcut) || length(clearcut)!=1L || is.na(clearcut)) stop("Invalid clearcut")
  if (!is.numeric(cells) || anyNA(cells) || any(!is.finite(cells)) || any(cells!=trunc(cells)) || any(cells<1 | cells>.Machine$integer.max) || anyDuplicated(cells)) stop("cells must be unique integer IDs")
  ids <- as.integer(simulation$ctx$ids)
  positions <- match(as.integer(cells),ids)
  if (anyNA(positions)) stop("Unknown or inactive cell")
  sp <- as.character(simulation$ctx$config$species$species)
  if (clearcut && !is.null(rules)) stop("clearcut and rules are mutually exclusive")
  rules <- .harvest_validate_rules(rules, sp)
  ans <- .Call("h03_apply",simulation$states,as.integer(simulation$cell_state),ids,as.integer(positions),sp,as.integer(match(rules$species,sp)-1L),as.numeric(rules$age_min),as.numeric(rules$age_max),as.numeric(rules$fraction),clearcut,PACKAGE="landisR.harvest")
  # Single commit, no expansion of unaffected cohorts; no RNG or time advance.
  if (length(ans$changed_cells)) {
    update <- landisR::landis_update(
      cohorts = ans$replacement,
      clear_cells = ans$clear_cells,
      preserve_cohort_order = preserve_cohort_order,
      metadata = list(source = "harvest_h03a", changed_cells = length(ans$changed_cells),
                      removed_g_m2_sum = sum(ans$sites$removed_g_m2))
    )
    simulation <- landisR::landis_apply(simulation, update)
  }
  ans$simulation <- simulation
  ans
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
