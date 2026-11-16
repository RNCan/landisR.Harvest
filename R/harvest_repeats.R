# Dates are rounded to the harvest interval, independent of the succession step.
.harvest_repeat_options <- function(p, rules, timestep) {
  type <- if ("repeat_type" %in% names(p)) as.character(p$repeat_type) else rep("", nrow(p))
  type[is.na(type)] <- ""
  if (!any(nzchar(type)) && !any(c("repeat_interval", "times_to_repeat") %in% names(p)))
    return(NULL)
  if (any(!type %in% c("", "SingleRepeat", "MultipleRepeat")))
    stop("repeat_type must be SingleRepeat, MultipleRepeat or empty", call. = FALSE)
  single <- type == "SingleRepeat"
  multiple <- type == "MultipleRepeat"
  read_integer <- function(column) {
    x <- p[[column]]
    if (is.null(x)) return(rep.int(0L, nrow(p)))
    if (all(is.na(x))) return(rep.int(0L, nrow(p)))
    x[is.na(x)] <- 0L
    .h04_integer(x, column, 0L)
  }
  requested <- read_integer("repeat_interval")
  count <- read_integer("times_to_repeat")
  if (any((requested > 0L) != (single | multiple)))
    stop("Provide a positive repeat_interval for each repeating prescription", call. = FALSE)
  if (any(count == 1L | (count > 0L & !multiple)))
    stop("times_to_repeat must be at least two for MultipleRepeat, or zero for unlimited",
         call. = FALSE)
  additional <- if (any(single)) paste0(p$prescription[single], "(1)") else character()
  if (any(grepl("[(][0-9]+[)]$", p$prescription)) ||
      any(!additional %in% names(rules)))
    stop("SingleRepeat needs a named treatment with a (1) suffix, and base names must not end with a repeat suffix",
         call. = FALSE)
  interval <- .h04_integer(ceiling(requested / timestep) * timestep, "effective repeat_interval", 0L)
  count[single] <- 1L
  count[multiple & count == 0L] <- .Machine$integer.max
  if (!any(single | multiple)) return(NULL)
  list(kind = as.integer(single) + 2L * as.integer(multiple), interval = interval, count = count)
}

.harvest_repeat_input <- function(history, p, topology, repeat_config) {
  repeat_state <- attr(history, "harvest_repeat_state", exact = TRUE)
  if (is.null(repeat_state)) return(list(event_id = 0L, rank = 0, jobs = list()))
  if (is.null(repeat_config))
    stop("Repeat history requires repeating prescriptions", call. = FALSE)
  jobs <- repeat_state$jobs
  if (!length(jobs)) return(repeat_state)
  prescription <- vapply(jobs, function(job) job$prescription, character(1L))
  areas <- vapply(jobs, function(job) job$management_area, integer(1L))
  rows <- match(paste(areas, prescription), paste(p$management_area, p$prescription))
  if (anyNA(rows))
    stop("Repeat history does not match the prescriptions", call. = FALSE)
  # Resolve all sparse job IDs in one batch. Rebuilding a lookup for each
  # stand or job would repeatedly walk the full active-cell vector.
  cell_lists <- do.call(c, lapply(jobs, function(job) job$cells))
  cell_counts <- lengths(cell_lists)
  positions <- match(as.integer(unlist(cell_lists, use.names = FALSE)), topology$ids) - 1L
  stand_counts <- lengths(lapply(jobs, function(job) job$stands))
  stands <- match(as.integer(unlist(lapply(jobs, function(job) job$stands), use.names = FALSE)),
                  topology$stands$stand)
  if (anyNA(positions) || anyNA(stands) ||
      !identical(topology$stands$management_area[stands], rep.int(areas, stand_counts)))
    stop("Repeat history does not match stand topology", call. = FALSE)
  cell_offsets <- c(0L, cumsum(cell_counts))
  stand_offsets <- c(0L, cumsum(stand_counts))
  repeat_state$jobs <- lapply(seq_along(jobs), function(i) {
    job <- jobs[[i]]
    j <- rows[i]
    if (repeat_config$kind[j] != job$kind || repeat_config$interval[j] != job$interval ||
        repeat_config$count[j] != job$count)
      stop("Repeat history does not match prescription: ", job$prescription, call. = FALSE)
    stands_in_job <- seq.int(stand_offsets[i] + 1L, length.out = stand_counts[i])
    cells <- lapply(stands_in_job, function(k) {
      positions[seq.int(cell_offsets[k] + 1L, length.out = cell_counts[k])]
    })
    list(prescription = j - 1L, year = job$year, iteration = job$iteration, event_id = job$event_id,
         rank = job$rank, stands = stands[stands_in_job] - 1L, cells = cells)
  })
  repeat_state
}

.harvest_repeat_output <- function(native, p, topology, config) {
  jobs <- lapply(native$repeat_jobs, function(job) {
    j <- job$prescription + 1L
    list(management_area = p$management_area[j], prescription = p$prescription[j],
      kind = config$kind[j], interval = config$interval[j], count = config$count[j],
      year = job$year, iteration = job$iteration, event_id = job$event_id, rank = job$rank,
      stands = topology$stands$stand[job$stands + 1L],
      cells = lapply(job$cells, function(x) topology$ids[x + 1L]))
  })
  list(event_id = native$last_event_id, rank = native$current_rank, jobs = jobs)
}
