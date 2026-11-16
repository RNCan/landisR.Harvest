# Map the sparse external or saved dates onto active cell positions.
.harvest_damage_positions <- function(records, ids, year) {
  if (is.null(records)) return(list(position = integer(), year = integer()))
  if (!is.data.frame(records) || !all(c("cell_id", "year") %in% names(records)))
    stop("damage_history requires cell_id and year", call. = FALSE)
  id <- .h04_integer(records$cell_id, "damage cell_id")
  time <- .h04_integer(records$year, "damage year", 0L)
  position <- match(id, ids) - 1L
  if (anyNA(position) || any(time > year))
    stop("Damage history must reference active cells and years up to the harvest year", call. = FALSE)
  list(position = position, year = time)
}
