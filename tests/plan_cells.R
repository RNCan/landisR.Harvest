index <- function(stands, partial, offsets, positions, ids) {
  .Call("h04e4d_plan_cells", stands, partial, offsets, positions, ids,
        PACKAGE = "landisR.harvest")
}
invisible(loadNamespace("landisR.harvest"))
reference <- function(stands, partial, offsets, positions, ids) {
  cells <- lapply(seq_along(stands), function(i) {
    selected <- partial[[i]]
    if (is.null(selected)) {
      stand <- stands[i]
      selected <- positions[seq.int(offsets[stand] + 1L, offsets[stand + 1L])]
    }
    as.integer(ids[selected + 1L])
  })
  list(by_event = cells, all_cells = as.integer(unlist(cells, use.names = FALSE)))
}
ids <- c(40L, 5L, 21L, 12L, 18L, 99L, 7L)
positions <- c(3L, 1L, 5L, 0L, 2L, 6L, 4L)
offsets <- c(0L, 2L, 5L, 7L)
cases <- list(
  list(integer(), list()),
  list(c(3L, 1L, 2L, 3L), vector("list", 4L)),
  list(c(3L, 1L, 2L, 3L, 1L),
       list(c(4L, 6L, 4L), NULL, integer(), c(6L, 4L), 1L))
)
for (x in cases) {
  stopifnot(identical(index(x[[1L]], x[[2L]], offsets, positions, ids),
                      reference(x[[1L]], x[[2L]], offsets, positions, ids)))
}
