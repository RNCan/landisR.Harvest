.harvest_elapsed <- function() unname(proc.time()[["elapsed"]])

.harvest_record_timing <- function(sink, level, year, times) {
  if (!length(times)) return(invisible(NULL))
  sink$phases[[length(sink$phases) + 1L]] <- data.frame(
    annee = as.integer(year), niveau = level, phase = names(times),
    secondes = unname(as.numeric(times))
  )
  invisible(NULL)
}

.harvest_timer <- function(level, year) {
  sink <- getOption("landisR.harvest.timing")
  if (!is.environment(sink)) return(function(phase = NULL) invisible(NULL))
  force(level)
  force(year)
  last <- .harvest_elapsed()
  times <- numeric()
  function(phase = NULL) {
    if (is.null(phase)) {
      .harvest_record_timing(sink, level, year, times)
    } else {
      elapsed <- .harvest_elapsed() - last
      previous <- times[phase]
      times[phase] <<- if (is.na(previous)) elapsed else previous + elapsed
      last <<- .harvest_elapsed()
    }
    invisible(NULL)
  }
}
