harvest_read <- function(path, project, scenario_directory = dirname(path),
                         raster_template = NULL) {
  if (!inherits(project, "landisR_project"))
    stop("Supply a landisR project with the landscape grid.", call. = FALSE)
  if (!requireNamespace("terra", quietly = TRUE))
    stop("Install terra to import LANDIS-II Harvest maps.", call. = FALSE)

  lines <- readLines(path, warn = FALSE)
  if (any(!validUTF8(lines))) Encoding(lines) <- "latin1"
  lines <- trimws(sub("<<.*$", "", sub("^\ufeff", "", enc2utf8(lines))))
  lines <- lines[nzchar(lines) & !startsWith(lines, ">>")]
  rows <- lapply(lines, function(line)
    scan(text = line, what = character(), quote = '"', comment.char = "", quiet = TRUE))
  keys <- vapply(rows, function(row) row[1L], character(1L))
  if (!length(rows) || !identical(rows[[1L]], c("LandisData", "Biomass Harvest")))
    stop("Expected a LANDIS-II Biomass Harvest parameter file.", call. = FALSE)

  starts <- which(keys == "Prescription")
  implementation <- which(keys == "HarvestImplementations")
  if (length(implementation) != 1L || !length(starts) || any(starts >= implementation))
    stop("Expected Prescription blocks followed by HarvestImplementations.", call. = FALSE)
  header <- rows[seq_len(starts[1L] - 1L)]
  value <- function(key) {
    index <- which(keys[seq_along(header)] == key)
    if (length(index) != 1L || length(header[[index]]) != 2L)
      stop("Expected one value for Harvest ", key, ".", call. = FALSE)
    header[[index]][2L]
  }
  timestep <- .harvest_read_number(value("Timestep"), "Timestep", whole = TRUE)
  if (timestep < 1L) stop("Harvest Timestep must be positive.", call. = FALSE)

  output_keys <- c("PrescriptionMaps", "BiomassMaps", "EventLog", "SummaryLog")
  end <- which(seq_along(rows) > implementation & keys %in% output_keys)
  end <- if (length(end)) end[1L] else length(rows) + 1L
  if (end <= implementation + 1L)
    stop("HarvestImplementations needs at least one assignment.", call. = FALSE)
  if (end <= length(rows) && any(!keys[seq.int(end, length(rows))] %in% output_keys))
    stop("Unexpected input after Harvest output settings.", call. = FALSE)

  blocks <- Map(function(first, last) {
    .harvest_read_prescription(rows[seq.int(first, last)], project$canonical$duration)
  }, starts, c(starts[-1L] - 1L, implementation - 1L))
  names(blocks) <- vapply(blocks, function(block) block$prescription$prescription, character(1L))
  if (anyDuplicated(names(blocks)))
    stop("Duplicate Harvest prescription names.", call. = FALSE)
  plans <- rows[seq.int(implementation + 1L, end - 1L)]
  if (any(lengths(plans) != 3L))
    stop("Expected management area, prescription and percent area in HarvestImplementations.",
         call. = FALSE)
  id <- .harvest_read_number(vapply(plans, `[[`, "", 1L), "management area", whole = TRUE)
  name <- vapply(plans, `[[`, "", 2L)
  fraction <- .harvest_read_percent(vapply(plans, `[[`, "", 3L), "area target")
  if (any(id == 0L))
    stop("HarvestImplementations requires positive management area codes.", call. = FALSE)
  if (any(!name %in% names(blocks)))
    stop("Unknown prescription in HarvestImplementations: ",
         paste(unique(name[!name %in% names(blocks)]), collapse = ", "), ".", call. = FALSE)
  prescriptions <- cbind(data.frame(management_area = id, fraction = fraction),
                          do.call(rbind, unname(lapply(blocks[name], `[[`, "prescription"))))
  rownames(prescriptions) <- NULL
  prescriptions$economic <- unname(lapply(blocks[name], `[[`, "economic"))

  map_path <- function(key) {
    name <- chartr("\\", "/", value(key))
    if (grepl("^(/|[A-Za-z]:/)", name)) return(name)
    file.path(scenario_directory, name)
  }
  areas <- terra::rast(map_path("ManagementAreas"))
  stands <- terra::rast(map_path("Stands"))
  if (is.null(raster_template) && is.list(project$source)) {
    source <- project$source
    raster_template <- if (is.list(source$landscape)) source$landscape$ecoregions_map
    if (is.null(raster_template)) raster_template <- source$ecoregions_map
  }
  if (is.null(raster_template)) raster_template <- areas
  if (is.character(raster_template)) raster_template <- terra::rast(raster_template)

  grid <- project$canonical$landscape
  if (!inherits(raster_template, "SpatRaster") || terra::nlyr(raster_template) != 1L ||
      terra::nrow(raster_template) != grid$nrow || terra::ncol(raster_template) != grid$ncol)
    stop("The Harvest raster template must have one layer and match the landscape dimensions.",
         call. = FALSE)
  maps <- list(ManagementAreas = areas, Stands = stands)
  for (key in names(maps)) {
    if (terra::nlyr(maps[[key]]) != 1L || !.harvest_read_aligned(raster_template, maps[[key]]))
      stop("Harvest ", key, " must have one layer and align with the reference grid.",
           call. = FALSE)
  }

  cells <- grid$cells[grid$ids, c("cell_id", "row", "col"), drop = FALSE]
  positions <- (cells$row - 1) * grid$ncol + cells$col
  area <- terra::values(areas, mat = FALSE)[positions]
  area[is.na(area)] <- 0
  area <- .harvest_read_number(area, "management area map", whole = TRUE)
  if (any(area == 0L)) {
    unmanaged <- max(c(area, prescriptions$management_area)) + 1
    if (unmanaged > .Machine$integer.max)
      stop("No available management area code for unmanaged sites.", call. = FALSE)
    area[area == 0L] <- as.integer(unmanaged)
  }
  stand <- .harvest_read_number(terra::values(stands, mat = FALSE)[positions],
                                "stand map", whole = TRUE)
  if (any(stand == 0L))
    stop("Each active landscape cell requires a positive stand code.", call. = FALSE)

  harvest_build(
    timestep = timestep,
    management = data.frame(cell_id = cells$cell_id, management_area = area, stand = stand),
    prescriptions = prescriptions,
    rules_by_prescription = lapply(blocks, `[[`, "rules")
  )
}

.harvest_read_number <- function(x, label, whole = FALSE) {
  value <- suppressWarnings(as.numeric(x))
  if (any(!is.finite(value) | value < 0) ||
      (whole && any(value != trunc(value) | value > .Machine$integer.max)))
    stop("Invalid Harvest ", label, ".", call. = FALSE)
  if (whole) as.integer(value) else value
}

.harvest_read_percent <- function(x, label) {
  if (any(!grepl("^[0-9]+(\\.[0-9]+)?%$", x)))
    stop("Expected a percentage for Harvest ", label, ".", call. = FALSE)
  values <- .harvest_read_number(sub("%$", "", x), label)
  if (any(values > 100)) stop("Harvest ", label, " must be at most 100%.", call. = FALSE)
  values / 100
}

.harvest_read_prescription <- function(rows, duration) {
  keys <- vapply(rows, function(row) row[1L], character(1L))
  name <- rows[[1L]][2L]
  required <- c("StandRanking", "MinimumAge", "SiteSelection", "CohortsRemoved")
  i <- match(required, keys)
  history <- which(keys == "MinimumTimeSinceLastHarvest")
  if (is.na(name) || length(rows[[1L]]) != 2L || anyNA(i) ||
      anyDuplicated(keys[keys %in% required]) || length(history) > 1L ||
      i[1L] != 2L || i[2L] <= i[1L] + 1L || i[4L] == length(rows) ||
      !identical(i[3:4], i[2L] + length(history) + 1:2) ||
      (length(history) && history != i[2L] + 1L) ||
      !identical(rows[[i[1L]]], c("StandRanking", "Economic")) ||
      !identical(rows[[i[4L]]], c("CohortsRemoved", "SpeciesList")) ||
      length(rows[[i[2L]]]) != 2L || length(rows[[i[3L]]]) != 4L ||
      rows[[i[3L]]][2L] != "CompleteStandSpread")
    stop("Unsupported Harvest prescription format: ", name,
         ". Expected Economic ranking, MinimumAge, optional MinimumTimeSinceLastHarvest, ",
         "CompleteStandSpread and SpeciesList.", call. = FALSE)

  economics <- rows[seq.int(i[1L] + 1L, i[2L] - 1L)]
  if (any(lengths(economics) != 3L))
    stop("Expected species, rank and merchantable age in ", name, ".", call. = FALSE)
  economic <- data.frame(species = vapply(economics, `[[`, "", 1L),
                         rank = .harvest_read_number(vapply(economics, `[[`, "", 2L),
                                                    "economic rank", whole = TRUE),
                         merchantable_age = .harvest_read_number(
                           vapply(economics, `[[`, "", 3L), "merchantable age", whole = TRUE))
  removal <- rows[seq.int(i[4L] + 1L, length(rows))]
  if (any(lengths(removal) != 2L))
    stop("Expected species and an age range with optional removal percent in ", name, ".",
         call. = FALSE)
  ranges <- vapply(removal, `[[`, "", 2L)
  bare <- grepl("^[0-9]+-[0-9]+$", ranges)
  ranges[bare] <- paste0(ranges[bare], "(100%)")
  pieces <- regmatches(ranges, regexec("^([0-9]+)-([0-9]+)\\(([0-9.]+%)\\)$", ranges))
  if (any(lengths(pieces) != 4L))
    stop("Invalid Harvest cohort removal range in ", name, ".", call. = FALSE)
  rules <- data.frame(species = vapply(removal, `[[`, "", 1L),
                      age_min = .harvest_read_number(vapply(pieces, `[[`, "", 2L),
                                                    "minimum removal age", whole = TRUE),
                      age_max = .harvest_read_number(vapply(pieces, `[[`, "", 3L),
                                                    "maximum removal age", whole = TRUE),
                      fraction = .harvest_read_percent(vapply(pieces, `[[`, "", 4L),
                                                       "cohort removal"))
  if (any(rules$age_min > rules$age_max))
    stop("Reversed Harvest cohort removal age range in ", name, ".", call. = FALSE)
  bounds <- .harvest_read_number(rows[[i[3L]]][3:4], "block area")
  if (bounds[2L] <= 0 || bounds[1L] > bounds[2L])
    stop("Invalid Harvest block area in ", name, ".", call. = FALSE)
  minimum <- .harvest_read_number(rows[[i[2L]]][2L], "minimum age", whole = TRUE)
  if (length(history) && length(rows[[history]]) != 2L)
    stop("Invalid MinimumTimeSinceLastHarvest in ", name, ".", call. = FALSE)
  waiting <- if (length(history)) .harvest_read_number(rows[[history]][2L],
                                                        "minimum time since harvest",
                                                        whole = TRUE) else 0L
  list(prescription = data.frame(
    prescription = name, ranking = "Economic", target_type = "area",
    minimum_age = minimum, maximum_age = 999L, begin_time = 0L,
    end_time = as.integer(duration), site_selection = "CompleteStandSpread",
    min_block_ha = bounds[1L], max_block_ha = bounds[2L],
    minimum_time_since_last_harvest = waiting
  ), economic = economic, rules = rules)
}

.harvest_read_aligned <- function(template, map) {
  if (!terra::compareGeom(template, map, res = TRUE, crs = FALSE, stopOnError = FALSE))
    return(FALSE)
  if (terra::compareGeom(template, map, res = TRUE, stopOnError = FALSE, messages = FALSE))
    return(TRUE)
  # The 6est Lambert grids encode equivalent standard parallels in different order.
  parameters <- function(x) strsplit(trimws(terra::crs(x, proj = TRUE)), "[[:space:]]+")[[1L]]
  a <- parameters(template)
  b <- parameters(map)
  lcc_key <- function(x) sort(sub("^[+]lat_[12]=", "+lat_parallel=", x))
  "+proj=lcc" %in% a && "+proj=lcc" %in% b && identical(lcc_key(a), lcc_key(b))
}
