# Accept the original removal table or a treatment with establishment directives.
.harvest_treatment <- function(treatment, species) {
  plant <- character()
  prevent <- FALSE
  rules <- treatment
  if (is.list(treatment) && !is.data.frame(treatment)) {
    if (is.null(names(treatment)) || anyDuplicated(names(treatment)) ||
        any(!names(treatment) %in% c("rules", "plant", "prevent_establishment")))
      stop("Treatment fields are rules, plant and prevent_establishment.", call. = FALSE)
    rules <- treatment$rules
    if (!is.null(treatment$plant)) plant <- treatment$plant
    if (!is.null(treatment$prevent_establishment)) prevent <- treatment$prevent_establishment
  }
  if (!is.character(plant) || anyNA(plant) || any(!plant %in% species))
    stop("plant must contain known species names.", call. = FALSE)
  if (!is.logical(prevent) || length(prevent) != 1L || is.na(prevent))
    stop("prevent_establishment must be TRUE or FALSE.", call. = FALSE)
  if (identical(rules, "PlantOnly")) rules <- NULL
  list(rules = .harvest_validate_rules(rules, species),
       plant = unique(plant), prevent_establishment = prevent)
}
