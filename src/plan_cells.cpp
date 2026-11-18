#include <Rcpp.h>
#include <vector>

namespace {
struct CellSpan {
  const int* positions = nullptr;
  R_xlen_t size = 0;
};
}

RcppExport SEXP h04e4d_plan_cells(SEXP stands_, SEXP partial_, SEXP offsets_,
                                  SEXP positions_, SEXP ids_) {
  BEGIN_RCPP
  const Rcpp::IntegerVector stands(stands_), offsets(offsets_),
    positions(positions_), ids(ids_);
  const Rcpp::List partial(partial_);
  if (partial.size() != stands.size())
    Rcpp::stop("Selected stand and explicit cell lists have different lengths.");

  // The planner stores zero-based positions in the active-cell vector.
  // A NULL selection takes the stand's full range. An empty integer vector
  // is an explicit selection of zero cells, including in a repeated entry.
  const int* stand_positions = INTEGER(positions);
  std::vector<CellSpan> spans(stands.size());
  R_xlen_t total = 0;
  for (R_xlen_t i = 0; i < stands.size(); ++i) {
    CellSpan& span = spans[i];
    SEXP explicit_cells = partial[i];
    if (Rf_isNull(explicit_cells)) {
      const int stand = stands[i];
      if (stand < 1 || stand >= offsets.size())
        Rcpp::stop("Invalid selected stand in Harvest plan.");
      const int start = offsets[stand - 1];
      const int end = offsets[stand];
      if (start < 0 || end < start || end > positions.size())
        Rcpp::stop("Invalid stand cell range in Harvest plan.");
      span = {stand_positions + start, end - start};
    } else {
      if (TYPEOF(explicit_cells) != INTSXP)
        Rcpp::stop("Explicit Harvest cell positions must be integers.");
      span = {INTEGER(explicit_cells), XLENGTH(explicit_cells)};
    }
    if (span.size > R_XLEN_T_MAX - total)
      Rcpp::stop("Too many cells in Harvest plan.");
    total += span.size;
  }

  Rcpp::List by_event(stands.size());
  Rcpp::IntegerVector all_cells(Rcpp::no_init(total));
  const int* active_ids = INTEGER(ids);
  int* flat = INTEGER(all_cells);
  for (R_xlen_t i = 0; i < stands.size(); ++i) {
    if (i % 4096 == 0) Rcpp::checkUserInterrupt();
    const CellSpan& span = spans[i];
    Rcpp::IntegerVector cells(Rcpp::no_init(span.size));
    int* local = INTEGER(cells);
    for (R_xlen_t j = 0; j < span.size; ++j) {
      const int pos = span.positions[j];
      if (pos < 0 || pos >= ids.size())
        Rcpp::stop("Invalid active cell position in Harvest plan.");
      *flat++ = local[j] = active_ids[pos];
    }
    by_event[i] = cells;
  }
  return Rcpp::List::create(Rcpp::Named("by_event") = by_event,
                            Rcpp::Named("all_cells") = all_cells);
  END_RCPP
}
