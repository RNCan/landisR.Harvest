#include <Rcpp.h>
#include <cmath>
#include <cstdint>
#include <vector>
#include "management_e2_math.hpp"
using namespace Rcpp;
RcppExport SEXP h04e2_plan(SEXP so_,SEXP mo_,SEXP first_,SEXP ages_,SEXP rank_,
  SEXP ma_,SEXP fraction_,SEXP type_,SEXP minimum_,SEXP maximum_,SEXP begin_,SEXP end_,
  SEXP year_,SEXP seed_,SEXP offset_,SEXP legacy_) {
 BEGIN_RCPP
 IntegerVector so(so_),mo(mo_),first(first_),ages(ages_),ma(ma_),type(type_),
               minimum(minimum_),maximum(maximum_),begin(begin_),end(end_);
 NumericVector rank(rank_),fraction(fraction_),yearv(year_),seedv(seed_),offsetv(offset_);
 LogicalVector legacy(legacy_);
 const int np=ma.size();
 if(np<1||fraction.size()!=np||type.size()!=np||minimum.size()!=np||maximum.size()!=np||
      begin.size()!=np||end.size()!=np||yearv.size()!=1||seedv.size()!=1||offsetv.size()!=1||
      legacy.size()!=1||legacy[0]==NA_LOGICAL||!std::isfinite(yearv[0])||
      yearv[0]!=std::floor(yearv[0])||yearv[0]<0||yearv[0]>1e9||
      !std::isfinite(seedv[0])||seedv[0]!=std::floor(seedv[0])||seedv[0]<0||seedv[0]>4294967295.0||
      !std::isfinite(offsetv[0])||offsetv[0]!=std::floor(offsetv[0])||offsetv[0]<0||offsetv[0]>1e9)
         stop("Invalid E2 planner inputs/MT seed");
 std::vector<h04e2::Prescription> ps;
 for(int j=0;j<np;++j)ps.push_back({ma[j],minimum[j],maximum[j],begin[j],end[j],type[j],fraction[j]});
 h04e2::Result out=h04e2::select(
   as<std::vector<int>>(so),as<std::vector<int>>(mo),
   as<std::vector<int>>(first),as<std::vector<int>>(ages),
   as<std::vector<double>>(rank),ps,static_cast<int>(yearv[0]),
   static_cast<std::uint32_t>(seedv[0]),static_cast<std::uint64_t>(offsetv[0]),
   static_cast<bool>(legacy[0]));
 IntegerVector selected(out.events.size()),prescriptions(out.events.size()),stands_age(out.events.size()),
               ma_events(out.events.size());
 NumericVector u(out.events.size()),prob(out.events.size()),drawindex(out.events.size());
 for(std::size_t k=0;k<out.events.size();++k){
   const auto&e=out.events[k];selected[k]=e.stand+1;prescriptions[k]=e.prescription+1;
   stands_age[k]=e.age;ma_events[k]=e.management_area+1;u[k]=e.selection_u;
   prob[k]=e.probability;drawindex[k]=static_cast<double>(e.rng_draw_index);
 }
 return List::create(_["selected_stands"]=selected,_["prescription_indices"]=prescriptions,
    _["ages"]=stands_age,_["management_indices"]=ma_events,
    _["draw_uniform"]=u,_["draw_probability"]=prob,_["draw_indices"]=drawindex,
    _["allocated_sites"]=wrap(out.sites),_["allocated_stands"]=wrap(out.stands),
    _["eligible_stands"]=wrap(out.eligible),_["total_sites"]=wrap(out.total_sites),
    _["rng_offset_end"]=static_cast<double>(out.next_offset),
    _["draws_consumed"]=static_cast<double>(out.next_offset-static_cast<std::uint64_t>(offsetv[0])));
 END_RCPP
}
