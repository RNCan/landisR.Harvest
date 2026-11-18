#include <Rcpp.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>
#include <unordered_map>
#include <vector>
#include <algorithm>
#include <memory>
#include "harvest_math.hpp"
using namespace Rcpp;
namespace {
struct Item { int sp; double age, before, take; };
struct Outcome { std::vector<Item> cohorts; std::vector<int> effects; double total=0; int full=0, partial=0; };
List run(List states, IntegerVector cell_state, IntegerVector ids,
         IntegerVector positions, CharacterVector names, IntegerVector rsp,
         NumericVector lo, NumericVector hi, NumericVector frac, bool clearcut) {
  IntegerVector off=as<IntegerVector>(states["offset"]), sp=as<IntegerVector>(states["species_id"]);
  NumericVector age=as<NumericVector>(states["age"]), bio=as<NumericVector>(states["biomass"]);
  const int ns=off.size()-1, nsp=names.size();
  if(ns<1 || off[0]!=0 || off[ns]!=sp.size() || age.size()!=sp.size() || bio.size()!=sp.size() || cell_state.size()!=ids.size()) stop("Invalid CSR layout");
  for(int s=0;s<ns;++s) if(off[s]<0 || off[s+1]<off[s] || off[s+1]>sp.size()) stop("Invalid CSR offsets");
  if(rsp.size()!=lo.size() || lo.size()!=hi.size() || hi.size()!=frac.size()) stop("Invalid rules");
  std::vector<std::vector<harvest_h03::Rule>> rules(nsp);
  for(int i=0;i<rsp.size();++i) {
    int s=rsp[i]; if(s<0 || s>=nsp) stop("Invalid rule species");
    rules[s].push_back({lo[i],hi[i],frac[i]});
  }
  std::unordered_map<int,Outcome> cache;
  std::vector<int> rc, ec, selected, changed, empty, nf, np;
  std::vector<std::string> rs, es;
  std::vector<double> ra, rb, ea, eb, er, ef, totals;
  for(int p=0;p<positions.size();++p) {
    if(p%4096==0) checkUserInterrupt();
    int pos=positions[p]-1;
    if(pos<0 || pos>=ids.size()) stop("Invalid active position");
    int s=cell_state[pos]-1, id=ids[pos];
    if(s<0 || s>=ns) stop("Invalid state index");
    auto it=cache.find(s);
    if(it==cache.end()) {
      Outcome out;
      for(int j=off[s];j<off[s+1];++j) {
        int k=sp[j]; double a=age[j], b=bio[j];
        if(k<0 || k>=nsp || !R_finite(a) || a<1 || a!=std::trunc(a)) stop("Invalid cohort");
        double f = clearcut ? 1.0 : harvest_h03::removal_fraction(a, rules[k]);
        double take=harvest_h03::removed(b,f);
        out.cohorts.push_back({k,a,b,take});
        if(take>0) {out.total+=take; if(take==b) ++out.full; else ++out.partial;}
      }
      // Preserve complete-removal phase before the partial-removal phase.
      for(int phase=0;phase<2;++phase) for(int j=0;j<(int)out.cohorts.size();++j) {
        const Item& x=out.cohorts[j];
        if(x.take>0 && ((phase==0 && x.take==x.before)||(phase==1 && x.take<x.before))) out.effects.push_back(j);
      }
      it=cache.emplace(s,std::move(out)).first;
    }
    const Outcome& out=it->second;
    selected.push_back(id); totals.push_back(out.total); nf.push_back(out.full); np.push_back(out.partial);
    if(out.total==0) continue;
    changed.push_back(id); bool survivors=false;
    for(const Item& x:out.cohorts) if(x.before>x.take) {
      survivors=true; rc.push_back(id); rs.push_back(as<std::string>(names[x.sp])); ra.push_back(x.age); rb.push_back(x.before-x.take);
    }
    if(!survivors) empty.push_back(id);
    for(int j:out.effects) {
      const Item& x=out.cohorts[j]; ec.push_back(id); es.push_back(as<std::string>(names[x.sp])); ea.push_back(x.age); eb.push_back(x.before); er.push_back(x.take); ef.push_back(x.before-x.take);
    }
  }
  return List::create(_["replacement"]=DataFrame::create(_["cell_id"]=rc,_["species"]=rs,_["age"]=ra,_["biomass"]=rb),
    _["effects"]=DataFrame::create(_["cell_id"]=ec,_["species"]=es,_["age"]=ea,_["before"]=eb,_["removed"]=er,_["after"]=ef),
    _["sites"]=DataFrame::create(_["cell_id"]=selected,_["removed_g_m2"]=totals,_["full"]=nf,_["partial"]=np),
    _["changed_cells"]=changed,_["clear_cells"]=empty,_["states_evaluated"]=static_cast<int>(cache.size()));
}
}
RcppExport SEXP h03_apply(SEXP a,SEXP b,SEXP c,SEXP d,SEXP e,SEXP f,SEXP g,SEXP h,SEXP i,SEXP j) {
  BEGIN_RCPP
  return run(as<List>(a),as<IntegerVector>(b),as<IntegerVector>(c),as<IntegerVector>(d),as<CharacterVector>(e),as<IntegerVector>(f),as<NumericVector>(g),as<NumericVector>(h),as<NumericVector>(i),as<bool>(j));
  END_RCPP
}
extern "C" SEXP h04_plan(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
extern "C" SEXP h04b2_plan(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
extern "C" SEXP h04c2_plan(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
extern "C" SEXP h04c3_random_plan(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
extern "C" SEXP h05_cell_topology(SEXP,SEXP,SEXP,SEXP);
extern "C" SEXP h04d2_topology(SEXP,SEXP,SEXP,SEXP);
extern "C" SEXP h04d2_plan(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
extern "C" SEXP h04e2_plan(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
extern "C" SEXP h04e3_plan(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
extern "C" SEXP h04e4d_plan(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
extern "C" SEXP h04e4d_plan_cells(SEXP,SEXP,SEXP,SEXP,SEXP);
static const R_CallMethodDef calls[]={{"h03_apply",(DL_FUNC)&h03_apply,10},{"h04_plan",(DL_FUNC)&h04_plan,7},{"h04b2_plan",(DL_FUNC)&h04b2_plan,12},{"h04c2_plan",(DL_FUNC)&h04c2_plan,23},{"h04c3_random_plan",(DL_FUNC)&h04c3_random_plan,11},{"h04d2_topology",(DL_FUNC)&h04d2_topology,4},{"h05_cell_topology",(DL_FUNC)&h05_cell_topology,4},{"h04d2_plan",(DL_FUNC)&h04d2_plan,16},{"h04e2_plan",(DL_FUNC)&h04e2_plan,16},{"h04e3_plan",(DL_FUNC)&h04e3_plan,13},{"h04e4d_plan",(DL_FUNC)&h04e4d_plan,13},{"h04e4d_plan_cells",(DL_FUNC)&h04e4d_plan_cells,5},{NULL,NULL,0}};
extern "C" void attribute_visible R_init_landisR_harvest(DllInfo* dll) {
 R_registerRoutines(dll,NULL,calls,NULL,NULL); R_useDynamicSymbols(dll,FALSE);
}
