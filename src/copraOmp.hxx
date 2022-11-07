#pragma once
#include <utility>
#include <vector>
#include <algorithm>
#include <omp.h>
#include "_main.hxx"
#include "vertices.hxx"
#include "edges.hxx"
#include "csr.hxx"
#include "copra.hxx"

using std::tuple;
using std::vector;
using std::make_pair;
using std::swap;




// COPRA-MOVE-ITERATION
// --------------------

/**
 * Move each vertex to its best community.
 * @param vcs communities vertex u is linked to (updated)
 * @param vcout total edge weight from vertex u to community C (updated)
 * @param vcom community set each vertex belongs to (updated)
 * @param x original graph
 * @param vtot total edge weight of each vertex
 * @param B belonging coefficient threshold
 * @returns number of changed vertices
 */
template <class G, class K, class V, size_t L, class FA, class FP>
K copraMoveIterationOmp(vector<vector<K>*>& vcs, vector<vector<V>*>& vcout, vector<Labelset<K, V, L>>& vcom, const G& x, const vector<V>& vtot, V B, FA fa, FP fp) {
  K a = K();
  K S = x.span();
  #pragma omp parallel for schedule(auto) reduction(+:a)
  for (K u=0; u<S; ++u) {
    int t = omp_get_thread_num();
    if (!x.hasVertex(u)) continue;
    if (!fa(u)) return;
    K d = vcom[u][0].first;
    copraClearScan(*vcs[t], *vcout[t]);
    copraScanCommunities(*vcs[t], *vcout[t], x, u, vcom);
    copraSortScan(*vcs[t], *vcout[t]);
    vcom[u] = copraChooseCommunity(x, u, vcom, *vcs[t], *vcout[t], B*vtot[u]);
    K c = vcom[u][0].first;
    if (c!=d) { ++a; fp(u); }
  }
  return a;
}




// COPRA-OMP
// ---------

template <size_t LABELS=COPRA_MAX_MEMBERSHIP, class G, class K, class FA, class FP>
CopraResult<K> copraSeq(const G& x, const vector<K>* q, const CopraOptions& o, FA fa, FP fp) {
  using V = typename G::edge_value_type;
  const size_t L = LABELS;
  int l = 0;
  int T = omp_get_max_threads();
  K S = x.span();
  K N = x.order();
  V B = V(1)/LABELS;
  vector<V> vtot(S);
  vector<vector<K>*> vcs(T);
  vector<vector<V>*> vcout(T);
  vector<Labelset<K, V, L>> vcom(S);
  for (int t=0; t<T; ++t) {
    vcs[t]   = new vector<K>();
    vcout[t] = new vector<V>(S);
  }
  float t = measureDuration([&]() {
    copraVertexWeights(vtot, x);
    copraInitialize(vcom, x);
    for (l=0; l<o.maxIterations;) {
      K n = copraMoveIterationOmp(vcs, vcout, vcom, x, vtot, B, fa, fp); ++l;
      PRINTFD("copraOmp(): l=%d, n=%d, N=%d, n/N=%f\n", l, n, N, float(n)/N);
      if (float(n)/N <= o.tolerance) break;
    }
  }, o.repeat);
  for (int t=0; t<T; ++t) {
    delete vcs[t];
    delete vcout[t];
  }
  return {copraBestCommunities(vcom), l, t};
}
template <size_t LABELS=COPRA_MAX_MEMBERSHIP, class G, class K, class FA>
inline CopraResult<K> copraOmp(const G& x, const vector<K>* q, const CopraOptions& o, FA fa) {
  auto fp = [](auto u) {};
  return copraOmp<LABELS>(x, q, o, fa, fp);
}
template <size_t LABELS=COPRA_MAX_MEMBERSHIP, class G, class K>
inline CopraResult<K> copraOmp(const G& x, const vector<K>* q, const CopraOptions& o) {
  auto fa = [](auto u) { return true; };
  return copraOmp<LABELS>(x, q, o, fa);
}




// COPRA-OMP-STATIC
// ----------------

template <size_t LABELS=COPRA_MAX_MEMBERSHIP, class G, class K>
inline CopraResult<K> copraOmpStatic(const G& x, const vector<K>* q=nullptr, const CopraOptions& o={}) {
  return copraOmp<LABELS>(x, q, o);
}




// COPRA-OMP-DYNAMIC-DELTA-SCREENING
// ---------------------------------

template <size_t LABELS=COPRA_MAX_MEMBERSHIP, class G, class K, class V>
inline CopraResult<K> copraOmpDynamicDeltaScreening(const G& x, const vector<tuple<K, K>>& deletions, const vector<tuple<K, K, V>>& insertions, const vector<K>* q, const CopraOptions& o={}) {
  const size_t L = LABELS;
  K S = x.span();
  const vector<Labelset<K, V, L>>& vcom = *q;
  auto vaff = copraAffectedVerticesDeltaScreening(x, deletions, insertions, vcom);
  auto fa   = [&](auto u) { return vaff[u]==true; };
  return copraOmp<LABELS>(x, q, o, fa);
}




// COPRA-OMP-DYNAMIC-FRONTIER
// --------------------------

template <size_t LABELS=COPRA_MAX_MEMBERSHIP, class G, class K, class V>
inline CopraResult<K> copraOmpDynamicFrontier(const G& x, const vector<tuple<K, K>>& deletions, const vector<tuple<K, K, V>>& insertions, const vector<K>* q, const CopraOptions& o={}) {
  const size_t L = LABELS;
  K S = x.span();
  const vector<Labelset<K, V, L>>& vcom = *q;
  auto vaff = copraAffectedVerticesFrontier(x, deletions, insertions, vcom);
  auto fa = [&](auto u) { return vaff[u]==true; };
  auto fp = [&](auto u) { x.forEachEdgeKey(u, [&](auto v) { vaff[v] = true; }); };
  return copraOmp<LABELS>(x, q, o, fa, fp);
}
