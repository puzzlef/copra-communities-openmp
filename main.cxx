#include <utility>
#include <vector>
#include <string>
#include <cstdio>
#include <iostream>
#include <omp.h>
#include "src/main.hxx"

using namespace std;




// You can define datatype with -DTYPE=...
#ifndef TYPE
#define TYPE float
#endif
// You can define number of threads with -DMAX_THREADS=...
#ifndef MAX_THREADS
#define MAX_THREADS 12
#endif




template <class G, class K, class V>
double getModularity(const G& x, const CopraResult<K>& a, V M) {
  auto fc = [&](auto u) { return a.membership[u]; };
  return modularityBy(x, fc, M, V(1));
}


template <class G>
void runExperiment(const G& x, int repeat) {
  using K = typename G::key_type;
  using V = typename G::edge_value_type;
  vector<K> *init = nullptr;
  auto M = edgeWeight(x)/2;
  auto Q = modularity(x, M, 1.0f);
  printf("[%01.6f modularity] noop\n", Q);
  CopraOptions o = {repeat};

  for (int i=0, f=10; f<=10000; f*=i&1? 5:2, ++i) {
    float tolerance = 1.0f / f;
    {
      // Find COPRA using a single thread (1 label).
      auto ak = copraSeqStatic<1>(x, init, {repeat, tolerance});
      printf("[%09.3f ms; %04d iters.; %01.9f modularity] copraSeqStatic {labels=%02d, tolerance=%.0e}\n", ak.time, ak.iterations, getModularity(x, ak, M), 1, tolerance);
      // Find COPRA using multiple threads (1 label).
      auto al = copraOmpStatic<1>(x, init, {repeat, tolerance});
      printf("[%09.3f ms; %04d iters.; %01.9f modularity] copraOmpStatic {labels=%02d, tolerance=%.0e}\n", al.time, al.iterations, getModularity(x, al, M), 1, tolerance);
    }
    {
      // Find COPRA using a single thread (2 labels).
      auto ak = copraSeqStatic<2>(x, init, {repeat, tolerance});
      printf("[%09.3f ms; %04d iters.; %01.9f modularity] copraSeqStatic {labels=%02d, tolerance=%.0e}\n", ak.time, ak.iterations, getModularity(x, ak, M), 2, tolerance);
      // Find COPRA using multiple threads (2 labels).
      auto al = copraOmpStatic<2>(x, init, {repeat, tolerance});
      printf("[%09.3f ms; %04d iters.; %01.9f modularity] copraOmpStatic {labels=%02d, tolerance=%.0e}\n", al.time, al.iterations, getModularity(x, al, M), 2, tolerance);
    }
    {
      // Find COPRA using a single thread (4 labels).
      auto ak = copraSeqStatic<4>(x, init, {repeat, tolerance});
      printf("[%09.3f ms; %04d iters.; %01.9f modularity] copraSeqStatic {labels=%02d, tolerance=%.0e}\n", ak.time, ak.iterations, getModularity(x, ak, M), 4, tolerance);
      // Find COPRA using multiple threads (4 labels).
      auto al = copraOmpStatic<4>(x, init, {repeat, tolerance});
      printf("[%09.3f ms; %04d iters.; %01.9f modularity] copraOmpStatic {labels=%02d, tolerance=%.0e}\n", al.time, al.iterations, getModularity(x, al, M), 4, tolerance);
    }
    {
      // Find COPRA using a single thread (8 labels).
      auto ak = copraSeqStatic<8>(x, init, {repeat, tolerance});
      printf("[%09.3f ms; %04d iters.; %01.9f modularity] copraSeqStatic {labels=%02d, tolerance=%.0e}\n", ak.time, ak.iterations, getModularity(x, ak, M), 8, tolerance);
      // Find COPRA using multiple threads (8 labels).
      auto al = copraOmpStatic<8>(x, init, {repeat, tolerance});
      printf("[%09.3f ms; %04d iters.; %01.9f modularity] copraOmpStatic {labels=%02d, tolerance=%.0e}\n", al.time, al.iterations, getModularity(x, al, M), 8, tolerance);
    }
    {
      // Find COPRA using a single thread (16 labels).
      auto ak = copraSeqStatic<16>(x, init, {repeat, tolerance});
      printf("[%09.3f ms; %04d iters.; %01.9f modularity] copraSeqStatic {labels=%02d, tolerance=%.0e}\n", ak.time, ak.iterations, getModularity(x, ak, M), 16, tolerance);
      // Find COPRA using multiple threads (16 labels).
      auto al = copraOmpStatic<16>(x, init, {repeat, tolerance});
      printf("[%09.3f ms; %04d iters.; %01.9f modularity] copraOmpStatic {labels=%02d, tolerance=%.0e}\n", al.time, al.iterations, getModularity(x, al, M), 16, tolerance);
    }
    {
      // Find COPRA using a single thread (32 labels).
      auto ak = copraSeqStatic<32>(x, init, {repeat, tolerance});
      printf("[%09.3f ms; %04d iters.; %01.9f modularity] copraSeqStatic {labels=%02d, tolerance=%.0e}\n", ak.time, ak.iterations, getModularity(x, ak, M), 32, tolerance);
      // Find COPRA using multiple threads (32 labels).
      auto al = copraOmpStatic<32>(x, init, {repeat, tolerance});
      printf("[%09.3f ms; %04d iters.; %01.9f modularity] copraOmpStatic {labels=%02d, tolerance=%.0e}\n", al.time, al.iterations, getModularity(x, al, M), 32, tolerance);
    }
  }
}


int main(int argc, char **argv) {
  using K = int;
  using V = TYPE;
  char *file = argv[1];
  int repeat = argc>2? stoi(argv[2]) : 5;
  OutDiGraph<K, None, V> x;  // V w = 1;
  printf("Loading graph %s ...\n", file);
  readMtxW<true>(x, file); println(x);
  auto y = symmetricize(x); print(y); printf(" (symmetricize)\n");
  // auto fl = [](auto u) { return true; };
  // selfLoopU(y, w, fl); print(y); printf(" (selfLoopAllVertices)\n");
  omp_set_num_threads(MAX_THREADS);
  printf("OMP_NUM_THREADS=%d\n", MAX_THREADS);
  runExperiment(y, repeat);
  printf("\n");
  return 0;
}
