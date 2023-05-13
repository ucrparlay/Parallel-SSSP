#pragma once
#include "graph.hpp"
#include "hashbag.h"
#include "parlay/internal/get_time.h"
using namespace std;
using namespace parlay;

char const *FILEPATH = nullptr;
constexpr int NUM_SRC = 10001;
constexpr int NUM_ROUND = 1;

// constexpr int BLOCK_SIZE = 1 << 12;
// constexpr size_t MIN_QUEUE = 1 << 14;
constexpr size_t LOCAL_QUEUE_SIZE = 1024;
constexpr size_t DEG_THLD = 20;
constexpr size_t SSSP_SAMPLES = 1000;

enum Algorithm { rho_stepping = 0, delta_stepping, bellman_ford };

class SSSP {
 private:
  const Graph &G;
  Algorithm algo;
  bool sparse;
  int sd_scale;
  EdgeTy delta;
  EdgeTy sample_dist[SSSP_SAMPLES];
  size_t sample_deg[SSSP_SAMPLES];
  size_t param;
  hashbag<NodeId> bag;
  sequence<EdgeTy> dist;
  sequence<NodeId> frontier;
  sequence<bool> in_frontier;
  sequence<bool> in_next_frontier;

  void degree_sampling(size_t sz);
  void sparse_sampling(size_t sz);
  size_t dense_sampling();
  size_t sparse_relax(size_t sz);
  size_t dense_relax();
  void sparse2dense(size_t sz);
  void dense2sparse();
  int pack();

 public:
  SSSP() = delete;
  SSSP(const Graph &_G, Algorithm _algo, size_t _param = 1 << 21)
      : G(_G), algo(_algo), param(_param), bag(G.n) {
    dist = sequence<EdgeTy>::uninitialized(G.n);
    frontier = sequence<NodeId>::uninitialized(G.n);
    in_frontier = sequence<bool>::uninitialized(G.n);
    in_next_frontier = sequence<bool>::uninitialized(G.n);
  }
  sequence<EdgeTy> sssp(int s);
  void set_sd_scale(int x) { sd_scale = x; }
};
