#pragma once
#include "graph.hpp"
#include "hashbag.h"
#include "parlay/internal/get_time.h"
using namespace std;
using namespace parlay;

char const *FILEPATH = nullptr;
constexpr int NUM_SRC = 100;
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
  void add_to_bag(NodeId v);
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

class Rho_Stepping {
  size_t rho;
  size_t seed;

  const bool &sparse;
  const size_t &size;
  const size_t &n;
  const sequence<EdgeTy> &dist;
  const sequence<bool> &in_frontier;
  const sequence<NodeId> &frontier;

  Rho_Stepping(size_t _rho, const bool &_sparse, const size_t &_size,
               const size_t &_n, const sequence<EdgeTy> &_dist,
               const sequence<bool> &_in_frontier,
               const sequence<NodeId> &_frontier)
      : rho(_rho),
        sparse(_sparse),
        size(_size),
        n(_n),
        dist(_dist),
        in_frontier(_in_frontier),
        frontier(_frontier) {}
  void init() { seed = 0; }
  EdgeTy get_threshold() {
    if (size <= rho) {
      return DIST_MAX;
    }
    EdgeTy sample_dist[SSSP_SAMPLES + 1];
    for (size_t i = 0; i <= SSSP_SAMPLES; i++) {
      if (sparse) {
        NodeId v = frontier[hash32(seed + i) % size];
        sample_dist[i] = dist[v];
      } else {
        NodeId v = hash32(seed + i) % n;
        if (in_frontier[v]) {
          sample_dist[i] = dist[v];
        } else {
          sample_dist[i] = DIST_MAX;
        }
      }
    }
    size_t id = 1.0 * rho / size * SSSP_SAMPLES;
    return sample_dist[id];
  }
};

class Delta_Stepping {
  EdgeTy delta;
  EdgeTy thres;

  Delta_Stepping(EdgeTy _delta) : delta(_delta) {}
  void init() { thres = 0; }
  EdgeTy get_threshold() {
    thres += delta;
    return thres;
  }
};

class Bellman_Ford {
  Bellman_Ford() {}
  void init() {}
  EdgeTy get_threshold() { return DIST_MAX; }
};
