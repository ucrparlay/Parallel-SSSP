#pragma once
#include "graph.hpp"
#include "pbbslib/get_time.h"
#include "pbbslib/parallel.h"
#include "pbbslib/sequence.h"
#include "pbbslib/utilities.h"
using namespace std;
using namespace pbbs;

char const *FILEPATH = nullptr;
constexpr int NUM_SRC = 1000;
constexpr int NUM_ROUND = 10;
constexpr uint32_t in_que = 1;
constexpr uint32_t to_add = 2;

constexpr int BLOCK_SIZE = 1 << 12;
constexpr size_t MIN_QUEUE = 1 << 14;
constexpr size_t DEG_THLD = 20;
constexpr size_t SSSP_SAMPLES = 1000;
constexpr size_t EXP_SAMPLES = 100;

enum Algorithm { rho_stepping = 0, delta_stepping, bellman_ford };

struct Information {
  EdgeTy dist;
  uint32_t fl;
  Information() : dist(INT_MAX), fl(0) {}
  Information(EdgeTy _dist, uint32_t _fl) : dist(_dist), fl(_fl) {}
};

class SSSP {
 private:
  const Graph &G;
  Algorithm algo;
  bool sparse;
  int cur, nxt;
  int doubling;
  int sd_scale = 1;
  EdgeTy delta;
  EdgeTy sample_dist[SSSP_SAMPLES];
  size_t sample_deg[SSSP_SAMPLES];
  size_t que_size;
  size_t param;
  size_t max_queue;
  sequence<Information> info;
  sequence<NodeId> que[2];
  sequence<NodeId> que_num;

  void degree_sampling(size_t sz);
  void sparse_sampling(size_t sz);
  size_t dense_sampling();
  void relax(size_t sz);
  int pack();

 public:
  SSSP() = delete;
  SSSP(const Graph &_G, Algorithm _algo, size_t _param = 1 << 21)
      : G(_G), algo(_algo), param(_param) {
    max_queue = 1ULL << static_cast<int>(ceil(log2(G.n)));
    doubling = ceil(log2(max_queue / MIN_QUEUE)) + 2;
    info = sequence<Information>(G.n);
    que[0] = que[1] = sequence<NodeId>(max_queue);
    que_num = sequence<NodeId>(max_queue);
  }
  void sssp(int s, EdgeTy *dist);
  void reset_timer();
  void set_sd_scale(int x) { sd_scale = x; }
  timer t_all;
};
