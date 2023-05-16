#pragma once
#include "graph.hpp"
#include "hashbag.h"
#include "parlay/internal/get_time.h"
using namespace std;
using namespace parlay;

constexpr int NUM_SRC = 1000;
constexpr int NUM_ROUND = 5;

constexpr size_t LOCAL_QUEUE_SIZE = 1024;
constexpr size_t DEG_THLD = 20;
constexpr size_t SSSP_SAMPLES = 1000;

enum Algorithm { rho_stepping = 0, delta_stepping, bellman_ford };

class SSSP {
 protected:
  const Graph &G;
  Algorithm algo;
  bool sparse;
  int sd_scale;
  size_t frontier_size;
  hashbag<NodeId> bag;
  sequence<EdgeTy> dist;
  sequence<NodeId> frontier;
  sequence<bool> in_frontier;
  sequence<bool> in_next_frontier;

  void add_to_bag(NodeId v) {
    // if (!in_next_frontier[v]) {
    if (compare_and_swap(&in_next_frontier[v], false, true)) {
      // in_next_frontier[v] = true;
      bag.insert(v);
    }
  }

  size_t estimate_size() {
    static uint32_t seed = 10086;
    size_t hits = 0;
    for (size_t i = 0; i < SSSP_SAMPLES; i++) {
      NodeId u = hash32(seed) % G.n;
      if (in_frontier[u]) {
        hits++;
      }
      seed++;
    }
    return 1.0 * hits / SSSP_SAMPLES * G.n;
  }

  size_t sparse_relax() {
    // static uint32_t seed = 353442899;
    // size_t sum_deg = 0;
    // for (size_t i = 0; i < SSSP_SAMPLES; i++) {
    // NodeId u = frontier[hash32(seed) % size];
    // sum_deg += G.offset[u + 1] - G.offset[u];
    // seed++;
    //}
    // size_t avg_deg = sum_deg / SSSP_SAMPLES;
    // bool super_sparse = (avg_deg <= DEG_THLD);
    bool super_sparse = true;

    printf(">>>frontier_size: %zu\n", frontier_size);
    EdgeTy th = get_threshold();
    printf("th: %u\n", th);

    parallel_for(0, frontier_size, [&](size_t i) {
      NodeId f = frontier[i];
      assert(in_frontier[f] == true);
      in_frontier[f] = false;
      if (dist[f] > th) {
        add_to_bag(f);
        // bag.insert(f);
      } else {
        size_t _n = G.offset[f + 1] - G.offset[f];
        if (super_sparse && _n < LOCAL_QUEUE_SIZE) {
          NodeId local_queue[LOCAL_QUEUE_SIZE];
          size_t front = 0, rear = 0;
          local_queue[rear++] = f;
          while (front < rear && rear < LOCAL_QUEUE_SIZE) {
            NodeId u = local_queue[front++];
            size_t deg = G.offset[u + 1] - G.offset[u];
            if (deg >= LOCAL_QUEUE_SIZE) {
              add_to_bag(u);
              //  bag.insert(u);
              continue;
            }
            if (dist[u] > th) {
              add_to_bag(u);
              // bag.insert(u);
              continue;
            }
            if (G.symmetrized) {
              EdgeTy temp_dis = dist[u];
              for (EdgeId es = G.offset[u]; es < G.offset[u + 1]; es++) {
                NodeId v = G.edge[es].v;
                EdgeTy w = G.edge[es].w;
                temp_dis = min(temp_dis, dist[v] + w);
              }
              write_min(&dist[u], temp_dis,
                        [](EdgeTy w1, EdgeTy w2) { return w1 < w2; });
            }
            for (EdgeId es = G.offset[u]; es < G.offset[u + 1]; es++) {
              NodeId v = G.edge[es].v;
              EdgeTy w = G.edge[es].w;
              if (write_min(&dist[v], dist[u] + w,
                            [](EdgeTy w1, EdgeTy w2) { return w1 < w2; })) {
                if (rear < LOCAL_QUEUE_SIZE) {
                  local_queue[rear++] = v;
                } else {
                  add_to_bag(v);
                  // bag.insert(v);
                }
              }
            }
          }
          for (size_t j = front; j < rear; j++) {
            add_to_bag(local_queue[j]);
            // bag.insert(local_queue[j]);
          }
        } else {
          blocked_for(
              G.offset[f], G.offset[f + 1], BLOCK_SIZE,
              [&](size_t, size_t _s, size_t _e) {
                if (G.symmetrized) {
                  EdgeTy temp_dist = dist[f];
                  for (EdgeId es = _s; es < _e; es++) {
                    NodeId v = G.edge[es].v;
                    EdgeTy w = G.edge[es].w;
                    temp_dist = min(temp_dist, dist[v] + w);
                  }
                  if (write_min(&dist[f], temp_dist,
                                [](EdgeTy w1, EdgeTy w2) { return w1 < w2; })) {
                    add_to_bag(f);
                    // bag.insert(f);
                  }
                }
                for (EdgeId es = _s; es < _e; es++) {
                  NodeId v = G.edge[es].v;
                  EdgeTy w = G.edge[es].w;
                  if (write_min(&dist[v], dist[f] + w,
                                [](EdgeTy w1, EdgeTy w2) { return w1 < w2; })) {
                    add_to_bag(v);
                    // bag.insert(v);
                  }
                }
              });
        }
      }
    });
    printf("done\n");
    swap(in_frontier, in_next_frontier);
    return bag.pack_into(make_slice(frontier));
  }

  size_t dense_relax() {
    int subround = 1;
    while (true) {
      size_t est_size = estimate_size();
      if (est_size == 0 || est_size < G.n / sd_scale) {
        break;
      }
      EdgeTy th = get_threshold();
      parallel_for(0, G.n, [&](NodeId u) {
        if (in_frontier[u]) {
          in_frontier[u] = false;
          if (dist[u] > th) {
            in_next_frontier[u] = true;
          } else {
            blocked_for(G.offset[u], G.offset[u + 1], BLOCK_SIZE,
                        [&](size_t, size_t _s, size_t _e) {
                          if (G.symmetrized) {
                            EdgeTy temp_dist = dist[u];
                            for (size_t es = _s; es < _e; es++) {
                              NodeId v = G.edge[es].v;
                              EdgeTy w = G.edge[es].w;
                              temp_dist = min(temp_dist, dist[v] + w);
                            }
                            if (write_min(&dist[u], temp_dist,
                                          [](EdgeTy w1, EdgeTy w2) {
                                            return w1 < w2;
                                          })) {
                              if (!in_next_frontier[u]) {
                                in_next_frontier[u] = true;
                              }
                            }
                          }
                          for (size_t es = _s; es < _e; es++) {
                            NodeId v = G.edge[es].v;
                            EdgeTy w = G.edge[es].w;
                            if (write_min(&dist[v], dist[u] + w,
                                          [](EdgeTy w1, EdgeTy w2) {
                                            return w1 < w2;
                                          })) {
                              if (!in_next_frontier[v]) {
                                in_next_frontier[v] = true;
                              }
                            }
                          }
                        });
          }
        }
      });
      subround++;
      swap(in_frontier, in_next_frontier);
    }
    return count(in_frontier, true);
  }

  void sparse2dense() {
    parallel_for(0, frontier_size, [&](size_t i) {
      NodeId u = frontier[i];
      assert(in_frontier[u] == true);
      // in_frontier[u] = true;
    });
  }

  void dense2sparse() {
    auto identity = delayed_seq<NodeId>(G.n, [&](NodeId i) { return i; });
    pack_into_uninitialized(identity, in_frontier, frontier);
  }

  function<void()> init;
  function<EdgeTy()> get_threshold;

 public:
  SSSP() = delete;
  SSSP(const Graph &_G) : G(_G), bag(G.n) {
    dist = sequence<EdgeTy>::uninitialized(G.n);
    frontier = sequence<NodeId>::uninitialized(G.n);
    in_frontier = sequence<bool>::uninitialized(G.n);
    in_next_frontier = sequence<bool>::uninitialized(G.n);
  }
  sequence<EdgeTy> sssp(NodeId s) {
    if (!G.weighted) {
      fprintf(stderr, "Error: Input graph is unweighted\n");
      exit(EXIT_FAILURE);
    }

    init();
    parallel_for(0, G.n, [&](NodeId i) {
      dist[i] = numeric_limits<EdgeTy>::max() / 2;
      in_frontier[i] = in_next_frontier[i] = false;
    });
    frontier_size = 1;
    dist[s] = 0;
    frontier[0] = s;
    in_frontier[s] = true;
    sparse = true;

    while (frontier_size) {
      printf("frontier_size: %zu, threshold: %zu, %s\n", frontier_size,
             G.n / sd_scale, sparse ? "sparse" : "dense");
      if (sparse) {
        frontier_size = sparse_relax();
      } else {
        frontier_size = dense_relax();
      }
      bool next_sparse = (frontier_size < G.n / sd_scale) ? true : false;
      if (sparse && !next_sparse) {
        sparse2dense();
      } else if (!sparse && next_sparse) {
        dense2sparse();
      }
      sparse = next_sparse;
    }
    return dist;
  }

  void set_sd_scale(int x) { sd_scale = x; }
};

class Rho_Stepping : public SSSP {
  size_t rho;
  uint32_t seed;

 public:
  Rho_Stepping(const Graph &_G, size_t _rho) : SSSP(_G), rho(_rho) {
    seed = 0;
    init = []() {};
    get_threshold = [&]() {
      printf("frontier_size: %zu, rho: %zu\n", frontier_size, rho);
      printf("G.n: %zu, G.m: %zu\n", G.n, G.m);
      printf("seed: %u\n", seed);
      frontier_size = 10;
      if (frontier_size <= rho) {
        return DIST_MAX;
      }
      EdgeTy sample_dist[SSSP_SAMPLES + 1];
      for (size_t i = 0; i <= SSSP_SAMPLES; i++) {
        if (sparse) {
          NodeId v = frontier[hash32(seed + i) % frontier_size];
          sample_dist[i] = dist[v];
        } else {
          NodeId v = hash32(seed + i) % G.n;
          if (in_frontier[v]) {
            sample_dist[i] = dist[v];
          } else {
            sample_dist[i] = DIST_MAX;
          }
        }
      }
      seed += SSSP_SAMPLES + 1;
      size_t id = 1.0 * rho / frontier_size * SSSP_SAMPLES;
      return sample_dist[id];
    };
  }
};

class Delta_Stepping : public SSSP {
  EdgeTy delta;
  EdgeTy thres;

 public:
  Delta_Stepping(const Graph &_G, EdgeTy _delta) : SSSP(_G), delta(_delta) {
    init = [&]() { thres = 0; };
    get_threshold = [&]() {
      thres += delta;
      return thres;
    };
  }
};

class Bellman_Ford : public SSSP {
 public:
  Bellman_Ford(const Graph &_G) : SSSP(_G) {
    init = []() {};
    get_threshold = []() { return DIST_MAX; };
  }
};
