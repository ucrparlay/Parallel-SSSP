#include "sssp.h"

#include <functional>
#include <numeric>

#include "dijkstra.hpp"
#include "parlay/parallel.h"

using namespace std;
using namespace parlay;

void SSSP::degree_sampling(size_t sz) {
  static uint32_t seed = 353442899;
  for (size_t i = 0; i < SSSP_SAMPLES; i++) {
    NodeId u = frontier[hash32(seed) % sz];
    sample_deg[i] = G.offset[u + 1] - G.offset[u];
    seed++;
  }
}

void SSSP::sparse_sampling(size_t sz) {
  static uint32_t seed = 998244353;
  for (size_t i = 0; i < SSSP_SAMPLES; i++) {
    NodeId u = frontier[hash32(seed) % sz];
    sample_dist[i] = dist[u];
    seed++;
  }
  sort(sample_dist, sample_dist + SSSP_SAMPLES);
}

size_t SSSP::dense_sampling() {
  static uint32_t seed = 10086;
  int hits = 0;
  for (size_t i = 0; i < SSSP_SAMPLES; i++) {
    NodeId u = hash32(seed) % G.n;
    if (in_frontier[u]) {
      sample_dist[i] = dist[u];
      hits++;
    } else {
      sample_dist[i] = DIST_MAX;
    }
    seed++;
  }
  sort(sample_dist, sample_dist + SSSP_SAMPLES);
  return 1.0 * hits / SSSP_SAMPLES * G.n;
}

size_t SSSP::sparse_relax(size_t sz) {
  degree_sampling(sz);
  size_t sum_deg = 0;
  for (size_t i = 0; i < SSSP_SAMPLES; i++) {
    sum_deg += sample_deg[i];
  }
  size_t avg_deg = sum_deg / SSSP_SAMPLES;
  bool super_sparse = (avg_deg <= DEG_THLD);
  EdgeTy th;
  if (algo == rho_stepping) {
    sparse_sampling(sz);
    int rate = min(SSSP_SAMPLES - 1, SSSP_SAMPLES * param / sz);
    th = sample_dist[rate];
  } else if (algo == delta_stepping) {
    th = delta;
    delta += param;
  } else {
    th = DIST_MAX;
  }
  // printf("th: %u\n", th);
  parallel_for(0, sz, [&](size_t i) {
    NodeId f = frontier[i];
    if (dist[f] > th) {
      bag.insert(f);
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
            bag.insert(u);
          }
          if (algo == delta_stepping && dist[u] > th) {
            bag.insert(u);
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
                bag.insert(v);
              }
            }
          }
        }
        for (size_t j = front; j < rear; j++) {
          bag.insert(local_queue[j]);
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
                  bag.insert(f);
                }
              }
              for (EdgeId es = _s; es < _e; es++) {
                NodeId v = G.edge[es].v;
                EdgeTy w = G.edge[es].w;
                if (write_min(&dist[v], dist[f] + w,
                              [](EdgeTy w1, EdgeTy w2) { return w1 < w2; })) {
                  bag.insert(v);
                }
              }
            });
      }
    }
  });
  return bag.pack_into(make_slice(frontier));
}

size_t SSSP::dense_relax() {
  int subround = 1;
  while (true) {
    size_t est_size = dense_sampling();
    if (est_size < G.n / sd_scale) {
      break;
    }
    // printf("subround %d\n", subround);
    EdgeTy th;
    if (algo == rho_stepping) {
      int rate;
      if (subround <= 2) {
        rate = min(SSSP_SAMPLES - 1, SSSP_SAMPLES * param / est_size / 10);
      } else {
        rate = min(SSSP_SAMPLES - 1, SSSP_SAMPLES * param / est_size);
      }
      th = sample_dist[rate];
      // printf("th: %u\n", th);
    } else if (algo == delta_stepping) {
      th = delta;
      delta += param;
    } else {
      th = DIST_MAX;
    }
    parallel_for(0, G.n, [&](size_t u) {
      if (dist[u] <= th && in_frontier[u]) {
        in_frontier[u] = false;
        blocked_for(
            G.offset[u], G.offset[u + 1], BLOCK_SIZE,
            [&](size_t, size_t _s, size_t _e) {
              if (G.symmetrized) {
                EdgeTy temp_dist = dist[u];
                for (size_t es = _s; es < _e; es++) {
                  NodeId v = G.edge[es].v;
                  EdgeTy w = G.edge[es].w;
                  temp_dist = min(temp_dist, dist[v] + w);
                }
                if (write_min(&dist[u], temp_dist,
                              [](EdgeTy w1, EdgeTy w2) { return w1 < w2; })) {
                  if (!in_next_frontier[u]) {
                    in_next_frontier[u] = true;
                  }
                }
              }
              for (size_t es = _s; es < _e; es++) {
                NodeId v = G.edge[es].v;
                EdgeTy w = G.edge[es].w;
                if (write_min(&dist[v], dist[u] + w,
                              [](EdgeTy w1, EdgeTy w2) { return w1 < w2; })) {
                  if (!in_next_frontier[v]) {
                    in_next_frontier[v] = true;
                  }
                }
              }
            });
      }
    });
    subround++;
    swap(in_frontier, in_next_frontier);
  }
  return count(in_frontier, true);
}

void SSSP::sparse2dense(size_t sz) {
  parallel_for(0, sz, [&](size_t i) {
    NodeId u = frontier[i];
    in_frontier[u] = true;
  });
}

void SSSP::dense2sparse() {
  auto identity = delayed_seq<NodeId>(G.n, [&](size_t i) { return i; });
  pack_into_uninitialized(identity, in_frontier, frontier);
}

int SSSP::pack() {
  size_t next_size = 0;
  bool next_sparse;
  if (sparse) {
    next_size = bag.pack_into(make_slice(frontier));
    next_sparse = (next_size < G.n / sd_scale);
    if (!next_sparse) {
    }
  } else {  // dense
    next_size = count(in_frontier, true);
    next_sparse = (next_size < G.n / sd_scale);
    if (next_sparse) {
    }
  }
  // printf("next size: %zu\n", next_size);
  return next_size;
}

sequence<EdgeTy> SSSP::sssp(int s) {
  if (!G.weighted) {
    fprintf(stderr, "Error: Input graph is unweighted\n");
    exit(EXIT_FAILURE);
  }
  if (algo == delta_stepping) {
    delta = param;
  }

  parallel_for(0, G.n, [&](size_t i) {
    dist[i] = numeric_limits<EdgeTy>::max() / 2;
    in_frontier[i] = in_next_frontier[i] = false;
  });
  size_t size = 1;
  frontier[0] = s;
  dist[s] = 0;
  sparse = true;

  while (size) {
    // printf("size: %zu, threshold: %zu, %s\n", size, G.n / sd_scale,
    // sparse ? "sparse" : "dense");
    if (sparse) {
      size = sparse_relax(size);
    } else {
      size = dense_relax();
    }
    bool next_sparse = (size < G.n / sd_scale) ? true : false;
    if (sparse && !next_sparse) {
      sparse2dense(size);
    } else if (!sparse && next_sparse) {
      dense2sparse();
    }
    sparse = next_sparse;
  }
  return dist;
}

int main(int argc, char *argv[]) {
  if (argc == 1) {
    fprintf(
        stderr,
        "Usage: %s [-i input_file] [-p parameter] [-w] [-s] [-v] [-a "
        "algorithm]\n"
        "Options:\n"
        "\t-i,\tinput file path\n"
        "\t-p,\tparameter(e.g. delta, rho)\n"
        "\t-w,\tweighted input graph\n"
        "\t-s,\tsymmetrized input graph\n"
        "\t-v,\tverify result\n"
        "\t-a,\talgorithm: [rho-stepping] [delta-stepping] [bellman-ford]\n",
        argv[0]);
    exit(EXIT_FAILURE);
  }
  char c;
  bool weighted = false;
  bool symmetrized = false;
  bool verify = false;
  size_t param = 1 << 21;
  Algorithm algo = rho_stepping;
  while ((c = getopt(argc, argv, "i:p:a:wsv")) != -1) {
    switch (c) {
      case 'i':
        FILEPATH = optarg;
        break;
      case 'p':
        param = atol(optarg);
        break;
      case 'a':
        if (!strcmp(optarg, "rho-stepping")) {
          algo = rho_stepping;
        } else if (!strcmp(optarg, "delta-stepping")) {
          algo = delta_stepping;
        } else if (!strcmp(optarg, "bellman-ford")) {
          algo = bellman_ford;
        } else {
          fprintf(stderr, "Error: Unknown algorithm %s\n", optarg);
          exit(EXIT_FAILURE);
        }
        break;
      case 'w':
        weighted = true;
        break;
      case 's':
        symmetrized = true;
        break;
      case 'v':
        verify = true;
        break;
      default:
        fprintf(stderr, "Error: Unknown option %c\n", optopt);
        exit(EXIT_FAILURE);
    }
  }
  Graph G(weighted, symmetrized);

  printf("Reading graph...\n");
  G.read_graph(FILEPATH);
  // G.generate_random_graph();
  if (!weighted) {
    printf("Generating edge weights...\n");
    G.generate_weight();
  }

  SSSP solver(G, algo, param);
  int sd_scale = G.m / G.n;
  solver.set_sd_scale(sd_scale);
  fprintf(stdout,
          "Running on %s: |V|=%zu, |E|=%zu, param=%zu, num_src=%d, "
          "num_round=%d\n",
          FILEPATH, G.n, G.m, param, NUM_SRC, NUM_ROUND);

  for (int v = 0; v < NUM_SRC; v++) {
    int s = hash32(v) % G.n;
    printf("source %d: %-10d\n", v, s);
    double total_time = 0;
    for (int i = 0; i <= NUM_ROUND; i++) {
      internal::timer t;
      solver.sssp(s);
      t.stop();
      if (i == 0) {
        printf("Warmup Round: %f\n", t.total_time());
      } else {
        printf("Round %d: %f\n", i, t.total_time());
        total_time += t.total_time();
      }
    }
    double average_time = total_time / NUM_ROUND;
    printf("Average time: %f\n", average_time);

    if (verify) {
      auto dist = solver.sssp(s);
      verifier(s, G, dist);
    }
    printf("\n");
  }
  return 0;
}
