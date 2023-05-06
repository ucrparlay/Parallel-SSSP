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
    NodeId u = que[cur][hash32(seed) % sz];
    sample_deg[i] = G.offset[u + 1] - G.offset[u];
    seed++;
  }
}

void SSSP::sparse_sampling(size_t sz) {
  static uint32_t seed = 998244353;
  for (size_t i = 0; i < SSSP_SAMPLES; i++) {
    NodeId u = que[cur][hash32(seed) % sz];
    sample_dist[i] = info[u].dist;
    seed++;
  }
  sort(sample_dist, sample_dist + SSSP_SAMPLES);
}

size_t SSSP::dense_sampling() {
  static uint32_t seed = 10086;
  int num_sample = 0;
  for (size_t i = 0; i < SSSP_SAMPLES;) {
    num_sample++;
    NodeId u = hash32(seed) % G.n;
    if (info[u].fl & in_que) {
      sample_dist[i] = info[u].dist;
      i++;
    }
    seed++;
  }
  sort(sample_dist, sample_dist + SSSP_SAMPLES);
  return 1.0 * SSSP_SAMPLES / num_sample * G.n;
}

void SSSP::relax(size_t sz) {
  if (sparse) {
    size_t qsize[doubling], cnt[doubling];
    qsize[0] = cnt[0] = cnt[1] = 0;
    qsize[1] = MIN_QUEUE;
    int db_len = 2;
    for (size_t s = MIN_QUEUE * 2; s <= max_queue; s *= 2) {
      qsize[db_len] = s;
      cnt[db_len] = 0;
      db_len++;
    }
    int pt = 1;
    auto add = [&](NodeId u) {
      if ((info[u].fl & to_add) ||
          !atomic_compare_and_swap(&info[u].fl, info[u].fl,
                                   info[u].fl | to_add)) {
        return;
      }
      int t_pt = pt;
      size_t pos =
          hash32(u) % (qsize[t_pt] - qsize[t_pt - 1]) + qsize[t_pt - 1];
      while (que[nxt][pos] != UINT_MAX ||
             !atomic_compare_and_swap(&que[nxt][pos], UINT_MAX, u)) {
        pos++;
        if (pos == qsize[t_pt]) {
          pos = qsize[t_pt - 1];
        }
      }
      // The queue should be half occupied when we have EXP_SAMPLES
      size_t len = qsize[t_pt] - qsize[t_pt - 1];
      size_t rate = len / (2 * EXP_SAMPLES);
      if (pos % rate == 0) {
        int ret = fetch_and_add(&cnt[t_pt], 1);
        if (ret + 1 == EXP_SAMPLES && t_pt + 1 < db_len) {
          atomic_compare_and_swap(&pt, t_pt, t_pt + 1);
        }
      }
    };

    auto relax_neighbors = [&](NodeId u, EdgeId _s, EdgeId _e) {
      _s += G.offset[u];
      _e += G.offset[u];
      if (G.symmetrized) {
        EdgeTy temp_dis = info[u].dist;
        for (EdgeId es = _s; es < _e; es++) {
          NodeId v = G.edge[es].v;
          EdgeTy w = G.edge[es].w;
          temp_dis = min(temp_dis, info[v].dist + w);
        }
        if (write_min(&info[u].dist, temp_dis,
                      [](EdgeTy w1, EdgeTy w2) { return w1 < w2; })) {
          add(u);
        }
      }
      for (EdgeId es = _s; es < _e; es++) {
        NodeId v = G.edge[es].v;
        EdgeTy w = G.edge[es].w;
        if (write_min(&info[v].dist, info[u].dist + w,
                      [](EdgeTy w1, EdgeTy w2) { return w1 < w2; })) {
          add(v);
        }
      }
    };
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
      th = UINT_MAX;
    }
    parallel_for(0, sz, [&](size_t i) {
      NodeId f = que[cur][i];
      que[cur][i] = UINT_MAX;
      if (info[f].dist > th) {
        add(f);
      } else {
        size_t _n = G.offset[f + 1] - G.offset[f];
        if (super_sparse && _n < BLOCK_SIZE) {
          sequence<NodeId> q(BLOCK_SIZE);
          int front = 0, rear = 0;
          q[rear++] = f;
          while (front < rear && rear < BLOCK_SIZE) {
            NodeId u = q[front];
            size_t deg = G.offset[u + 1] - G.offset[u];
            if (deg >= BLOCK_SIZE) {
              break;
            }
            front++;
            if (algo == delta_stepping && info[u].dist > th) {
              add(u);
              continue;
            }
            if (G.symmetrized) {
              EdgeTy temp_dis = info[u].dist;
              for (EdgeId es = G.offset[u]; es < G.offset[u + 1]; es++) {
                NodeId v = G.edge[es].v;
                EdgeTy w = G.edge[es].w;
                temp_dis = min(temp_dis, info[v].dist + w);
              }
              write_min(&info[u].dist, temp_dis,
                        [](EdgeTy w1, EdgeTy w2) { return w1 < w2; });
            }
            for (EdgeId es = G.offset[u]; es < G.offset[u + 1]; es++) {
              NodeId v = G.edge[es].v;
              EdgeTy w = G.edge[es].w;
              if (write_min(&info[v].dist, info[u].dist + w,
                            [](EdgeTy w1, EdgeTy w2) { return w1 < w2; })) {
                if (rear < BLOCK_SIZE) {
                  q[rear++] = v;
                } else {
                  add(v);
                }
              }
            }
          }
          while (front < rear) {
            NodeId u = q[front++];
            add(u);
          }
        } else {
          blocked_for(0, _n, BLOCK_SIZE, [&](size_t, size_t _s, size_t _e) {
            relax_neighbors(f, _s, _e);
          });
        }
      }
    });
    que_size = qsize[pt];
  } else {  // dense
    auto relax_neighbors = [&](NodeId u, EdgeId _s, EdgeId _e) {
      _s += G.offset[u];
      _e += G.offset[u];
      if (G.symmetrized) {
        EdgeTy temp_dis = info[u].dist;
        for (size_t es = _s; es < _e; es++) {
          NodeId v = G.edge[es].v;
          EdgeTy w = G.edge[es].w;
          temp_dis = min(temp_dis, info[v].dist + w);
        }
        if (write_min(&info[u].dist, temp_dis,
                      [](EdgeTy w1, EdgeTy w2) { return w1 < w2; })) {
          if (!(info[u].fl & in_que)) {
            info[u].fl |= in_que;
          }
        }
      }
      for (size_t es = _s; es < _e; es++) {
        NodeId v = G.edge[es].v;
        EdgeTy w = G.edge[es].w;
        if (write_min(&info[v].dist, info[u].dist + w,
                      [](EdgeTy w1, EdgeTy w2) { return w1 < w2; })) {
          if (!(info[v].fl & in_que)) {
            info[v].fl |= in_que;
          }
        }
      }
    };

    int subround = 1;
    while (true) {
      size_t est_size = dense_sampling();
      if (est_size < G.n / sd_scale) {
        break;
      }
      EdgeTy th;
      if (algo == rho_stepping) {
        int rate;
        if (subround <= 2) {
          rate = min(SSSP_SAMPLES - 1, SSSP_SAMPLES * param / est_size / 10);
        } else {
          rate = min(SSSP_SAMPLES - 1, SSSP_SAMPLES * param / est_size);
        }
        th = sample_dist[rate];
      } else if (algo == delta_stepping) {
        th = delta;
        delta += param;
      } else {
        th = UINT_MAX;
      }
      parallel_for(0, G.n, [&](size_t u) {
        if (info[u].dist <= th && (info[u].fl & in_que)) {
          info[u].fl &= ~in_que;
          size_t _n = G.offset[u + 1] - G.offset[u];
          blocked_for(0, _n, BLOCK_SIZE,
                      [&]([[maybe_unused]] size_t j, size_t _s, size_t _e) {
                        relax_neighbors(u, _s, _e);
                      });
        }
      });
      subround++;
    }
  }
}

int SSSP::pack() {
  size_t nxt_sz = 0;
  bool next_sparse;
  if (sparse) {
    parallel_for(0, que_size,
                 [&](size_t i) { que_num[i] = (que[nxt][i] != UINT_MAX); });
    nxt_sz = scan_inplace(que_num.cut(0, que_size),
                          monoid([](NodeId a, NodeId b) { return a + b; }, 0));
    next_sparse = (nxt_sz < G.n / sd_scale);
    if (next_sparse) {
      sequence<NodeId> tmp(nxt_sz);
      parallel_for(0, que_size, [&](size_t i) {
        if (que[nxt][i] != UINT_MAX) {
          NodeId u = que[nxt][i];
          que[nxt][i] = UINT_MAX;
          info[u].fl ^= to_add;
          tmp[que_num[i]] = u;
        }
      });
      parallel_for(0, nxt_sz, [&](size_t i) { que[nxt][i] = tmp[i]; });
    } else {
      parallel_for(0, que_size, [&](size_t i) {
        if (que[nxt][i] != UINT_MAX) {
          NodeId u = que[nxt][i];
          que[nxt][i] = UINT_MAX;
          info[u].fl |= in_que;
          info[u].fl ^= to_add;
        }
      });
    }
  } else {  // dense
    auto que_num0 = dseq(
        G.n, [&](size_t i) -> NodeId { return (info[i].fl & in_que) ? 1 : 0; });
    nxt_sz = internal::scan_(
        que_num0, make_slice(que_num),
        monoid([](NodeId a, NodeId b) { return a + b; }, 0), no_flag);
    next_sparse = (nxt_sz < G.n / sd_scale);
    if (next_sparse) {
      parallel_for(0, G.n, [&](size_t i) {
        if (info[i].fl & in_que) {
          info[i].fl &= ~in_que;
          que[nxt][que_num[i]] = i;
        }
      });
    }
  }
  swap(cur, nxt);
  return nxt_sz;
}

void SSSP::reset_timer() { t_all.reset(); }

void SSSP::sssp(int s, EdgeTy *_dist) {
  if (!G.weighted) {
    fprintf(stderr, "Error: Input graph is unweighted\n");
    exit(EXIT_FAILURE);
  }
  t_all.start();
  cur = 0, nxt = 1;
  if (algo == delta_stepping) {
    delta = param;
  }
  parallel_for(0, que[0].size(), [&](size_t i) { que[0][i] = UINT_MAX; });
  parallel_for(0, que[1].size(), [&](size_t i) { que[1][i] = UINT_MAX; });
  parallel_for(0, info.size(),
               [&](size_t i) { info[i] = Information(INT_MAX / 2, 0); });

  size_t sz = 1;
  que[cur][0] = s;
  info[s].dist = 0;
  sparse = true;

  while (sz) {
    relax(sz);
    sz = pack();
    if (sz >= G.n / sd_scale) {
      sparse = false;
    } else {
      sparse = true;
    };
  }
  t_all.stop();
  parallel_for(0, G.n, [&](size_t i) { _dist[i] = info[i].dist; });
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

  printf("Info: Reading graph\n");
  G.read_graph(FILEPATH);
  if (!weighted) {
    printf("Info: Generating edge weights\n");
    G.generate_weight();
  }

  SSSP solver(G, algo, param);
  int sd_scale = G.m / G.n;
  solver.set_sd_scale(sd_scale);
  fprintf(stdout,
          "Running on %s: |V|=%zu, |E|=%zu, param=%zu, num_src=%d, "
          "num_round=%d\n",
          FILEPATH, G.n, G.m, param, NUM_SRC, NUM_ROUND);
  EdgeTy *dijkstra_dist = new EdgeTy[G.n];
  EdgeTy *my_dist = new EdgeTy[G.n];

  for (int v = 0; v < NUM_SRC; v++) {
    int s = hash32(v) % G.n;
    printf("source: %-10d\n", s);
    vector<double> sssp_time;
    // first time warmup
    solver.reset_timer();
    solver.sssp(s, my_dist);
    printf("warmup round (not counted): %f\n", solver.t_all.total_time());

    for (int i = 0; i < NUM_ROUND; i++) {
      solver.reset_timer();
      solver.sssp(s, my_dist);
      sssp_time.push_back(solver.t_all.total_time());

      printf("round %d: %f\n", i + 1, solver.t_all.total_time());
    }
    sort(begin(sssp_time), end(sssp_time));
    printf("median running time: %f\n", sssp_time[(sssp_time.size() - 1) / 2]);
    printf("average running time: %f\n",
           accumulate(begin(sssp_time), end(sssp_time), 0.0) / NUM_ROUND);

    if (verify) {
      printf("Info: Running verifier\n");
      verifier(s, G, my_dist);
    }
  }
  delete[] dijkstra_dist;
  delete[] my_dist;

  return 0;
}
