#include "sssp.h"

#include <fstream>
#include <functional>
#include <numeric>

#include "dijkstra.hpp"
#include "parlay/parallel.h"

using namespace std;
using namespace parlay;

template <class Algo>
void run(Algo &algo, const Graph &G, bool verify) {
  for (int v = 0; v < NUM_SRC; v++) {
    NodeId s = hash32(v) % G.n;
    printf("source %d: %-10d\n", v, s);
    double total_time = 0;
    for (int i = 0; i <= NUM_ROUND; i++) {
      internal::timer t;
      algo.sssp(s);
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

    ofstream ofs("sssp.tsv", ios_base::app);
    ofs << average_time << '\t';
    ofs.close();

    if (verify) {
      printf("Running verifier...\n");
      internal::timer t;
      auto dist = algo.sssp(s);
      t.stop();
      printf("Our running time: %f\n", t.total_time());
      verifier(s, G, dist);
    }
    printf("\n");
  }
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
  size_t param = ULLONG_MAX;
  int algo = rho_stepping;
  char const *FILEPATH = nullptr;
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
  if (!weighted) {
    printf("Generating edge weights...\n");
    G.generate_weight();
  }

  fprintf(stdout,
          "Running on %s: |V|=%zu, |E|=%zu, param=%zu, num_src=%d, "
          "num_round=%d\n",
          FILEPATH, G.n, G.m, param, NUM_SRC, NUM_ROUND);

  int sd_scale = G.m / G.n;
  if (algo == rho_stepping) {
    Rho_Stepping solver(G);
    solver.set_sd_scale(sd_scale);
    run(solver, G, verify);
  } else if (algo == delta_stepping) {
    Delta_Stepping solver(G);
    solver.set_sd_scale(sd_scale);
    run(solver, G, verify);
  } else if (algo == bellman_ford) {
    Bellman_Ford solver(G);
    solver.set_sd_scale(sd_scale);
    run(solver, G, verify);
  }
  return 0;
}
