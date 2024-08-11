#include <fstream>
#include <functional>
#include <numeric>

#include "graph.h"
using namespace std;
using namespace parlay;

int main(int argc, char *argv[]) {
  if (argc == 1) {
    fprintf(stderr,
            "Usage: %s [-i input_file] [-o output_file] [-w] [-s]\n"
            "Options:\n"
            "\t-i,\tinput file path\n"
            "\t-o,\toutput file path\n"
            "\t-w,\tweighted input graph\n"
            "\t-s,\tsymmetrized input graph\n",
            argv[0]);
    exit(EXIT_FAILURE);
  }
  char c;
  bool weighted = false;
  bool symmetrized = false;
  char const *input_path = nullptr;
  char const *output_path = nullptr;
  while ((c = getopt(argc, argv, "i:o:ws")) != -1) {
    switch (c) {
      case 'i':
        input_path = optarg;
        break;
      case 'o':
        output_path = optarg;
        break;
      case 'w':
        weighted = true;
        break;
      case 's':
        symmetrized = true;
        break;
      default:
        fprintf(stderr, "Error: Unknown option %c\n", optopt);
        exit(EXIT_FAILURE);
    }
  }
  Graph G(weighted, symmetrized);

  printf("Reading graph...\n");
  G.read_graph(input_path);
  printf("n: %zu, m: %zu\n", G.n, G.m);
  auto G2 = G.symmetrize_weighted_graph();
  G2.write_pbbs_format(output_path);
  return 0;
}
