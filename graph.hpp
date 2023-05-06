#pragma once
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#include <cassert>
#include <cinttypes>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "parlay/delayed_sequence.h"
#include "parlay/primitives.h"
#include "parlay/sequence.h"
#include "parlay/utilities.h"
#include "utils.h"
using namespace std;
using namespace parlay;

using NodeId = uint32_t;
using EdgeId = uint64_t;
using EdgeTy = uint32_t;
constexpr int LOG2_WEIGHT = 18;
constexpr int WEIGHT = 1 << LOG2_WEIGHT;

struct Edge {
  NodeId v;
  EdgeTy w;
  Edge() : v(0), w(0){};
  Edge(NodeId _v, EdgeTy _w) : v(_v), w(_w) {}
  bool operator<(const Edge& rhs) const {
    if (v != rhs.v) {
      return v < rhs.v;
    }
    return w < rhs.w;
  }
  bool operator!=(const Edge& rhs) const { return v != rhs.v || w != rhs.w; }
};

class Graph {
 public:
  uint64_t n, m;
  sequence<Edge> edge;
  sequence<EdgeId> offset;
  bool weighted;
  bool symmetrized;

  Graph() = delete;
  Graph(bool _weighted = false, bool _symmetrized = false)
      : weighted(_weighted), symmetrized(_symmetrized) {}
  void generate_weight() {
    if (weighted) {
      fprintf(stderr, "Warning: Overwrite original weight\n");
    } else {
      weighted = true;
    }
    parallel_for(0, n, [&](size_t i) {
      for (size_t j = offset[i]; j < offset[i + 1]; j++) {
        edge[j].w = ((hash32(i) ^ hash32(edge[j].v)) & (WEIGHT - 1)) + 1;
      }
    });
  }
  void read_pbbs_format(char const* filename) {
    FILE* fp = fopen(filename, "r");
    if (fp == nullptr) {
      fprintf(stderr, "Error: file %s does not exist\n", filename);
      exit(EXIT_FAILURE);
    }
    fseek(fp, 0, SEEK_END);
    size_t size = ftell(fp);
    rewind(fp);
    sequence<char> buf(size);
    if (fread(buf.begin(), 1, size, fp) != size) {
      fprintf(stderr, "Error: Read failed\n");
      exit(EXIT_FAILURE);
    }
    sequence<bool> digit(size);
    auto idx = delayed_seq<size_t>(size, [](size_t i) { return i; });
    auto st = filter(idx, [&](size_t i) {
      return isdigit(buf[i]) && (i == 0 || !isdigit(buf[i - 1]));
    });
    auto ed = filter(idx, [&](size_t i) {
      return isdigit(buf[i]) && (i == size - 1 || !isdigit(buf[i + 1]));
    });
    assert(st.size() == ed.size());
    size_t num_sum = st.size();
    auto num = delayed_seq<size_t>(num_sum, [&](size_t i) {
      return stol(string(buf.begin() + st[i], buf.begin() + ed[i] + 1));
    });
    n = num[0], m = num[1];
    if (weighted) {
      assert(num.size() == n + m + m + 2);
    } else {
      assert(num.size() == n + m + 2);
    }
    offset = sequence<EdgeId>(n + 1);
    edge = sequence<Edge>(m);
    parallel_for(0, n, [&](size_t i) { offset[i] = num[i + 2]; });
    offset[n] = m;
    parallel_for(0, m, [&](size_t i) { edge[i].v = num[i + n + 2]; });
    if (weighted) {
      parallel_for(0, m, [&](size_t i) { edge[i].w = num[i + n + m + 2]; });
    }
    fclose(fp);
  }
  void read_gapbs_format(char const* filename) {
    ifstream ifs(filename);
    if (!ifs.is_open()) {
      fprintf(stderr, "Error: file %s does not exist\n", filename);
      exit(EXIT_FAILURE);
    }
    bool directed;
    ifs.read(reinterpret_cast<char*>(&directed), sizeof(bool));
    assert(directed == !symmetrized);
    ifs.read(reinterpret_cast<char*>(&m), sizeof(size_t));
    ifs.read(reinterpret_cast<char*>(&n), sizeof(size_t));
    offset = sequence<EdgeId>(n + 1);
    edge = sequence<Edge>(m);
    ifs.read(reinterpret_cast<char*>(offset.begin()), (n + 1) * sizeof(EdgeId));
    ifs.read(reinterpret_cast<char*>(edge.begin()), m * sizeof(Edge));
    if (directed) {
      sequence<EdgeId> inv_offset(n + 1);
      sequence<Edge> inv_edge(m);
      ifs.read(reinterpret_cast<char*>(inv_offset.begin()),
               (n + 1) * sizeof(EdgeId));
      ifs.read(reinterpret_cast<char*>(inv_edge.begin()), m * sizeof(Edge));
    }
    if (ifs.peek() != EOF) {
      fprintf(stderr, "Error: Bad data\n");
      exit(EXIT_FAILURE);
    }
    ifs.close();
  }
  void read_galois_format(char const* filename) {
    FILE* fp = fopen(filename, "r");
    if (fp == nullptr) {
      fprintf(stderr, "Error: file %s does not exist\n", filename);
      exit(EXIT_FAILURE);
    }
    fseek(fp, 0, SEEK_END);
    size_t size = ftell(fp);
    rewind(fp);
    vector<char> buf(size);
    if (fread(buf.data(), 1, size, fp) != size) {
      fprintf(stderr, "Error: Read failed\n");
      exit(EXIT_FAILURE);
    }
    uint64_t* fptr = (uint64_t*)buf.data();
    size_t version = *fptr++;
    size_t sizeof_edge_data = *fptr++;
    assert(version == 1);
    assert(sizeof_edge_data == sizeof(EdgeTy));
    n = *fptr++;
    m = *fptr++;
    offset = sequence<EdgeId>(n + 1);
    edge = sequence<Edge>(m);
    offset[0] = 0;
    for (size_t i = 1; i <= n; i++) {
      offset[i] = *fptr++;
    }
    uint32_t* fptr32 = (uint32_t*)fptr;
    for (size_t i = 0; i < m; i++) {
      edge[i].v = *fptr32++;
    }
    if (m % 2) fptr32++;  // padding
    for (size_t i = 0; i < m; i++) {
      edge[i].w = *fptr32++;
    }
    assert((void*)fptr32 == buf.data() + size);
    fclose(fp);
  }
  void read_binary_format(char const* filename) {
    // use mmap by default
    if (weighted == true) {
      fprintf(stderr, "Error: Binary format does not support weighted input\n");
      exit(EXIT_FAILURE);
    }
    struct stat sb;
    int fd = open(filename, O_RDONLY);
    if (fd == -1) {
      fprintf(stderr, "Error: Cannot open file %s\n", filename);
      exit(EXIT_FAILURE);
    }
    if (fstat(fd, &sb) == -1) {
      fprintf(stderr, "Error: Unable to acquire file stat\n");
      exit(EXIT_FAILURE);
    }
    char* data =
        static_cast<char*>(mmap(0, sb.st_size, PROT_READ, MAP_PRIVATE, fd, 0));
    size_t len = sb.st_size;
    n = reinterpret_cast<uint64_t*>(data)[0];
    m = reinterpret_cast<uint64_t*>(data)[1];
    size_t sizes = reinterpret_cast<uint64_t*>(data)[2];
    assert(sizes == (n + 1) * 8 + m * 4 + 3 * 8);
    this->n = n, this->m = m;
    offset = sequence<EdgeId>(n + 1);
    edge = sequence<Edge>(m);
    parallel_for(0, n + 1, [&](size_t i) {
      offset[i] = reinterpret_cast<uint64_t*>(data + 3 * 8)[i];
    });
    parallel_for(0, m, [&](size_t i) {
      edge[i].v = reinterpret_cast<uint32_t*>(data + 3 * 8 + (n + 1) * 8)[i];
    });

    if (data) {
      const void* b = data;
      munmap(const_cast<void*>(b), len);
    }
  }
  void read_graph(char const* filename) {
    size_t idx = string(filename).find_last_of('.');
    if (idx == string::npos) {
      fprintf(stderr, "Error: No file extension provided\n");
      exit(EXIT_FAILURE);
    }
    string subfix = string(filename).substr(idx + 1);
    if (subfix == "adj") {
      printf("Info: Reading pbbs format\n");
      read_pbbs_format(filename);
    } else if (subfix == "wsg") {
      printf("Info: Reading gapbs format\n");
      read_gapbs_format(filename);
    } else if (subfix == "gr") {
      printf("Info: Reading galois format\n");
      read_galois_format(filename);
    } else if (subfix == "bin") {
      read_binary_format(filename);
    } else {
      fprintf(stderr, "Error: Unrecognized file extension\n");
      exit(EXIT_FAILURE);
    }
  }
  void write_pbbs_format(char const* filename) {
    printf("Info: Writing pbbs format\n");
    FILE* fp = fopen(filename, "w");
    if (weighted) {
      fprintf(fp, "WeightedAdjacencyGraph\n");
    } else {
      fprintf(fp, "AdjacencyGraph\n");
    }
    fprintf(fp, "%zu\n", n);
    fprintf(fp, "%zu\n", m);
    for (size_t i = 0; i < n; i++) {
      fprintf(fp, "%" PRIu64 "\n", offset[i]);
    }
    for (size_t i = 0; i < m; i++) {
      fprintf(fp, "%" PRIu32 "\n", edge[i].v);
    }
    if (weighted) {
      for (size_t i = 0; i < m; i++) {
        fprintf(fp, "%" PRIu32 "\n", edge[i].w);
      }
    }
    fclose(fp);
  }
  void write_gapbs_format(char const* filename) {
    printf("Info: Writing gapbs format\n");
    sequence<EdgeId> inv_offset(n + 1);
    sequence<Edge> inv_edge(m);
    parallel_for(0, n + 1, [&](size_t i) { inv_offset[i] = 0; });
    parallel_for(0, n, [&](size_t i) {
      for (size_t j = offset[i]; j < offset[i + 1]; j++) {
        write_add(&inv_offset[edge[j].v], 1);
      }
    });
    scan_inplace(make_slice(inv_offset),
                 monoid([](size_t a, size_t b) { return a + b; }, 0));
    sequence<EdgeId> tmp_offset = inv_offset;
    parallel_for(0, n, [&](size_t i) {
      parallel_for(offset[i], offset[i + 1], [&](size_t j) {
        size_t pos = fetch_and_add(&tmp_offset[edge[j].v], 1);
        inv_edge[pos] = Edge(i, edge[j].w);
      });
    });
    parallel_for(0, n, [&](size_t i) {
      sort_inplace(inv_edge.cut(inv_offset[i], inv_offset[i + 1]),
                   [](Edge a, Edge b) { return a < b; });
    });
    ofstream ofs(filename);
    if (!ofs.is_open()) {
      fprintf(stderr, "Error: Open %s failed\n", filename);
      exit(EXIT_FAILURE);
    }
    bool directed = !symmetrized;
    ofs.write(reinterpret_cast<char*>(&directed), sizeof(bool));
    ofs.write(reinterpret_cast<char*>(&m), sizeof(size_t));
    ofs.write(reinterpret_cast<char*>(&n), sizeof(size_t));
    ofs.write(reinterpret_cast<char*>(offset.begin()),
              (n + 1) * sizeof(EdgeId));
    ofs.write(reinterpret_cast<char*>(edge.begin()), m * sizeof(Edge));
    if (directed) {
      ofs.write(reinterpret_cast<char*>(inv_offset.begin()),
                (n + 1) * sizeof(EdgeId));
      ofs.write(reinterpret_cast<char*>(inv_edge.begin()), m * sizeof(Edge));
    }
    ofs.close();
  }
  void write_galois_format(char const* filename) {
    printf("Info: Writing galois format\n");
    FILE* fp = fopen(filename, "w");
    uint64_t graph_version = 1, sizeof_edge_data = sizeof(EdgeTy);
    fwrite(&graph_version, sizeof(uint64_t), 1, fp);
    fwrite(&sizeof_edge_data, sizeof(uint64_t), 1, fp);
    fwrite(&n, sizeof(uint64_t), 1, fp);
    fwrite(&m, sizeof(uint64_t), 1, fp);
    for (size_t i = 1; i <= n; i++) {
      fwrite(&offset[i], sizeof(uint64_t), 1, fp);
    }
    for (size_t i = 0; i < m; i++) {
      fwrite(&edge[i].v, sizeof(uint32_t), 1, fp);
    }
    if (m % 2) {
      uint32_t padding = 0;
      fwrite(&padding, sizeof(uint32_t), 1, fp);
    }
    for (size_t i = 0; i < m; i++) {
      fwrite(&edge[i].w, sizeof(EdgeTy), 1, fp);
    }
    fclose(fp);
  }
  void check_order() {
    bool ordered = true;
    parallel_for(0, n, [&](size_t i) {
      parallel_for(offset[i], offset[i + 1], [&](size_t j) {
        if (j + 1 < offset[i + 1]) {
          if (edge[j].v >= edge[j + 1].v) {
            ordered = false;
          }
        }
      });
    });
    if (!ordered) {
      fprintf(stderr, "Warning: Graph is not ordered, reordering\n");
      parallel_for(0, n, [&](size_t i) {
        sort_inplace(edge.cut(offset[i], offset[i + 1]),
                  [](Edge a, Edge b) { return a < b; });
      });
      check_order();
    }
  }
  void check_symmetrized() {
    bool _symmetrized = true;
    parallel_for(0, n, [&](size_t i) {
      parallel_for(offset[i], offset[i + 1], [&](size_t j) {
        Edge es = Edge(i, edge[j].w);
        NodeId v = edge[j].v;
        if (*lower_bound(edge.begin() + offset[v], edge.begin() + offset[v + 1],
                         es) != es) {
          _symmetrized = false;
        }
      });
    });
    if (_symmetrized != symmetrized) {
      auto type = [](bool s) { return s ? "symmetrized" : "unsymmetrized"; };
      fprintf(stderr, "Warning: Graph is specified as %s but is %s\n",
              type(symmetrized), type(_symmetrized));
    }
  }
  void check_self_loop() {
    size_t self_loop = 0, duplicate_edge = 0;
    parallel_for(0, n, [&](size_t i) {
      size_t pre = ULONG_MAX;
      for (size_t j = offset[i]; j < offset[i + 1]; j++) {
        NodeId v = edge[j].v;
        if (i == v) {
          write_add(&self_loop, 1);
        }
        if (v == pre) {
          write_add(&duplicate_edge, 1);
        }
        pre = v;
      }
    });
    printf("num of self-loop: %zu\n", self_loop);
    printf("num of duplicate edge: %zu\n", duplicate_edge);
  }
  void degree_distribution() {
    constexpr int LOG2_MAX_DEG = 30;
    NodeId deg[LOG2_MAX_DEG] = {};
    parallel_for(0, n, [&](size_t i) {
      NodeId cur = ceil(log2(offset[i + 1] - offset[i] + 1));
      write_add(&deg[cur], 1);
    });
    printf("Degree distribution:\n");
    printf("degree between [%10d, %10d]: %u\n", 0, 0, deg[0]);
    for (int i = 1; i < LOG2_MAX_DEG; i++) {
      printf("degree between [%10d, %10d): %u\n", 1 << (i - 1), 1 << i, deg[i]);
    }
  }
  void weight_distribution() {
    constexpr int LOG2_MAX = 30;
    EdgeTy weight[LOG2_MAX] = {};
    parallel_for(0, m, [&](size_t i) {
      int v = ceil(log2(edge[i].v + 1));
      write_add(&weight[v], 1);
    });
    printf("Weight distribution:\n");
    printf("weight between [%10d, %10d]: %u\n", 0, 0, weight[0]);
    for (int i = 1; i < LOG2_MAX; i++) {
      printf("weight between [%10d, %10d): %u\n", 1 << (i - 1), 1 << i,
             weight[i]);
    }
  }
};
