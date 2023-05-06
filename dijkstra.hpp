#pragma once
#include <queue>

#include "graph.hpp"
#include "parlay/internal/get_time.h"
void dijkstra(size_t s, const Graph &G, EdgeTy *dist) {
  fill(dist, dist + G.n, INT_MAX / 2);
  dist[s] = 0;
  priority_queue<pair<EdgeTy, NodeId>, vector<pair<EdgeTy, NodeId>>,
                 greater<pair<EdgeTy, NodeId>>>
      pq;
  pq.push(make_pair(dist[s], s));
  while (!pq.empty()) {
    pair<EdgeTy, NodeId> dist_and_node = pq.top();
    pq.pop();
    EdgeTy d = dist_and_node.first;
    NodeId u = dist_and_node.second;
    if (dist[u] < d) continue;
    for (size_t j = G.offset[u]; j < G.offset[u + 1]; j++) {
      NodeId v = G.edge[j].v;
      EdgeTy w = G.edge[j].w;
      if (dist[v] > dist[u] + w) {
        dist[v] = dist[u] + w;
        pq.push(make_pair(dist[v], v));
      }
    }
  }
}

void verifier(size_t s, const Graph &G, EdgeTy *ch_dist) {
  EdgeTy *cor_dist = new EdgeTy[G.n];
  internal::timer tm;
  dijkstra(s, G, cor_dist);
  tm.stop();
  printf("dijkstra running time: %-10f\n", tm.total_time());
  parallel_for(0, G.n, [&](size_t i) {
    if (cor_dist[i] != ch_dist[i]) {
      printf("dijkstra_dist[%zu]=%d, my_dist[%zu]=%d\n", i, cor_dist[i], i,
             ch_dist[i]);
    }
    assert(cor_dist[i] == ch_dist[i]);
  });
  delete[] cor_dist;
}
