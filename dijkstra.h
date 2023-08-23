#pragma once
#include <queue>

#include "graph.h"
#include "parlay/internal/get_time.h"

sequence<EdgeTy> dijkstra(size_t s, const Graph &G) {
  sequence<EdgeTy> dist(G.n, DIST_MAX);
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
  return dist;
}

void verifier(size_t s, const Graph &G, const sequence<EdgeTy> &act_dist) {
  internal::timer tm;
  auto exp_dist = dijkstra(s, G);
  tm.stop();
  printf("dijkstra running time: %-10f\n", tm.total_time());
  parallel_for(0, G.n, [&](size_t i) {
    if (exp_dist[i] != act_dist[i]) {
      string out =
          "exp_dist[" + to_string(i) + "]: " + to_string(exp_dist[i]) + " vs. ";
      out += "act_dist[" + to_string(i) + "]: " + to_string(act_dist[i]) + "\n";
      cout << out;
    }
    assert(exp_dist[i] == act_dist[i]);
  });
}
