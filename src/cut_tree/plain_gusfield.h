#pragma once
#include <base/base.h>
#include <graph/graph.h>
#include <vector>
#include <queue>
#include <iostream>
#include <limits>

namespace agl {
namespace cut_tree_internal {
template<class max_flow_t>
class plain_gusfield {
  int query_dfs(V v, V t, int cost, V par = -1) const {
    if (v == t) return cost;
    for (const auto& to_cost : binary_tree_edges_[v]) {
      V to; int edge_cost; std::tie(to, edge_cost) = to_cost;
      if (to == par) continue;
      int ncost = query_dfs(to, t, std::min(cost, edge_cost), v);
      if (ncost != -1) return ncost;
    }

    return -1;
  }

public:

  void add_edge(V s, V t, int cost) {
    binary_tree_edges_[s].emplace_back(t, cost);
    binary_tree_edges_[t].emplace_back(s, cost);
  }

  plain_gusfield(G& g) : num_vertices_(g.num_vertices()), binary_tree_edges_(g.num_vertices()) {
    union_find uf(num_vertices_);

    for (int v = 0; v < num_vertices_; v++) for (auto& e : g.edges(v)) {
      uf.unite(v, to(e));
    }
    std::vector<int> p(num_vertices_);
    for (int v = 0; v < num_vertices_; v++) {
      if (v == uf.root(v)) p[v] = -1;
      else p[v] = uf.root(v);
    }
    std::vector<int> root_vtxs;
    for (int v = 0; v < num_vertices_; v++) if (p[v] == -1) root_vtxs.push_back(v);
    for (int i = 0; i < int(root_vtxs.size()) - 1; i++) {
      add_edge(root_vtxs[i], root_vtxs[i + 1], 0);
    }

    max_flow_t dc(g);

    for (int s = 0; s < num_vertices_; s++) {
      if (p[s] == -1) continue;
      V t = p[s];

      int cost = dc.max_flow(s, t);
      // fprintf(stderr, "(%d,%d) cost = %d\n", s, t, cost);
      add_edge(s, t, cost);


      std::vector<char> used(num_vertices_);
      std::queue<int> q;
      q.push(s);
      used[s] = true;
      while (!q.empty()) {
        V v = q.front(); q.pop();
        for (auto& e : dc.edges(v)) {
          if (dc.cap(e) == 0 || used[dc.to(e)]) continue;
          used[dc.to(e)] = true;
          q.push(dc.to(e));
          if (p[dc.to(e)] == t) p[dc.to(e)] = s;
        }
      }
    }
  }

  int query(V u, V v) const {
    CHECK(u != v);
    CHECK(u < num_vertices_ && v < num_vertices_);
    int ans = query_dfs(u, v, std::numeric_limits<int>::max());
    if (ans == -1) {
      return 0; // 到達できなかった
    }
    return ans;
  }

  void print_gomory_hu_tree_dfs(V v, V par, std::ostream& os) {
    for (auto& to : binary_tree_edges_[v]) {
      if (to.first == par) continue;
      os << v << " " << to.first << " " << to.second << std::endl;
      print_gomory_hu_tree_dfs(to.first, v, os);
    }
  }

  void print_gomory_hu_tree(std::ostream& os) {
    if(num_vertices_ > 0) print_gomory_hu_tree_dfs(0, -1, os);
  }

private:
  const int num_vertices_;
  std::vector<std::vector<std::pair<V, int>>> binary_tree_edges_;
};
} // namespace cut_tree_internal
} // namespace agl