#pragma once
#include <unordered_set>
#include "disjoint_cut_set.h"

namespace agl {
namespace cut_tree_internal {
namespace plain_gomory_hu {
class gomory_hu_tree_builder {
  void dfs(V v, V par = -1) {
    int dep = (par == -1) ? 0 : depth_[par] + 1;
    depth_[v] = dep;
    for (auto& to : edges_[v]) {
      if (std::get<0>(to) == par) continue;
      parent_cost_[std::get<0>(to)] = std::make_pair(v, std::get<1>(to));
      dfs(std::get<0>(to), v);
    }
  }

public:
  gomory_hu_tree_builder(int n) : n_(n), edges_(n), depth_(n), parent_cost_(n) {
    add_edge_count_ = 0;
  }

  void add_degree2_edge(V u, V v) {
    add_edge_count_++;
    degree2_edges_.emplace_back(u, v);
  }

  void contraction(V s, V t, V sside_new_vtx, V tside_new_vtx) {
    edges_.resize(edges_.size() + 2);
    CHECK(std::get<0>(edges_[s].back()) == t);
    CHECK(std::get<0>(edges_[t].back()) == s);
    int f = std::get<1>(edges_[s].back());

    std::get<0>(edges_[s].back()) = sside_new_vtx;
    std::get<2>(edges_[s].back()) = 0;
    edges_[sside_new_vtx].emplace_back(s, f, edges_[s].size() - 1);

    std::get<0>(edges_[t].back()) = tside_new_vtx;
    std::get<2>(edges_[t].back()) = 0;
    edges_[tside_new_vtx].emplace_back(t, f, edges_[t].size() - 1);

    edges_[sside_new_vtx].emplace_back(tside_new_vtx, f, 1);
    edges_[tside_new_vtx].emplace_back(sside_new_vtx, f, 1);
  }

  void add_edge(V u, V v, int cost, const std::vector<V>& vs, const disjoint_cut_set& dcs) {
    CHECK(u != v);
    add_edge_count_++;
    for (V w : vs) {
      if (!dcs.is_same_group(u, w)) continue;
      for (size_t i = 0; i < edges_[w].size(); ++i) {
        V t = std::get<0>(edges_[w][i]);
        int r = std::get<2>(edges_[w][i]);
        edges_[v].emplace_back(edges_[w][i]);
        std::get<0>(edges_[t][r]) = v;
        std::get<2>(edges_[t][r]) = edges_[v].size() - 1;
      }
      edges_[w].clear();
    }
    for (V w : vs) {
      if (w == u) continue;
      for (size_t i = 0; i < edges_[w].size(); ++i) {
        V t = std::get<0>(edges_[w][i]);
        if (!dcs.is_same_group(t, v)) continue;
        int r = std::get<2>(edges_[w][i]);
        edges_[t][r] = edges_[t].back();
        std::get<2>(edges_[std::get<0>(edges_[t][r])][std::get<2>(edges_[t][r])]) = r;
        edges_[t].pop_back();
        std::get<0>(edges_[w][i]) = u;
        std::get<2>(edges_[w][i]) = edges_[u].size();
        edges_[u].emplace_back(w, std::get<1>(edges_[w][i]), i);
      }
    }
    edges_[u].emplace_back(v, cost, edges_[v].size());
    edges_[v].emplace_back(u, cost, edges_[u].size() - 1);
  }

  void build() {
    for (auto e : degree2_edges_) {
      edges_[e.first].emplace_back(e.second, 2, edges_[e.second].size());
      edges_[e.second].emplace_back(e.first, 2, edges_[e.first].size() - 1);
    }
    degree2_edges_.clear();
    for (V v = n_; v < (V)(edges_.size()); ++v) {
      CHECK(edges_[v].size() == 2);
      for (auto& e : edges_[std::get<0>(edges_[v][0])]) {
        if (std::get<0>(e) == v) {
          std::get<0>(e) = std::get<0>(edges_[v][1]);
        }
      }
      for (auto& e : edges_[std::get<0>(edges_[v][1])]) {
        if (std::get<0>(e) == v) {
          std::get<0>(e) = std::get<0>(edges_[v][0]);
        }
      }
    }
    CHECK(add_edge_count_ == n_ - 1);
    parent_cost_[0] = std::make_pair(-1, 0);
    dfs(0);
    edges_.clear(); edges_.shrink_to_fit();
  }

  int query(V u, V v) const {
    CHECK(u != v);
    CHECK(u < n_ && v < n_);
    int ans = std::numeric_limits<int>::max();
    while (u != v) {
      if (depth_[u] > depth_[v]) {
        ans = std::min(ans, parent_cost_[u].second);
        u = parent_cost_[u].first;
      } else {
        ans = std::min(ans, parent_cost_[v].second);
        v = parent_cost_[v].first;
      }
    }
    return ans;
  }

  const std::vector<std::pair<V, int>>& parent_weight() const {
    return parent_cost_;
  }

  int debug_add_edge_count() const {
    return add_edge_count_;
  }

  void test(const G& g) {
    if (n_ == 1) {
      return;
    }
    std::queue<V> que;
    std::vector<std::vector<V>> edge(n_), children(n_);
    std::vector<V> cnt(n_);
    std::vector<std::unordered_set<V>> d(n_);
    for (V v : irange<V>(n_)) {
      for (auto e : g.neighbors(v)) {
        edge[v].emplace_back(to(e));
        edge[to(e)].emplace_back(v);
      }
    }
    for (V v : irange<V>(n_)) {
      if (parent_cost_[v].first != -1) {
        children[parent_cost_[v].first].emplace_back(v);
        ++cnt[parent_cost_[v].first];
      }
    }
    for (V v : irange<V>(n_)) {
      if (cnt[v] == 0) {
        que.emplace(v);
      }
    }
    while (!que.empty()) {
      V v = que.front();
      que.pop();
      for (V u : children[v]) {
        for (V t : d[u]) {
          d[v].emplace(t);
        }
      }
      d[v].emplace(v);
      V cut = 0;
      for (V u : d[v]) {
        for (V t : edge[u]) {
          if (d[v].count(t) == 0) {
            ++cut;
          }
        }
      }
      CHECK(cut == parent_cost_[v].second);
      if (parent_cost_[v].first >= 0) {
        --cnt[parent_cost_[v].first];
        if (cnt[parent_cost_[v].first] == 0) {
          que.emplace(parent_cost_[v].first);
        }
      }
    }
  }

private:
  int add_edge_count_;
  int n_;
  std::vector<std::vector<std::tuple<V, int, int>>> edges_;
  std::vector<int> depth_;
  std::vector<std::pair<V, int>> parent_cost_;
  std::vector<std::pair<V, V>> degree2_edges_;
};
} // namespace plain_gomory_hu
} // namespace cut_tree_internal
} // namespace agl
